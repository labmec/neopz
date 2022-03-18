/**
 * @file
 * @brief Contains implementations of the TPZMGAnalysis methods.
 */

#include <stdlib.h>
#include "pzmganalysis.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzgeoel.h"
#include "pztransfer.h"
#include "pzadmchunk.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "pzskylmat.h"
#include "pzskylstrmatrix.h"

#include "TPZFrontSym.h"
#include "TPZFrontNonSym.h"
#include "TPZFrontStructMatrix.h"
#include "pzmgsolver.h"
#include "pzseqsolver.h"
#include "pzstepsolver.h"

#include "pzquad.h"
#include "TPZMaterial.h"
#include "pztransfer.h"

using namespace std;

TPZMGAnalysis::TPZMGAnalysis(TPZCompMesh *cmesh) : TPZLinearAnalysis(cmesh) {
	fMeshes.Push(cmesh);
}

TPZMGAnalysis::~TPZMGAnalysis() {
	while (fMeshes.NElements()) delete fMeshes.Pop();
	while(fSolutions.NElements()) delete fSolutions.Pop();
	while(fSolvers.NElements()) delete fSolvers.Pop();
	while(fPrecondition.NElements()) delete fPrecondition.Pop();
}

void TPZMGAnalysis::AppendMesh(TPZCompMesh * mesh){
	if(fMeshes.NElements() != fSolvers.NElements() || fMeshes.NElements() != fSolutions.NElements() ||
	   fPrecondition.NElements() != fMeshes.NElements()) {
		cout << "TPZMGAnalysis::AppendMesh can only be called after solving the coarse mesh\n";
		return;
	}
	fMeshes.Push(mesh);
	TPZCompMesh *coarse = fCompMesh;
	
	fCompMesh = mesh;
	OptimizeBandwidth();
	TPZSkylineStructMatrix<STATE> skstr(mesh);
	SetStructuralMatrix(skstr);
	int nmeshes = fSolvers.NElements();
	TPZTransfer<STATE> *tr = new TPZTransfer<STATE>;
	TPZAutoPointer<TPZMatrix<STATE> > trauto(tr);
	mesh->BuildTransferMatrix(*coarse,*tr);
	TPZFMatrix<STATE> sol(mesh->NEquations(),1);
	tr->TransferSolution(fSolution,sol);
	fSolution = sol;
	//TODOCOMPLEX
	mesh->LoadSolution(fSolution);
	TPZBlockDiagonalStructMatrix<STATE> bdstr(mesh);
	TPZBlockDiagonal<STATE> *bd = (TPZBlockDiagonal<STATE> *) bdstr.Create();
	bdstr.AssembleBlockDiagonal(*bd);
	TPZAutoPointer<TPZMatrix<STATE> > bdauto(bd);
	TPZStepSolver<STATE> s1(bdauto);
	s1.SetDirect(ELU);
	int nvar = mesh->MaterialVec().begin()->second->NStateVariables();
	
	TPZMatrixSolver<STATE> *prec;
	prec = fPrecondition[nmeshes-1];
	if(!prec) prec = fSolvers[nmeshes-1];
	TPZAutoPointer<TPZGuiInterface> guiInterface = new TPZGuiInterface;
	TPZAutoPointer<TPZMatrix<STATE> > skauto =
		dynamic_cast<TPZMatrix<STATE>*>(skstr.CreateAssemble(fRhs, guiInterface));
	TPZMGSolver<STATE> s2(trauto,*prec,nvar,skauto);
	TPZSequenceSolver<STATE> s3;
	s3.ShareMatrix(s2);
	s3.AppendSolver(s1);
	s3.AppendSolver(s2);
	s3.AppendSolver(s1);
	
	TPZStepSolver<STATE> s4;
	s4.ShareMatrix(s3);
	fPrecondition.Push((TPZMatrixSolver<STATE> *)s3.Clone());
	s4.SetCG(200,s3,1.e-6,1);
	SetSolver(s4);
	fSolvers.Push((TPZMatrixSolver<STATE> *) s4.Clone());
}

TPZCompMesh *TPZMGAnalysis::PopMesh() {
	if(fMeshes.NElements() == 1) {
		cout << "TPZMGAnalysis cannot delete the root mesh, sorry\n";
		return 0;
	}
	if(fSolutions.NElements() == fMeshes.NElements()) delete fSolutions.Pop();
	
	delete fSolvers.Pop();
	delete fPrecondition.Pop();
	
	SetSolver(*fSolvers[fSolvers.NElements()-1]);
	fCompMesh = fMeshes[fMeshes.NElements()-2];
	fSolution = *fSolutions[fSolutions.NElements()-1];
	return fMeshes.Pop();
}

void TPZMGAnalysis::MeshError(TPZCompMesh *fine, TPZCompMesh *coarse, TPZVec<REAL> &ervec,
							  void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),TPZVec<REAL> &truervec){
	coarse->Reference()->ResetReference();
	coarse->LoadReferences();
	int dim = fine->MaterialVec().begin()->second->Dimension();
	ervec.Resize(coarse->NElements());
	if(f) {
		truervec.Resize(coarse->NElements());
		truervec.Fill(0.,0);
	}
	ervec.Fill(0.,0);
	TPZCompEl *cel;
	TPZAdmChunkVector<TPZCompEl *> &elementvec = fine->ElementVec();
	int numel = elementvec.NElements();
	int el;
	for(el=0; el<numel; el++) {
		cel = elementvec[el];
		if (!cel) continue;
		TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
		if(!cint) continue;
		int ncon = cint->NConnects();
		TPZGeoElSide gelside(cint->Reference(),ncon-1);
		if(gelside.Dimension() != dim) continue;
		TPZGeoElSide gellarge(gelside);
		while(!gellarge.Reference().Exists() && gellarge.Father2().Exists()) gellarge = gellarge.Father2();
		if(!gellarge.Reference().Exists()) {
			cout << "TPZMGAnalsysis::BuildTransferMatrix element " << el << " found no corresponding element\n";
			continue;
		}
		TPZCompElSide cellargeside = gellarge.Reference();
		TPZCompEl *cellarge = cellargeside.Element();
		TPZInterpolatedElement *cintlarge = (TPZInterpolatedElement *) cellarge;
		TPZTransform<> transform(gelside.Dimension(),gellarge.Dimension());
		gelside.SideTransform3(gellarge,transform);
		int index = cellarge->Index();
		REAL truerror = 0.;
		ervec[index] += ElementError(cint,cintlarge,transform,f,truerror);
		if(f) truervec[index]  += truerror;
	}
}


REAL TPZMGAnalysis::ElementError(TPZInterpolatedElement *fine, TPZInterpolatedElement *coarse, TPZTransform<> &tr,
								 void (*f)(const TPZVec<REAL> &loc, TPZVec<STATE> &val, TPZFMatrix<STATE> &deriv),REAL &truerror){
	// accumulates the transfer coefficients between the current element and the
	// coarse element into the transfer matrix, using the transformation t
	int locnod = fine->NConnects();
	int cornod = coarse->NConnects();
	int locmatsize = fine->NShapeF();
	int cormatsize = coarse->NShapeF();
	REAL error = 0.;
	truerror = 0.;
	
	STATE loclocmatstore[500] = {0.},loccormatstore[500] = {0.};
	TPZFMatrix<STATE> loclocmat(locmatsize,locmatsize,loclocmatstore,500);
	TPZFMatrix<STATE> loccormat(locmatsize,cormatsize,loccormatstore,500);
	
	TPZAutoPointer<TPZIntPoints> intrule = fine->GetIntegrationRule().Clone();
	int dimension = fine->Dimension();
	int numdof = fine->Material()->NStateVariables();
	TPZBlock &locblock = fine->Mesh()->Block();
	TPZFMatrix<STATE> &locsolmesh = fine->Mesh()->Solution();
	
	TPZBlock &corblock = coarse->Mesh()->Block();
	TPZFMatrix<STATE> &corsolmesh = coarse->Mesh()->Solution();
	
	TPZVec<STATE> locsol(numdof);
	TPZFMatrix<STATE> locdsol(dimension,numdof);
	
	TPZVec<STATE> corsol(numdof);
	TPZFMatrix<STATE> cordsol(dimension,numdof);
	
	TPZManVector<int> prevorder(dimension),
    order(dimension);
	intrule->GetOrder(prevorder);
	
	TPZManVector<int> interpolation(dimension);
	fine->GetInterpolationOrder(interpolation);
	
	// compute the interpolation order of the shapefunctions squared
	int dim;
	int maxorder = interpolation[0];
	for(dim=0; dim<interpolation.NElements(); dim++) {
		maxorder = interpolation[dim] < maxorder ? maxorder : interpolation[dim];
	}
	for(dim=0; dim<dimension; dim++) {
		order[dim] = 20;
	}
	intrule->SetOrder(order);
	
	
	REAL locphistore[50]={0.},locdphistore[150]={0.};
	TPZFMatrix<REAL> locphi(locmatsize,1,locphistore,50);
	TPZFMatrix<REAL> locdphi(dimension,locmatsize,locdphistore,150),locdphix(dimension,locmatsize);
	// derivative of the shape function
	// in the master domain
	
	REAL corphistore[50]={0.},cordphistore[150]={0.};
	TPZFMatrix<REAL> corphi(cormatsize,1,corphistore,50);
	TPZFMatrix<REAL> cordphi(dimension,cormatsize,cordphistore,150), cordphix(dimension,cormatsize);
	// derivative of the shape function
	// in the master domain
	
	REAL jacobianstore[9],
    axesstore[9];
	TPZManVector<REAL> int_point(dimension),
    coarse_int_point(dimension);
	TPZFMatrix<REAL> jacfine(dimension,dimension,jacobianstore,9),jacinvfine(dimension,dimension);
	TPZFMatrix<REAL> axesfine(3,3,axesstore,9);
	TPZManVector<REAL> xfine(3);
	TPZFMatrix<REAL> jaccoarse(dimension,dimension),jacinvcoarse(dimension,dimension);
	TPZFMatrix<REAL> axescoarse(3,3), axesinner(dimension,dimension);
	TPZManVector<REAL> xcoarse(3);
	
	REAL jacdetcoarse;
	int numintpoints = intrule->NPoints();
	REAL weight;
	int i,j,k;
	
	TPZVec<STATE> truesol(numdof);
	TPZFMatrix<STATE> truedsol(dimension,numdof);
	for(int int_ind = 0; int_ind < numintpoints; ++int_ind) {
		
		intrule->Point(int_ind,int_point,weight);
		REAL jacdetfine;
		fine->Reference()->Jacobian( int_point, jacfine , axesfine, jacdetfine, jacinvfine);
		fine->Reference()->X(int_point, xfine);
		if(f) f(xfine,truesol,truedsol);
		fine->Shape(int_point,locphi,locdphi);
		tr.Apply(int_point,coarse_int_point);
		coarse->Shape(coarse_int_point,corphi,cordphi);
		coarse->Reference()->Jacobian( coarse_int_point,jaccoarse,axescoarse, jacdetcoarse, jacinvcoarse);
		coarse->Reference()->X(coarse_int_point,xcoarse);
		REAL dist = (xfine[0]-xcoarse[0])*(xfine[0]-xcoarse[0])+(xfine[1]-xcoarse[1])*(xfine[1]-xcoarse[1])+(xfine[2]-xcoarse[2])*(xfine[2]-xcoarse[2]);
		if(dist > 1.e-6) cout << "TPZMGAnalysis::ElementError transformation between fine and coarse is wrong\n";
		for(i=0; i<dimension; i++) {
			for(j=0; j<dimension; j++) {
				axesinner(i,j) = 0.;
				for(k=0; k<3; k++) axesinner(i,j) += axesfine(i,k)*axescoarse(j,k);
			}
		}
		if(fabs(axesinner(0,0)-1.) > 1.e-6 || fabs(axesinner(1,1)-1.) > 1.e-6 || fabs(axesinner(0,1)) > 1.e-6 || fabs(axesinner(1,0)) > 1.e-6) {
			cout << "TPZMGAnalysis axesinner is not identify?\n";
		}
		weight *= fabs(jacdetfine);
		
		locdphix.Zero();
		
		int ieq,d;
		switch(dim) {
			case 0:
				break;
			case 1:
				for(d=0; d<dimension; d++) {
					for(ieq=0; ieq<locmatsize; ieq++) locdphix(d,ieq) = locdphi(d,ieq)*(1./jacdetfine);
					for(ieq=0; ieq<cormatsize; ieq++) cordphix(d,ieq) = cordphi(d,ieq)*(axesinner(0,0)/jacdetcoarse);
				}
				break;
			case 2:
				for(ieq = 0; ieq < locmatsize; ieq++) {
					locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(1,0)*locdphi(1,ieq);
					locdphix(1,ieq) = jacinvfine(0,1)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq);
					REAL tmp[2];
					tmp[0] = locdphix(0,ieq)*axesfine(0,0)+locdphix(1,ieq)*axesfine(1,0);
					tmp[1] = locdphix(0,ieq)*axesfine(0,1)+locdphix(1,ieq)*axesfine(1,1);
					locdphix(0,ieq) = tmp[0];
					locdphix(1,ieq) = tmp[1];
				}
				for(ieq = 0; ieq < cormatsize; ieq++) {
					cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(1,0)*cordphi(1,ieq);
					cordphix(1,ieq) = jacinvcoarse(0,1)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq);
					REAL tmp[2];
					tmp[0] = cordphix(0,ieq)*axescoarse(0,0)+cordphix(1,ieq)*axescoarse(1,0);
					tmp[1] = cordphix(0,ieq)*axescoarse(0,1)+cordphix(1,ieq)*axescoarse(1,1);
					cordphix(0,ieq) = tmp[0];
					cordphix(1,ieq) = tmp[1];
				}
				break;
			case 3:
				for(ieq = 0; ieq < locmatsize; ieq++) {
					locdphix(0,ieq) = jacinvfine(0,0)*locdphi(0,ieq) + jacinvfine(0,1)*locdphi(1,ieq) + jacinvfine(0,2)*locdphi(2,ieq);
					locdphix(1,ieq) = jacinvfine(1,0)*locdphi(0,ieq) + jacinvfine(1,1)*locdphi(1,ieq) + jacinvfine(1,2)*locdphi(2,ieq);
					locdphix(2,ieq) = jacinvfine(2,0)*locdphi(0,ieq) + jacinvfine(2,1)*locdphi(1,ieq) + jacinvfine(2,2)*locdphi(2,ieq);
				}
				for(ieq = 0; ieq < cormatsize; ieq++) {
					cordphix(0,ieq) = jacinvcoarse(0,0)*cordphi(0,ieq) + jacinvcoarse(0,1)*cordphi(1,ieq) + jacinvcoarse(0,2)*cordphi(2,ieq);
					cordphix(1,ieq) = jacinvcoarse(1,0)*cordphi(0,ieq) + jacinvcoarse(1,1)*cordphi(1,ieq) + jacinvcoarse(1,2)*cordphi(2,ieq);
					cordphix(2,ieq) = jacinvcoarse(2,0)*cordphi(0,ieq) + jacinvcoarse(2,1)*cordphi(1,ieq) + jacinvcoarse(2,2)*cordphi(2,ieq);
					REAL tmp[3];
					tmp[0] = cordphix(0,ieq)*axesinner(0,0)+cordphix(1,ieq)*axesinner(0,1)+cordphix(2,ieq)*axesinner(0,2);
					tmp[1] = cordphix(0,ieq)*axesinner(1,0)+cordphix(1,ieq)*axesinner(1,1)+cordphix(2,ieq)*axesinner(1,2);
					tmp[2] = cordphix(0,ieq)*axesinner(2,0)+cordphix(1,ieq)*axesinner(2,1)+cordphix(2,ieq)*axesinner(2,2);
					cordphix(0,ieq) = tmp[0];
					cordphix(1,ieq) = tmp[1];
					cordphix(2,ieq) = tmp[2];
				}
				break;
			default:
				PZError << "pzintel.c please implement the " << dim << "d Jacobian and inverse\n";
				PZError.flush();
		}
		int iv=0;
		locsol.Fill(0.);
		locdsol.Zero();
		iv=0;
		int in;
		for(in=0; in<locnod; in++) {
			TPZConnect *df = &fine->Connect(in);
			int dfseq = df->SequenceNumber();
			int dfvar = locblock.Size(dfseq);
			int pos = locblock.Position(dfseq);
			
			for(int jn=0; jn<dfvar; jn++) {
				locsol[iv%numdof] += (STATE)locphi(iv/numdof,0)*locsolmesh(pos+jn,0);
				for(d=0; d<dim; d++)
					locdsol(d,iv%numdof) += (STATE)locdphix(d,iv/numdof)*locsolmesh(pos+jn,0);
				iv++;
			}
		}
		corsol.Fill(0.);
		cordsol.Zero();
		iv=0;
		for(in=0; in<cornod; in++) {
			TPZConnect *df = &coarse->Connect(in);
			int dfseq = df->SequenceNumber();
			int dfvar = corblock.Size(dfseq);
			int pos = corblock.Position(dfseq);
			for(int jn=0; jn<dfvar; jn++) {
				corsol[iv%numdof] += (STATE)corphi(iv/numdof,0)*corsolmesh(pos+jn,0);
				for(d=0; d<dim; d++)
					cordsol(d,iv%numdof) += (STATE)cordphix(d,iv/numdof)*corsolmesh(pos+jn,0);
				iv++;
			}
		}
		int jn;
		for(jn=0; jn<numdof; jn++) {
			//       error += (locsol[jn]-corsol[jn])*(locsol[jn]-corsol[jn])*weight;
			//       if(f) truerror += (corsol[jn]-truesol[jn])*(corsol[jn]-truesol[jn])*weight;
			for(d=0; d<dim; d++) {
				error += fabs((locdsol(d,jn)-cordsol(d,jn))*(locdsol(d,jn)-cordsol(d,jn))*(STATE)weight);
				if(f) truerror += fabs((cordsol(d,jn)-truedsol(d,jn))*(cordsol(d,jn)-truedsol(d,jn))*(STATE)weight);
			}
		}
	}
	intrule->SetOrder(prevorder);
	return error;
}


template<class TVar>
void TPZMGAnalysis::SolveInternal()
{
  if(fMeshes.NElements() == 1) {
		TPZLinearAnalysis::Solve();
		if(fSolvers.NElements() == 0) {
			fSolvers.Push((TPZMatrixSolver<TVar> *) fSolver->Clone());
		}
		if(fPrecondition.NElements() == 0) {
			fPrecondition.Push(0);
		}
		if(fSolutions.NElements() == 0) {
			fSolutions.Push(new TPZSolutionMatrix(fSolution));
		} else {
			int nsol = fSolutions.NElements();
			*(fSolutions[nsol-1]) = fSolution;
		}
		return;
	}
	int numeq = fCompMesh->NEquations();
	if(fRhs.Rows() != numeq ) return;
	int nsolvers = fSolvers.NElements();
	
	TPZFMatrix<TVar> residual(fRhs);
	TPZFMatrix<TVar> delu(numeq,1,0.);
	TPZMatrixSolver<TVar> *solve = dynamic_cast<TPZMatrixSolver<TVar> *> (fSolvers[nsolvers-1]);
	if(fSolution.Rows() != numeq) {
		fSolution.Redim(numeq,1);
	} else {
		solve->Matrix()->Residual(fSolution,fRhs,residual);
	}
	
	REAL normrhs = Norm(fRhs);
	REAL normres  = Norm(residual);
	if(normrhs*1.e-6 >= normres) {
		cout << "TPZMGAnalysis::Solve no need for iterations normrhs = " << normrhs << " normres = " << normres << endl;
		if(fSolutions.NElements() < fMeshes.NElements()) {
			fSolutions.Push(new TPZSolutionMatrix(fSolution));
		} else {
			int nsol = fSolutions.NElements();
			*(fSolutions[nsol-1]) = fSolution;
		}
		return ;
	}
	
	TPZStepSolver<TVar> *stepsolve = dynamic_cast<TPZStepSolver<TVar> *> (solve);
	if(stepsolve) stepsolve->SetTolerance(1.e-6*normrhs/normres);
	cout << "TPZMGAnalysis::Run res : " << Norm(residual) << " neq " << numeq << endl;
	solve->Solve(residual, delu);
	TPZFMatrix<TVar> &mysol = fSolution;
	mysol += delu;

	fCompMesh->LoadSolution(fSolution);
	if(fSolutions.NElements() < fMeshes.NElements()) {
		fSolutions.Push(new TPZSolutionMatrix(fSolution));
	} else {
		int nsol = fSolutions.NElements();
		*(fSolutions[nsol-1]) = fSolution;
	}
}
void TPZMGAnalysis::Solve() {
	return SolveInternal<STATE>();
}

TPZCompMesh  *TPZMGAnalysis::UniformlyRefineMesh(TPZCompMesh *mesh, bool withP) {
	
	TPZGeoMesh *gmesh = mesh->Reference();
	if(!gmesh) {
		cout << "TPZMGAnalysis::UniformlyRefineMesh encountered no geometric mesh\n";
		return 0;
	}
	gmesh->ResetReference();
	TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
	mesh->CopyMaterials(*cmesh);
	TPZAdmChunkVector<TPZCompEl *> &elementvec = mesh->ElementVec();
	int el,nelem = elementvec.NElements();
	for(el=0; el<nelem; el++) {
		TPZCompEl *cel = elementvec[el];
		if(!cel) continue;
		TPZInterpolatedElement *cint = dynamic_cast<TPZInterpolatedElement *> (cel);
		if(!cint) {
			cout << "TPZMGAnalysis::UniformlyRefineMesh encountered a non interpolated element\n";
			continue;
		}
		int ncon = cint->NConnects();
		int porder = cint->PreferredSideOrder(ncon-1);
		
		TPZGeoEl *gel = cint->Reference();
		if(!gel) {
			cout << "TPZMGAnalysis::UniformlyRefineMesh encountered an element without geometric reference\n";
			continue;
		}
		TPZVec<TPZGeoEl *> sub;
		gel->Divide(sub);
		int nsub = sub.NElements();
		int isub;
		for(isub=0; isub<nsub; isub++) {
			TPZInterpolatedElement *csint;
			csint = (TPZInterpolatedElement *) cmesh->CreateCompEl(sub[isub]);
			if(withP) csint->PRefine(porder+1);
			else csint->PRefine(porder);
		}
	}
	return cmesh;
}

void TPZMGAnalysis::ComputeError(TPZVec<REAL> &error) {
	int nmesh = fMeshes.NElements();
	int nsol = fSolutions.NElements();
	if(nsol != nmesh || nsol < 2) {
		cout << "TPZMGAnalysis cannot compute the errors\n";
		return;
	}
	fMeshes[nsol-2]->LoadSolution(*fSolutions[nsol-2]);
	fMeshes[nsol-1]->LoadSolution(*fSolutions[nsol-1]);
	TPZVec<REAL> truerror;
	MeshError(fMeshes[nsol-1],fMeshes[nsol-2],error,0,truerror);
}
