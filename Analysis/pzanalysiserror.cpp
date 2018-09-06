/**
 * @file
 * @brief Contains implementations of the TPZAnalysisError methods.
 */

#include "pzanalysiserror.h"
#include <stdlib.h>           // for exit
#include <cmath>              // for log, pow, sqrt, fabs
#include <iostream>           // for cout
#include "pzadmchunk.h"       // for TPZAdmChunkVector
#include "pzblock.h"          // for TPZBlock
#include "TPZChunkVector.h"          // for TPZChunkVector
#include "pzcmesh.h"          // for TPZCompMesh
#include "pzconnect.h"        // for TPZConnect
#include "pzerror.h"          // for PZError
#include "pzgeoel.h"          // for TPZGeoEl
#include "pzgeoelside.h"      // for TPZGeoElSide
#include "pzgmesh.h"          // for TPZGeoMesh
#include "pzgnode.h"          // for TPZGeoNode
#include "pzintel.h"          // for TPZInterpolatedElement
#include "TPZMaterial.h"       // for TPZMaterial
#include "pzmatrix.h"         // for TPZFMatrix, TPZMatrix, DecomposeType::E...
#include "pzskylstrmatrix.h"  // for TPZSkylineStructMatrix
#include "pzsolve.h"          // for TPZMatrixSolver<>::MSolver, TPZSolver
#include "pzstepsolver.h"     // for TPZStepSolver
#include "pzvec.h"            // for TPZVec
#include <functional>

using namespace std;

TPZAnalysisError::TPZAnalysisError(TPZCompMesh *mesh,std::ostream &out) : TPZAnalysis(mesh,true,out),fElIndexes(0),fElErrors(0),
fSingular(),fTotalError(0.),fAdmissibleError(0.0),fEtaAdmissible(0.05),fNIterations(4) {}

void TPZAnalysisError::SetAdaptivityParameters(REAL EtaAdmissible, int64_t NIterations) {
	fEtaAdmissible = EtaAdmissible;
	fNIterations = NIterations;
}
/** @brief Output file with number of iteration made. */
std::ofstream arq("Param.dat");
void TPZAnalysisError::hp_Adaptive_Mesh_Design(std::ostream &out,REAL &CurrentEtaAdmissible) {
	int64_t iter = 0;//iteracao atual
	cout << "\n\nIteration  1\n";
	out << "\n   Iteration  1\n";
	Run(out);//solucao malha inicial
	TPZManVector<REAL,3> errors(3);
	errors.Fill(0.0);
	TPZVec<REAL> flux(0);
    bool store_error = false;
	fNIterations--;
	while(iter++ < fNIterations) {
		arq << "\n iter = " << iter << endl;
		CurrentEtaAdmissible = pow(fEtaAdmissible,sqrt((REAL)(iter*1./(fNIterations*1.))));//error admissivel corrigido
		// Print data
		PlotLocal(iter,CurrentEtaAdmissible,out);
		//if more norms than 3 are available, the pzvec is resized in the material error method
		Mesh()->EvaluateError(fExact,store_error, errors);
		if (errors.NElements() < 3) {
			PZError << endl << "TPZAnalysisError::hp_Adaptive_Mesh_Design - At least 3 norms are expected." << endl;
			exit (-1);
		}
		out << " - Error  Energy Norm :  " << errors[0] << endl;
		out << " - FE Sol Energy Norm :  " << errors[1] << endl;
		out << " - Ex Sol Energy Norm :  " << errors[2] << endl;
		
		for(int ier = 3; ier < errors.NElements(); ier++)
			out << "Other norms : " << errors[ier] << endl;
		
		out << " - Eta Admissible     :  " << CurrentEtaAdmissible << endl;
		//      out << " - Eta Reached        :  " << true_error/exactnorm << endl;
		out << " - Eta Reached        :  " << errors[0]/errors[2] << endl;
		//Code isn't place to chat
		//#warning Philippe, nao entendo nada!!!!! //<!>
		//#warning De fato Thiago, voce tem razao
		
		out << " - Number of D.O.F.   :  " << fCompMesh->NEquations() << endl;
		out.flush();	  
		HPAdapt(CurrentEtaAdmissible,out);//processa os restantes elementos ; (nadmerror)
		Mesh()->AdjustBoundaryElements();
		OptimizeBandwidth();
		TPZSkylineStructMatrix skystr(fCompMesh);
		SetStructuralMatrix(skystr);
		TPZStepSolver<STATE> sol;
		sol.ShareMatrix(Solver());
		
		sol.SetDirect(ECholesky);//ECholesky
		SetSolver(sol);
		cout << "\n\nIteration " << (iter+1) << endl;
		out << "\n   Iteration " << (iter+1) << endl;
		Run(out);
	}
	//Code isn't place to chat
	//#warning Philippe, nao parece igual acima ?? //<!>
	//#warning Olhar aqui //<!>
	errors.Resize(3);
	errors.Fill(0.0);
	
	PlotLocal(iter,CurrentEtaAdmissible,out);
	Mesh()->EvaluateError(fExact,store_error, errors);
	
	if (errors.NElements() < 3) {
		PZError << endl << "TPZAnalysisError::hp_Adaptive_Mesh_Design - At least 3 norms are expected." << endl;
		exit (-1);
	}
	
	out << " - Error  Energy Norm :  " << errors[0] << endl;
	out << " - FE Sol Energy Norm :  " << errors[1] << endl;
	out << " - Ex Sol Energy Norm :  " << errors[2] << endl;
	
	for(int ier = 3; ier < errors.NElements(); ier++)
		out << "Other norms : " << errors[ier] << endl;
	
	out << " - Eta Admissible     :  " << CurrentEtaAdmissible << endl;
	//   out << " - Eta Reached        :  " << true_error/exactnorm << endl; <!>
	out << " - Eta Reached        :  " << errors[0]/errors[2] << endl;
	out << " - Number of D.O.F.   :  " << fCompMesh->NEquations() << endl;
	out.flush();	 
}

void TPZAnalysisError::PlotLocal(int64_t iter, REAL CurrentEtaAdmissible, std::ostream &out) {
    bool store_error = false;
	EvaluateError(CurrentEtaAdmissible,store_error, out);
	TPZManVector<std::string> solution(1);
	solution[0] = "Solution";
	TPZVec<std::string> vecnames(0);
	{
		std::stringstream plotfile;
		plotfile << "plot";
		plotfile <<  '0' << iter;
		plotfile << '\0';
		plotfile << ".plt";
		DefineGraphMesh(2, solution, vecnames, plotfile.str());
	}
	PostProcess(0);
	solution.Resize(2);
	solution[0] = "POrder";
	solution[1] = "Error";
	{
		std::stringstream plotfile;
		plotfile << "plop";
		plotfile <<  '0' << iter;
		plotfile << '\0';
		plotfile << ".plt";
		DefineGraphMesh(2, solution, vecnames, plotfile.str());
	}
	PostProcess(3);
}

void TPZAnalysisError::GetSingularElements(TPZStack<TPZCompElSide> &elist) {
	//point deve ser interior a algum elemento da malha
}
//SingularElements(..)
void TPZAnalysisError::ZoomInSingularity(REAL csi, TPZCompElSide elside, REAL singularity_strength) {

	REAL hn = 1./pow(csi,1./singularity_strength);
	REAL Q=2.;
	REAL NcReal = log( 1.+(1./hn - 1.)*(Q - 1.) )/log(Q);
	int Nc = 0;
	while(REAL(Nc) < (NcReal+0.5)) Nc++;
	int minporder = 2;
	
	TPZStack<TPZCompElSide> ElToRefine;
	TPZStack<int> POrder;
	TPZStack<TPZGeoElSide> subelements;
	TPZStack<int64_t> csubindex;
	ElToRefine.Push(elside);
	POrder.Push(Nc);
	while(ElToRefine.NElements()) {
		/** Take the next element and its interpolation order from the stack*/
		TPZCompElSide curelside = ElToRefine.Pop();
		int curporder = POrder.Pop();
		if(!curelside.Exists()) continue;
		int64_t cindex = curelside.Element()->Index();
		if(cindex < 0) continue;
		
		/** Cast the element to an interpolated element if possible*/
		TPZCompEl *cel = curelside.Element();
		TPZInterpolatedElement *cintel = 0;
		cintel = dynamic_cast<TPZInterpolatedElement *> (cel);
		/** If the element is not interpolated, nothing to do */
		if(!cintel) continue;
		/** Set the interpolation order of the current element to curporder*/
		if(curporder == minporder) {
			cintel->PRefine(Nc);
			fSingular.Push(curelside);
		} else {
			cintel->PRefine(curporder);
			cintel->Divide(cindex,csubindex,1);
			/** Identify the subelements along the side and push them on the stack*/
		}
		TPZGeoElSide gelside = curelside.Reference();
		if(!gelside.Exists()) continue;
		gelside.GetSubElements2(subelements);
		int64_t ns = subelements.NElements();
		curporder--;
		int64_t is;
		for(is=0; is<ns; is++) {
			TPZGeoElSide sub = subelements[is];
			TPZCompElSide csub = sub.Reference();
			if(csub.Exists()) {
				ElToRefine.Push(csub);
				POrder.Push(curporder);
			}
		}
	}
	ExpandConnected(fSingular);
	
	/*
	 REAL H1_error,L2_error,estimate;
	 TPZBlock *flux=0;
	 int64_t nel = fElIndexes.NElements();
	 for(int64_t elloc=0;elloc<nel;elloc++) {
	 int64_t el = fElIndexes[elloc];
	 estimate = fElErrors[elloc];
	 REAL csi = estimate / fAdmissibleError;
	 REAL h = h_Parameter(intellist[el]);
	 REAL hn = h/pow(csi,1./.9);
	 REAL Nc = log( 1.+(h/hn - 1.)*(Q - 1.) )/log(Q);
	 if(hn > 1.3*h) hn = 2.0*h*hn / (h + hn);
	 REAL hsub = h;//100.0;//pode ser = h ; Cedric
	 TPZCompEl *locel = intellist[el];
	 //obter um subelemento que contem o ponto singular e tem tamanho <= hn
	 TPZAdmChunkVector<TPZCompEl *> sublist;
	 while(hsub > hn) {
	 TPZVec<int64_t> indexsubs;
	 int64_t index = locel->Index();
	 locel->Divide(index,indexsubs,1);
	 int64_t nsub = indexsubs.NElements();
	 TPZAdmChunkVector<TPZCompEl *> listsub(0);
	 for(int64_t k=0;k<nsub;k++) {
	 index = listsub.AllocateNewElement();
	 listsub[index] = Mesh()->ElementVec()[indexsubs[k]];
	 }
	 //existe um unico filho que contem o ponto singular
	 SingularElement(point,listsub,sublist);
	 hsub = h_Parameter(sublist[0]);
	 }
	 TPZInterpolatedElement *intel = (TPZInterpolatedElement *) locel;
	 intel->PRefine(Nc+1);
	 indexlist.Push(intel->Index());
	 //os elemento viz devem ter ordens menores a cel quanto mais longe de point
	 TPZInterpolatedElement *neighkeep,*neigh;
	 //feito s�para o caso 1d , extender para o caso geral
	 int dim = intel->Dimension();
	 if(dim != 1) {
	 cout << "TPZAnalysisError::Step3 not dimension implemented , dimension = " << intellist[el]->Dimension() << endl;
	 return ;//exit(1);
	 }
	 for(int side=0;side<2;side++) {
	 int ly = 1;
	 TPZGeoElSide neighside = intel->Reference()->Neighbour(side);
	 TPZGeoElSide neighsidekeep = neighside;
	 TPZCompElSide neighsidecomp(0,0);
	 TPZStack<TPZCompElSide> elvec(0);
	 TPZCompElSide thisside(intel,side);
	 if(!neighsidekeep.Exists()) thisside.HigherLevelElementList(elvec,1,1);
	 if(!neighsidekeep.Exists() && elvec.NElements() == 0) {
	 neighsidekeep = thisside.LowerLevelElementList(1).Reference();
	 } else if(elvec.NElements() != 0) {
	 neighsidekeep = elvec[0].Reference();
	 }
	 while(ly < (Nc+1) && neighsidekeep.Exists() && neighsidekeep.Element()->Reference()->Material()->Id() > -1) {
	 neigh = (TPZInterpolatedElement *) neighsidekeep.Element()->Reference();
	 if(neigh) {
	 neigh->PRefine(ly);
	 int otherside = (neighsidekeep.Side()+1)%2;
	 neighsidekeep.SetSide(otherside);
	 indexlist.Push(neighsidekeep.Reference().Element()->Index());
	 }
	 neighside = neighsidekeep.Neighbour();
	 while(!neighside.Exists()) {
	 neighsidecomp = neighsidekeep.Reference();
	 neighsidecomp.HigherLevelElementList(elvec,1,1);
	 if(elvec.NElements()) {
	 neighside = elvec[0].Reference();
	 break;
	 }
	 neighside = neighsidecomp.LowerLevelElementList(1).Reference();
	 if(!neighside.Exists()) break;
	 }
	 neighsidekeep = neighside;
	 ly++;
	 }
	 }
	 }//for
	 Mesh()->InitializeBlock();
	 */
}

//void DivideRecursive(TPZCompEl *locel,int index,TPZVec<int> indexsubs,int hn);
void TPZAnalysisError::HPAdapt(REAL CurrentEtaAdmissible, std::ostream &out) {
	
	arq << "CurrentEtaAdmissible "  << CurrentEtaAdmissible << endl;
	
	TPZAdmChunkVector<TPZCompEl *>&listel = Mesh()->ElementVec();
    bool store_error = false;
	EvaluateError(CurrentEtaAdmissible,store_error, out);
	TPZVec<TPZCompElSide> SingLocal(fSingular);
	fSingular.Resize(0);
	
	int64_t nel = fElIndexes.NElements();
	for(int64_t ielloc=0;ielloc<nel;ielloc++) {
		int64_t iel = fElIndexes[ielloc];
		// if the element has already been treated (e.g. singularity) skip the process
		if(iel == -1) continue;
		TPZInterpolatedElement *elem = (TPZInterpolatedElement *) listel[iel];
		if(!elem || elem->Material()->Id() < 0) {
			PZError << "TPZAnalysisError::HPAdapt boundary element with error?\n";
			PZError.flush();
			continue;
		}
		REAL csi = fElErrors[ielloc] / fAdmissibleError;
		// Verificar se o element atual e um elemento singular
		int64_t ising, nsing = SingLocal.NElements();
		for(ising=0; ising<nsing; ising++) {
			if(elem == SingLocal[ising].Element()) {
				ZoomInSingularity(csi,SingLocal[ising]);
				break;
			}
		}
		// Go to the end of the loop if the element was handled by the singularity
		if(ising < nsing) continue;
		//calculo da ordem pn do elemento atual
		int nsides = elem->Reference()->NSides();
		REAL pFo = double(elem->EffectiveSideOrder(nsides-1));//quadrilatero
		
		// Newton's Method -> compute pNew    
		REAL pFn = pFo, res = 10.0, phi, del, dph, tol = 0.001;
		int64_t MaxIter = 100; int64_t iter=0;
		while (iter < MaxIter && res > tol) {
			phi = pFo+log(csi)-pFo*log(pFn/pFo);
			dph = pFn/(pFo+pFn);
			del = dph*(phi-pFn);
			if (del+pFn <= 0.0) // caiu fora do intervalo!
				del = 0.5 - pFn;
			res = fabs(del);
			pFn = del+pFn;
			iter++;
		} // end of Newton's Method
		
		if (iter == MaxIter)
			PZError << "\n - Newton's Method Failed at Element = " << elem->Reference()->Id() << endl;
		if(pFn < 1.) pFn = 1.;
		REAL factor = pow(csi,-1.0/pFn)*(pow(pFn/pFo,pFo/pFn));
		
		if(factor <= 0.75) 
			factor = 0.5;   // refine h
		else factor = 1.0; // don't refine h
		
		// Newton's Method -> compute once again pNew    
		pFn = pFo; iter = 0; res = 10.0;
		while (iter < MaxIter && res > tol) {
			phi = -pFn*log(factor)-log(csi)+pFo*log(pFn/pFo);
			dph = pFn/(pFo-pFn*log(factor));
			del = -dph*phi;
			if (del+pFn <= 0.0) // caiu fora do intervalo!
				del = 0.5 - pFn;
			res = fabs(del);
			pFn = del+pFn;
			iter++;
		} // end of Newton's Method
		
		if (iter == MaxIter)
			PZError << "\n - Newton's Method Failed at Element = " << elem->Reference()->Id() << endl;
		if(pFn < 1.) pFn = 1.;
		int pNew = 0;//(int) floor(pFn + 0.5);  // get the integer
		while(REAL(pNew) < (pFn+0.5)) pNew++;
		TPZVec<REAL> x(3),cs(2,0.);
		elem->Reference()->X(cs,x);
		
		elem->PRefine(pNew);
		TPZCompEl *locel = elem;
		//Divide elements
		if(factor == 0.5 ) {
			TPZVec<int64_t> indexsubs;
			int64_t index = locel->Index();
			elem->Divide(index,indexsubs,1);
		} 
	}
}


REAL TPZAnalysisError::h_Parameter(TPZCompEl *cel) {
	
 	REAL h = 0.,cicjdist;
	TPZGeoEl *gel = cel->Reference();
	int nconn = gel->NCornerNodes();
	for(int conni=0;conni<nconn;conni++) {
		for(int connj=conni;connj<nconn;connj++) {
			cicjdist = 0.;
			for(int coordi=0;coordi<3;coordi++) {
				REAL coor1 = gel->NodePtr(conni)->Coord(coordi);
				REAL coor2 = gel->NodePtr(connj)->Coord(coordi);
				cicjdist += pow( coor1 - coor2, (REAL)2.0);
			}
			cicjdist = sqrt(cicjdist);
			if(h < cicjdist) h = cicjdist;
		}
	}
	return h;
}

/** @brief Function to zeroes data */
void NullFunction(const TPZVec<REAL> &point,TPZVec<STATE>&val,TPZFMatrix<STATE> &deriv);

void NullFunction(const TPZVec<REAL> &point,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv) {
	
    val[0] = 0.*point[0];
    deriv(0,0) = 0.;
    deriv(1,0) = 0.;
}

void TPZAnalysisError::MathematicaPlot() {
	
	TPZGeoMesh *gmesh = fCompMesh->Reference();
	TPZAdmChunkVector<TPZGeoNode> &listnodes = gmesh->NodeVec();
	int64_t nnodes = gmesh->NNodes();
	TPZAdmChunkVector<TPZGeoNode> nodes(listnodes);
	TPZVec<int64_t> nodeindex(0);
	int64_t keepindex;
	int64_t in;
    TPZGeoNode *nodei = 0;   /// jorge 2017 - I think it's unnecessary, because if node is created its exists always
	for(in=0;in<nnodes;in++) {
		nodei = &nodes[in];
		REAL xi;
		if(nodei) xi = nodei->Coord(0);
		else continue;
		keepindex = in;
		for(int64_t jn=in;jn<nnodes;jn++) {
			TPZGeoNode *nodej = &nodes[jn];
			REAL xj;
			if(nodej) xj = nodej->Coord(0); else continue;
			if(xj < xi) {
				keepindex = jn;
				xi = xj;
			}
		}
		if(keepindex!=in) {
			TPZGeoNode commut = nodes[in];
			nodes[in] = nodes[keepindex];
			nodes[keepindex] = commut;
		}
	}
	ofstream mesh("Malha.dat");
	ofstream graph("Graphic.nb");
	mesh << "\nDistribuicao de nos\n\n";
	int64_t i;
	for(i=0;i<nnodes;i++) {
        nodei = &nodes[i];
		if(nodei) mesh << nodei->Coord(0) << endl;
	}
	//2a parte
	TPZVec<int64_t> locnodid(nnodes,0);
	TPZVec<TPZGeoEl *> gelptr(nnodes);
	int64_t nel = gmesh->NElements();
	int64_t count = 0;
	for(in=0;in<nnodes;in++) {
		nodei = &nodes[in];
		if(!nodei) continue;
		for(int64_t iel=0;iel<nel;iel++) {
			TPZGeoEl *gel = gmesh->ElementVec()[iel];
			if(!gel || !gel->Reference() || gel->MaterialId() < 0) continue;
			//int numnodes = gel->NNodes();
			int ic=0;//for(int ic=0;ic<numnodes;ic++)
			if(in==nnodes-1) ic = 1;
			if(nodei->Id() == gel->NodePtr(ic)->Id()) {
				gelptr[count] = gel;
				locnodid[count++] = ic;
				break;
            }
		}
	}
	TPZVec<int64_t> connects(nnodes);
	int64_t iel;
	for(iel=0;iel<nnodes;iel++) {
		//if(!gelptr[iel] || !(gelptr[iel]->Reference())) continue;
		connects[iel] = gelptr[iel]->Reference()->ConnectIndex(locnodid[iel]);
	}
	TPZVec<STATE> sol(nnodes);
	for(i=0; i<nnodes; i++) {
        TPZConnect *df = &fCompMesh->ConnectVec()[connects[i]];
        if(!df) {
         	cout << "\nError in structure of dates\n";
            mesh << "\nError in structure of dates\n";
            return;//exit(-1);
        }
        int64_t seqnum = df->SequenceNumber();
        int64_t pos = fCompMesh->Block().Position(seqnum);
        sol[i] = fCompMesh->Solution()(pos,0);
	}
	mesh << "\nSolucao nodal\n\n";
	for(i=0;i<nnodes;i++) {
		mesh << sol[i] << endl;
	}
	// expanding solution
	int64_t numsols = 5*(nnodes-1)+1;   // 5 values by element: 3 interpolated (interior) + 2 over corners
	TPZVec<STATE> expand_sol(numsols);
	TPZVec<REAL> expand_nodes(numsols);
 	TPZVec<REAL> qsi(1);
	int64_t exp_iel = -1;
	for(iel=0;iel<nnodes-1;iel++) {
		TPZManVector<STATE> locsol(1);
		exp_iel++;
		expand_sol[exp_iel] = sol[iel];
		qsi[0] = -1.;
		expand_nodes[exp_iel] = nodes[iel].Coord(0);
		REAL h = h_Parameter(gelptr[iel]->Reference());
		for(int i=1;i<5;i++) {
			exp_iel++;
			expand_nodes[exp_iel] = i*h/5. + nodes[iel].Coord(0);
			qsi[0] += .4;
			gelptr[iel]->Reference()->Solution(qsi,0,locsol);
			expand_sol[exp_iel] = locsol[0];
		}
	}
	expand_nodes[numsols-1] = nodes[nnodes-1].Coord(0);
	expand_sol[numsols-1] = sol[nnodes-1];
	// Mathematica output format
	graph << "list = {" << endl;
	for(int64_t isol=0;isol<numsols;isol++) {
		if(isol > 0) graph << ",";
		STATE expandsol = expand_sol[isol];
		if(fabs(expandsol) < 1.e-10) expandsol = 0.;
		graph << "{" << expand_nodes[isol] << "," << expandsol << "}";
		if( !((isol+1)%5) ) graph << endl;
	}
	graph << "\n};\n";
	graph << "ListPlot[list, PlotJoined->True, PlotRange->All];" << endl;
}
void TPZAnalysisError::EvaluateError(REAL CurrentEtaAdmissible, bool store_error, std::ostream &out) {
	//Code isn�t place to chat
	//#warning Philippe, tambem nao entendo aqui //<!>
	
	TPZManVector<REAL,3> elerror(3);
	elerror.Fill(0.);
	TPZManVector<REAL,3> errorSum(3);
	errorSum.Fill(0.0);
	
	TPZBlock<REAL> *flux = 0;
	int64_t elcounter=0;
	int64_t numel = Mesh()->ElementVec().NElements();
	fElErrors.Resize(numel);
	fElIndexes.Resize(numel);
	Mesh()->ElementSolution().Redim(numel,1);
	// Sum of the errors over all computational elements
	int64_t el;
	for(el=0;el< numel;el++) {
		TPZCompEl *elptr = Mesh()->ElementVec()[el];
		if(elptr && !(elptr->Material()->Id() < 0)) {
			elptr->EvaluateError(fExact,elerror, flux);
			int nerrors = elerror.NElements();
			errorSum.Resize(nerrors, 0.);
			for(int ii = 0; ii < nerrors; ii++)
				errorSum[ii] += elerror[ii]*elerror[ii];
			
			fElErrors[elcounter] = elerror[0];
			(Mesh()->ElementSolution())(el,0) = elerror[0];
			fElIndexes[elcounter++] = el;
		} else {
			(Mesh()->ElementSolution())(el,0) = 0.;
		}
	}
	fElErrors.Resize(elcounter);
	fElIndexes.Resize(elcounter);
	fTotalError = sqrt(errorSum[0]);
//    void NullFunction(TPZVec<REAL> &point,TPZVec<STATE>&val,TPZFMatrix<STATE> &deriv);

    std::function<void(const TPZVec<REAL> &,TPZVec<STATE>&,TPZFMatrix<STATE> &)> func(NullFunction);
	Mesh()->EvaluateError(func,store_error, elerror);
	fAdmissibleError = CurrentEtaAdmissible*sqrt(elerror[0]*elerror[0] + fTotalError*fTotalError) / sqrt(1.*elcounter);
}

void TPZAnalysisError::ExpandConnected(TPZStack<TPZCompElSide> &singel){
	int64_t nelem = singel.NElements();
	TPZStack<TPZGeoElSide> gelstack;
	TPZStack<TPZCompElSide> celstack;
	int64_t iel;
	for(iel=0; iel<nelem; iel++) {
		TPZCompElSide celside = singel[iel];
		TPZGeoElSide gelside;
		gelside = celside.Reference();
		if(!gelside.Exists()) continue;
		gelstack.Resize(0);
		cout << "This part needs to be fixed\n";
		//		gelside.Element()->LowerDimensionSides(gelside.Side(),gelstack);
		while(gelstack.NElements()) {
			TPZGeoElSide gelsideloc;
			gelsideloc = gelstack.Pop();
			gelsideloc.EqualLevelCompElementList(celstack,1,0);
			while(celstack.NElements()) {
				TPZCompElSide celsideloc = celstack.Pop();
				if(! celsideloc.Exists()) continue;
				int64_t nelsing = singel.NElements();
				int64_t smel;
				for(smel=0; smel<nelsing; smel++) if(singel[smel].Element() == celsideloc.Element()) break;
				if(smel != nelsing) singel.Push(celsideloc);
			}
		}
	}
}
