// //$Id: pzexplfinvolanal.cpp,v 1.11 2009-11-24 17:13:09 fortiago Exp $

#include "pzexplfinvolanal.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzseqsolver.h"
#include "checkconv.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"
#include "pzeuler.h"
#include "adapt.h"

using namespace std;


TPZExplFinVolAnal::TPZExplFinVolAnal(TPZCompMesh *mesh, std::ostream &out):TPZAnalysis(mesh,out){
  this->fTimeStep = -1.;
  this->fNMaxIter = -1.;
  this->fSimulationTime = 0.;
  this->Set(-1.,0.,0.);
  this->SetInitialSolutionAsZero();
  this->SetSaveFrequency(0,0);
  this->CleanAuxiliarVariables();
}

TPZExplFinVolAnal::~TPZExplFinVolAnal(){
//nothing to be done here
}

void TPZExplFinVolAnal::SetInitialSolution(TPZFMatrix<REAL> & InitialSol){
  const int nrows = this->Mesh()->Solution().Rows();
  const int ncols = this->Mesh()->Solution().Cols();
  if ( (InitialSol.Rows() != nrows) || (InitialSol.Cols() != ncols) ){
    PZError << "ERROR! " << __PRETTY_FUNCTION__ << " at line " << __LINE__ << std::endl;
  }
  else{
    this->fSolution = InitialSol;
  }
  TPZAnalysis::LoadSolution();
}

void TPZExplFinVolAnal::SetInitialSolutionAsZero(){
  TPZFMatrix<REAL> & MeshSol = this->Mesh()->Solution();
  this->fSolution.Redim( MeshSol.Rows(), MeshSol.Cols() );
  this->fSolution.Zero();
}

void TPZExplFinVolAnal::MultiResolution(double Epsl, std::ostream &out){

   this->DX(0, "testeinicial_");
   cout << "Starting adapt .... \n";cout.flush();
   TPZCompMesh * cmesh = this->Mesh();
   GetAdaptedMesh(cmesh, Epsl);
   cout << "Adapt finished .... \n";cout.flush();
    this->SetCompMesh(cmesh);   
    this->fSolution = this->Mesh()->Solution();
    TPZAnalysis::LoadSolution();
    this->Mesh()->LoadReferences();
    this->Mesh()->Reference()->RestoreReference(this->Mesh());


  fSimulationTime = 0.;

  this->DX(0, "testeinicialadapt_");
  TPZFMatrix<REAL> LastSol, NextSol;
  LastSol = fSolution;

  for(int iter = 0; iter < this->fNMaxIter; iter++){

    this->InitializeAuxiliarVariables();
    this->TimeEvolution(LastSol,NextSol);

      //checking steady state
      double steadynorm = 0.;
      for(int i = 0; i < LastSol.Rows(); i++){
        double val = (NextSol(i,0)-LastSol(i,0));
        steadynorm += val*val;
      }
      steadynorm = sqrt(steadynorm);
   std::cout << "*********** Steady state error at iteration " << iter << ", " << iter*fTimeStep << " = " << steadynorm << "\n\n";
   cout << "Norm(LastSol) = " << Norm(LastSol) << " , Norm(NextSol) = " << Norm(NextSol) << endl;
   cout << "LastSol.Rows = " << LastSol.Rows() << " , NextSol.Rows = " << NextSol.Rows() << endl;
   cout << "cmesh->NEquations() = " << cmesh->NEquations() << endl;

    this->CleanAuxiliarVariables();
    LastSol = NextSol;
    fSolution = NextSol;
    TPZAnalysis::LoadSolution();

  if( iter % 1  == 0){
//     this->DX(iter, "antes_");
    cout << "Starting adapt .... \n";cout.flush();
    GetAdaptedMesh(this->Mesh(), Epsl);
    cout << "Adapt finished .... \n";cout.flush();
    this->SetCompMesh(this->Mesh());   
//    this->DX(iter, "depois_");
    this->fSolution = this->Mesh()->Solution();
    LastSol = fSolution;
    TPZAnalysis::LoadSolution();
    this->Mesh()->LoadReferences();
    this->Mesh()->Reference()->RestoreReference(this->Mesh());
  }

    if (this->fSaveFrequency){
      if ((!(iter % fSaveFrequency)) || (iter == (this->fNMaxIter-1))){
        this->DX(iter,"teste_");
      }
    }
	//  this->Mesh()->Print();
	//  this->Mesh()->Reference()->Print();

    std::cout << "iter = " << iter << "\n";
    std::cout.flush();

  }//time step iterations

}//void

void TPZExplFinVolAnal::DX(int iter, string filename){
  TPZVec<string> scal(3-2),vec(0);
  scal[0] = "density";
  //  scal[1] = "energy";
  //  scal[2] = "Mach";
  stringstream nome; nome << filename << iter <<".dx";
  this->DefineGraphMesh(3,scal,vec,nome.str());  
  this->fSimulationTime = (iter+1.)*fTimeStep;
  this->PostProcess(0);
  this->CloseGraphMesh();
}

void TPZExplFinVolAnal::Run(std::ostream &out){

  cout << "  PASSO DE TEMPO = " << this->fTimeStep << endl;

  this->InitializeAuxiliarVariables();

  fSimulationTime = 0.;

  this->PostProcess(this->fDXResolution);

  TPZFMatrix<REAL> LastSol, NextSol;
  LastSol = fSolution;

  REAL steadynorm;

  for(int iter = 0; iter < this->fNMaxIter; iter++){

      this->TimeEvolution(LastSol,NextSol);

      if(iter == (fNMaxIter -1)){
        ofstream solfile("solucaoOlivier.txt");
        solfile.precision(12);
        for(int is = 0; is < NextSol.Rows(); is++){
          solfile << NextSol(is,0) << "\t";
          if((is+1)%20 == 0) solfile << "\n";
        }
        solfile << "\n";
      }

      //checking steady state
      steadynorm = 0.;
      for(int i = 0; i < LastSol.Rows(); i++){
        double val = (NextSol(i,0)-LastSol(i,0));
        steadynorm += val*val;
      }
      steadynorm = sqrt(steadynorm);

      LastSol = NextSol;
      fSolution = NextSol;
      TPZAnalysis::LoadSolution();

   if (this->fSaveFrequency){
    if ((!(iter % fSaveFrequency)) || (iter == (this->fNMaxIter-1))){
      this->fSimulationTime = (iter+1.)*fTimeStep;
      this->PostProcess(this->fDXResolution);
    }
   }

   std::cout << "*********** Steady state error at iteration " << iter << ", " << iter*fTimeStep << " = " << steadynorm << "\n\n";
   if (!fForceAllSteps){
    if (steadynorm < this->fSteadyTol){
      std::cout << "Steady state solution achieved\n\n";
      break;
    }
   }
   std::cout.flush();

  }//time step iterations

  this->CleanAuxiliarVariables();

}//method

void TPZExplFinVolAnal::PostProcess(int resolution, int dimension){
    REAL T = this->fSimulationTime;
    this->fTime = T;
    TPZAnalysis::PostProcess(resolution, dimension);
}//method


void TPZExplFinVolAnal::Set(REAL timestep, int niter, REAL eps, bool ForceAllSteps){
  this->fTimeStep = timestep;
  this->fNMaxIter = niter;
  this->fSteadyTol = eps;
  this->fForceAllSteps = ForceAllSteps;
}

void TPZExplFinVolAnal::SetSaveFrequency(int SaveFrequency, int resolution){
  this->fSaveFrequency = SaveFrequency;
  this->fDXResolution = resolution;
}

REAL TPZExplFinVolAnal::TimeStep(){
  return this->fTimeStep;
}


void TPZExplFinVolAnal::FromConservativeToPrimitiveAndLoad(const TPZFMatrix<REAL> & Solution){
  this->fSolution.Redim(this->Mesh()->NEquations(),1);//fSolution is now zeroed
  const int nv = this->NStateVariables();
  const int dim = this->Dimension();
  const int locSize = nv+nv*dim;//vars + gradients
  const int nvols = Solution.Rows()/locSize;
  TPZVec<REAL> sol(nv);
  for(int ivol = 0; ivol < nvols; ivol++){
    for(int i = 0; i < nv; i++){
      sol[i] = Solution.GetVal(ivol*locSize+i,0);
    }//i
    TPZEulerEquation::FromConservativeToPrimitive(sol,TPZEulerEquation::Gamma());
    for(int i = 0; i < nv; i++){
      fSolution(ivol*locSize+i,0) = sol[i];
    }//i
  }
  this->LoadSolution();
}//void

void TPZExplFinVolAnal::GetSol(TPZCompElDisc * disc, TPZVec<REAL> &sol){
  if(!disc || disc->NConnects() == 0){
    sol.Resize(0);
    return;
  }

  TPZCompMesh *cmesh = disc->Mesh();
  int bl = disc->Connect(0).SequenceNumber();
  int blpos = cmesh->Block().Position(bl);
  int blocksize = cmesh->Block().Size(bl);

  int nstate = this->NStateVariables();
  int nvars = nstate + nstate*this->Dimension();//variables + gradients
  sol.Resize(nvars);
  for(int i = 0; i < nvars; i++){
    sol[i] = this->fSolution(blpos+blocksize-nvars+i,0);
  }

}//void

void TPZExplFinVolAnal::GetNeighbourSolution(TPZInterfaceElement *face, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol){
  TPZCompElDisc * left = dynamic_cast<TPZCompElDisc*>(face->LeftElement());
  this->GetSol(left,LeftSol);
  TPZCompElDisc * right = dynamic_cast<TPZCompElDisc*>(face->RightElement());
  this->GetSol(right,RightSol);
}

void TPZExplFinVolAnal::ComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList){
  TPZElementMatrix ef(this->Mesh(), TPZElementMatrix::EF);

  TPZManVector<REAL,10> solL, solR;
  std::list< TPZInterfaceElement* >::iterator w;
  int iPira;
  for(iPira = 0, w = FacePtrList.begin(); w != FacePtrList.end(); w++, iPira++){
    TPZInterfaceElement * face = *w;
    if(!face){
      PZError << "\nFatal error in " << __PRETTY_FUNCTION__ << "\n";
      DebugStop();
    }

    this->GetNeighbourSolution(face,solL,solR);
    this->CalcResidualFiniteVolumeMethod(face,ef,solL,solR);
    ef.ComputeDestinationIndices();
    this->fRhs.AddFel(ef.fMat,ef.fSourceIndex,ef.fDestinationIndex);

  }//for iel = TPZInterfaceElement
}//void

void TPZExplFinVolAnal::DivideByVolume(TPZFMatrix<REAL> &vec, double alpha){
  std::map< TPZInterpolationSpace*, std::pair< REAL, TPZVec<int> > >:: iterator wVol;
  for(wVol = this->fVolumeData.begin(); wVol != this->fVolumeData.end(); wVol++){
    TPZInterpolationSpace *sp = wVol->first;
    if(!sp){
      PZError << "\nFatal error in " << __PRETTY_FUNCTION__ << "\n";
      DebugStop();
    }

    const REAL volume = wVol->second.first;
    TPZVec<int> destIndices = wVol->second.second;
    const int n = destIndices.NElements();
    for(int i = 0; i < n; i++){
      int pos = destIndices[i];
      vec(pos,0) *= alpha/volume;
    }

  }//for iel = TPZInterpolationSpace
}//void

void TPZExplFinVolAnal::ComputeGradient(const TPZFMatrix<REAL> & SolutionConsVars){
  this->FromConservativeToPrimitiveAndLoad(SolutionConsVars);

  //compute gradient into fRhs
  //fSolution must have zeros in gradient positions
  this->fRhs.Zero();
  TPZEulerEquation::SetComputeGradient();
  this->ParallelComputeFlux( this->fFacePtrList );
  this->DivideByVolume(fRhs,1.);
  fSolution += fRhs;//fRhs has zeros in state variables position and fSolution has zeros in gradient positions
}

void TPZExplFinVolAnal::AssembleFluxes2ndOrder(const TPZFMatrix<REAL> & Solution){
  this->FromConservativeToPrimitiveAndLoad(Solution);

  //compute gradient into fRhs
  //fSolution must have zeros in gradient positions
  this->fRhs.Zero();
  TPZEulerEquation::SetComputeGradient();
  this->ComputeFlux( this->fFacePtrList );
  this->DivideByVolume(fRhs,1.);
  fSolution += fRhs;//fRhs has zeros in state variables position and fSolution has zeros in gradient positions

  fRhs.Zero();
  TPZEulerEquation::SetComputeFlux();
  this->ParallelComputeFlux( this->fFacePtrList );
  this->DivideByVolume(fRhs,fTimeStep);
}//void

void * TPZExplFinVolAnal::ExecuteParallelComputeFlux(void * ExtData){
  TMTFaceData * data = static_cast< TMTFaceData* > (ExtData);
  data->fAn->ComputeFlux(data->fFaces);
  return NULL;
}

const int nthreads = 2;
void TPZExplFinVolAnal::ParallelComputeFlux(std::list< TPZInterfaceElement* > &FacePtrList){

  pthread_t allthreads[nthreads];
  for(int ithread = 0; ithread < nthreads; ithread++){
    allthreads[ithread] = NULL;
    pthread_create(&allthreads[ithread],NULL,ExecuteParallelComputeFlux, fVecFaces[ithread]);
  }//threads

  for(int i=0;i<nthreads;i++){
    if(!allthreads[i]) continue;
    pthread_join(allthreads[i], NULL);
  }

//   for(int i = 0; i < nthreads; i++){
//     delete vecFaces[i];
//     vecFaces[i] = NULL;
//   }
//   delete []allthreads;

}//void

void TPZExplFinVolAnal::TimeEvolution(TPZFMatrix<REAL> &LastSol, TPZFMatrix<REAL> &NextSol){

  const int order = 3;
  if(order == 1){
    //Euler explicit:
    int sz = fCompMesh->NEquations();
    fRhs.Redim(sz,1);
    this->AssembleFluxes(LastSol);
    NextSol = LastSol;
    NextSol -= this->fRhs;
      //un+1 = un -rhs(un)
  }
  if(order == 2){
    //RK2
    int sz = fCompMesh->NEquations();
    this->fRhs.Redim(sz,1);
    //uEtoile = un - rhs(un) //uEtoile is stored in NextSol to avoid another vector
    this->fRhs.Zero();
    this->AssembleFluxes(LastSol);//D = -rhs
    NextSol = LastSol;
    NextSol -= this->fRhs;

    //un+1 = 0.5 * ( un + uEtoile -rhs(uEtoile) )
    this->fRhs.Zero();
    this->AssembleFluxes(NextSol);
    NextSol += LastSol;
    NextSol -= this->fRhs;
    NextSol *= 0.5;
  }
  if(order == 3){
    //RK3: un+1 = ( un + 2 u** - 2 rhs(u**) ) (1/3)
    int sz = fCompMesh->NEquations();
    this->fRhs.Redim(sz,1);
    this->AssembleFluxes(LastSol);//rhs(un)
    NextSol = LastSol;
    NextSol -= this->fRhs;//u*
    this->fRhs.Zero();
    this->AssembleFluxes(NextSol);//rhs(u*)
    NextSol.ZAXPY(3.,LastSol); //NextSol += 3*LastSol
    NextSol -= this->fRhs;
    NextSol *= (1./4.);//u**
    this->fRhs.Zero();
    this->AssembleFluxes(NextSol);//rhs(u**)
    NextSol *= 2.;
    NextSol += LastSol;
    NextSol.ZAXPY(-2.,this->fRhs);
    NextSol *= (1./3.);//un+1 = ( un + 2 u** - 2 rhs(u**) ) (1/3)
  }

}//void

void TPZExplFinVolAnal::InitializeAuxiliarVariables(){
  //cleaning data structure
  this->CleanAuxiliarVariables();

  const int nelem = this->Mesh()->NElements();

  //finding interface elements
  for(int iel = 0; iel < nelem; iel++){
    TPZCompEl *el = this->Mesh()->ElementVec()[iel];
    if(!el) continue;
    TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(el);
    if(face){
      this->fFacePtrList.push_back(face);
    }
  }//for iel = TPZInterfaceElement

  //finding volume elements
  TPZElementMatrix ef(this->Mesh(), TPZElementMatrix::EF);
  for(int iel=0; iel < nelem; iel++){
    TPZCompEl *el = this->Mesh()->ElementVec()[iel];
    if(!el) continue;
    TPZInterpolationSpace* sp = dynamic_cast<TPZInterpolationSpace*>(el);
    if(!sp) continue;
    if(el->Reference()->Dimension() != 3) continue;

    const REAL volume = el->Reference()->Volume();
    sp->InitializeElementMatrix(ef);
    ef.ComputeDestinationIndices();

    std::pair< REAL, TPZVec<int> > myPair;
    myPair.first = volume;
    myPair.second = ef.fDestinationIndex;

    fVolumeData[sp] = myPair;
  }//for iel = TPZInterpolationSpace

  //splitting the list
  fVecFaces.Resize(nthreads);
  for(int it = 0; it < nthreads; it++){
    fVecFaces[it] = new TMTFaceData;
    fVecFaces[it]->fAn = this;
  }
  std::list< TPZInterfaceElement* >::iterator w;
  int it;
  for(w = fFacePtrList.begin(), it = 0; w != fFacePtrList.end(); w++){
    fVecFaces[it]->fFaces.push_back( *w );
    it++;
    if(it == nthreads) it = 0;
  }

}//void


void TPZExplFinVolAnal::CleanAuxiliarVariables(){
  this->fFacePtrList.clear();
  this->fVolumeData.clear();
}//void

void TPZExplFinVolAnal::CalcResidualFiniteVolumeMethod(TPZInterfaceElement *face, TPZElementMatrix &ef, TPZVec<REAL> &LeftSol, TPZVec<REAL> &RightSol){
  TPZDiscontinuousGalerkin *mat = dynamic_cast<TPZDiscontinuousGalerkin *>(face->Material().operator ->());
  if(!mat || mat->Name() == "no_name"){
      PZError << "TPZInterfaceElement::CalcResidual interface material null, do nothing\n";
      ef.Reset();
      return;
   }

   TPZCompElDisc * left = dynamic_cast<TPZCompElDisc*>(face->LeftElement());
   TPZCompElDisc * right = dynamic_cast<TPZCompElDisc*>(face->RightElement());

   if (!left || !right){
     PZError << "\nError at " << __PRETTY_FUNCTION__ << " null neighbour\n";
     ef.Reset();
     return;
   }
   if(!left->Material() || !right->Material()){
      PZError << "\n Error at " << __PRETTY_FUNCTION__ << " null material\n";
      ef.Reset();
      return;
   }

  TPZMaterialData data,dataleft,dataright;
  dataleft.sol[0] = LeftSol;
  dataright.sol[0] = RightSol;
  //neighbour centers
  {
    TPZManVector<REAL,3> qsi(3);
    dataleft.XCenter.Resize(3);
    dataright.XCenter.Resize(3);
    TPZGeoEl * gel = face->LeftElement()->Reference();
    gel->CenterPoint(gel->NSides()-1,qsi);
      gel->X(qsi,dataleft.XCenter);
    gel = face->RightElement()->Reference();
    gel->CenterPoint(gel->NSides()-1,qsi);
    gel->X(qsi,dataright.XCenter);
  }

  //Init ef parameter
   TPZManVector<TPZConnect*> ConnectL, ConnectR;
   TPZManVector<int> ConnectIndexL, ConnectIndexR;

   face->GetConnects( face->LeftElementSide(),  ConnectL, ConnectIndexL );
   face->GetConnects( face->RightElementSide(), ConnectR, ConnectIndexR );

   const int dim = face->Dimension();
   int nshapel = left->NShapeF();
   int nshaper = right->NShapeF();

   const int nstatel = left->Material()->NStateVariables();
   const int nstater = right->Material()->NStateVariables();
   const int ncon = ConnectL.NElements() + ConnectR.NElements();
   const int neql = nshapel * nstatel;
   const int neqr = nshaper * nstater;
   const int neq = neql + neqr;
   ef.fMat.Redim(neq,1);
   ef.fBlock.SetNBlocks(ncon);
   ef.fConnect.Resize(ncon);

   int ic = 0;
   int n = ConnectL.NElements();
   for(int i = 0; i < n; i++) {
    const int nshape = left->NConnectShapeF(i);
    const int con_neq = nstatel * nshape;
    ef.fBlock.Set(ic,con_neq);
    (ef.fConnect)[ic] = ConnectIndexL[i];
    ic++;
   }
   n = ConnectR.NElements();
   for(int i = 0; i < n; i++) {
    const int nshape = right->NConnectShapeF(i);
    const int con_neq = nstater * nshape;
    ef.fBlock.Set(ic,con_neq);
    (ef.fConnect)[ic] = ConnectIndexR[i];
    ic++;
   }
   ef.fBlock.Resequence();
   //ef till here

   TPZGeoEl *ref = face->Reference();
   TPZIntPoints *intrule = ref->CreateSideIntegrationRule(ref->NSides()-1, 0);
   const int npoints = intrule->NPoints();

   TPZManVector<REAL,3> intpoint(dim);
   data.x.Resize(3);
   REAL weight;

   //phil must have a size. I copy the solution
   dataleft.phi.Redim(dataleft.sol.NElements(),1);
   dataright.phi.Redim(dataright.sol.NElements(),1);

   //LOOP OVER INTEGRATION POINTS
   for(int ip = 0; ip < npoints; ip++){
      intrule->Point(ip,intpoint,weight);
      ref->Jacobian( intpoint, data.jacobian, data.axes, data.detjac, data.jacinv);
      ref->X(intpoint,data.x);
      weight *= fabs(data.detjac);
      face->Normal(data.axes,data.normal);
      mat->ContributeInterface(data, dataleft, dataright, weight, ef.fMat);
   }//loop over integration points

   delete intrule;
}

void TPZExplFinVolAnal::ComputeGradientForDetails(const TPZFMatrix<REAL> & PrimitiveSolution, TPZFMatrix<REAL> & SolutionWithGrad){
	this->InitializeAuxiliarVariables();
	fSolution = PrimitiveSolution;
	this->LoadSolution();	
	//compute gradient into fRhs
	//fSolution must have zeros in gradient positions
  this->fRhs.Redim(this->Mesh()->NEquations(),1);
	this->fRhs.Zero();
 	TPZEulerEquation::SetComputeGradient();
	this->ComputeFlux( this->fFacePtrList );
	this->DivideByVolume(fRhs,1.);
	fSolution += fRhs;//fRhs has zeros in state variables position and fSolution has zeros in gradient positions
	this->CleanAuxiliarVariables();
	SolutionWithGrad = this->fSolution;
}//void

