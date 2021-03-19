
//
//  toolstransienttime.cpp
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>

#include "toolstransienttime.h"


#include "pzmat1dlin.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZJIntegral2D.h"
#include "pzreducedspace.h"
#include "pzbndcond.h"
#include "pzl2projection.h"
#include "tpzmathtools.cpp"
#include "TPZVTKGeoMesh.h"
#include "pzvtkmesh.h"
#include "TPZSkylineNSymStructMatrix.h"

//Plasticidade
#include "pzelastoplasticanalysis.h"
#include "TPZSandlerDimaggio.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "TPZTensor.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZMatElastoPlastic2D.h"
//Plasticidade

//Teste CohesiveBC
#include "pznlelasmat.h"
#include "pznonlinanalysis.h"
#include "TPZCohesiveBC.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "pzcompelwithmem.h"
#include "CohesiveTests.h"
//Teste CohesiveBC

#include "pzgeoquad.h"
#include "pzlog.h"

#ifdef PZ_LOG
static PZLogger logger("pz.reducedspace.data");
#endif

REAL uglyGlobTol = 1.e-10;

REAL mypow(REAL a, int n)
{
	if (n == 0) return 1.;
	return (a * mypow(a,n-1));
}

ToolsTransient::ToolsTransient() : flastElastSol(){
  fMustStop = true;
  fCouplingMaterial1 = NULL;
	fCouplingMaterialH1 = NULL;
	fCohesiveMaterial = NULL;
  fCohesiveMaterialFirst = NULL;
  fgmesh = NULL;
  fmeshvec.Resize(0);
  fmphysics = NULL;
	fSetRunH1 = false;
  fPostProcessNumber = 0;
  fwhichPropag = 0;
  DebugStop();//Nao deveria utilizar este construtor
}

ToolsTransient::ToolsTransient(int pOrder) : fPostprocess(), flastElastSol()
{
  fpOrder = pOrder;
  fPostProcessNumber = 0;
  fwhichPropag = 0;
	fSetRunH1 = false;
  fMustStop = false;
  
  int dim = 2;
  //fCouplingMaterial1 = new TPZPlasticFrac2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> >(globMultiFisicMatId1,dim, globFractInputData.E1(), globFractInputData.Poisson1(), globFractInputData.Visc());
  
#ifdef PlasticMC
  fCouplingMaterial1 = new TPZPlasticFrac2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>(globMultiFisicMatId1,dim, globFractInputData.E1(), globFractInputData.Poisson1(), globFractInputData.Visc());
#endif
  
#ifdef PlasticSDi
  fCouplingMaterial1 = new TPZPlasticFrac2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>(globMultiFisicMatId1,dim, globFractInputData.E1(), globFractInputData.Poisson1(), globFractInputData.Visc());
#endif
  
  fCouplingMaterialH1 = new TPZH1PlasticFrac2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>(globMultiFisicMatId1,dim, globFractInputData.E1(), globFractInputData.Poisson1(), globFractInputData.Visc());

  
  if (globFractInputData.IsMohrCoulomb()) {
    this->SetMohrCoulombParameters(globFractInputData.Poisson1(), globFractInputData.E1(), globFractInputData.Cohesion(), globFractInputData.PhiMC(), globFractInputData.PhiMC());
  }
  if (globFractInputData.IsSandler()){
    this->SetSandlerParameters();
  }
  
	fCohesiveMaterial = new TPZCohesiveBC(globCohesiveMatId);
  fCohesiveMaterialFirst = new TPZCohesiveBC(globCohesiveMatIdHalf);
	
  int planestrain = 0;
  fCouplingMaterial1->SetfPlaneProblem(planestrain);
	fCouplingMaterialH1->SetfPlaneProblem(planestrain);
  
  fgmesh = NULL;
  fmeshvec.Resize(2);
  fmphysics = NULL;
}

void ToolsTransient::SetMohrCoulombParameters(REAL poisson, REAL elast, REAL cohesion, REAL Phi, REAL Psi)
{
#ifdef PlasticMC
  TPZElasticResponse ER;
  ER.SetEngineeringData(elast,poisson);
  fPlasticStepPV.fYC.SetUp(Phi, Psi, cohesion, ER);
  fPlasticStepPV.fER.SetEngineeringData(elast,poisson);
#endif
}

void ToolsTransient::SetSandlerParameters()
{
#ifdef PlasticSDi
  
  REAL poisson = 0.2;
  REAL elast = 25.e3;
  REAL A = 155.54; //152.54
  REAL B = 0.0015489;
  REAL C = 146.29;
  REAL R = 0.91969;
  REAL D = 0.018768;
  REAL W = 0.006605;
  TPZElasticResponse ER;
  ER.SetUp(elast,poisson);
  STATE G=elast/(2.*(1.+poisson));
  STATE K=elast/(3.*(1.-2*poisson));
  STATE phi=0,psi=1.,N=0.;
  fPlasticStepPV.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
  fPlasticStepPV.fER.SetUp(elast,poisson);
#endif
}

ToolsTransient::~ToolsTransient(){
	delete fCouplingMaterial1;
	delete fCouplingMaterialH1;
	delete fCohesiveMaterial;
  delete fCohesiveMaterialFirst;
}

void ToolsTransient::SetRunH1(bool RunH1)
{
	fSetRunH1 = RunH1;
}


void ToolsTransient::Run()
{
  //this->SetRunH1();
  
  int anCount = 0;
  if (globFractInputData.IsMohrCoulomb() || globFractInputData.IsSandler()) {
    fCouplingMaterial1->SetRunPlasticity(); // AQUINATHAN PLASTICIDADE
		fCouplingMaterialH1->SetRunPlasticity();
    uglyGlobTol = 5.e-6;
  }
  this->Mesh2D();  //Principal Geometric Mesh (Lf initial)
  this->InitializeUncoupledMeshesAttributes();
  this->CMeshMultiphysics();
  std::cout << "elastic nel = " << fmeshvec[0]->NElements() << std::endl;
  std::cout << "elastic neq = " << fmeshvec[0]->NEquations() << std::endl;
  std::cout << "fluid nel = " << fmeshvec[1]->NElements() << std::endl;
  std::cout << "fluid neq = " << fmeshvec[1]->NEquations() << std::endl;
  std::cout << "multphysics nel = " << fmphysics->NElements() << std::endl;
  std::cout << "multphysics neq = " << fmphysics->NEquations() << std::endl;
  
	bool OptimizeBandwidth = false;
  TPZAnalysis * an = new TPZAnalysis(fmphysics,OptimizeBandwidth);
  if (globPlotVTK) {
    if (globFractInputData.IsElastic()){
      globFractOutputData.PlotElasticVTK(an, fPostProcessNumber++);
    }
    else{
      CreatePostProcessingMesh();
      this->PostProcess();
    }
  }

  PostprocessPressure();
  PostProcessAcumVolW();
  PostProcessVolLeakoff();
  bool initialElasticKickIsNeeded = true;
  
  
	// CUIDADO ESSE CHECKCONV FERRA O RESULTADO FINAL! E ESTA FEITO PARA O CASO DE PRESTRESS
	bool IWantCheckConv = false;
	if(IWantCheckConv){
		CheckConv(OptimizeBandwidth);
		DebugStop();
	}

  while(fMustStop == false)
  {
    bool propagate = this->SolveSistTransient(an, initialElasticKickIsNeeded);
    initialElasticKickIsNeeded = false;
    anCount = an->GetStep();
    
    if(propagate == true)///Setting new Lfrac and tranferring solution
    {
      
      REAL newLfrac = globFractInputData.Lf() + globFractInputData.Lmax_edge();
      globFractInputData.SetLf(newLfrac);
      std::cout << "Nel propaga = " << globFractInputData.GetnElPropag() << std::endl;
      TPZCompMeshReferred *fmeshref = dynamic_cast<TPZCompMeshReferred *>(fmeshvec[0]);
      TPZCompMeshReferred *lastElastReferredCMesh = NULL;
      if (fmeshref) {
        lastElastReferredCMesh = this->CMeshReduced(fmeshref->ReferredMesh());
      }
      else{
        DebugStop();
      }
      
      lastElastReferredCMesh->LoadSolution(flastElastSol);
      
      fwhichPropag++;
      if (fwhichPropag >= globMaxNumberOfPropagations){
        break;
      }
      globFractInputData.GetnElPropag() = 0;
      globFractInputData.SetNotPropagated(); // so it can run the next time steps
      int lastElIndex = -6378;
      if(!this->FindElementAfterFracture(lastElIndex)){
        DebugStop(); // This element must exist!!!
      }
      std::cout << "Index of Quad Element After Fracture = " << lastElIndex << std::endl;
      TPZCompEl *cel = fmphysics->ElementVec()[lastElIndex];
      TPZGeoEl *gel = cel->Reference();
      const int lastmatid = globFractInputData.GetPressureMatIds_StripeId_ElastId().rbegin()->first;
      const int newmatid = lastmatid+1;
      
      TPZGeoEl *newGeoEl = gel->CreateBCGeoEl(4,newmatid);
      int64_t index;
      
      //fmeshvec[0]->CreateCompEl(newGeoEl, index);
      
      TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
      TPZMaterial *mat1 = fmphysics->FindMaterial(globMultiFisicMatId1);
      TPZMaterial *NewBCPressure = mat1->CreateBC(fCouplingMaterial1, newmatid, typeNeumann, val1, val2);
      fmphysics->InsertMaterialObject(NewBCPressure);
      TPZCompEl *celnew = fmphysics->CreateCompEl(newGeoEl, index);
      if (globFractInputData.NStripes() == 1) {
        globFractInputData.InsertBCId_StripeId_ElastId(newmatid, 0, globMultiFisicMatId1);
      }
      else{
        PZError << "not prepared to handle more than 1 stripe!!" << std::endl;
        DebugStop();
      }
      
#ifdef PZDEBUG
      std::ofstream outg("GeoMeshAfterPropagation.vtk");
      TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, outg, true);
#endif
      
      fgmesh->ResetReference();
      fmphysics->LoadReferences();
      //delete fmeshvec[0];
      //delete fmeshvec[1];
      std::cout << "*********************** FRACTURE PROPAGATED ************************" << std::endl;
      std::cout << "Number of propagations: " << fwhichPropag << std::endl;
      this->InitializeUncoupledMeshesAttributes();
      // Creating multiphysic elements into mphysics computational mesh
      TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fmphysics);
      TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fmphysics);
      celnew->PrepareIntPtIndices(); // soh posso fazer isso depois do multifisico enxergar seus compels

      // Setando o material como metade da forca
      fmphysics->LoadReferences();
      TPZGeoElSide side(gel,1);
      TPZGeoElSide neigh = side.Neighbour();
      while (side != neigh) {
        if (neigh.Element()->MaterialId() == globCohesiveMatId) {
          break;
        }
        neigh = neigh.Neighbour();
      }
      if (side == neigh) {
        DebugStop(); // should find the cohesive element
      }
      neigh.Element()->SetMaterialId(globCohesiveMatIdHalf);
      TPZCompEl *celpoint = neigh.Element()->Reference();
      celpoint->ForcePrepareIntPtIndices();
//      TPZMultiphysicsElement *mcel = dynamic_cast<TPZMultiphysicsElement *>(celpoint);
//      if (!mcel) {
//        DebugStop(); // Reference set wrong!!!
//      }
      
      //bool usingLeakOff = globFractInputData.GetIfUsingLeakOff();
      //globFractInputData.SetUsingLeakOff(false);
      
      an = new TPZAnalysis(fmphysics,OptimizeBandwidth);
      TransferElasticSolution(an,lastElastReferredCMesh);
      
      //globFractInputData.SetUsingLeakOff(usingLeakOff);
      
      //TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
      
      //fgmesh->ResetReference();
      fmphysics->LoadReferences();
      /*
      if (globFractInputData.IsElastic()) {
        globFractOutputData.PlotElasticVTK(an,fPostProcessNumber++);
      }
      else{
        CreatePostProcessingMesh();
        this->PostProcess();
      }
       */
      
    }
  }
  
  globFractOutputData.SetQinj1WingAndLfracmax(globFractInputData.Qinj(), globFractInputData.Lf());
  std::ofstream outMath("OutputMathematica.txt");
  globFractOutputData.PrintMathematica(outMath);
  
  if (fwhichPropag >= globMaxNumberOfPropagations){
    DebugStop(); // You didnt ran all the time steps!!!
  }
}


//------------------------------------------------------------------------------------


void ToolsTransient::InitializeUncoupledMeshesAttributes()
{
  
  if (fSetRunH1) {
    fmeshvec[0] = this->CMeshElasticH1();
  }
  else{
    TPZCompMesh * elastReference = ElastCMeshReferenceProcessed();
    std::cout << "elast nel = " << elastReference->NElements() << std::endl;
    std::cout << "elast neq = " << elastReference->NEquations() << std::endl;
    fmeshvec[0] = this->CMeshReduced(elastReference);
  }

  fmeshvec[1] = this->CMeshPressure();
  //fgmesh->ResetReference();

	
	bool SeeSol = false;
	if (SeeSol) {
    REAL oldval0 = fmeshvec[0]->Solution()(0,0);
    REAL oldval1 = fmeshvec[0]->Solution()(1,0);
    REAL oldval2 = fmeshvec[0]->Solution()(2,0);
    fmeshvec[0]->Solution()(0,0) = 12.9274;
    fmeshvec[0]->Solution()(1,0) = 3.46963;
    fmeshvec[0]->Solution()(2,0) = 3.64654;
		//fmeshvec[0]->Solution()(2,0) = 10.;
    fmeshvec[0]->Solution().Print("Solu");
		TPZMaterial * mat = fmeshvec[0]->FindMaterial(globReservMatId1);
		TPZManVector<std::string> scalnames(2),vecnames(1),tennames(0);
		scalnames[0]= "SigmaX";
		scalnames[1]= "SigmaY";
		vecnames[0] = "Displacement";
		std::string postprocessname("ElasticAfterMultPhysics.vtk");
		int dim = 2;
        if(!mat){
            DebugStop();
        }
        std::set<int> matids;
        int matid = mat->Id();
        matids.insert(matid);
		TPZVTKGraphMesh vtkmesh(fmeshvec[0],dim,matids,scalnames,vecnames,tennames);
		vtkmesh.SetFileName(postprocessname);
		vtkmesh.SetResolution(0);
		int numcases = 1;
		
		// Iteracoes de tempo
		int istep = 0;
		vtkmesh.DrawMesh(numcases);
		vtkmesh.DrawSolution(istep, 1.);
    
    fmeshvec[0]->Solution()(0,0) = oldval0;
    fmeshvec[0]->Solution()(1,0) = oldval1;
    fmeshvec[0]->Solution()(2,0) = oldval2;
		//fmeshvec[0]->Solution()(2,0) = 0.;
	}
  
}

TPZCompMesh * ToolsTransient::ElastCMeshReferenceProcessed()
{
  
  if(0) // test using h1 and plane stress linear elasticity
	{
    REAL pressure = 1.;
		TPZCompMesh *cmesh_h1 = this->CMeshCohesiveH1(pressure);
		TPZNonLinearAnalysis *nlan = new TPZNonLinearAnalysis(cmesh_h1,std::cout);
		TPZSkylineStructMatrix skyl(cmesh_h1);
		nlan->SetStructuralMatrix(skyl);
		TPZStepSolver<REAL> step;
		step.SetDirect(ELDLt);
		nlan->SetSolver(step);
		
		this->SolveNLElasticity(cmesh_h1,*nlan);
		
		int dim = 2;
		TPZStack<std::string> scalnames,vecnames;
		vecnames.push_back("Strain");
		vecnames.push_back("Displacement");
		scalnames.push_back("SigmaX");
		scalnames.push_back("SigmaY");
		
		nlan->DefineGraphMesh(dim, scalnames, vecnames, "CMeshH1.vtk");
		
		nlan->PostProcess(0);
    this->ShowDisplacementSigmaYCohesive(cmesh_h1);
		this->ShowDisplacementSigmaYBottom(cmesh_h1);
	}
	
  TPZCompMesh * cmesh_elast = this->CMeshElastic();
  TPZAnalysis * an = new TPZAnalysis(cmesh_elast,false);
	TPZSkylineStructMatrix full(cmesh_elast); //caso simetrico
	if (globFractInputData.NthreadsForAssemble() > 0) {
		full.SetNumThreads(globFractInputData.NthreadsForAssemble());
	}
	an->SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	an->SetSolver(step);
	this->SolveInitialElasticity(*an, cmesh_elast); // Resolvendo todos os problemas AO MESMO TEMPO! UAU!!!!
  TPZFMatrix<STATE> solutions = cmesh_elast->Solution();
	
#ifdef PZDEBUG
	bool IWantToSeeElastSol = true;
	if (IWantToSeeElastSol) {
		TPZStack<std::string> scalnames,vecnames;
		vecnames.Push("Displacement");
		scalnames.Push("SigmaX");
		scalnames.Push("SigmaY");
		std::string filename = "ElasticUncoupledSolutions.vtk";
		int dim = 2;
		an->DefineGraphMesh(dim, scalnames, vecnames, filename);
    TPZMaterial *mat = cmesh_elast->MaterialVec()[globReservMatId1];
    if (!mat) {
      DebugStop();
    }
    int NStripes = globFractInputData.NStripes();
    int NHats = globNHat;
    NHats = 1;
    for (int i=0; i<NStripes+NHats; i++) {
      mat->SetPostProcessIndex(i);
      an->PostProcess(0);
    }
	}
#endif
   int NStripes = globFractInputData.NStripes();
  /*
   // colocando as solucoes com as stripes
   solutions.Resize(solutions.Rows(), NStripes);
   for(int stripe = 0; stripe < NStripes; stripe++)
   {
   SetSigmaNStripeNum(cmesh_elast,stripe);
   if (stripe == 0) {
   an->Assemble();
   }
   else {
   an->AssembleResidual();
   }
   an->Solve();
   if (IWantToSeeElastSol) {
   an->PostProcess(0);
   }
   int oldsz = solutions.Cols();
   for(int r = 0; r < cmesh_elast->Solution().Rows(); r++)
   {
   solutions(r,stripe) = cmesh_elast->Solution()(r,0);// - solutions(r,0);
   }
   }
   */
   // colocando as funcoes hat
  int nhat = globNHat;
  nhat-=1;
  int porder = fpOrder;
  int oldsz = solutions.Cols();
  solutions.Resize(solutions.Rows(), oldsz+nhat);
  int dirid = globDirichletRecElastMatId1Cohe-2*(fwhichPropag+1);
  bool IWantToSeeHat = true;
  for (int ihat = 0; ihat < nhat; ihat++) {
    TPZCompMesh *cmesh_hat = CMeshHat(dirid, porder);
    an->SetCompMesh(cmesh_hat, false);
    this->SolveInitialElasticity(*an, cmesh_hat);
    if (IWantToSeeHat) {
      TPZStack<std::string> scalnames,vecnames;
      vecnames.Push("Displacement");
      scalnames.Push("SigmaX");
      scalnames.Push("SigmaY");
      std::string filename = "ElasticUncoupledSolutions.vtk";
      int dim = 2;
      an->DefineGraphMesh(dim, scalnames, vecnames, filename);
      an->SetStep(NStripes+1+ihat);
      if (globPlotVTK){
        an->PostProcess(0);
      }
    }
    int rowshat = cmesh_hat->Solution().Rows();
    if (rowshat != solutions.Rows()) {
      DebugStop(); //they have to be equal
    }
    for(int r = 0; r < rowshat; r++)
    {
      solutions(r,oldsz+ihat) = cmesh_hat->Solution()(r,0);// - solutions(r,0);
    }
    
    dirid-=2;
  }
  
  cmesh_elast->LoadSolution(solutions);
  
  return cmesh_elast;
}

// Used when was using prestress
void ToolsTransient::ApplyEquationFilter(TPZAnalysis * an)
{
	std::set<int64_t> eqOut;
	
	int neq = this->fmphysics->NEquations();
	int64_t blockAlphaElast = this->fmphysics->ConnectVec()[0].SequenceNumber();// Eu sei que o prestress eh na primeira posicao
	int64_t posBlock = this->fmphysics->Block().Position(blockAlphaElast);
	eqOut.insert(posBlock);
	
	TPZVec<int64_t> actEquations(neq-eqOut.size());
	int p = 0;
	for (int64_t eq = 0; eq < neq; eq++) {
		if (eqOut.find(eq) == eqOut.end()) {
			actEquations[p] = eq;
			p++;
		}
	}
	TPZEquationFilter eqF(neq);
	eqF.SetActiveEquations(actEquations);
	an->StructMatrix()->EquationFilter() = eqF;
	
}

// equation filter in one hat function
void ToolsTransient::ApplyEquationFilterInAllHats(TPZAnalysis * an)
{
	std::set<int64_t> eqOut;
	
	int neq = this->fmphysics->NEquations();
	int64_t posblock = globFractInputData.NStripes();
  for (int i = 0 ; i < globNHat; i++) {
    eqOut.insert(posblock+i);
  }

	TPZVec<int64_t> actEquations(neq-eqOut.size());
	int p = 0;
	for (int64_t eq = 0; eq < neq; eq++) {
		if (eqOut.find(eq) == eqOut.end()) {
			actEquations[p] = eq;
			p++;
		}
	}
	TPZEquationFilter eqF(neq);
	eqF.SetActiveEquations(actEquations);
	an->StructMatrix()->EquationFilter() = eqF;
	
}


// Used for tests. MustOptimizeBandwitch must be false
void ToolsTransient::ApplyEquationFilterOnPressure(TPZAnalysis * an)
{
	std::set<int64_t> eqOut;
	
	int neq = this->fmphysics->NEquations();
	int neqel = this->fmeshvec[0]->NEquations();
	
	for (int i = neqel; i < neq; i++) {
		eqOut.insert(i);
	}
	
	TPZVec<int64_t> actEquations(neq-eqOut.size());
	int p = 0;
	for (int64_t eq = 0; eq < neq; eq++) {
		if (eqOut.find(eq) == eqOut.end()) {
			actEquations[p] = eq;
			p++;
		}
	}
	TPZEquationFilter eqF(neq);
	eqF.SetActiveEquations(actEquations);
	an->StructMatrix()->EquationFilter() = eqF;
}

void ToolsTransient::Mesh2D()
{
  fgmesh = new TPZGeoMesh;
	bool IsPG = false;
	
	// PG Values mesh
	REAL q = globFractInputData.GMeshq();
	int ndivV = globFractInputData.NdivV();
	int ndivH = globFractInputData.NdivH();
	REAL a1V = globFractInputData.Lx() * (q - 1.)/(mypow(q,ndivV) - 1.);
	REAL a1H = globFractInputData.Ly() * (q - 1.)/(mypow(q,ndivH) - 1.);
	REAL posV = 0., posH = 0., acumV = 0., acumH = 0.;
	
	// Normal mesh values
	if (!IsPG) {
		ndivV = int(globFractInputData.Lx()/globFractInputData.Lmax_edge() + 0.5);
		ndivH = int(globFractInputData.Ly()/globFractInputData.Lmax_edge() + 0.5);
	}
	REAL deltadivV = globFractInputData.Lx()/ndivV;
	REAL deltandivH = globFractInputData.Ly()/ndivH;
  
  
  int64_t ncols = ndivV + 1;
  int64_t nrows = ndivH + 1;
  int64_t nnodes = nrows*ncols;
  
  fgmesh->NodeVec().Resize(nnodes);
  
  // Creating nodes
  int64_t nid = 0;
  REAL cracktipDist = globFractInputData.Lf();
  int colCracktip = -1;
	
  for(int64_t r = 0; r < nrows; r++)
  {
		if (IsPG) {
			posH += acumH*a1H;
		}
		else {
			posH = r*deltandivH;
		}
    
    
    for(int64_t c = 0; c < ncols; c++)
    {
			if (IsPG) {
				posV += acumV*a1V;
			}
			else {
				posV = c*deltadivV;
			}
      
      //REAL x = c*deltadivV;
      //REAL y = r*deltandivH;
      
      REAL dist = fabs(globFractInputData.Lf()-posV);
      if(r == 0 && dist < cracktipDist)
      {
        cracktipDist = dist;
        colCracktip = c;
      }
      
      TPZVec<REAL> coord(3,0.);
      coord[0] = posV;
      coord[1] = posH;
      fgmesh->NodeVec()[r*ncols + c].SetCoord(coord);
      fgmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
      nid++;
			if (IsPG) {
				if (c == 0) {
					acumV = 1.;
				}
				else {
					acumV *= q;
				}
			}
    }
		if (IsPG) {
			if (r == 0) {
				acumH = 1.;
			}
			else {
				acumH *= q;
			}
			posV = 0.;
			acumV = 0.;
		}
  }
	
	
  if(colCracktip == 0)
  {
    colCracktip = 1;//fratura minima corresponde aa distancia entre coluna 0 e coluna 1
  }
  
  // Creating elements
  TPZGeoEl * gel = NULL;
  TPZVec<int64_t> topol(4);
  int64_t indx = 0;
  for(int64_t r = 0; r < nrows-1; r++)
  {
    for(int64_t c = 0; c < ncols-1; c++)
    {
      topol[0] = r*(ncols) + c;
      topol[1] = r*(ncols) + c + 1;
      topol[2] = r*(ncols) + c + 1 + ncols;
      topol[3] = r*(ncols) + c + ncols;
      
      gel = fgmesh->CreateGeoElement(EQuadrilateral, topol, globReservMatId1, indx);
      gel->SetId(indx);
      indx++;
    }
  }
  
  fgmesh->BuildConnectivity();
  
  // Creating the BCs
  REAL stripeWidth = globFractInputData.Lf() / globFractInputData.NStripes();
  int64_t nelem = fgmesh->NElements();
  int bcId = globPressureMatId;
  for(int64_t el = 0; el < nelem; el++)
  {
    TPZGeoEl * gel = fgmesh->ElementVec()[el];
    
    //south BC
    TPZGeoElSide sideS(gel,4);
    TPZGeoElSide neighS(sideS.Neighbour());
    if(sideS == neighS)
    {
      if(el < colCracktip)
      {
        TPZGeoEl * bcFrac = gel->CreateBCGeoEl(4,bcId);
        
        //Increasing materialId number with respect with what stripe contains bcFrac
        TPZVec<REAL> centerQsi(bcFrac->Dimension(),0.);
        TPZVec<REAL> centerX(3,0.);
        bcFrac->CenterPoint(bcFrac->NSides()-1, centerQsi);
        bcFrac->X(centerQsi, centerX);
        REAL xCoord = centerX[0];
        int stripeNumber = (int)(xCoord/stripeWidth);
        
        globFractInputData.InsertBCId_StripeId_ElastId(bcId, stripeNumber, gel->MaterialId());
        
        bcId++;
        ////////

        
        if (el+1 == colCracktip) {
          gel->CreateBCGeoEl(1, globCohesiveMatIdHalf);
        }
      }
      else
      {
        if(gel->MaterialId() == globReservMatId1)
        {
          gel->CreateBCGeoEl(4, globDirichletBottom);
          gel->CreateBCGeoEl(1, globCohesiveMatId);
        }
        else
        {
          DebugStop(); // only using one mat, should never enter here
          gel->CreateBCGeoEl(4, globDirichletElastMatId2);
        }
      }
    }
    
    //east BC
    TPZGeoElSide sideE(gel,5);
    TPZGeoElSide neighE(sideE.Neighbour());
    if(sideE == neighE)
    {
      if(gel->MaterialId() == globReservMatId1)
      {
        gel->CreateBCGeoEl(5, globDirichletElastMatId1);
      }
      else
      {
				DebugStop(); // only using one mat, should never enter here
        gel->CreateBCGeoEl(5, globDirichletElastMatId2);
      }
      
    }
    
    //north BC
    TPZGeoElSide sideN(gel,6);
    TPZGeoElSide neighN(sideN.Neighbour());
    if(sideN == neighN)
    {
      if(gel->MaterialId() == globReservMatId1)
      {
        gel->CreateBCGeoEl(6, globDirichletElastMatId1);
      }
      else
      {
				DebugStop();
        gel->CreateBCGeoEl(6, globDirichletElastMatId2);
      }
    }
    
    //west BC
    TPZGeoElSide sideW(gel,7);
    TPZGeoElSide neighW(sideW.Neighbour());
    if(sideW == neighW)
    {
      gel->CreateBCGeoEl(7, globMixedElastMatId);
    }
  }
  
  
  globFractInputData.SetLastFracMatId(bcId);
  topol.Resize(1);
  for(int64_t p = 0; p < ncols; p++)
  {
    topol[0] = p;
    if(p == 0)
    {
      gel = fgmesh->CreateGeoElement(EPoint, topol, globBCfluxIn, indx);
    }
    else if(p == colCracktip)
    {
      gel = fgmesh->CreateGeoElement(EPoint, topol, globCracktip, indx);
    }
    indx++;
  }
  
  
  int diridhat = globDirichletRecElastMatId1Cohe;
  int recidhat = diridhat - 1;
  int ihat = 0;
  int iAllowThisNumberOfPropagations = globMaxNumberOfPropagations;
  int nhat = globNHat+iAllowThisNumberOfPropagations;
  
  for (ihat = 0 ; ihat < nhat ; ihat++){
    int ieltohat = 0;
    int whathat = 0;
    std::pair<int,int> twogeoelsindex;
    for(int64_t el = 0; el < nelem; el++)
    {
      TPZGeoEl * gel = fgmesh->ElementVec()[el];
      
      if (gel->MaterialId() == globReservMatId1){
        //south BC
        //TPZGeoElSide sideS(gel,4);
        //TPZGeoElSide neighS(sideS.Neighbour());
        if(el < ndivV)
        {
          if(el < colCracktip-1)
          {
            //gel->CreateBCGeoEl(4, diridhat); AQUINATHAN ideia do phil de mudar o espaco das hat
            gel->CreateBCGeoEl(0,globPointZeroDisplacement);
          }
          else if (ieltohat == 0 && whathat == ihat)
          {
            gel->CreateBCGeoEl(1, recidhat);
            ieltohat++;
            
            int gelindex = gel->Index();
            if (gelindex != el) {
              DebugStop(); // Id is diferent for position in the vector?
            }
            twogeoelsindex.first = gelindex;
            gel->CreateBCGeoEl(0,globPointZeroDisplacement);
          }
          else if(ieltohat == 1){
            ieltohat++;
            
            int gelindex = gel->Index();
            if (gelindex != el) {
              DebugStop(); // Id is diferent for position in the vector?
            }
            twogeoelsindex.second = gelindex;
            gel->CreateBCGeoEl(1, diridhat);
          }
          else{
						if (ieltohat > 1) {
							gel->CreateBCGeoEl(1, diridhat);
            }
            gel->CreateBCGeoEl(0,globPointZeroDisplacement);
            whathat++;
          }
        }
        
        //east BC
        if((el+1)%(ndivV) == 0 && el != 0)
        {
          gel->CreateBCGeoEl(5, diridhat);
        }
        
        //north BC
        if(el > ((ndivH-1)*ndivV)-1)
        {
          gel->CreateBCGeoEl(6, diridhat);
        }
        
        //west BC
        if(el%ndivV == 0)
        {
          //gel->CreateBCGeoEl(7, globMixedElastMatId); ja existe!
        }
      }
    }
    
    globFractInputData.GetfMatID_Rec_GeoEl()[diridhat] = twogeoelsindex;
    diridhat-=2;
    recidhat-=2;
  }
	
  fgmesh->BuildConnectivity();
	
#ifdef PZDEBUG
  std::map<int,std::pair<int, int> >::iterator it = globFractInputData.GetfMatID_Rec_GeoEl().begin();
  for (; it != globFractInputData.GetfMatID_Rec_GeoEl().end(); it++) {
    std::cout << "diriid = " << it->first << std::endl;
    std::cout << "Geoel 1 = " << it->second.first << "\tGeoel 2 = " << it->second.second << std::endl;
    if (it->second.second - it->second.first != 1) {
      DebugStop(); // Why arent the elements side by side? is there any new refinement?
    }
		if (it->second.first == 0) {
			DebugStop(); // There should be previous elements of the fracture
		}
  }
#endif
  
  
#ifdef PZDEBUG
  std::ofstream outg("GeoMeshForAll.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(fgmesh, outg, true);
#endif
  
}

TPZCompMesh * ToolsTransient::CMeshHat(int &dirid, int &porder)
{
  /// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestrain = 0;
  
  TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                globFractInputData.E1(),
                                                                globFractInputData.Poisson1(),
                                                                globFractInputData.Fx(),
                                                                globFractInputData.Fy(),
                                                                planestrain);
  
  
  TPZMaterial * mat1(material1);
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
  cmesh->SetDefaultOrder(fpOrder);
	cmesh->SetDimModel(dim);
  cmesh->InsertMaterialObject(mat1);
  
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  // Dirichlet em volta
  val1.Redim(2,2);
  val2.Redim(2,1);
  TPZMaterial * BCond11 = material1->CreateBC(mat1, dirid, typeDirichlet, val1, val2);
  
  //Dirichlet na fratura
  TPZMaterial *BCDiriFrac = material1->CreateBC(mat1, globPointZeroDisplacement, typeDirichlet, val1, val2);
  cmesh->InsertMaterialObject(BCDiriFrac);
  
  // mista na esquerda
  val1(0,0) = big;
  TPZMaterial * BCond21 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
  
	//Recalque unitario para fazer a funcao chapeu
  val1.Redim(2,2);
  val2.Redim(2,1);
  val2(1,0) = 0.0001;
  int recid = dirid - 1;
  TPZMaterial * BCond31 = material1->CreateBC(mat1, recid, typeDirichlet, val1, val2);
  
  cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(BCond11);
	cmesh->InsertMaterialObject(BCond21);
	cmesh->InsertMaterialObject(BCond31);
  
  cmesh->SetDefaultOrder(porder);
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
  // Setando order 1 nos sides sul dos elementos ao lado do recalque
  std::map< int,std::pair<int,int> >::iterator it = globFractInputData.GetfMatID_Rec_GeoEl().find(dirid);
  if (it == globFractInputData.GetfMatID_Rec_GeoEl().end()) {
    DebugStop(); // essa BC nao eXISTE!
  }
  
  
  std::pair<int, int> mygeoels;
  mygeoels = it->second;
  int geoel1 = mygeoels.first;
  int geoel2 = mygeoels.second;
  
  TPZGeoEl *gel1 = cmesh->Reference()->ElementVec()[geoel1];
  TPZGeoEl *gel2 = cmesh->Reference()->ElementVec()[geoel2];
  TPZCompEl *cel1 = gel1->Reference();
  TPZCompEl *cel2 = gel2->Reference();
  if(!cel1 || !cel2){
    DebugStop(); // Porque esse compels nao existem?
  }
  TPZInterpolatedElement *intel1 = dynamic_cast<TPZInterpolatedElement*>(cel1);
  TPZInterpolatedElement *intel2 = dynamic_cast<TPZInterpolatedElement*>(cel2);
  if (!intel1 || !intel2) {
    DebugStop(); // nao existe o interpoaletedelement?
  }
  
  //int side = 4;
  //int neworder = 1;
  //intel1->SetSideOrder(side,neworder); //AQUINATHAN, devo mudar a ordem?
  //intel2->SetSideOrder(side,neworder);
  
	return cmesh;
}

TPZCompMesh * ToolsTransient::CMeshElastic()
{
  /// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestrain = 0;
  
  TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                globFractInputData.E1(),
                                                                globFractInputData.Poisson1(),
                                                                globFractInputData.Fx(),
                                                                globFractInputData.Fy(),
                                                                planestrain);
	
	
  
  //material1->SetPreStress(globFractInputData.PreStressXX(), globFractInputData.PreStressYY(), globFractInputData.PreStressXY(), 0.);
  
  
  unsigned int nstripes = globFractInputData.NStripes();
  int nhats = globNHat;
  nhats = 1; // isso porque nao consigo rodar as outras hat aqui
  const int nloadcases = nstripes + nhats;
  material1->SetNumLoadCases(nloadcases);
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
  cmesh->SetDefaultOrder(fpOrder);
	cmesh->SetDimModel(dim);
  
	// material coesivo
	int cohesiveid = globCohesiveMatId;
  TPZCohesiveBC * material3 = new TPZCohesiveBC(cohesiveid);
  const REAL SigmaT = 0., DeltaC = -1.;
	const REAL DeltaT = 0. * DeltaC;
  material3->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *CoheMat(material3);
	cmesh->InsertMaterialObject(CoheMat);
  
  TPZCohesiveBC *CoheMatFirst = new TPZCohesiveBC(globCohesiveMatIdHalf);
  CoheMatFirst->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  cmesh->InsertMaterialObject(CoheMatFirst);
  
  TPZMaterial * mat1(material1);
  cmesh->InsertMaterialObject(mat1);
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  // tototototototo
  val2(1,0) = 1.;
  
  std::map< int,std::pair<int,int> >::iterator it;
  TPZVec<TPZFMatrix<STATE> > val2vec(nloadcases);
  val1.Redim(2,2);
  val2.Redim(2,1);
  for (int i=0; i<nloadcases; i++) {
    val2vec[i] = val2;
  }
  
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    int stripeindex = it->second.first;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
      TPZBndCond* BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2); //AQUINATHAN typeNeumann
      val2vec[stripeindex](1,0) = 1.;
      BCond11->SetLoadCases(val2vec);
      val2vec[stripeindex](1,0) = 0.;
      cmesh->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
      DebugStop(); // Nunca deveria entra nesse caso pois soh uso um material
    }
  }
  
  // Dirichlet no topo e direita
  TPZBndCond * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
 	cmesh->InsertMaterialObject(BCond21);
  
  if (globNHat == 0){
    TPZMaterial * BCondDiri = material1->CreateBC(mat1, globDirichletBottom, typeDirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCondDiri);
  }
  else{
    TPZBndCond * BCond22 = material1->CreateBC(mat1, globDirichletRecElastMatId1Cohe-(fwhichPropag+nhats-1)*2, typeDirichlet, val1, val2);
    cmesh->InsertMaterialObject(BCond22);    
  }
  
  for (int i=0; i<nloadcases; i++) {
    val2vec[i] = val2;
  }
  val1(0,0) = big;
  val1(1,1) = big;
  //val1.Zero();
  for (int ihat = 0; ihat<nhats; ihat++)
  {
    TPZBndCond * BCond23 = material1->CreateBC(mat1, globDirichletRecElastMatId1Cohe-1-2*(ihat+fwhichPropag), typeMixed, val1, val2);
    BCond23->SetNumLoadCases(nloadcases);
    
    val2vec[nstripes+ihat](1,0) = big*0.0001;
    BCond23->SetLoadCases(val2vec);
    val2vec[nstripes+ihat](1,0) = 0.;
    /*
    for(int j=ihat; j<nhats; j++)
    {
      val2vec[nstripes+j](1,0) = big*0.0001;
    }
    BCond23->SetLoadCases(val2vec);
    for(int j=ihat; j<nhats; j++)
    {
      val2vec[nstripes+j](1,0) = 0.;
    }
     */
    
    cmesh->InsertMaterialObject(BCond23);
  }
  
  
  val1(0,0) = big;
  val1(1,1) = 0.;
  
  TPZBndCond * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
	cmesh->InsertMaterialObject(BCond31);
  
  cmesh->SetAllCreateFunctionsContinuous();
  
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
	return cmesh;
}

void ToolsTransient::SetSigmaNStripeNum(TPZCompMesh * cmeshref, int actStripe)
{
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int matId = it->first;
    int stripe = it->second.first;
    TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
    if(!mat)
    {
      DebugStop();
    }
    TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
    if(!bcmat)
    {
      DebugStop();
    }
    if(stripe < actStripe)
    {
      bcmat->Val2()(1,0) = 0.;
    }
    else if(stripe == actStripe)
    {
      bcmat->Val2()(1,0) = globFractInputData.SigN();
    }
    else if(stripe > actStripe)
    {
      return;
    }
  }
}

TPZCompMeshReferred * ToolsTransient::CMeshReduced(TPZCompMesh *cmeshref){
  /// criar materiais
	int dim = 2;
  
  TPZVec<REAL> force(dim,0.);
  int planestrain = 0;
  
  TPZCompMeshReferred * cmeshreferred = new TPZCompMeshReferred(fgmesh);
  
  TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                globFractInputData.E1(),
                                                                globFractInputData.Poisson1(),
                                                                globFractInputData.Fx(),
                                                                globFractInputData.Fy(),
                                                                planestrain);
  
	// material coesivo
	const REAL SigmaT = globFractInputData.SigmaT();
	const REAL DeltaC = globFractInputData.DeltaC();
	const REAL DeltaT = globFractInputData.DeltaT();
  fCohesiveMaterial->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *CoheMat(fCohesiveMaterial);
	cmeshreferred->InsertMaterialObject(CoheMat);
  
  fCohesiveMaterialFirst->SetCohesiveData(SigmaT/2., DeltaC, DeltaT);
  TPZMaterial *CoheMatFirst(fCohesiveMaterialFirst);
	cmeshreferred->InsertMaterialObject(CoheMatFirst);
  
  //material1->SetPreStress(globFractInputData.PreStressXX(), globFractInputData.PreStressYY(), globFractInputData.PreStressXY(), 0.); NOT USING
  
  cmeshreferred->SetDimModel(dim);
  TPZMaterial * mat1(material1);
  cmeshreferred->InsertMaterialObject(mat1);
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
      TPZMaterial * BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2);
      cmeshreferred->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
			DebugStop(); // Soh tenho um!
    }
  }
  
  val1.Redim(2,2);
  val2.Redim(2,1);
  TPZMaterial * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
  cmeshreferred->InsertMaterialObject(BCond21);
  
  val1(0,0) = big;
  TPZMaterial * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
  cmeshreferred->InsertMaterialObject(BCond31);
  
  int numsol = cmeshref->Solution().Cols();
  cmeshreferred->AllocateNewConnect(numsol, 1, 1);
  
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
  
	cmeshreferred->SetDefaultOrder(fpOrder);
  cmeshreferred->SetDimModel(dim);
	
  fgmesh->ResetReference();
	cmeshreferred->AutoBuild();
  cmeshref->AdjustBoundaryElements();
	cmeshref->CleanUpUnconnectedNodes();
  cmeshreferred->LoadReferred(cmeshref);
	
	std::ofstream out("CMeshReduced.txt");
	cmeshreferred->Print(out);
  
  return cmeshreferred;
}

TPZCompMesh * ToolsTransient::CMeshElasticH1(){
  /// criar materiais
	int dim = 2;
  
  TPZVec<REAL> force(dim,0.);
  int planestrain = 0;
  
  TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);

  cmesh->SetDimModel(dim);
  
  // Material do reservatorio
  TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                globFractInputData.E1(),
                                                                globFractInputData.Poisson1(),
                                                                globFractInputData.Fx(),
                                                                globFractInputData.Fy(),
                                                                planestrain);
  

  TPZMaterial * mat1(material1);
  cmesh->InsertMaterialObject(mat1);

  // material coesivo
	const REAL SigmaT = globFractInputData.SigmaT();
	const REAL DeltaC = globFractInputData.DeltaC();
	const REAL DeltaT = globFractInputData.DeltaT();
  fCohesiveMaterial->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *CoheMat(fCohesiveMaterial);
	cmesh->InsertMaterialObject(CoheMat);
  
  fCohesiveMaterialFirst->SetCohesiveData(SigmaT/2., DeltaC, DeltaT);
  TPZMaterial *CoheMatFirst(fCohesiveMaterialFirst);
	cmesh->InsertMaterialObject(CoheMatFirst);
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  // Pressao
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
      TPZMaterial * BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2);
      cmesh->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
			DebugStop(); // Soh tenho um!
    }
  }
  
  
  val1.Redim(2,2);
  val2.Redim(2,1);
  
  // Dirichlet no topo e esquerda
  TPZMaterial * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
  cmesh->InsertMaterialObject(BCond21);
  
  // Dirichlet na esquerda soh em x
  val1(0,0) = big;
  TPZMaterial * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
  cmesh->InsertMaterialObject(BCond31);
  
	cmesh->SetAllCreateFunctionsContinuousWithMem();
  
	cmesh->SetDefaultOrder(fpOrder);
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}

TPZCompMesh * ToolsTransient::CMeshElasticH1ForPostProc(){
  /// criar materiais
	int dim = 2;
  
  TPZVec<REAL> force(dim,0.);
  int planestrain = 0;
  
  TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
  
  cmesh->SetDimModel(dim);
  
  // Material do reservatorio
#ifdef PlasticMC
  TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse>, TPZElastoPlasticMem> *material1 = new TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem>(globReservMatId1,1);
  material1->SetPlasticityModel(fPlasticStepPV);
#endif
#ifdef PlasticSDi
  TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse>, TPZElastoPlasticMem> *material1 = new TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem>(globReservMatId1,1);
  material1->SetPlasticityModel(fPlasticStepPV);
#endif
  
  TPZMaterial * mat1(material1);
  cmesh->InsertMaterialObject(mat1);
  
  // material coesivo
	const REAL SigmaT = globFractInputData.SigmaT();
	const REAL DeltaC = globFractInputData.DeltaC();
	const REAL DeltaT = globFractInputData.DeltaT();
  fCohesiveMaterial->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *CoheMat(fCohesiveMaterial);
	cmesh->InsertMaterialObject(CoheMat);
  
  fCohesiveMaterialFirst->SetCohesiveData(SigmaT/2., DeltaC, DeltaT);
  TPZMaterial *CoheMatFirst(fCohesiveMaterialFirst);
	cmesh->InsertMaterialObject(CoheMatFirst);
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  
  // Pressao
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
      TPZMaterial * BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2);
      cmesh->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
			DebugStop(); // Soh tenho um!
    }
  }
  
  
  val1.Redim(2,2);
  val2.Redim(2,1);
  
  // Dirichlet no topo e esquerda
  TPZMaterial * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
  cmesh->InsertMaterialObject(BCond21);
  
  // Dirichlet na esquerda soh em x
  val1(0,0) = big;
  TPZMaterial * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
  cmesh->InsertMaterialObject(BCond31);
  
	cmesh->SetAllCreateFunctionsContinuousWithMem();
  
	cmesh->SetDefaultOrder(fpOrder);
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
#ifdef PZDEBUG
  std::ofstream out("malhaPostProc");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * ToolsTransient::CMeshPressure(){
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
  cmesh->SetDefaultOrder(fpOrder);
  int dim = 1;
	cmesh->SetDimModel(dim);
  
  /// criar materiais
  TPZFMatrix<REAL> xk(1,1,1.);
  TPZFMatrix<REAL> xc(1,1,0.);
  TPZFMatrix<REAL> xb(1,1,0.);
  TPZFMatrix<REAL> xf(1,1,-2.);
  
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    TPZMat1dLin * material = new TPZMat1dLin(bcId);
    material->SetMaterial(xk,xc,xb,xf);
    
    TPZMaterial * mat(material);
    
    cmesh->InsertMaterialObject(mat);
    
    if(it == globFractInputData.GetPressureMatIds_StripeId_ElastId().begin())
    {
      ///Inserir condicao de contorno
      TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
      val2(0,0) = globFractInputData.Qinj();
      TPZMaterial * BCond1 = material->CreateBC(mat, globBCfluxIn, typeNeumann, val1, val2);
      cmesh->InsertMaterialObject(BCond1);
    }
    
    cmesh->InsertMaterialObject(mat);
  }
  cmesh->SetAllCreateFunctionsContinuous();
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
  
	return cmesh;
}



void ToolsTransient::SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	an.SetSolver(step);
	an.Run();
}

void ToolsTransient::CMeshMultiphysics()
{
  //Creating computational mesh for multiphysic elements
	fgmesh->ResetReference();
	fmphysics = new TPZCompMesh(fgmesh);
  fmphysics->SetDefaultOrder(fpOrder);
  
  fCouplingMaterial1->SetPlasticityModel(fPlasticStepPV);
  if (fSetRunH1){
#ifdef PlasticMC
    fCouplingMaterialH1->SetPlasticityModel(fPlasticStepPV);
#endif
  }
	
	TPZMaterial *mat1 = NULL;
	if (fSetRunH1) {
		mat1 = fCouplingMaterialH1;
	}
	else {
	  mat1 = fCouplingMaterial1;
	}
	if (!mat1) {
		DebugStop();
	}
  
  fmphysics->InsertMaterialObject(mat1);
  
  ///Inserir condicao de contorno
  REAL big = TPZMaterial::gBigNumber;
  
  TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
  
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
			TPZMaterial * BCond11 = NULL;
			if (fSetRunH1) {
				BCond11 = mat1->CreateBC(fCouplingMaterialH1, bcId, typeNeumann, val1, val2);
			}
			else {
				BCond11 = mat1->CreateBC(fCouplingMaterial1, bcId, typeNeumann, val1, val2);
			}
			if (!BCond11) {
				DebugStop();
			}
      
      fmphysics->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
      //      TPZMaterial * BCond12 = mat2->CreateBC(fCouplingMaterial2, bcId, typeNeumann, val1, val2);
      //      fmphysics->InsertMaterialObject(BCond12);
			DebugStop(); // not using matid2
    }
  }
  
	// dirichlet na direita e topo
  val2.Redim(3,1);
  val1.Redim(3,2);
  TPZMaterial * BCond21 = mat1->CreateBC(mat1, globDirichletElastMatId1, typeDir_elast, val1, val2);
  fmphysics->InsertMaterialObject(BCond21);
  
 	// dirichlet em x na esquerda
  val1(0,0) = big;
  TPZMaterial * BCond31 = mat1->CreateBC(mat1, globMixedElastMatId, typeMix_elast, val1, val2);
  fmphysics->InsertMaterialObject(BCond31);
  
	// vazao de entrada
  val2.Redim(3,1);
  val1.Redim(3,2);
  val2(0,0) = globFractInputData.Qinj();
  TPZMaterial * BCond41 = mat1->CreateBC(mat1, globBCfluxIn, typeNeum_pressure, val1, val2);
	fmphysics->InsertMaterialObject(BCond41);
  
	// material coesivo
	const REAL SigmaT = globFractInputData.SigmaT();
	const REAL DeltaC = globFractInputData.DeltaC();
	const REAL DeltaT = globFractInputData.DeltaT();
  fCohesiveMaterial->SetCohesiveData(SigmaT, DeltaC, DeltaT);
  TPZMaterial *CoheMat(fCohesiveMaterial);
	fmphysics->InsertMaterialObject(CoheMat);
  
  fCohesiveMaterialFirst->SetCohesiveData(SigmaT/2., DeltaC, DeltaT);
  TPZMaterial *CoheMatFirst(fCohesiveMaterialFirst);
  fmphysics->InsertMaterialObject(CoheMatFirst);
  
  // Setando o espaco
  fmphysics->SetAllCreateFunctionsMultiphysicElemWithMem();
  
  fmphysics->AutoBuild();
	fmphysics->AdjustBoundaryElements();
	fmphysics->CleanUpUnconnectedNodes();
  
  // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fmphysics);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fmphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
  
  // Preparando os index dos pontos de integracao.
	int64_t nel = fmphysics->NElements();
	for (int64_t iel = 0; iel < nel; iel++) {
		TPZCompEl *cel = fmphysics->ElementVec()[iel];
		cel->PrepareIntPtIndices();
	}
	
#ifdef PZDEBUG
   std::ofstream out("cmeshMult.txt");
   fmphysics->Print(out);
#endif
}

//------------------------------------------------------------------------------------

TPZCompMesh* ToolsTransient::CMeshCohesiveH1(REAL pressure)
{
  //Creating computational mesh for multiphysic elements
	fgmesh->ResetReference();
	TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  cmesh->SetDefaultOrder(fpOrder);
  
  /// criar materiais
	int dim = 2;
	
  TPZVec<REAL> force(dim,0.);
  
  int planestress = 1; // AQUI ESTOU USANDO PLANESTRESS
  //int planestrain = 0; // AQUI ESTOU USANDO PLANESTRAIN
  
  TPZNLElasticityMaterial * material1 = new TPZNLElasticityMaterial(globReservMatId1,
                                                                    globFractInputData.E1(),
                                                                    globFractInputData.Poisson1(),
                                                                    globFractInputData.Fx(),
                                                                    globFractInputData.Fy(),
                                                                    planestress);
	
	
  cmesh->InsertMaterialObject(material1);
  
  ///Inserir condicao de contorno
  REAL big = material1->gBigNumber;
  
  TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
  
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    int elastId = it->second.second;
    if(elastId == globReservMatId1)
    {//estou no globReservMatId1
			val2(1,0) = pressure;
      TPZMaterial * BCond11 = material1->CreateBC(material1, bcId, 1, val1, val2);
      cmesh->InsertMaterialObject(BCond11);
    }
    else
    {//estou no globReservMatId2
			DebugStop(); // not using matid2
    }
  }
  
	// dirichlet na direita e topo
  val2.Redim(3,1);
  val1.Redim(3,2);
  TPZMaterial * BCond21 = material1->CreateBC(material1, globDirichletElastMatId1, 0, val1, val2);
  cmesh->InsertMaterialObject(BCond21);
	
 	// dirichlet em x na esquerda
  val1(0,0) = big;
  TPZMaterial * BCond31 = material1->CreateBC(material1, globMixedElastMatId, 2, val1, val2);
  cmesh->InsertMaterialObject(BCond31);
  
	
	// material coesivo
	const REAL SigmaT = globFractInputData.SigmaT();
	const REAL DeltaC = globFractInputData.DeltaC();
	const REAL DeltaT = globFractInputData.DeltaT();
	TPZCohesiveBC *CoheMat = new TPZCohesiveBC(globCohesiveMatId);
  CoheMat->SetCohesiveData(SigmaT, DeltaC, DeltaT);
	cmesh->InsertMaterialObject(CoheMat);
  
  // primeiro coesivo
  TPZCohesiveBC *CoheMatFirst = new TPZCohesiveBC(globCohesiveMatIdHalf);
  CoheMatFirst->SetCohesiveData(SigmaT/2., DeltaC, DeltaT);
	cmesh->InsertMaterialObject(CoheMatFirst);
  
  //	std::ofstream outg("GeoMeshH1.vtk");
  //  TPZVTKGeoMesh::PrintCMeshVTK(cmesh, outg, true);
	
	cmesh->SetAllCreateFunctionsContinuousWithMem();
  cmesh->AutoBuild();
	cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
	/*
   const int neworder = 4;
   const int side = 4;
   for (int i = 0 ; i < 6 ; i++){
   TPZCompEl *cel = cmesh->ElementVec()[i];
   TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
   intel->SetSideOrder(side,neworder);
   }
	 */
	
	return cmesh;
}

//------------------------------------------------------------------------------------

void ToolsTransient::SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an)
{
	int iter = 0;
	REAL error = 1.e10;
	int numeq = cmesh->NEquations();
	
	TPZFMatrix<STATE> prevsol(an.Solution());
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
	
  REAL tol = uglyGlobTol;
  int numiter = 20;
	while(error > tol && iter < numiter) {
		
		an.Solution().Redim(0,0);
		an.Assemble();
    
		an.Solve();
    an.Solution() += prevsol;
    
		prevsol -= an.Solution();
		REAL normDeltaSol = Norm(prevsol);
		prevsol = an.Solution();
		an.LoadSolution(prevsol);
		an.AssembleResidual();
		double NormResLambda = Norm(an.Rhs());
		double norm = normDeltaSol;
    std::cout << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
		
		if(norm < tol) {
      std::cout << "\nTolerancia atingida na iteracao : " << (iter+1) << std::endl;
      std::cout << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;
			
		} else
			if( (norm - error) > 1.e-9 ) {
        std::cout << "\nDivergent Method\n";
			}
		error = norm;
		iter++;
    std::cout.flush();
	}
}

//------------------------------------------------------------------------------------

void ToolsTransient::StiffMatrixLoadVec(TPZAnalysis *an, TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec, bool IsEqFilter)
{
	if (fSetRunH1) {
		fCouplingMaterialH1->SetCurrentState();
	}
	else {
		fCouplingMaterial1->SetCurrentState();
	}
	
	fCohesiveMaterial->SetCurrentState();
  fCohesiveMaterialFirst->SetCurrentState();
  
  //TPZSkylineNSymStructMatrix matsk(fmphysics);
  TPZFStructMatrix matsk(fmphysics);
	if (globFractInputData.NthreadsForAssemble() > 0) {
		matsk.SetNumThreads(globFractInputData.NthreadsForAssemble());
	}
  
	an->SetStructuralMatrix(matsk);
	
	if (IsEqFilter) {
    //this->ApplyEquationFilterOnPressure(an);
    this->ApplyEquationFilterInAllHats(an);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec,fmphysics);
  }
	
	
	TPZStepSolver<REAL> step;
  
	step.SetDirect(ELU);
	an->SetSolver(step);
  
  an->Assemble();
	
  matK1 = an->Solver().Matrix();
  
	fvec = an->Rhs();
  
  if (fSetRunH1){
    AdjustTangentMatrix(matK1);
    AdjustResidual(fvec);
  }
}

void ToolsTransient::TransferSolutions(TPZCompMesh * lastMPhysicsCMesh, TPZCompMesh * lastElastReferredCMesh)
{
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
  //TransferElasticSolution(lastElastReferredCMesh);
  TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
  
  TransferLeakoff(lastMPhysicsCMesh);
}

void ToolsTransient::TransferElasticSolution(TPZAnalysis *an, TPZCompMesh * cmeshFrom)
{
#ifdef PZDEBUG
  if(!cmeshFrom)
  {
    DebugStop();
  }
#endif
  
  std::cout << "\n******************** TRANSFERINDO ********************\n\n";
  
  TPZMaterial *mat = fmphysics->FindMaterial(globBCfluxIn);
  if (!mat) {
    DebugStop();
  }
  /*
  TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
  TPZFMatrix<> oldval2 = bnd->Val2();
  bnd->Val2().Zero();
   */
  fCouplingMaterial1->SetLastElastCMesh(cmeshFrom);
  cmeshFrom->Solution().Print("cmeshfrom");
  this->RunOneStep(an);
  fCouplingMaterial1->SetLastElastCMesh(NULL);
  //bnd->Val2() = oldval2;
  
  /*
  REAL AreaFrom = IntegrateSolution(cmeshFrom, 0);
  fCouplingMaterial1->SetLastElastMesh(cmeshFrom);
  
  TPZAutoPointer< TPZFunction<STATE> > func = new TElastSolFunction<STATE>(cmeshFrom);

  
  //Setting old solution as forcing function on new cmesh (fmeshvec[0])
  TPZMaterial * ReservMat = NULL;
  std::map<int,TPZMaterial*>::iterator it = fmeshvec[0]->MaterialVec().find(globReservMatId1);
  if(it != fmeshvec[0]->MaterialVec().end())
  {
    ReservMat = it->second;
    ReservMat->SetForcingFunction(func);
  }
  else
  {
		DebugStop(); //only using one mat
    it = fmeshvec[0]->MaterialVec().find(globReservMatId2);
    if(it != fmeshvec[0]->MaterialVec().end())
    {
      ReservMat = it->second;
      ReservMat->SetForcingFunction(func);
    }
    else
    {
      DebugStop();//cade o mardito material???
    }
  }
  
  //////Solving
  cmeshFrom->LoadReferences();
  TPZAnalysis anTo(fmeshvec[0]);
  TPZSkylineStructMatrix full(fmeshvec[0]);
  anTo.SetStructuralMatrix(full);
  TPZStepSolver<REAL> step;
  step.SetDirect(ELDLt);
  anTo.SetSolver(step);
  anTo.Run();
  anTo.LoadSolution();
  fmeshvec[0]->Solution().Print("sol");
  anTo.Solution().Print("ansol");
  
  //Restoring original state
  ReservMat->SetForcingFunction(NULL);
  
  ///Integral correction
  REAL AreaTo = IntegrateSolution(fmeshvec[0], 0);
  
  if(fabs(AreaTo) > 1.E-18)
  {
    REAL alpha = AreaFrom/AreaTo;
    
    TPZFMatrix<REAL> solutionTo = fmeshvec[0]->Solution();
    for(int r = 0; r < solutionTo.Rows(); r++)
    {
      for(int c = 0; c < solutionTo.Cols(); c++)
      {
        solutionTo(r,c) *= alpha;
      }
    }
    
    fmeshvec[0]->LoadSolution(solutionTo);
  }
   */
}

REAL ToolsTransient::IntegrateSolution(TPZCompMesh * cmesh, int variable)
{
  REAL integral = 0.;
  for(int c = 0; c < cmesh->NElements(); c++)
  {
    TPZCompEl * cel = cmesh->ElementVec()[c];
    if(!cel || globFractInputData.IsBC(cel->Reference()->MaterialId()) == false)
    {
      continue;
    }
    TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(cel);
    if(!intel)
    {
      DebugStop();
    }
    TPZVec<REAL> value;
    intel->Integrate(variable, value);
    
    if(value.NElements() == 1)//integrando Leakoff (vlAcum)
    {
      integral += value[0];
    }
    else if(value.NElements() == 2)//Integrando uy (w)
    {
      integral += value[1];
    }
  }
  
  return integral;
}


void ToolsTransient::TransferLeakoff(TPZCompMesh * oldMphysicsCMesh)
{
  //    {
  //        std::ofstream outAntes("LeakoffANTES.txt");
  //        std::map<int,REAL>::iterator it;
  //        for(it = globFractInputData.GetLeakoffmap().begin(); it != globFractInputData.GetLeakoffmap().end(); it++)
  //        {
  //            outAntes << it->second << "\n";
  //        }
  //    }
  
  TPZCompMesh * cmeshTemp = new TPZCompMesh(fmeshvec[1]->Reference());
  cmeshTemp->SetAllCreateFunctionsDiscontinuous();
  cmeshTemp->AdjustBoundaryElements();
  
  //////L2Projection material
  int dim = 1;
  int pOrder = 0;
  int nsol = 1;
  TPZVec<REAL> solini(nsol,0.);
  
  TPZAutoPointer< TPZFunction<STATE> > func = new TLeakoffFunction<STATE>(oldMphysicsCMesh);
  std::map< int,std::pair<int,int> >::iterator it;
  for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
      it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
      it++)
  {
    int bcId = it->first;
    TPZL2Projection * materialL2 = new TPZL2Projection(bcId, dim, nsol, solini, pOrder);
    materialL2->SetForcingFunction(func);
    
    //////Inserindo na malha 1D
    cmeshTemp->InsertMaterialObject(materialL2);
  }
  
	cmeshTemp->AutoBuild();
  
  ///////Solving
  TPZAnalysis anTemp(cmeshTemp);
  TPZSkylineStructMatrix full(cmeshTemp);
  anTemp.SetStructuralMatrix(full);
  TPZStepSolver<REAL> step;
  step.SetDirect(ELDLt);
  anTemp.SetSolver(step);
  anTemp.Run();
  
  anTemp.LoadSolution();
  
  std::map<int,REAL> newLeakoff;
  for(int cel = 0; cel < cmeshTemp->NElements(); cel++)
  {
    TPZCompEl * compEl = cmeshTemp->ElementVec()[cel];
    TPZGeoEl * geoEl = compEl->Reference();
    
    int gelId = geoEl->Id();
    TPZVec<REAL> center(geoEl->Dimension());
    geoEl->CenterPoint(geoEl->NSides()-1, center);
    
    TPZInterpolationSpace * intpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
    TPZMaterialData data;
    intpEl->InitMaterialData(data);
    intpEl->ComputeShape(center, data);
    intpEl->ComputeSolution(center, data);
    TPZL2Projection * L2proj = dynamic_cast<TPZL2Projection *>(compEl->Material());
    TPZVec<REAL> Solout(1);
    int var = 1;//TPZL2Projection::ESolution
    L2proj->Solution(data, var, Solout);
    
    REAL vl = Solout[0];
    if(vl < 0.)
    {
      vl = 0.;
    }
    if(geoEl->Neighbour(2).Element()->MaterialId() == globReservMatId1 ||
       geoEl->Neighbour(2).Element()->MaterialId() == globReservMatId2)
    {
      newLeakoff[gelId] = vl;
    }
    else
    {
      DebugStop();
    }
  }
  globFractInputData.GetLeakoffmap() = newLeakoff;
  
  //    {
  //        std::ofstream outDepois("LeakoffDEPOIS.txt");
  //        std::map<int,REAL>::iterator it;
  //        for(it = globFractInputData.GetLeakoffmap().begin(); it != globFractInputData.GetLeakoffmap().end(); it++)
  //        {
  //            outDepois << it->second << "\n";
  //        }
  //    }
}

void ToolsTransient::MassMatrix(TPZFMatrix<REAL> & Un, bool IsFirstTime)
{
	if (fSetRunH1) {
		fCouplingMaterialH1->SetLastState();
	}
	else {
		fCouplingMaterial1->SetLastState();
	}
	
	fCohesiveMaterial->SetLastState();
  fCohesiveMaterialFirst->SetLastState();
  TPZFMatrix<> Oldsol = fmphysics->Solution();
  if (IsFirstTime) {
    TPZFMatrix<> ZeroSol(Oldsol.Rows(),Oldsol.Cols(),0.);
    fmphysics->LoadSolution(ZeroSol);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
  }
	TPZSpStructMatrix matsp(fmphysics);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
  matsp.CreateAssemble(Un,guiInterface);
  fmphysics->LoadSolution(Oldsol);
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
}

bool ToolsTransient::SolveSistTransient(TPZAnalysis *an, bool initialElasticKickIsNeeded)
{
  
  
  if(initialElasticKickIsNeeded)
  {
    //TPZFMatrix<REAL> chutenewton(fmeshvec[0]->Solution().Rows(), fmeshvec[0]->Solution().Cols(), 0.);
    //fmeshvec[0]->LoadSolution(chutenewton);
    //TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
  }
  
  REAL tol = uglyGlobTol;
  int maxit = 50;
  bool IsFirstTime = false;
  if(initialElasticKickIsNeeded)
  {
    REAL res = 0.;
    if (fSetRunH1) {
      res = this->IterativeProcess(IsFirstTime,an,maxit,tol,true); // Newton Method for one time step!
    }
    else{
      bool IWantEqFilter = true;
      bool usingLeakOff = globFractInputData.GetIfUsingLeakOff();
      globFractInputData.SetUsingLeakOff(false);
      res = this->IterativeProcess(IsFirstTime,an,maxit,tol,false,IWantEqFilter); // Newton Method for one time step!
      globFractInputData.SetUsingLeakOff(usingLeakOff);
    }
    IsFirstTime = true;    
    an->StructMatrix()->EquationFilter().Reset();
  }
  
  bool propagate = false;

	while( fMustStop == false && propagate == false )
	{
		
    if (fSetRunH1) {
      this->AddNoPenetration(globDirichletElastMatId1, 0);
      this->AddNoPenetration(globDirichletElastMatId1, 1);
      this->AddNoPenetration(globMixedElastMatId, 0);
      this->IdentifyEquationsToZero();
    }
    
    REAL res = 0.;
    
    flastElastSol = fmeshvec[0]->Solution();
    
    flastElastSol.Print("elassolll");
    
    if (fSetRunH1) {
      res = this->IterativeProcess(IsFirstTime,an,maxit,tol,true); // Newton Method for one time step!
    }
    else{
      res = this->IterativeProcess(IsFirstTime,an,maxit,tol); // Newton Method for one time step!
    }
    IsFirstTime = false;
    
    flastElastSol.Print("elassolll");
    
    if(res >= tol)
    {
      std::cout << "\nAtingido o numero maximo de iteracoes, nao convergindo portanto!!!\n";
      std::cout << "||Res|| = " << res << std::endl;
      DebugStop();
    }
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
		
		this->AcceptSolution(an); // this will update the memory
    
    if(globFractInputData.IsPropagated())
    {//propagou!!!
      globFractInputData.SetMinDeltaT();
      propagate = true;
    }
    else
    {//nao propagou!!!
      globFractInputData.UpdateLeakoff(fmeshvec[1]);
      globFractInputData.UpdateActTime();
      PostprocessPressure();
      
      globFractInputData.SetNextDeltaT();
      //fmat.Zero();
      //MassMatrix(fmat);
      if (globPlotVTK) {
        if (globFractInputData.IsElastic()) {
          globFractOutputData.PlotElasticVTK(an,fPostProcessNumber++);
        }
        else{
          CreatePostProcessingMesh();
          this->PostProcess();
        }
      }

      PostProcessAcumVolW();
      PostProcessVolLeakoff();
			//ShowDisplacementSigmaYCohesive(fmphysics);
    }
    REAL peteleco = 1.E-8;
    if( globFractInputData.actTime() > (globFractInputData.Ttot() - peteleco) )
    {
      fMustStop = true;
    }
    if (!globFractInputData.IsElastic())
    {
      this->SetUpdateToUseFullU(an);
    }
  }
  return propagate;
}

bool ToolsTransient::RunOneStep(TPZAnalysis *an)
{
  
  REAL tol = uglyGlobTol;
  int maxit = 50;
  bool propagate = false;
  
  if (fSetRunH1) {
    this->AddNoPenetration(globDirichletElastMatId1, 0);
    this->AddNoPenetration(globDirichletElastMatId1, 1);
    this->AddNoPenetration(globMixedElastMatId, 0);
    this->IdentifyEquationsToZero();
  }
  
  bool IsFirstTime = false;
  
  REAL res = 0.;
  if (fSetRunH1) {
    res = this->IterativeProcess(IsFirstTime,an,maxit,tol,true); // Newton Method for one time step!
  }
  else{
    res = this->IterativeProcess(IsFirstTime,an,maxit,tol); // Newton Method for one time step!
  }
  IsFirstTime = false;
  
  if(res >= tol)
  {
    std::cout << "\nAtingido o numero maximo de iteracoes, nao convergindo portanto!!!\n";
    std::cout << "||Res|| = " << res << std::endl;
    DebugStop();
  }
  
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
  
  this->AcceptSolution(an); // this will update the memory
  
  if(globFractInputData.IsPropagated())
  {//propagou!!!
    DebugStop(); // not handling this case, please lower deltaT
    globFractInputData.SetMinDeltaT();
    propagate = true;
  }
  else
  {//nao propagou!!!
    globFractInputData.UpdateLeakoff(fmeshvec[1]);
    globFractInputData.UpdateActTime();
    
    PostprocessPressure();
    globFractInputData.SetNextDeltaT();
    //fmat.Zero();
    //MassMatrix(fmat);
    if (globPlotVTK) {
      if (globFractInputData.IsElastic()) {
        globFractOutputData.PlotElasticVTK(an,fPostProcessNumber++);
      }
      else{
        CreatePostProcessingMesh();
        this->PostProcess();
      }
    }
    
    PostProcessAcumVolW();
    PostProcessVolLeakoff();
    //ShowDisplacementSigmaYCohesive(fmphysics);
  }
  REAL peteleco = 1.E-8;
  if( globFractInputData.actTime() > (globFractInputData.Ttot() - peteleco) )
  {
    fMustStop = true;
  }
  if (!globFractInputData.IsElastic())
  {
    this->SetUpdateToUseFullU(an);
  }

  return propagate;
}

REAL ToolsTransient::IterativeProcess(bool IsFirstTime, TPZAnalysis *an, int maxit, REAL tol, bool linesearch, bool IsEqFilter)
{
	int nrows = an->Solution().Rows();
	
	TPZFMatrix<REAL> res_total(nrows,1,0.);
  TPZFMatrix<REAL> SolIterK = fmphysics->Solution();
	SolIterK.Print("SolIterK");
	
	TPZAutoPointer< TPZMatrix<REAL> > matK;
	TPZFMatrix<REAL> fres(fmphysics->NEquations(),1);
  TPZFMatrix<REAL> fmat(fmphysics->NEquations(),1);
  fres.Zero();
  fmat.Zero();
  
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
  
	MassMatrix(fmat,IsFirstTime);
  if (fCouplingMaterial1->GetLastElastCMesh()) {
    fmat.Zero();
    for (int i = 0; i < flastMass.Rows(); i++) {
      fmat(i,0) = flastMass(i,0);
    }
  }
  else{
    flastMass = fmat;
  }

  fmat.Print("fmat");
  
  
	fres.Zero();
	StiffMatrixLoadVec(an, matK, fres, IsEqFilter);
	
	res_total = fres + fmat;
	REAL lastres = Norm(res_total);
	REAL res = lastres;
	int nit = 0;
	an->Solution().Print("sol");
	std::cout << "iter = " << nit << "\t||res|| = " << res << std::endl;
	
#ifdef PZ_LOG
	if (logger.isDebugEnabled()) {
		std::stringstream totout;
		totout << "************** Comecando as iteracoes do metodo de Newton *************** " << std::endl;
		totout << "iter = " << nit << "\t||res|| = " << res << std::endl;
		an->Solution().Print("Solution",totout);
		res_total.Print("restotal",totout);
		fres.Print("fres",totout);
		fmat.Print("fmat",totout);
		//matK->Print("matk",totout);
		LOGPZ_DEBUG(logger,totout.str())
	}
#endif
	
  bool linesearchconv = true;
	while(res > tol && nit < maxit) //itercao de Newton
	{
		an->Rhs() = res_total;
		an->Solve();
    
		if (linesearch){
			TPZFMatrix<REAL> nextSol;
			const int niter = 10;
      REAL RhsNormResult = 0.;
			this->LineSearch(an,SolIterK, an->Solution(), fmat, nextSol, res, RhsNormResult, niter,linesearchconv);
			SolIterK = nextSol;
		}
		else{
      SolIterK += an->Solution();
		}

    an->LoadSolution(SolIterK);
    
  
		TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
		
		SolIterK = an->Solution();
		
		fres.Zero();
		StiffMatrixLoadVec(an, matK, fres, IsEqFilter);
    if (fSetRunH1){
      AdjustResidual(fres);
    }

		
		res_total = fres + fmat;
		
#ifdef PZ_LOG
		if (logger.isDebugEnabled()) {
			std::stringstream totout;
      /*
			totout << "iter = " << nit+1 << "\t||res|| = " << res << std::endl;
			SolIterK.Print("Solution",totout);
			res_total.Print("restotal",totout);
			fres.Print("fres",totout);
			fmat.Print("fmat",totout);
			matK->Print("matk",totout);
      */
			LOGPZ_DEBUG(logger,totout.str())
		}
#endif
		
		res = Norm(res_total);
		SolIterK.Print("sol");
		std::cout << "iter = " << nit+1 << "\t||res|| = " << res << std::endl;
		if (res <= tol) {
			std::cout << "Convergencia atingida!!" << std::endl;
		}
    
		nit++;
	}
	return res;
}


void ToolsTransient::SetUpdateMem(bool update)
{
	if(!fmphysics) return;
	
	std::map<int, TPZMaterial *> & refMatVec = fmphysics->MaterialVec();
  
  std::map<int, TPZMaterial * >::iterator mit;
	
	TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem; // defined in file pzelastoplastic.h
  
  TPZCohesiveBC *CoheBC;
  
  for(mit=refMatVec.begin(); mit!= refMatVec.end(); mit++)
  {
    pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *>( mit->second );
		if(pMatWithMem != NULL)
    {
      pMatWithMem->SetUpdateMem(update);
    }
    CoheBC = dynamic_cast<TPZCohesiveBC *>(mit->second);
    if (CoheBC) {
      CoheBC->SetUpdateMem(update);
    }
  }
}

void ToolsTransient::AcceptSolution(TPZAnalysis *an)
{
	this->SetUpdateMem(true);
	
	an->Rhs().Zero();
	
  an->AssembleResidual();
  //AdjustResidual(fRhs);
	//REAL norm = Norm(fRhs);
	
	this->SetUpdateMem(false);
	
	//an->Solution().Zero();
	
	//an->LoadSolution();
  SetOpening();
}

void ToolsTransient::SetUpdateToUseFullU(TPZAnalysis *an)
{
  fCouplingMaterial1->SetUpdateToUseFullU(true);
	
	an->Rhs().Zero();
	
  an->AssembleResidual();
  //AdjustResidual(fRhs);
	//REAL norm = Norm(fRhs);
	
  fCouplingMaterial1->SetUpdateToUseFullU(false);

}

void ToolsTransient::PostprocessPressure()
{
  std::map<REAL,REAL> pos_pressure;
  
  for(int i = 0;  i < fmeshvec[1]->ElementVec().NElements(); i++)
  {
    TPZCompEl * cel = fmeshvec[1]->ElementVec()[i];
    TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*>(cel);
    if(!sp)
    {
      continue;
    }
    TPZMaterialData data;
    sp->InitMaterialData(data);
    
    for(int qsiPt = -1; qsiPt <= +1; qsiPt++)
    {
      TPZVec<REAL> qsi(1,0.), Xqsi(3,0.);
      qsi[0] = qsiPt;
      cel->Reference()->X(qsi,Xqsi);
      
      sp->ComputeShape(qsi, data);
      sp->ComputeSolution(qsi, data);
      
      REAL pos = Xqsi[0];
      REAL press = data.sol[0][0];
      pos_pressure[pos] = press;
    }
  }
  
  globFractOutputData.InsertTposP(globFractInputData.actTimeStep(), pos_pressure);
}

void ToolsTransient::SetOpening()
{
  TPZCompEl *cel = fmphysics->Element(0);
  TPZManVector<REAL,3> qsi(3,-1.), sol(3,0.);
  int varu = 7;
  cel->Solution(qsi, varu, sol);
  globFractInputData.SetOpening(sol[1]);  
}

void ToolsTransient::PostProcessAcumVolW()
{
  REAL wAcum = 2. * IntegrateSolution(fmeshvec[0], 0);//aqui jah sao consideradas as 2 faces da fratura
  int time = globFractInputData.actTimeStep();
  globFractOutputData.InsertTAcumVolW(time, wAcum);
}

void ToolsTransient::PostProcessVolLeakoff()
{
  //volume por trecho
  TPZAutoPointer< TPZFunction<STATE> > func = new TLeakoffFunction<STATE>(fmphysics);
  
  REAL vlAcum = 0.;
  
  for(int i = 0;  i < fmeshvec[1]->ElementVec().NElements(); i++)
  {
    TPZCompEl * cel = fmeshvec[1]->ElementVec()[i];
    int matid = cel->Reference()->MaterialId();
    if(globFractInputData.IsBC(matid) == false)
    {
      continue;
    }
    
#ifdef PZDEBUG
    if(!cel || cel->Reference()->Dimension() != 1)
    {
      DebugStop();
    }
#endif
    cel->Material()->SetForcingFunction(func);
    TPZGeoEl * fluidGel = cel->Reference();
    
#ifdef PZDEBUG
    if(!fluidGel)
    {
      DebugStop();
    }
#endif
    
    TPZVec<REAL> qsi(1,0.), xLeft(3,0.), xMiddle(3,0.), xRight(3,0.), leakoff(1,0);
    
    qsi[0] = -1.;
    fluidGel->X(qsi,xLeft);
    REAL posLeft = ((int)(xLeft[0]*1000.))/1000.;
    
    qsi[0] = +1.;
    fluidGel->X(qsi,xRight);
    REAL posRight = ((int)(xRight[0]*1000.))/1000.;
    
    REAL lengthElement = xRight[0] - xLeft[0];
    
    xMiddle[0] = (xLeft[0] + xRight[0])/2.;
    func->Execute(xMiddle, leakoff);
    
    //soh para nao coincidir pontos no mapa
    xLeft[0] = xMiddle[0] - 0.99*lengthElement/2.;
    xRight[0] = xMiddle[0] + 0.99*lengthElement/2.;
    
    REAL lengthLeakoff = leakoff[0];
    globFractOutputData.InsertTposVolLeakoff(globFractInputData.actTimeStep(), posRight, lengthLeakoff);
    globFractOutputData.InsertTposVolLeakoff(globFractInputData.actTimeStep(), posLeft, lengthLeakoff);
    
    cel->Material()->SetForcingFunction(NULL);
    
    //volume Acumulado
    vlAcum += 2. * (lengthLeakoff * lengthElement);//sao duas faces da fratura, lembra?
  }
  
  globFractOutputData.InsertTAcumVolLeakoff(globFractInputData.actTimeStep(), vlAcum);
}

void ToolsTransient::CheckConv(bool OptimizeBandwidth)
{
	
  std::cout.precision(15);
	int64_t neq = fmphysics->NEquations();
  int nsteps = 4;
	 
  TPZFMatrix<REAL> xIni(neq,1);
  for(int64_t i = 0; i < xIni.Rows(); i++)
  {
    REAL val = (double)(rand())*(1.e-8);
    xIni(i,0) = val;
  }
	
  TPZAnalysis *an = new TPZAnalysis(fmphysics,OptimizeBandwidth);
  an->LoadSolution(xIni);
  TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
	{
		fmphysics->Solution().Print(std::cout);
		fmeshvec[0]->Solution().Print(std::cout);
		fmeshvec[1]->Solution().Print(std::cout);
	}
  
  TPZFMatrix<REAL> actX = xIni;
  
  TPZAutoPointer< TPZMatrix<REAL> > fL_xIni;
  TPZFMatrix<REAL> f_xIni(neq,1);
  
  StiffMatrixLoadVec(an, fL_xIni, f_xIni);
	
  if(fL_xIni->Rows() != neq || fL_xIni->Cols() != neq || fL_xIni->IsDecomposed())
  {
    DebugStop();
  }
  
  TPZFMatrix<REAL> fAprox_x(neq,1);
  TPZFMatrix<REAL> fExato_x(neq,1);
  
  TPZFMatrix<REAL> errorVec(neq,1,0.);
	TPZFMatrix<REAL> errorNorm(nsteps,1,0.);
  
  
  TPZAutoPointer< TPZMatrix<REAL> > fLtemp;
  TPZFMatrix<REAL> dFx(neq,1);
  
  TPZVec<REAL> deltaX(neq,0.001), alphas(nsteps);
	for (int64_t i = 0; i < neq; i++) {
		REAL valdx = rand() % 10 + 1;
		valdx *= 0.000001;
		deltaX[i] = valdx;
	}
	
  double alpha;
  
  std::stringstream exatoSS, aproxSS;
  exatoSS << "exato={";
  aproxSS << "aprox={";
  for(int i = 0; i < nsteps; i++)
  {
    alpha = (i+1)/10.;
    alphas[i] = alpha;
    
    ///Fx aproximado
    dFx.Zero();
    for(int64_t r = 0; r < neq; r++)
    {
      for(int64_t c = 0; c < neq; c++)
      {
        dFx(r,0) +=  (-1.) * fL_xIni->GetVal(r,c) * (alpha * deltaX[c]); // (-1) porque dFx = -D[res,sol]*alpha*deltaX
      }
    }
    fAprox_x = f_xIni + dFx;
    
    int wantToSeeRow = 5;
    REAL aproxSol = fAprox_x(wantToSeeRow,0);
    std::cout << "Aprox : " << aproxSol << std::endl;
    {
      aproxSS << aproxSol;
      if(i < nsteps-1)
      {
        aproxSS << ",";
      }
    }
    
    ///Fx exato
    for(int64_t r = 0; r < neq; r++)
    {
      actX(r,0) = xIni(r,0) + (alpha * deltaX[r]);
    }
    an->LoadSolution(actX);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
    
    fExato_x.Zero();
    if(fLtemp) fLtemp->Zero();
    StiffMatrixLoadVec(an, fLtemp, fExato_x);
    
    REAL exatoSol = fExato_x(wantToSeeRow,0);
    std::cout << "Exato : " << exatoSol << std::endl;
    {
      exatoSS << exatoSol;
      if(i < nsteps-1)
      {
        exatoSS << ",";
      }
    }
    
    ///Erro
    errorVec.Zero();
    for(int64_t r = 0; r < neq; r++)
    {
      errorVec(r,0) = fExato_x(r,0) - fAprox_x(r,0);
    }
    
    ///Norma do erro
    double XDiffNorm = Norm(errorVec);
    errorNorm(i,0) = XDiffNorm;
  }
  aproxSS << "};";
  exatoSS << "};";
  std::cout << aproxSS.str() << std::endl;
  std::cout << exatoSS.str() << std::endl;
  std::cout << "Show[ListPlot[aprox, Joined -> True, PlotStyle -> Red],ListPlot[exato, Joined -> True]]\n";
  
  std::cout << "Convergence Order:\n";
  for(int j = 1; j < nsteps; j++)
  {
    std::cout << ( log(errorNorm(j,0)) - log(errorNorm(j-1,0)) )/( log(alphas[j]) - log(alphas[j-1]) ) << "\n";
  }
}


//---------------------------------------------------------------


template<class TVar>
TElastSolFunction<TVar>::TElastSolFunction()
{
  fIniElIndex = 0;
  fcmesh = NULL;
}

template<class TVar>
TElastSolFunction<TVar>::TElastSolFunction(TPZCompMesh * cmesh)
{
  fIniElIndex = 0;
  this->fcmesh = cmesh;
}

template<class TVar>
TElastSolFunction<TVar>::~TElastSolFunction()
{
  fIniElIndex = 0;
  fcmesh = NULL;
}

template<class TVar>
void TElastSolFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f)
{
  TPZVec<REAL> xcp(x);
  TPZVec<REAL> qsi2D(2,0.);
  TPZGeoEl * gel = fcmesh->Reference()->FindElement(xcp, qsi2D, fIniElIndex, 2);
  //TPZGeoEl *gel = fcmesh->ElementVec()[0]->Reference();
  
  if(!gel)
  {
    DebugStop();
  }
  
  TPZElasticityMaterial * elast = dynamic_cast<TPZElasticityMaterial *>(gel->Reference()->Material());
  if(!elast)
  {
    DebugStop();
  }
  
  f.Resize(2);
  
  TPZReducedSpace * intpEl = dynamic_cast<TPZReducedSpace *>(gel->Reference());
  if(!intpEl)
  {
    DebugStop();
  }
  TPZMaterialData data;
  intpEl->InitMaterialData(data);
  
  intpEl->ComputeShape(qsi2D, data);
  intpEl->ComputeSolution(qsi2D, data);
  
  int var = 0;
  elast->Solution(data, var, f);
  f.Resize(3);
  
}

template<class TVar>
void TElastSolFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
{
  DebugStop();
}

template<class TVar>
int TElastSolFunction<TVar>::NFunctions() const
{
  return 2;
}

template<class TVar>
int TElastSolFunction<TVar>::PolynomialOrder() const
{
  return fcmesh->GetDefaultOrder();
}


//---------------------------------------------------------------


template<class TVar>
TLeakoffFunction<TVar>::TLeakoffFunction()
{
  DebugStop();//use o outro construtor
}

template<class TVar>
TLeakoffFunction<TVar>::TLeakoffFunction(TPZCompMesh * cmesh)
{
  this->fIniElIndex = 0;
  this->fcmesh = cmesh;
}

template<class TVar>
TLeakoffFunction<TVar>::~TLeakoffFunction()
{
  fIniElIndex = 0;
  fcmesh = NULL;
}

template<class TVar>
void TLeakoffFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f)
{
  TPZVec<REAL> xcp(x);
  TPZVec<REAL> qsi1D(1,0.);
  TPZGeoEl * gel = fcmesh->Reference()->FindElement(xcp, qsi1D, fIniElIndex, 1);
  if(!gel)
  {
    DebugStop();
  }
  
  if(globFractInputData.IsBC(gel->MaterialId()) == false)
  {
    TPZGeoElSide gelSide(gel,gel->NSides()-1);
    TPZGeoElSide neighSide(gelSide.Neighbour());
    while(neighSide != gelSide)
    {
      if(globFractInputData.IsBC(neighSide.Element()->MaterialId()))
      {
        gel = neighSide.Element();
        break;
      }
      neighSide = neighSide.Neighbour();
    }
  }
  if(globFractInputData.IsBC(gel->MaterialId()) == false)
  {
    f[0] = 0.;
    return;
  }
  
  f.Resize(1);
  int elId = gel->Id();
  std::map<int,REAL>::iterator it = globFractInputData.GetLeakoffmap().find(elId);
  if(it == globFractInputData.GetLeakoffmap().end())
  {
    f[0] = 0.;
  }
  else
  {
    f[0] = it->second;
  }
}

template<class TVar>
void TLeakoffFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
{
  DebugStop();
}

template<class TVar>
int TLeakoffFunction<TVar>::NFunctions() const
{
  return 1;
}

template<class TVar>
int TLeakoffFunction<TVar>::PolynomialOrder() const 
{
  return 0;
}

//--------------------------------------------------- Tests and Prints -------------------------------------------------------------

void ToolsTransient::RunPlasticity()
{
  //std::map<int,REAL> leakoffMap;
  
  this->Mesh2D();
  
  TPZCompMesh * cmesh = CMeshElastoPlastic(fgmesh, globFractInputData.SigN());
  TPZElastoPlasticAnalysis an(cmesh,std::cout);
  
  this->SolveInitialElastoPlasticity(an, cmesh);
  
  std::string vtkFile = "pocoplastico.vtk";
  TPZPostProcAnalysis ppanalysis(cmesh);
  ppanalysis.SetStep(0);
  TPZFStructMatrix structmatrix(ppanalysis.Mesh());
  structmatrix.SetNumThreads(8);
  ppanalysis.SetStructuralMatrix(structmatrix);
  
  TPZVec<int> PostProcMatIds(1,1);
  TPZStack<std::string> PostProcVars, scalNames, vecNames;
  scalNames.Push("Alpha");
  scalNames.Push("PlasticSqJ2");
  scalNames.Push("PlasticSqJ2El");
  scalNames.Push("POrder");
  
  scalNames.Push("I1Stress");
  scalNames.Push("J2Stress");
  
  vecNames.Push("Displacement");
  vecNames.Push("YieldSurface");
  vecNames.Push("NormalStress");
  vecNames.Push("ShearStress");
  vecNames.Push("NormalStrain");
  vecNames.Push("ShearStrain");
  vecNames.Push("DisplacementMem");
  for (int i=0; i<scalNames.size(); i++)
  {
    PostProcVars.Push(scalNames[i]);
  }
  for (int i=0; i<vecNames.size(); i++)
  {
    PostProcVars.Push(vecNames[i]);
  }
  //
  ppanalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
  //
  ppanalysis.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
  //
  ppanalysis.TransferSolution();
  ppanalysis.PostProcess(0);// pOrder
}


void ToolsTransient::SolveInitialElastoPlasticity(TPZElastoPlasticAnalysis &analysis, TPZCompMesh *Cmesh)
{
	TPZSkylineStructMatrix full(Cmesh);
  full.SetNumThreads(0);
	analysis.SetStructuralMatrix(full);
  
	TPZStepSolver<REAL> step;
  step.SetDirect(ELDLt);
	analysis.SetSolver(step);
  
  int NumIter = 2;
  bool linesearch = true;
  bool checkconv = false;
  
  analysis.IterativeProcess(cout, 1.e-6, NumIter, linesearch, checkconv);
  
  analysis.AcceptSolution();
}

TPZCompMesh * ToolsTransient::CMeshElastoPlastic(TPZGeoMesh *gmesh, REAL SigmaN)
{
  /// criar materiais
  int dim = 2;
  
  TPZVec<REAL> force(dim,0.);
  
  //int planestress = 1;
  int planestrain = 1;
  
  TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> SD;
  REAL poisson = 0.203;
  REAL elast = 29269.;
  REAL A = 152.54;
  REAL B = 0.0015489;
  REAL C = 146.29;
  REAL R = 0.91969;
  REAL D = 0.018768;
  REAL W = 0.006605;
  SD.SetUp(poisson, elast, A, B, C, R, D, W);
  SD.SetResidualTolerance(1.e-10);
  SD.fIntegrTol = 10.;
  
  TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > * PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > (globReservMatId1,planestrain);
  
  TPZMaterial * mat(PlasticSD);
  PlasticSD->SetPlasticityModel(SD);
  
  TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (mat);
  
  ///criar malha computacional
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDefaultOrder(fpOrder);
  cmesh->SetDimModel(dim);
  cmesh->InsertMaterialObject(pMatWithMem);
  
  TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(cmesh);
  
  
  ///Inserir condicao de contorno
  REAL big = mat->gBigNumber;
  
  TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
  val2(1,0) = SigmaN;
  TPZMaterial * BCond1 = PlasticSD->CreateBC(pMatWithMem, globPressureMatId, typeNeumann, val1, val2);///NATHAN, falta incrementar stripenumber!!!
  
  val1.Redim(2,2);
  val2.Redim(2,1);
  TPZMaterial * BCond2 = PlasticSD->CreateBC(pMatWithMem, globDirichletElastMatId1, typeDirichlet, val1, val2);///NATHAN, agora temos 2 faixas de mat elastico!!!
  
  val1(0,0) = big;
  TPZMaterial * BCond3 = PlasticSD->CreateBC(pMatWithMem, globMixedElastMatId, typeMixed, val1, val2);
  
  //cmesh->SetAllCreateFunctionsContinuous();
  cmesh->InsertMaterialObject(BCond1);
  cmesh->InsertMaterialObject(BCond2);
  cmesh->InsertMaterialObject(BCond3);
  
  //Ajuste da estrutura de dados computacional
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}


void ToolsTransient::ShowDisplacementSigmaYCohesive(TPZCompMesh *cmesh)
{
	int nel = cmesh->NElements();
	
	TPZStack<std::pair<REAL,REAL> > cohepoints;
	int nelcomputed = 0;
	for (int iel = 0; iel < nel; iel++) {
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		int matid = gel->MaterialId();
		if (matid != globMultiFisicMatId1) {
			continue;
		}
		
		TPZGeoElSide gelside(gel,1);
		TPZGeoElSide neighbour = gelside.Neighbour();
		while (gelside != neighbour) {
			
			if (neighbour.Element()->Dimension() == 0) {
				if (neighbour.Element()->MaterialId() == globCohesiveMatId) {
					break;
				}
			}
			neighbour = neighbour.Neighbour();
		}
		if (gelside == neighbour) {
			continue;
		}
    if (neighbour.Element()->MaterialId() != globCohesiveMatId) {
      DebugStop();
    }
    
		//TPZMultiphysicsCompEl *mpcel = dynamic_cast<TPZMultiphysicsCompEl*> (cel);
		int var = -6378;
		
		int npoints = 3;
		TPZVec <REAL> qsi(3,-1.), Solout(3,0.);
		REAL delta =  2./ (REAL) (npoints-1);
		REAL displa,sigma;
    
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
    int var1 = -1, var2 = -1;
    if (intel) {
      var1 = 9;
      var2 = 6;
    }
    else{
      var1 = 9;
      var2 = 6;
    }
    
		for (int i = 0 ; i < npoints; i++) {
			qsi[0] = -1 + i * delta;
			var = var1;
			cel->Solution(qsi, var, Solout);
			displa = Solout[1];
			
			var = var2;
			cel->Solution(qsi, var, Solout);
			sigma = Solout[0];
			
			std::pair<REAL,REAL> mypair(displa,sigma);
			cohepoints.Push(mypair);
		}
		
		nelcomputed++;
		if (nelcomputed == 20) {
			break;
		}
	}
	
	std::cout << "cohepoints = " << "{"
	<< "{" << cohepoints[0].first << "," << cohepoints[0].second << "}";
	for (int i = 1; i < cohepoints.NElements(); i++) {
		std::cout << ",{" << cohepoints[i].first << "," << cohepoints[i].second << "}";
	}
	std::cout << "};" << std::endl;
	std::cout << "ListPlot[cohepoints,Joined->True]" << std::endl;
	
}

void ToolsTransient::ShowDisplacementSigmaYBottom(TPZCompMesh *cmesh)
{
	int nel = cmesh->NElements();
	
	TPZStack<std::pair<REAL,REAL> > cohepoints;
	int nelcomputed = 0;
	for (int iel = 0; iel < nel; iel++) {
		TPZCompEl *cel = cmesh->ElementVec()[iel];
		if (!cel) {
			continue;
		}
		TPZGeoEl *gel = cel->Reference();
		int matid = gel->MaterialId();
		if (matid != globMultiFisicMatId1) {
			continue;
		}
		
		TPZGeoElSide gelside(gel,1);
		TPZGeoElSide neighbour = gelside.Neighbour();
		while (gelside != neighbour) {
			
			if (neighbour.Element()->Dimension() == 0) {
				if (neighbour.Element()->MaterialId() == globCohesiveMatId) {
					break;
				}
			}
			neighbour = neighbour.Neighbour();
		}
		if (gelside == neighbour) {
			TPZGeoElSide gelside(gel,4);
			TPZGeoElSide neighbour = gelside.Neighbour();
			while (gelside != neighbour) {
				
				if (neighbour.Element()->Dimension() == 1) {
					if (globFractInputData.IsBC(neighbour.Element()->MaterialId())) {
						break;
					}
				}
				neighbour = neighbour.Neighbour();
			}
			if (gelside == neighbour) {
				continue;
			}
		}
		
		//TPZMultiphysicsCompEl *mpcel = dynamic_cast<TPZMultiphysicsCompEl*> (cel);
		int var = -6378;
		
		int npoints = 3;
		TPZVec <REAL> qsi(3,-1.), Solout(3,0.);
		REAL delta =  2./ (REAL) (npoints-1);
		REAL displa,sigma;
    
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(cel);
    int var1 = -1, var2 = -1;
    if (intel) {
      var1 = 9;
      var2 = 6;
    }
    else{
      var1 = 9;
      var2 = 6;
    }
		
    
		qsi[0] = -1;
		var = var1;
		cel->Solution(qsi, var, Solout);
		displa = Solout[1];
		
		var = var2;
		cel->Solution(qsi, var, Solout);
		sigma = Solout[0];
		
		std::pair<REAL,REAL> mypair(displa,sigma);
		cohepoints.Push(mypair);
		
		if (nelcomputed == 49) {
			qsi[0] = 1;
			var = var1;
			cel->Solution(qsi, var, Solout);
			displa = Solout[1];
			
			var = var2;
			cel->Solution(qsi, var, Solout);
			sigma = Solout[0];
			
			std::pair<REAL,REAL> mypair(displa,sigma);
			cohepoints.Push(mypair);
		}
    
		nelcomputed++;
		if (nelcomputed > 50) {
			DebugStop();
		}
		
	}
	
	std::cout << "cohepoints = " << "{"
	<< "{" << cohepoints[0].first << "," << cohepoints[0].second << "}";
	for (int i = 1; i < cohepoints.NElements(); i++) {
		std::cout << ",{" << cohepoints[i].first << "," << cohepoints[i].second << "}";
	}
	std::cout << "};" << std::endl;
	std::cout << "ListPlot[cohepoints,Joined->True]" << std::endl;
	
}

void ToolsTransient::IdentifyEquationsToZero()
{
  fEquationstoZero.clear();
  int64_t nel = fmphysics->NElements();
  for (int64_t iel=0; iel<nel; iel++) {
    TPZCompEl *cel = fmphysics->ElementVec()[iel];
    if (!cel) {
      continue;
    }
    TPZMaterial *mat = cel->Material();
    if (!mat) {
      continue;
    }
    int matid = mat->Id();
    if (fMaterialIds.find(matid) == fMaterialIds.end()) {
      continue;
    }
    int direction = fMaterialIds[matid];
    int64_t nc = cel->NConnects();
    for (int64_t ic=0; ic<nc; ic++) {
      TPZConnect &c = cel->Connect(ic);
      int64_t seqnum = c.SequenceNumber();
      int64_t pos = fmphysics->Block().Position(seqnum);
      int blsize = fmphysics->Block().Size(seqnum);
      for (int64_t i=pos+direction; i<pos+blsize; i+=2) {
        fEquationstoZero.insert(i);
      }
    }
  }
#ifdef PZ_LOG
  {
    std::stringstream sout;
    sout << "Equations to zero ";
    std::set<int64_t>::iterator it;
    for (it=fEquationstoZero.begin(); it!= fEquationstoZero.end(); it++) {
      sout << *it << " ";
    }
    LOGPZ_DEBUG(logger, sout.str())
  }
#endif
}

/// Apply zero on the lines and columns of the Dirichlet boundary conditions
void ToolsTransient::AdjustTangentMatrix(TPZMatrix<STATE> &matrix)
{
  std::set<int64_t>::iterator it;
  int64_t size = matrix.Rows();
  for (it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
    int eq = *it;
    for (int i=0; i<size; i++) {
      matrix.Put(eq, i, 0.);
      matrix.Put(i, eq, 0.);
    }
    matrix.Put(eq, eq, 1.);
  }
  
}

/// Apply zero to the equations of the Dirichlet boundary conditions
void ToolsTransient::AdjustResidual(TPZFMatrix<STATE> &rhs)
{
  std::set<int64_t>::iterator it;
  int64_t size = rhs.Rows();
  int64_t cols = rhs.Cols();
  for (it = fEquationstoZero.begin(); it != fEquationstoZero.end(); it++) {
    int64_t eq = *it;
    for (int64_t i=0; i<size; i++) {
      for (int64_t j=0; j<cols; j++)
      {
        rhs.Put(eq, j, 0.);
      }
    }
  }
}

//this->LineSearch(prevsol, fSolution, nextSol, RhsNormPrev, RhsNormResult, niter,linesearchconv);
REAL ToolsTransient::LineSearch(TPZAnalysis *an, const TPZFMatrix<REAL> &Wn, const TPZFMatrix<REAL> &DeltaW, TPZFMatrix<REAL> &MassVec, TPZFMatrix<REAL> &NextW, REAL RhsNormPrev, REAL &RhsNormResult, int niter, bool & converging){
  
  TPZFMatrix<REAL> Interval = DeltaW;
  
  int neq = Wn.Rows();
  TPZAutoPointer<TPZMatrix<> > matK = new TPZFMatrix<>(neq,neq);
  TPZFMatrix<> rhs(neq,0);
  REAL scalefactor = 1.;
  int iter = 0;
  do {
    Interval *= scalefactor;
    NextW = Wn;
    NextW += Interval;
    an->LoadSolution(NextW);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
    this->StiffMatrixLoadVec(an, matK, rhs);
    this->AdjustResidual(rhs);
    rhs += MassVec;
    RhsNormResult = Norm(rhs);
    std::cout << "Scale factor " << scalefactor << " resnorm " << RhsNormResult << std::endl;
    scalefactor *= 0.5;
    iter++;
  } while (RhsNormResult > RhsNormPrev && iter < 20);
  if(fabs(RhsNormResult - RhsNormPrev)<1.e-6 )
  {
    converging=false;
  }
  else
  {
    converging=true;
  }
  scalefactor *= 2.;
	return scalefactor;
	
}//void

bool ToolsTransient::FindElementAfterFracture(int &index)
{
  bool found = false;
  int nel = fmphysics->NElements();
  int iel = 0;
  for (iel = 0; iel < nel; iel++) {
    TPZCompEl *cel = fmphysics->ElementVec()[iel];
    if (!cel) {
      continue;
    }
    TPZGeoEl *gel = cel->Reference();
    if (gel->MaterialId() != globMultiFisicMatId1) {
      continue;
    }
    TPZGeoElSide side(gel,0);
    TPZGeoElSide neigh = side.Neighbour();
    while (neigh != side) {
      if (neigh.Element()->MaterialId() == globCohesiveMatIdHalf) {
        break;
      }
      neigh = neigh.Neighbour();
    }
    if (side == neigh) {
      continue;
    }
    TPZGeoElSide sidebot(gel,4);
    neigh = sidebot.Neighbour();
    while (neigh != sidebot) {
      if (globFractInputData.IsBC(neigh.Element()->MaterialId()) ) {
        break;
      }
      neigh = neigh.Neighbour();
    }
    if (sidebot != neigh) {
      continue;
    }
    
    found = true;
    break;
  }
  index = iel;
  return found;
}



void ToolsTransient::CreatePostProcessingMesh()
{
  TPZCompMesh *cmeshJustForPostProc = CMeshElasticH1ForPostProc();
  TPZMaterial * mat = cmeshJustForPostProc->FindMaterial(1);
#ifdef PlasticMC
  TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> *plasticmat = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> , TPZElastoPlasticMem> *>(mat);
#endif
#ifdef PlasticSDi
  TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> *plasticmat = dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> , TPZElastoPlasticMem> *>(mat);
  
#endif
  if (!plasticmat) {
    DebugStop();
  }
  plasticmat->GetMemory() = fCouplingMaterial1->GetMemory();
  
  
  
  //cmeshJustForPostProc->LoadReferred(fmphysics);
  //fPostprocess.SetMultiPhysicsMesh(fmphysics);
  if (fPostprocess.ReferenceCompMesh() != cmeshJustForPostProc)
  {
    
    fPostprocess.SetCompMesh(cmeshJustForPostProc);
    TPZFStructMatrix structmatrix(fPostprocess.Mesh());
    structmatrix.SetNumThreads(0);
    fPostprocess.SetStructuralMatrix(structmatrix);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZStack<std::string> PostProcVars, scalNames, vecNames;
    this->PostProcessVariables(scalNames, vecNames);
    
    for (int i=0; i<scalNames.size(); i++) {
      PostProcVars.Push(scalNames[i]);
    }
    for (int i=0; i<vecNames.size(); i++) {
      PostProcVars.Push(vecNames[i]);
    }
    //
    fPostprocess.SetPostProcessVariables(PostProcMatIds, PostProcVars);
  }
  //
  fPostprocess.TransferSolution();
}

/// Get the post processing variables
void ToolsTransient::PostProcessVariables(TPZStack<std::string> &scalNames, TPZStack<std::string> &vecNames)
{
  scalNames.Resize(0);
  vecNames.Resize(0);
  
  vecNames.Push("Displacement"); //1 2 3
  vecNames.Push("DisplacementMem"); // x y z
  vecNames.Push("YieldSurface");
  vecNames.Push("TotalPlasticStrain");

  /*
  vecNames.Push("NormalStress");
  vecNames.Push("ShearStress");
  vecNames.Push("NormalStrain");
  vecNames.Push("ShearStrain");
   */
}

void ToolsTransient::PostProcess(int resolution)
{

  std::string vtkFile = "TransientSolutionPlasticPapa.vtk";
  TPZStack<std::string> scalNames,vecNames;
  PostProcessVariables(scalNames,vecNames);

  fPostprocess.DefineGraphMesh(2,scalNames,vecNames,vtkFile);

  fPostprocess.SetStep(fPostProcessNumber);
  fPostprocess.PostProcess(resolution);

  fPostProcessNumber++;
}

// DRAFT ********************* DRAFT *********************** DRAFT *********************** DRAFT ******************** DRAFT

//AQUINATHAN EquationFilter

//    int neqpr = this->fmeshvec[1]->NEquations();
//    for (int i = 0; i < neqpr; i++) {
//      fmeshvec[1]->Solution()(i,0) = pressure*;
//    }

//AQUINATHAN EquationFilter da hat

//fmeshvec[0]->Solution()(1,0) = 0.;
//TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);

//AQUINATHAN

//TPZFMatrix<REAL> fmat(fmphysics->NEquations(),1);


/*
 REAL delta = 1.;
 REAL aj = -6378;
 TPZStack<std::pair<REAL,REAL> > mypoints;
 
 TPZFMatrix<REAL> fmat(fmphysics->NEquations(),1);
 fmat.Zero();
 TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
 MassMatrix(fmat);
 for (int i = 0; i < 30; i++) {
 aj = i * delta;
 fmat.Zero();
 an->Solution().Zero();
 an->LoadSolution();
 TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
 fmeshvec[0]->Solution()(1,0) = aj;
 TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
 an->Solution().Print("sol");
 fmphysics->Solution().Print("solmphysics");
 MassMatrix(fmat);
 REAL res = this->IterativeProcess(an,maxit,tol); // Newton Method for one time step!
 fmat.Print("fmat");
 TPZFMatrix<> thisrhs = an->Rhs();
 thisrhs += fmat;
 thisrhs.Print("rhsEqF:");
 an->StructMatrix()->EquationFilter().Reset();
 an->AssembleResidual();
 thisrhs = an->Rhs();
 thisrhs.Print("rhsbefore:");
 thisrhs += fmat;
 thisrhs.Print("rhs:");
 std::pair<REAL,REAL> mypair(aj,thisrhs(1,0));
 mypoints.Push(mypair);
 }
 std::cout << "RhsForEachAlphaJ = " << "{"
 << "{" << mypoints[0].first << "," << mypoints[0].second << "}";
 for (int i = 1; i < mypoints.NElements(); i++) {
 std::cout << ",{" << mypoints[i].first << "," << mypoints[i].second << "}";
 }
 std::cout << "};" << std::endl;
 std::cout << "ListPlot[RhsForEachAlphaJ,Joined->True]" << std::endl;
 */

// DRAFT ********************* DRAFT *********************** DRAFT *********************** DRAFT ******************** DRAFT
