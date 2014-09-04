#include "pzlog.h"
#include "TPZFracAnalysis.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzanalysis.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

TPZFracAnalysis::TPZFracAnalysis(TPZAutoPointer<TPZFracData> Data)
{
  fData = Data;
  fmeshvec.Resize(2);
  fgmesh = NULL;
  fcmeshMixed = NULL;
  for (int i = 0; i < 2; i++) {
    fmeshvec[i] = NULL;
  }
  fLastStepRhs.Redim(0, 0);
}


TPZFracAnalysis::~TPZFracAnalysis()
{
  for (int i = 0; i < 2; i++) {
    delete fmeshvec[i];
  }
  delete fcmeshMixed;
  delete fgmesh;
}

void TPZFracAnalysis::Run()
{
  
  // Malha geometrica
  fgmesh = CreateGMesh();
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  fcmeshMixed = CreateCMeshMixed();
  
  // Analysis
  bool mustOptimizeBandwidth = true;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  // Plot of first solution
  const int dim = 1;
  TPZStack<std::string> scalnames, vecnames;
  scalnames.Push("Pressure");
  scalnames.Push("Flow");
  scalnames.Push("Opening");
  an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
  an->PostProcess(0);
  
  // Solving transiente system
  SolveSistTransient(an);
  
  delete an;
}

TPZGeoMesh * TPZFracAnalysis::CreateGMesh()
{
  const int nel = fData->NelFrac();
  const REAL lfrac = fData->Lfrac();
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  
  // Nos
  const int nnodes = nel+1;
  const REAL elsize = lfrac / nel;
  gmesh->NodeVec().Resize(nnodes);
  TPZVec<REAL> coord(3,0.);
  for (int in = 0; in < nnodes; in++) {
    coord[0] = in * elsize;
    gmesh->NodeVec()[in].SetNodeId(in);
    gmesh->NodeVec()[in].SetCoord(coord);
  }
  
  // Elementos
  TPZVec<long> TopolLine(2,0);
  long index = 0;
  for (int iel = 0; iel < nel; iel++) {
    TopolLine[0] = iel;
    TopolLine[1] = iel+1;
    gmesh->CreateGeoElement(EOned, TopolLine, matIdFrac, index);
  }
  
  gmesh->BuildConnectivity();
  
  // Left
  TPZVec<long> TopolPoint(1,0);
  gmesh->CreateGeoElement(EPoint, TopolPoint, bcLeftId, index);
  
  // Right
  TopolPoint[0] = nnodes-1;
  gmesh->CreateGeoElement(EPoint,TopolPoint,bcRightId,index);
  
  gmesh->BuildConnectivity();
  
#ifdef DEBUG
  std::ofstream out("GeoMesh.vtk");
  TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
#endif
  
  return gmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshFluxH1()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int fluxorder = fData->PorderFlow();
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);

  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando H1
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(fluxorder);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  
#ifdef DEBUG
  std::ofstream out("CMeshFluxH1.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshPressureL2()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int pressureorder = fData->PorderPressure();
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando L2
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(pressureorder);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  std::set<int> materialIds,materialIdsBC;
  materialIds.insert(matIdFrac);
  materialIdsBC.insert(bcLeftId);
  materialIdsBC.insert(bcRightId);

  cmesh->AutoBuild();

  int ncon = cmesh->NConnects();
  for(int i=0; i<ncon; i++)
  {
    TPZConnect &newnod = cmesh->ConnectVec()[i];
    newnod.SetLagrangeMultiplier(1);
  }
  
  for (int i = 0; i < cmesh->Solution().Rows(); i++) {
    cmesh->Solution()(i,0) = 1.01 * fData->SigmaConf();
  }
  
#ifdef DEBUG
  std::ofstream out("CMeshPressureL2.txt");
  cmesh->Print(out);
#endif
  
  return cmesh;
}

TPZCompMesh * TPZFracAnalysis::CreateCMeshMixed()
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(fgmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  mat->SetSimulationData(fData);
  cmesh->InsertMaterialObject(mat);
  
  // Condicao de contorno na esquerda
  val2(0,0) = fData->Q();
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  val2(0,0) = fData->SigmaConf();
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando Multifisico
  cmesh->SetDimModel(1);
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->AutoBuild();
  
  // Transferindo para a multifisica
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, cmesh);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, cmesh);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, cmesh);
  
#ifdef DEBUG
  std::ofstream out("CMeshMultiPhysics.txt");
  cmesh->Print(out);
#endif

  return cmesh;
}

void TPZFracAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out)
{
	int iter = 0;
	REAL error = 1.e10;
  const REAL tol = 1.e-8;
  const int numiter = 50;
  
  fData->SetCurrentState();
	int numeq = an->Mesh()->NEquations();
	
	TPZFMatrix<STATE> prevsol(an->Solution());
  TPZFMatrix<STATE> SoliterK(prevsol);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
  
  an->Assemble();
  an->Rhs() += fLastStepRhs;
  TPZAutoPointer< TPZMatrix<REAL> > matK;
  
	while(error > tol && iter < numiter) {
		
		an->Solve(); // o an->Solution() eh o deltaU aqui
    SoliterK = prevsol - an->Solution();
		REAL normDeltaSol = Norm(an->Solution());
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      matK=an->Solver().Matrix();
      matK->Print("matK = ", sout,EMathematicaInput);
      an->Solution().Print("DeltaX = ", sout,EMathematicaInput);
      SoliterK.Print("Xk = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    an->LoadSolution(SoliterK); // Aqui o an->Solution() eh o U
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
    an->Assemble();
    an->Rhs() += fLastStepRhs;
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
      std::stringstream sout;
      an->Rhs().Print("Res = ", sout,EMathematicaInput);
      LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    double NormResLambda = Norm(an->Rhs());
		double norm = normDeltaSol;
		out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << std::endl;
    
		if(norm < tol) {
			out << "\nTolerancia atingida na iteracao : " << (iter+1) << std::endl;
			out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << std::endl << std::endl;

		}
    else if( (norm - error) > 1.e-9 ) {
      out << "\nDivergent Method\n";
    }
    
		error = norm;
		iter++;
    prevsol = SoliterK;
		out.flush();
	}
  
  if (error > tol) {
    DebugStop(); // Metodo nao convergiu!!
  }
  
}

void TPZFracAnalysis::AssembleLastStep(TPZAnalysis *an)
{
  fData->SetLastState();
  an->Assemble();
  fLastStepRhs = an->Rhs();
}

void TPZFracAnalysis::SolveSistTransient(TPZAnalysis *an)
{
  
  bool propagate = false, mustStop = false;
  while (mustStop == false && propagate == false) {
    
    AssembleLastStep(an);
    IterativeProcess(an, std::cout);
    
    fData->SetNextTime();
    
    const int dim = 1;
    TPZStack<std::string> scalnames, vecnames;
    scalnames.Push("Pressure");
    scalnames.Push("Flow");
    scalnames.Push("Opening");
    an->DefineGraphMesh(dim, scalnames, vecnames, fData->PostProcessFileName());
    an->PostProcess(0);
    
    REAL peteleco = 1.E-8;
    if( fData->Time() > (fData->TotalTime() - peteleco) )
    {
      mustStop = true;
    }
  }
}