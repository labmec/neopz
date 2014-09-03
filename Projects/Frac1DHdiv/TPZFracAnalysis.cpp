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
}


TPZFracAnalysis::~TPZFracAnalysis()
{
  
}

void TPZFracAnalysis::Run()
{
  // Parametros
  const int nel = 20;
  
  // Malha geometrica
  fgmesh = CreateGMesh(nel);
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  fmeshvec[0] = CreateCMeshFluxH1();
  fmeshvec[1] = CreateCMeshPressureL2();
  fcmeshMixed = CreateCMeshMixed();
  
  // Transferindo para a multifisica
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fcmeshMixed);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fcmeshMixed);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fcmeshMixed);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(fcmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(fcmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  IterativeProcess(an, std::cout);
  
  const int dim = 1;
  TPZStack<std::string> scalnames, vecnames;
  std::string plotfile = "1DMixedResults.vtk";
  scalnames.Push("Pressure");
  scalnames.Push("Flow");
  an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
  an->PostProcess(0);
}

TPZGeoMesh * TPZFracAnalysis::CreateGMesh(const int nel)
{
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
  val2(0,0) = 1.;
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando H1
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(fluxorder);
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  
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
  val2(0,0) = 1.;
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando L2
  cmesh->SetDimModel(1);
  cmesh->SetDefaultOrder(pressureorder);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  std::set<int> materialIds;
  materialIds.insert(matIdFrac);
  cmesh->AutoBuild(materialIds);
  
  int ncon = cmesh->NConnects();
  for(int i=0; i<ncon; i++)
  {
    TPZConnect &newnod = cmesh->ConnectVec()[i];
    newnod.SetLagrangeMultiplier(1);
  }
  
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
  val2(0,0) = 1.;
  TPZBndCond * bcLeft = mat->CreateBC(mat, bcLeftId, typeFlux, val1, val2);
 	cmesh->InsertMaterialObject(bcLeft);
  
  val2(0,0) = 10;
  TPZBndCond * bcRight = mat->CreateBC(mat, bcRightId, typePressure, val1, val2);
 	cmesh->InsertMaterialObject(bcRight);
  
  // Setando Multifisico
  cmesh->SetDimModel(1);
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->AutoBuild();
  return cmesh;
}

void TPZFracAnalysis::IterativeProcess(TPZAnalysis *an, std::ostream &out)
{
	int iter = 0;
	REAL error = 1.e10;
  const REAL tol = 1.e-8;
  const int numiter = 20;
  
	int numeq = an->Mesh()->NEquations();
	
	TPZFMatrix<STATE> prevsol(an->Solution());
  TPZFMatrix<STATE> SoliterK(prevsol);
	if(prevsol.Rows() != numeq) prevsol.Redim(numeq,1);
  
  an->Assemble();
  
  TPZAutoPointer< TPZMatrix<REAL> > matK;
  
	while(error > tol && iter < numiter) {
		
		an->Solve();
    
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
    
    an->LoadSolution(SoliterK);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fcmeshMixed);
    an->Assemble();
    
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
			
		} else
			if( (norm - error) > 1.e-9 ) {
				out << "\nDivergent Method\n";
			}
		error = norm;
		iter++;
		out.flush();
    
    prevsol = SoliterK;
	}
  
  if (error > tol) {
    DebugStop(); // Metodo nao convergiu!!
  }
  
}