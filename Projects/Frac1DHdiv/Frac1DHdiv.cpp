#include "pzlog.h"
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
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif

TPZGeoMesh * CreateGMesh(const REAL lfrac, const int nel);
TPZCompMesh * CreateCMeshFluxH1(TPZGeoMesh *gmesh);
TPZCompMesh * CreateCMeshPressureL2(TPZGeoMesh *gmesh);
TPZCompMesh * CreateCMeshMixed(TPZGeoMesh *gmesh);
void IterativeProcess(TPZVec<TPZCompMesh*> &meshvec,TPZCompMesh *cmeshMixed, TPZAnalysis *an, std::ostream &out);

int main()
{
#ifdef LOG4CXX
  InitializePZLOG();
  
#endif
  
  // Parametros
  const REAL lfrac = 1.;
  const int nel = 20;
  
  // Malha geometrica
  TPZGeoMesh *gmesh = CreateGMesh(lfrac, nel);
  
  // Malhas computacionais - FluxoH1, PressaoL2, Multifisica para sistema misto quente
  TPZCompMesh *cmeshFlux = CreateCMeshFluxH1(gmesh);
  TPZCompMesh *cmeshPressure = CreateCMeshPressureL2(gmesh);
  TPZCompMesh *cmeshMixed = CreateCMeshMixed(gmesh);
  
  // Criando vetor de malhas
  TPZVec<TPZCompMesh *> meshvec(2);
  meshvec[0] = cmeshFlux;
  meshvec[1] = cmeshPressure;
  
  // Transferindo para a multifisica
	TPZBuildMultiphysicsMesh::AddElements(meshvec, cmeshMixed);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmeshMixed);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmeshMixed);
  
  // Analysis
  bool mustOptimizeBandwidth = false;
  TPZAnalysis *an = new TPZAnalysis(cmeshMixed,mustOptimizeBandwidth);
  TPZSkylineNSymStructMatrix skyl(cmeshMixed);
  TPZStepSolver<STATE> step;
  step.SetDirect(ELU);
  an->SetSolver(step);
  an->SetStructuralMatrix(skyl);
  
  IterativeProcess(meshvec, cmeshMixed, an, std::cout);
  
  const int dim = 1;
  TPZStack<std::string> scalnames, vecnames;
  std::string plotfile = "1DMixedResults.vtk";
  scalnames.Push("Pressure");
  scalnames.Push("Flow");
  an->DefineGraphMesh(dim, scalnames, vecnames, plotfile);
  an->PostProcess(0);
  
  return 0;
}

TPZGeoMesh * CreateGMesh(const REAL lfrac, const int nel)
{
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

TPZCompMesh * CreateCMeshFluxH1(TPZGeoMesh *gmesh)
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int fluxorder = 4;
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  mat->SetTimeStep(TPZMaterial::gBigNumber); // bignumber para sumir com o transiente
  mat->SetViscosity(1./12.);
  mat->SetTime(0.);
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

TPZCompMesh * CreateCMeshPressureL2(TPZGeoMesh *gmesh)
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  const int pressureorder = 3;
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  mat->SetTimeStep(TPZMaterial::gBigNumber); // bignumber para sumir com o transiente
  mat->SetViscosity(1./12.);
  mat->SetTime(0.);
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

TPZCompMesh * CreateCMeshMixed(TPZGeoMesh *gmesh)
{
  // Definicao de ids e tipos
  const int matIdFrac = 1, bcLeftId = -1, bcRightId = -2;
  const int typeFlux = 0, typePressure = 1;
  TPZFMatrix<> val1(1,1,0.), val2(1,1,0.);
  
  // Malha computacional
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  
  // Material da fratura
  TPZMatfrac1dhdiv *mat = new TPZMatfrac1dhdiv(matIdFrac);
  mat->SetTimeStep(TPZMaterial::gBigNumber); // bignumber para sumir com o transiente
  mat->SetViscosity(1./12.);
  mat->SetTime(0.);
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

void IterativeProcess(TPZVec<TPZCompMesh*> &meshvec,TPZCompMesh *cmeshMixed, TPZAnalysis *an, std::ostream &out)
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
    
    SoliterK.Print("SoliterK");
    
    SoliterK = prevsol - an->Solution();
    
    SoliterK.Print("SoliterKafter");
    
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
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, cmeshMixed);
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
		out << "Iteracao n : " << (iter+1) << " : normas |Delta(Un)| e |Delta(rhs)| : " << normDeltaSol << " / " << NormResLambda << endl;
    
    
		if(norm < tol) {
			out << "\nTolerancia atingida na iteracao : " << (iter+1) << endl;
			out << "\n\nNorma da solucao |Delta(Un)|  : " << norm << endl << endl;
			
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