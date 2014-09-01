#include "pzlog.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzbuildmultiphysicsmesh.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

TPZGeoMesh * CreateGMesh(const REAL lfrac, const int nel);
TPZCompMesh * CreateCMeshFluxH1(TPZGeoMesh *gmesh);
TPZCompMesh * CreateCMeshPressureL2(TPZGeoMesh *gmesh);
TPZCompMesh * CreateCMeshMixed(TPZGeoMesh *gmesh);

int main()
{
#ifdef LOG4CXX
  InitializePZLOG();
#endif
  
  // Parametros
  const REAL lfrac = 1.;
  const int nel = 10;
  
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
  
  // AQUI FAZER O SOLVESISTTRANSIENT!
  
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

  // Setando H1
  cmesh->SetAllCreateFunctionsContinuous();
  cmesh->AutoBuild();
  
  return cmesh;
}

TPZCompMesh * CreateCMeshPressureL2(TPZGeoMesh *gmesh)
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
  
  // Setando L2
  cmesh->SetAllCreateFunctionsDiscontinuous();
  cmesh->AutoBuild();
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
  
  // Setando Multifisico
  cmesh->SetAllCreateFunctionsMultiphysicElem();
  cmesh->AutoBuild();
  return cmesh;
}