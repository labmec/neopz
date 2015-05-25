
#include "tpzautopointer.h"
#include "pzlog.h"
#include <time.h>

#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternTools.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"

#include "tpzhierarquicalgrid.h"
#include "TPZVTKGeoMesh.h"
#include "TPZMonoPhaseWell.h"
#include "pzfstrmatrix.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzbndcond.h"

#include "TPZVTKGeoMesh.h"
#include "TPZReadGIDGrid.h"

#include "pznonlinanalysis.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"



#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.WellFlow"));
#endif

void ParametricfunctionS(const TPZVec<STATE> &par, TPZVec<STATE> &X);
TPZGeoMesh * WellMesh(REAL s, REAL ds, int nelements);
void PrintGeoMesh(TPZGeoMesh * gmesh);

TPZCompMesh * CmeshFlux(int qorder, TPZGeoMesh * gmesh);
TPZCompMesh * CmeshPressure(int porder, TPZGeoMesh * gmesh);
TPZCompMesh * CmeshMixed(TPZGeoMesh * gmesh);

int main()
{
    std::string dirname = PZSOURCEDIR;
    gRefDBase.InitializeUniformRefPattern(EOned);
    
    // Geometry of well
    REAL s= 0.0;
    REAL ds = 2.0;
    int nelements = 100;
    TPZGeoMesh * gmesh = WellMesh(s, ds, nelements);
    PrintGeoMesh(gmesh);
    
    // Computational Mesh

    int q = 1;
    int p = 1;
    
    TPZManVector<TPZCompMesh *>meshvec(2);
    
    meshvec[0] = CmeshFlux(q,gmesh);
    meshvec[1] = CmeshPressure(p,gmesh);
    
    TPZCompMesh * cmeshwell = CmeshMixed(gmesh);
    
    // Transferindo para a multifisica
    TPZBuildMultiphysicsMesh::AddElements(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::AddConnects(meshvec, cmeshwell);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, cmeshwell);
    
#ifdef DEBUG
    std::ofstream dumpfile("ComputationaMeshMultiphysic.txt");
    cmeshwell->Print(dumpfile);
#endif
    
    
    // Setting up the analysis configuration
    int numofThreads = 0;
    TPZNonLinearAnalysis * an = new TPZNonLinearAnalysis(cmeshwell,std::cout);
    TPZSkylineNSymStructMatrix skylnsym(cmeshwell);
    TPZStepSolver<STATE> step;
    skylnsym.SetNumThreads(numofThreads);
    step.SetDirect(ELU);
    an->SetStructuralMatrix(skylnsym);
    an->SetSolver(step);

    
    // Initialize solution
    
    
    
    // Time-foward procedure
    
    
    
    
  
    std::cout << " Process complete normally." << std::endl;
    return 0;
}

TPZCompMesh * CmeshMixed(TPZGeoMesh * gmesh)
{
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    // Material medio poroso
    TPZMonoPhaseWell * mat = new TPZMonoPhaseWell(wellId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    val2(0,0) = 0.0;
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    
    // Bc Top
    val2(0,0) = 0.0;
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typePressure, val1, val2);

    cmesh->InsertMaterialObject(bcBottom);
    cmesh->InsertMaterialObject(bcTop);
    
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsMultiphysicElem();
    cmesh->AutoBuild();
    
    return cmesh;
}

TPZCompMesh * CmeshFlux(int qorder, TPZGeoMesh * gmesh)
{
    
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Malha computacional
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMonoPhaseWell * mat = new TPZMonoPhaseWell(wellId);
    cmesh->InsertMaterialObject(mat);
    
    // Bc Bottom
    TPZBndCond * bcBottom = mat->CreateBC(mat, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Top
    TPZBndCond * bcTop = mat->CreateBC(mat, topId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    
    // Setando Hdiv
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(qorder);
    cmesh->SetAllCreateFunctionsContinuous();
    
    
    cmesh->AutoBuild();
    
    
#ifdef DEBUG
    std::ofstream out("cmeshFlux.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
}

TPZCompMesh * CmeshPressure(int porder, TPZGeoMesh * gmesh)
{
    
    int dim = 1;
    int wellId = 1;
    int bottomId = 2;
    int topId = 3;
    
    const int typeFlux = 0, typePressure = 1;
    TPZFMatrix<STATE> val1(3,2,0.), val2(3,1,0.);
    
    // Computational mesh
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    
    TPZMonoPhaseWell * material = new TPZMonoPhaseWell(wellId);
    cmesh->InsertMaterialObject(material);
    
    // Bc Bottom
    TPZMaterial * bcBottom = material->CreateBC(material, bottomId, typeFlux, val1, val2);
    cmesh->InsertMaterialObject(bcBottom);
    
    // Bc Top
    TPZBndCond * bcTop = material->CreateBC(material, topId, typePressure, val1, val2);
    cmesh->InsertMaterialObject(bcTop);
    
    // Setando L2
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild();
    
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
#ifdef DEBUG
    std::ofstream out("cmeshPress.txt");
    cmesh->Print(out);
#endif
    
    return cmesh;
    
}






void ParametricfunctionS(const TPZVec<STATE> &par, TPZVec<STATE> &X)
{
    X[0] = 0.0;
    X[1] = 0.0;
    X[2] = par[0];
}

void PrintGeoMesh(TPZGeoMesh * gmesh)
{
    
    
#ifdef DEBUG
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    gmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh,Dummyfile, true);
    
#endif
}

TPZGeoMesh * WellMesh(REAL s, REAL ds,int nelements)
{
    
    // Creating a 0D element to be extruded
    TPZGeoMesh * GeoMesh0D = new TPZGeoMesh;
    GeoMesh0D->NodeVec().Resize(1);
    TPZGeoNode Node;
    TPZVec<REAL> coors(3,0.0);
    Node.SetCoord(coors);
    Node.SetNodeId(0);
    GeoMesh0D->NodeVec()[0]=Node;
    
    TPZVec<long> Topology(1,0);
    int elid=0;
    int matid=1;
    
    new TPZGeoElRefPattern < pzgeom::TPZGeoPoint >(elid,Topology,matid,*GeoMesh0D);
    GeoMesh0D->BuildConnectivity();
    GeoMesh0D->SetDimension(0);
    
    TPZHierarquicalGrid CreateGridFrom(GeoMesh0D);
    TPZAutoPointer<TPZFunction<STATE> > ParFuncX = new TPZDummyFunction<STATE>(ParametricfunctionS);
    CreateGridFrom.SetParametricFunction(ParFuncX);
    CreateGridFrom.SetFrontBackMatId(2,3);
    
    // Computing Mesh extruded along the parametric curve -> ParametricfunctionS
    TPZGeoMesh * GeoMesh1D = CreateGridFrom.ComputeExtrusion(s, ds, nelements);
    
    return GeoMesh1D;
    
}