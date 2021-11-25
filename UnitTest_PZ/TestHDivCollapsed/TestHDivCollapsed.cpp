// ----- PZ includes -----
#include <TPZGenGrid3D.h>
#include <pzgmesh.h>
#include <TPZVTKGeoMesh.h>
#include <pzgeoquad.h>
#include <tpzgeoelrefpattern.h>

#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <pzbuildmultiphysicsmesh.h>
#include "pzshapequad.h"
#include "TPZCompElHDivCollapsed.h"

#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzstepsolver.h>
#include "pzfstrmatrix.h"

#include <TPZNullMaterial.h>
#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "DarcyFlow/TPZMixedDarcyFractureFlow.h"

#include <pzlog.h>

// ----- Unit test includes -----
#include <catch2/catch.hpp>

// ----- Functions -----
void TestHdivCollapsed(const bool& isRefMesh);
TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh);
void SplitConnectsAtInterface(TPZCompElSide& compside);
TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh);
TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh);
void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh);
void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh);

using namespace std;

enum EMatid {ENone, EDomainNotUSE, EInlet, EOutlet, ENoflux, EPressure, EIntersection, EIntersectionEnd, EVolume, EFaceBCPressure};
int globFracID = 10;

//-------------------------------------------------------------------------------------------------
//   __  __      _      _   _   _     
//  |  \/  |    / \    | | | \ | |
//  | |\/| |   / _ \   | | |  \| |
//  | |  | |  / ___ \  | | | |\  |
//  |_|  |_| /_/   \_\ |_| |_| \_|
//-------------------------------------------------------------------------------------------------
// ----- Test cases -----
// ---- Test 0 ----
TEST_CASE("1_frac_element","[hdivcollapsed]"){
    TestHdivCollapsed(false);
}
TEST_CASE("uniformly_refined_mesh","[hdivcollapsed]"){
    TestHdivCollapsed(true);
}
/* Left in case needs some serious debugging. Catch2 does not stop in debugstops
int main(int argc, char* argv[]){
    const bool isRefMesh = true;
    TestHdivCollapsed(isRefMesh);
    
    return 0;
}
 */

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void TestHdivCollapsed(const bool& isRefMesh){
    
    // ----- dimension of the problem -----
    constexpr int dim{3};
    constexpr int pOrder{1};

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif

    // ----- Create Geo Mesh -----
    const TPZVec<REAL> minX = {-1.,-1.,-1.};
    const TPZVec<REAL> maxX = {1.,1.,1.};
    const TPZVec<int> nelDiv = {1,1,2};
    const MMeshType elType = MMeshType::EHexahedral;
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,elType);
    TPZGeoMesh *gmesh = gen3d.BuildVolumetricElements(EVolume);
    
    // ----- Fracture element -----
    int64_t index;
    TPZManVector<int64_t,2> nodesId = {4,5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
    nodesId = {5,7};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
    nodesId = {7,6};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
    nodesId = {6,4};
    new TPZGeoElRefPattern<pzgeom::TPZGeoLinear>(nodesId,EPressure,*gmesh,index);
    TPZManVector<int64_t,4> nodesIdVec = {4,6,7,5};
    new TPZGeoElRefPattern<pzgeom::TPZGeoQuad>(nodesIdVec,globFracID,*gmesh,index);
    gmesh = gen3d.BuildBoundaryElements(EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure, EFaceBCPressure);
    gmesh->BuildConnectivity();
    
    if (isRefMesh) {
        gRefDBase.InitializeUniformRefPattern(ECube);
        gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
        gRefDBase.InitializeUniformRefPattern(EOned);
        for (auto gel : gmesh->ElementVec()){
            TPZManVector<TPZGeoEl*,10> children;
            gel->Divide(children);
        }
        gmesh->BuildConnectivity();
    }
    
    // ----- Print gmesh -----
    std::ofstream outgmesh("GeoMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, outgmesh);
    
    // ----- Create Flux mesh -----
    TPZCompMesh * cmeshflux = FluxCMesh(dim,pOrder,gmesh);

    // ----- Pressure mesh -----
    TPZCompMesh * cmeshpressure = PressureCMesh(dim,pOrder,gmesh);

    // ----- Multiphysics mesh -----
    TPZManVector< TPZCompMesh *, 2> meshvector(2);
    meshvector[0] = cmeshflux;
    meshvector[1] = cmeshpressure;
    TPZCompMesh * cmesh = MultiphysicCMesh(dim,pOrder,meshvector,gmesh);
    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;

    // ----- Solve system -----
    TPZLinearAnalysis an(cmesh,false);
    SolveProblemDirect(an,cmesh);

    // ----- Print results -----
    PrintResultsMultiphysic(dim,meshvector,an,cmesh);

    // ----- Compute integral of pressure over domain and compare with analytical solution -----
    const std::string varname = "Pressure";
    std::set<int> matids;
    matids.insert(EVolume);
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences(); // compute integral in the multiphysics mesh
    TPZVec<STATE> vecint = cmesh->Integrate(varname, matids);
    if (vecint.size() != 1){
        DebugStop();
    }
    const STATE integratedpressure = vecint[0];
    std::cout << "\nintegral of pressure  = " << integratedpressure << std::endl;
    
    // Since pressure is unit constant and the domain volume is 2*2*2=8, the integral of pressure
    // should be 1*8 = 8.
    REQUIRE( integratedpressure == Approx( 8.0 ) ); // Approx is from catch2 lib
    
    delete gmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *FluxCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    // ===> Creating compmesh
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);

    // ===> Domain material
    const int nstate = 1;
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EVolume,dim,nstate);
    cmesh->InsertMaterialObject(mat);

    // ===> Boundary conditions
    TPZNullMaterial<> *matbc = new TPZNullMaterial<>(EFaceBCPressure,dim-1,nstate);
    cmesh->InsertMaterialObject(matbc);
    
    // ===> Fracture material
    TPZNullMaterial<> *matfrac = new TPZNullMaterial<>(globFracID,dim-1,nstate);
    cmesh->InsertMaterialObject(matfrac);
    
    // ===> Fracture boundary conditions
    TPZNullMaterial<> *matbcfrac = new TPZNullMaterial<>(EPressure,dim-2,nstate);
    cmesh->InsertMaterialObject(matbcfrac);
    
    // ===> Creating space for 3D elements
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false); // we need to disconnect by hand at the fracture location later
    cmesh->SetDefaultOrder(pOrder);
    std::set<int> buildmatids;
    buildmatids.insert(EVolume);
    cmesh->AutoBuild(buildmatids);
    
    // Print flux mesh
    std::ofstream myfile("FluxMeshOnlyVol.txt");
    cmesh->Print(myfile);
    
    // ===> Creating fracture element
    CreateFractureHDivCollapsedEl(cmesh);
    cmesh->CleanUpUnconnectedNodes();
    
    // ===> Create BCs for fracture
    buildmatids.clear();
    buildmatids.insert(EPressure);
    gmesh->ResetReference();
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        if (cel->Reference()->MaterialId() == globFracID) {
            cel->LoadElementReference();
        }
    }
    cmesh->SetDimModel(dim-1);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(buildmatids);
    
    
    // ===> Create BCs for 3D domain
    buildmatids.clear();
    buildmatids.insert(EFaceBCPressure);
    gmesh->ResetReference();
    for (auto cel : cmesh->ElementVec()) {
        if (!cel) {
            continue;
        }
        if (cel->Reference()->MaterialId() == EVolume) {
            cel->LoadElementReference();
        }
    }
    cmesh->SetDimModel(dim);
    cmesh->SetAllCreateFunctionsHDiv();
    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    cmesh->AutoBuild(buildmatids);
    cmesh->CleanUpUnconnectedNodes();
    
    // ===> Fix blocks
    cmesh->SetDimModel(dim);
    cmesh->InitializeBlock();

    // Print flux mesh
    std::ofstream myfile2("FluxMeshEnd.txt");
    cmesh->Print(myfile2);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void CreateFractureHDivCollapsedEl(TPZCompMesh* cmesh) {
    
    // ===> Create the fracture hdivcollapsed elements
    // These should be unconnected with the 3D elements so we reset references in the gmesh
    // They are, however, connected (with connects) between themselves
    TPZGeoMesh* gmesh = cmesh->Reference();
    gmesh->ResetReference(); // So it does not try to use connects of neighbor elements
    const int gmeshdim = gmesh->Dimension();
    cmesh->SetDimModel(gmeshdim-1);
    for(auto gel : gmesh->ElementVec()) {
        const int matid = gel->MaterialId();
        if (matid != globFracID) continue;
        const int hassubel = gel->HasSubElement();
        if (hassubel) { // the mesh can be uniformly refined
            continue;
        }
        int64_t index;
        TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* hdivcollapsed = new TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>(*cmesh,gel,index);
    }
    
    // ===> Set HdivCollapsed elements connection with adjacent 3d elements
    // This is done through top and bottom connects of the hdivcollapsed el.
    // Since bottom connect's default sideorient is negative one (-1),
    // we should connect it with the 3D element side that has sideorient positive one (+1)
    cmesh->LoadReferences(); // So the hdivcollapsed element can see the 3D the flux mesh adjacent elements to it
    for(auto cel : cmesh->ElementVec()) {
        if (!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        const int matid = gel->MaterialId();
        if (matid != globFracID) continue;
        int64_t index;
        TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* hdivcollapsed = dynamic_cast<TPZCompElHDivCollapsed<pzshape::TPZShapeQuad>* >(cel);
        if (!hdivcollapsed) {
            DebugStop();
        }
        
        int nconnects = hdivcollapsed->NConnects();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neigh = gelside.Neighbour();
        int icon = 0;
        while(neigh != gelside){
            TPZGeoEl* gelnei = neigh.Element();
            if (!gelnei || gelnei->Dimension() != 3 || gelnei->MaterialId() != EVolume) {
                neigh++;
                continue;
            }
            TPZCompElSide compside = TPZCompElSide(gelnei->Reference(), neigh.Side());
            if (icon == 0){
                // split the connect shared by two adjacent 3d elements if it is the first time
                SplitConnectsAtInterface(compside);
            }
            
            const int conindex = compside.ConnectIndex();
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement*>(compside.Element());
            const int sideorient = intel->GetSideOrient(compside.Side());
            if (sideorient == 1) { // this is a bottom connection (bottom is -1)
                hdivcollapsed->SetConnectIndex(nconnects, conindex);
                icon++;
            }
            else if (sideorient == -1){
                hdivcollapsed->SetConnectIndex(nconnects+1, conindex);
                icon++;
            }
            else{
                DebugStop();
            }
            neigh++;
        }
        if (icon != 2) {
            DebugStop();
        }
    }
    
    cmesh->ExpandSolution();
//    std::ofstream out("cmeshaftercollapsed.txt");
//    cmesh->Print(out);
    
    cmesh->SetDimModel(gmeshdim);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void SplitConnectsAtInterface(TPZCompElSide& compside) {
    
    TPZCompMesh* fluxmesh = compside.Element()->Mesh();
    
    // ===> Find 3D element neighbor
    TPZCompElSide compsideright;
    TPZGeoElSide gleft(compside.Reference());
    TPZGeoElSide neigh = gleft.Neighbour();
    int icon = 0;
    while(neigh != gleft){
        TPZGeoEl* gelnei = neigh.Element();
        if (!gelnei || gelnei->Dimension() != 3 || gelnei->MaterialId() != EVolume) {
            neigh++;
            continue;
        }
        compsideright = TPZCompElSide(gelnei->Reference(), neigh.Side());
        break;
    }
    
    compside.SplitConnect(compsideright);
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZCompMesh *PressureCMesh(int dim, int pOrder, TPZGeoMesh *gmesh)
{
    // ===> Creating compmesh
    gmesh->ResetReference();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);

    // ===> Domain material
    const int nstate = 1;
    TPZNullMaterial<> *mat = new TPZNullMaterial<>(EVolume,dim,nstate);
    cmesh->InsertMaterialObject(mat);

    // ===> Fracture material
    TPZNullMaterial<> *matfrac = new TPZNullMaterial<>(globFracID,dim-1,nstate);
    cmesh->InsertMaterialObject(matfrac);
    
    // ===> P order distinction
    if(pOrder == 0)
    {
        cmesh->SetAllCreateFunctionsDiscontinuous();
    }
    else
    {
        cmesh->SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }
    
    // ===> Creating space for 3D elements
    std::set<int> buildmatids;
    buildmatids.insert(EVolume);
    buildmatids.insert(globFracID);
    cmesh->AutoBuild(buildmatids);
        
    // ===> Set lagrange multiplier for order of assembly
    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++){
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(1);
    }
    
    // Print pressure mesh
    std::ofstream myfile("PressureMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

TPZMultiphysicsCompMesh *MultiphysicCMesh(int dim, int pOrder, TPZVec<TPZCompMesh *> meshvector,TPZGeoMesh * gmesh)
{
    // ===> Creating mp cmesh
    gmesh->ResetReference();
    auto cmesh = new TPZMultiphysicsCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
    cmesh->SetDimModel(dim);
    
    // ===> Domain mat
    auto mat = new TPZMixedDarcyFlow(EVolume, dim);
    mat->SetConstantPermeability(1.);
    cmesh->InsertMaterialObject(mat);
    
    // ===> Fracture mat
    auto matfrac = new TPZMixedDarcyFractureFlow(globFracID, dim);
    matfrac->SetConstantPermeability(1.);
    cmesh->InsertMaterialObject(matfrac);
        
    // ===> Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,0.);
    TPZManVector<STATE> val2(1,1.);
    // domain bcs
    auto * BCond0 = mat->CreateBC(mat, EFaceBCPressure, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond0);
    // frac bcs
    auto * BCond1 = mat->CreateBC(mat, EPressure, 0, val1, val2);
    cmesh->InsertMaterialObject(BCond1);
    
    // ===> Build space
    TPZManVector<int> active(2,1);
    cmesh->BuildMultiphysicsSpace(active, meshvector);

//    Prints Multiphysics mesh
    std::ofstream myfile("MultiPhysicsMesh.txt");
    cmesh->Print(myfile);

    return cmesh;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void SolveProblemDirect(TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    constexpr int nThreads{0};
//    TPZSkylineStructMatrix<STATE> matskl(cmesh);
    TPZFStructMatrix<STATE> matskl(cmesh);
    matskl.SetNumThreads(nThreads);
    an.SetStructuralMatrix(matskl);
    
    ///Setting a direct solver
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);//ELU //ECholesky // ELDLt
    an.SetSolver(step);
    
    //assembles the system
    an.Assemble();
    
  {
    std::ofstream outmat("mat.txt");
    TPZMatrixSolver<STATE>* matsol = dynamic_cast<TPZMatrixSolver<STATE>*>(an.Solver());
    matsol->Matrix()->Print("singmat=",outmat,EMathematicaInput);
    std::ofstream outrhs("rhs.txt");
    an.Rhs().Print("rhs=",outrhs,EMathematicaInput);
  }
    
    ///solves the system
    an.Solve();
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------

void PrintResultsMultiphysic(int dim, TPZVec<TPZCompMesh *> meshvector, TPZLinearAnalysis &an, TPZCompMesh *cmesh)
{
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvector, cmesh);
    TPZManVector<std::string,2> scalnames(1), vecnames(1);
    
    scalnames[0] = "Pressure";
    vecnames[0]= "Flux";
    
    int div = 0;
    std::string plotfile = "PostProcess1Frac.vtk";
    an.DefineGraphMesh(dim,scalnames,vecnames,plotfile);
    an.PostProcess(div,dim);
    
    return;
}

// ---------------------------------------------------------------------
// ---------------------------------------------------------------------
