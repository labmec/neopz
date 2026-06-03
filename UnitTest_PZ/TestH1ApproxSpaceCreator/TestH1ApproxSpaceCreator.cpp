//
// Created by victor on 26/09/2022.
//

// ----- PZ includes -----
#include <TPZGenGrid3D.h>
#include <TPZGenGrid2D.h>
#include <pzgmesh.h>
#include <TPZVTKGeoMesh.h>

#include <pzcmesh.h>
#include <TPZMultiphysicsCompMesh.h>
#include <TPZLinearAnalysis.h>
#include <pzskylstrmatrix.h>
#include <pzfstrmatrix.h>

#include <DarcyFlow/TPZHybridDarcyFlow.h>
#include "TPZH1ApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZAnalyticSolution.h"
#include "pzstepsolver.h"

#include <pzlog.h>
#include "TPZSimpleTimer.h"
#ifdef PZ_USING_MKL
#include "TPZSSpStructMatrix.h"
#endif

// ----- Unit test includes -----
// #define USE_MAIN

static const std::map<ProblemType, const char *> ProblemNames({{ProblemType::EElastic, "Elasticity"},
                                                               {ProblemType::EDarcy, "Darcy"}});

static const std::map<HybridizationType, const char *> HybridNames({{HybridizationType::ENone, "None"},
                                                                    {HybridizationType::EStandard, "Standard"},
                                                                    {HybridizationType::EStandardSquared, "StandardSquared"}});
#ifndef USE_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#include <catch2/catch_approx.hpp>
#endif

TPZGeoMesh *Create2DGeoMesh();

TPZGeoMesh *Create2DGeoMesh(ProblemType& pType, MMeshType &mType);

TPZGeoMesh *Create3DGeoMesh(ProblemType& pType, MMeshType &mType);

void InsertMaterials(TPZH1ApproxCreator &approxCreator);

void Refinement(TPZGeoMesh *gmesh);

void CheckNEqCondensedProb(TPZCompMesh* cmesh, TPZH1ApproxCreator& h1creator, MMeshType& elType);

void SolveSystem(TPZCompMesh* cmesh, const bool isTestKnownSol);

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZCompMesh* cmesh);

void PostProcessVTK(TPZCompMesh* cmesh, ProblemType probType);

void Test2DRigidBodyMotion(H1Family h1Fam, HybridizationType hybType, ProblemType probType, int pOrder, int plusOrder, bool IsRigidBodySpaces, bool shouldCondense);

void TestH1ApproxSpaceCreator(H1Family h1Fam, ProblemType probType, int pOrder, bool isRigidBodySpaces, MMeshType mType, int extrapOrder, bool isCondensed, HybridizationType hType, bool isRef);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, H1Family h1fam);

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann,EBCDisplacementLeft,EBCDisplacementRight};

constexpr const char* H1FamilyToChar(H1Family h1fam) {
    switch (h1fam){
        case H1Family::EH1Standard: return "EH1Standard";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}
constexpr const char* HybridizationTypeToChar(HybridizationType hType) {
    switch (hType){
        case HybridizationType::ENone: return "ENone";
        case HybridizationType::EStandard: return "EStandard";
        case HybridizationType::EStandardSquared: return "EStandardSquared";
        case HybridizationType::ESemi: return "ESemi";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}
constexpr const char* ProblemTypeToChar(ProblemType ptype) {
    switch (ptype){
        case ProblemType::EDarcy: return "EDarcy";
        case ProblemType::EElastic: return "EElastic";
        case ProblemType::EStokes: return "EStokes";
        case ProblemType::ENone: return "ENone";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item";//silences compiler warning on gcc
}

constexpr const char* MeshTypeToChar(MMeshType mtype)
{
    switch (mtype){
        case MMeshType::EQuadrilateral: return "EQuadrilateral";
        case MMeshType::ETriangular: return "ETriangular";
        case MMeshType::EHexahedral: return "EHexahedral";
        case MMeshType::ETetrahedral: return "ETetrahedral";
        default: std::invalid_argument("Unimplemented item");
    }
    return "Unimplemented item"; // silences compiler warning on gcc
}

auto exactSolDarcy = [](const TPZVec<REAL> &loc,
                        TPZVec<STATE>&u,
                        TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    
    u[0] = 1.;
    gradU.Zero();
};

auto exactSolElastic = [](const TPZVec<REAL> &loc,
                          TPZVec<STATE>&disp,
                          TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    
    disp.Fill(0.);
    disp[0] = 1.;
    disp[1] = 1.;
    
    gradU.Zero();
    gradU(0,0) = 0.;
};

TElasticity2DAnalytic gElast2d;
TElasticity3DAnalytic gElast3d;
TLaplaceExample1 gDarcy;

#ifndef USE_MAIN
#define TESTALL
TEST_CASE("2D and 3D linear solution", "[h1_space_creator_test]")
{
#ifdef TESTALL
    H1Family sType = GENERATE(H1Family::EH1Standard);
#else
    H1Family sType = GENERATE(H1Family::EH1Standard);
#endif
    SECTION("H1 family type: " + std::string(H1FamilyToChar(sType))) {
#ifdef TESTALL
        ProblemType pType = GENERATE(ProblemType::EDarcy, ProblemType::EElastic);
#else
        ProblemType pType = GENERATE(ProblemType::EElastic);
#endif
        SECTION("Problem type: " + std::string(ProblemTypeToChar(pType))) {
#ifdef TESTALL
            MMeshType mType = GENERATE(MMeshType::EQuadrilateral, MMeshType::ETriangular, MMeshType::ETetrahedral, MMeshType::EHexahedral);
#else
            MMeshType mType = GENERATE(MMeshType::EHexahedral);
#endif
            SECTION("Mesh type: " + std::string(MeshTypeToChar(mType))) {
#ifdef TESTALL
                int pOrder = GENERATE(1, 2);
#else
                int pOrder = GENERATE(1);
#endif
                SECTION("pOrder=" + std::to_string(pOrder)) {
#ifdef TESTALL
                    int extraporder = GENERATE(0, 1);
#else
                    int extraporder = GENERATE(0);
#endif
                    SECTION("extraporder=" + std::to_string(extraporder)) {
#ifdef TESTALL
                        HybridizationType hType = GENERATE(HybridizationType::ENone, HybridizationType::EStandard, HybridizationType::EStandardSquared);
#else
                        HybridizationType hType = GENERATE(HybridizationType::ENone);
#endif
                        SECTION("Hybridization type: " + std::string(HybridizationTypeToChar(hType))) {
                            bool isRBSpaces = GENERATE(false, true);
                            SECTION("isRigidBodySpaces=" + std::to_string(isRBSpaces)) {
#ifdef TESTALL
                                bool isCondensed = GENERATE(false, true);
#else
                                bool isCondensed = GENERATE(false);
#endif
                                SECTION("isCondensed=" + std::to_string(isCondensed)) {
#ifdef TESTALL
                                    bool isRef = GENERATE(false, true);
#else
                                    bool isRef = GENERATE(true);
#endif
                                    SECTION("isRef=" + std::to_string(isRef)) {
                                        TestH1ApproxSpaceCreator(sType, pType, pOrder, isRBSpaces, mType, extraporder, isCondensed, hType, isRef);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

TEST_CASE("2D Constant solution", "[h1_space_creator_test]") {
    H1Family sType = H1Family::EH1Standard;
    HybridizationType hybtype = GENERATE(HybridizationType::ENone, HybridizationType::EStandard, HybridizationType::EStandardSquared);
    SECTION("HybridizationType = " + std::string(HybridNames.at(hybtype))) {
        ProblemType pType = GENERATE(ProblemType::EElastic, ProblemType::EDarcy);
        SECTION("ProblemType = " + std::string(ProblemNames.at(pType))) {
            int pOrder = GENERATE(1, 2);
            SECTION("pOrder=" + std::to_string(pOrder)) {
                int plusOrder = GENERATE(2, 3);
                SECTION("plusOrder=" + std::to_string(plusOrder)) {
                    bool isRigidBody = GENERATE(true, false);
                    SECTION("isRigidBody=" + std::to_string(isRigidBody)) {
                        bool shouldCondense = GENERATE(true, false);
                        SECTION("shouldCondense=" + std::to_string(shouldCondense)) {
                            Test2DRigidBodyMotion(sType, hybtype, pType, pOrder, plusOrder, isRigidBody, shouldCondense);
                        }
                    }
                }
            }
        }
    }
}
 
#else
int main() {

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    HybridizationType hybtype = HybridizationType::ENone;
    Test2DRigidBodyMotion(H1Family::EH1Standard, hybtype, ProblemType::EDarcy, 2, 1, true, true);
    return 0;
}
#endif


TPZGeoMesh *Create2DGeoMesh() {

    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,2> minX = {-1.,-1.};
    const TPZManVector<REAL,2> maxX = {1.,1.};
    const TPZManVector<int,2> nelDiv = {2,2};
    TPZGenGrid2D gen2d(nelDiv,minX,maxX);
    const MMeshType elType = MMeshType::EQuadrilateral;
    gen2d.SetElementType(elType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen2d.Read(gmesh,EDomain);

    gen2d.SetBC(gmesh, 4, EBCDirichlet);
    gen2d.SetBC(gmesh, 5, EBCDirichlet);
    gen2d.SetBC(gmesh, 6, EBCDirichlet);
    gen2d.SetBC(gmesh, 7, EBCDirichlet);

    return gmesh;
}

TPZGeoMesh *Create2DGeoMesh(ProblemType& pType, MMeshType &mType) {
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,2> minX = {-1.,-1.};
    const TPZManVector<REAL,2> maxX = {1.,1.};
    const TPZManVector<int,2> nelDiv = {2,1};
    TPZGenGrid2D gen2d(nelDiv,minX,maxX);
    gen2d.SetElementType(mType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gen2d.Read(gmesh,EDomain);
    
    if(pType == ProblemType::EDarcy){
        gen2d.SetBC(gmesh, 4, EBCDirichlet);
        gen2d.SetBC(gmesh, 5, EBCDirichlet);
        gen2d.SetBC(gmesh, 6, EBCDirichlet);
        gen2d.SetBC(gmesh, 7, EBCDirichlet);
    }
    else if (pType == ProblemType::EElastic){
        const bool isrbmProb = false;
        if(isrbmProb){
            gen2d.SetBC(gmesh, 4, EBCDirichlet);
            gen2d.SetBC(gmesh, 5, EBCDirichlet);
            gen2d.SetBC(gmesh, 6, EBCDirichlet);
            gen2d.SetBC(gmesh, 7, EBCDirichlet);
        }
        else{
            gen2d.SetBC(gmesh, 4, EBCNeumann);
            gen2d.SetBC(gmesh, 5, EBCDisplacementRight);
            gen2d.SetBC(gmesh, 6, EBCNeumann);
            gen2d.SetBC(gmesh, 7, EBCDisplacementLeft);
        }
    }
    
    return gmesh;
}

TPZGeoMesh *Create3DGeoMesh(ProblemType& pType, MMeshType &mType) {
    
    // ----- Create Geo Mesh -----
    const TPZManVector<REAL,3> minX = {-1.,-1.,-1.};
    const TPZManVector<REAL,3> maxX = {1.,1.,1.};
    const TPZManVector<int,3> nelDiv = {2,1,1};
    TPZGenGrid3D gen3d(minX,maxX,nelDiv,mType);
    TPZGeoMesh* gmesh = new TPZGeoMesh;
    gmesh = gen3d.BuildVolumetricElements(EDomain);
    
    if(pType == ProblemType::EDarcy){
        gmesh = gen3d.BuildBoundaryElements(EBCNeumann,EBCNeumann,EBCNeumann,EBCDirichlet,EBCNeumann,EBCNeumann);
    }
    else if (pType == ProblemType::EElastic){
        const bool isrbmProb = false;
        if(isrbmProb){
            gmesh = gen3d.BuildBoundaryElements(EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet,EBCDirichlet);
        }
        else {
            gmesh = gen3d.BuildBoundaryElements(EBCNeumann,EBCDisplacementLeft,EBCNeumann,EBCDisplacementRight,EBCNeumann,EBCNeumann);
        }
    }
    
    return gmesh;
}

#include <Elasticity/TPZHybridElasticity2D.h>
#include <Elasticity/TPZElasticity2D.h>
#include <Elasticity/TPZHybridElasticity3D.h>
#include <Elasticity/TPZElasticity3D.h>

void InsertMaterials(TPZH1ApproxCreator &approxCreator){

    bool isdarcy = approxCreator.ProbType() == ProblemType::EDarcy;
    bool iselast = approxCreator.ProbType() == ProblemType::EElastic;

    int dim = approxCreator.GeoMesh()->Dimension();
    HybridizationType hybType = approxCreator.HybridType();
    TPZMaterialT<STATE> *mat = 0;
    int nstate = 0;
    if(hybType == HybridizationType::ENone){
        if(isdarcy) {
            auto matdarcy = new TPZDarcyFlow(EDomain, dim);
            matdarcy->SetConstantPermeability(1.);
            matdarcy->SetForcingFunction(gDarcy.ForceFunc(), 2);
            matdarcy->SetExactSol(gDarcy.ExactSolution(), 2);
            mat = matdarcy;
            nstate = 1;
        }
        else if(iselast) {
            REAL E = 1.;
            REAL nu = 0.3;
            REAL fx = 0., fy = 0.;
            if(dim == 2) {
                auto *matelast = new TPZElasticity2D(EDomain,E,nu,fx,fy);
                matelast->SetForcingFunction(gElast2d.ForceFunc(), 2);
                matelast->SetExactSol(gElast2d.ExactSolution(), 2);
                mat = matelast;
                nstate = 2;
            } else {
                TPZManVector<STATE,3> force(3,0.);
                auto *matelast = new TPZElasticity3D(EDomain,E,nu,force);
                matelast->SetForcingFunction(gElast3d.ForceFunc(), 2);
                matelast->SetExactSol(gElast3d.ExactSolution(), 2);
                mat = matelast;
                nstate = 3;
            }
        }
    }
    else {
        if(isdarcy) {
            auto matdarcy = new TPZHybridDarcyFlow(EDomain, dim);
            matdarcy->SetForcingFunction(gDarcy.ForceFunc(), 2);
            matdarcy->SetExactSol(gDarcy.ExactSolution(), 2);
            matdarcy->SetConstantPermeability(1.);
            mat = matdarcy;
            nstate = 1;
        }
        else if(iselast && dim == 2) {
            REAL E = 1.;
            REAL nu = 0.3;
            REAL fx = 0., fy = 0.;
            auto *matelast = new TPZHybridElasticity2D(EDomain,E,nu,fx,fy);
            matelast->SetForcingFunction(gElast2d.ForceFunc(), 2);
            matelast->SetExactSol(gElast2d.ExactSolution(), 2);
            mat = matelast;
            nstate = 2;
        } else if(iselast && dim == 3){
            REAL E = 1.;
            REAL nu = 0.3;
            TPZManVector<REAL,3> f(3,0.);
            auto *matelast = new TPZHybridElasticity3D(EDomain,E,nu,f);
            matelast->SetForcingFunction(gElast3d.ForceFunc(), 2);
            matelast->SetExactSol(gElast3d.ExactSolution(), 2);
            mat = matelast;
            nstate = 3;

        } else {
            DebugStop();
        }
    }
    // mat->SetForcingFunction()
    approxCreator.InsertMaterialObject(mat);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(nstate, nstate, 1.);
    TPZManVector<STATE> val2(nstate, 1.);

    TPZManVector<TPZBndCondT<STATE> *,4> BCond(4,0);
    //Dirichlet Boundary Conditions
    for(int i=0; i<4; i++) {
        BCond[i] = mat->CreateBC(mat, i+1, 0, val1, val2);
        approxCreator.InsertMaterialObject(BCond[i]);
    }

    for(int i=0; i<4; i++) {
        if(nstate == 1) {
            BCond[i]->SetForcingFunctionBC(gDarcy.ExactSolution(), 2);
        } else if (nstate == 2) {
            BCond[i]->SetForcingFunctionBC(gElast2d.ExactSolution(), 2);
        } else if (nstate == 3) {
            BCond[i]->SetForcingFunctionBC(gElast3d.ExactSolution(), 2);
        }
    }
}

//2D tests:
// Classic H1 (OK)
// Hybrid H1:
//    fCondensed = false:
//        ExtraInternalOrder = 0:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//        ExtraInternalOrder = 1:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//    fCondensed = true:
//        ExtraInternalOrder = 0:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//        ExtraInternalOrder = 1:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
// Hybrid2 H1:
//    fCondensed = false:
//        ExtraInternalOrder = 0:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//        ExtraInternalOrder = 1:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//    fCondensed = true:
//        ExtraInternalOrder = 0:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
//        ExtraInternalOrder = 1:
//            RBS = false (NOT TESTED)
//            RBS = true  (NOT TESTED)
void Test2DRigidBodyMotion(H1Family h1Fam, HybridizationType hybtype ,ProblemType probType, int pOrder, int plusOrder, bool IsRigidBodySpaces, bool shouldCondense){

    gDarcy.fExact = TLaplaceExample1::EConst;
    gElast2d.fProblemType = TElasticity2DAnalytic::EDispxy;
    gElast3d.fProblemType = TElasticity3DAnalytic::EDispxyz;
    TPZGeoMesh *gmesh = Create2DGeoMesh();

    TPZH1ApproxCreator h1Creator(gmesh);

    h1Creator.SetShouldCondense(shouldCondense);
    h1Creator.SetHybridType(hybtype);
    h1Creator.SetProbType(probType);
    h1Creator.IsRigidBodySpaces() = IsRigidBodySpaces;
    h1Creator.SetDefaultOrder(pOrder);
    h1Creator.SetExtraInternalOrder(plusOrder);
    InsertMaterials(h1Creator);

    TPZCompMesh *cmesh;
    if(h1Creator.HybridType() == HybridizationType::ENone) {
        cmesh = h1Creator.CreateClassicH1ApproximationSpace();
    } else {
        h1Creator.SetHybridizeBoundary();
        cmesh = h1Creator.CreateApproximationSpace();
    }

    {
        std::ofstream out("cmesh.txt");
        cmesh->Print(out);
    }
#ifdef USE_MAIN
    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    std::string geovtk = "geomesh.vtk";
    std::ofstream mygeomesh(geovtk);
    TPZVTKGeoMesh::PrintGMeshVTK(cmesh->Reference(),mygeomesh,true);
    TPZMultiphysicsCompMesh *mcmesh = dynamic_cast<TPZMultiphysicsCompMesh *>(cmesh);
    if(mcmesh) {
        mcmesh->Print(myfile);
        for(int i =0 ; i <mcmesh->MeshVector().size(); i++) {
            std::stringstream ss; ss << "atomicmesh" << i << ".txt";
            std::ofstream ofs(ss.str());
            mcmesh->MeshVector()[i]->Print(ofs);
        }
    }
    else
        cmesh->Print(myfile);
#endif


    std::cout <<"porder:\t" << pOrder <<  "\tplusOrder:\t" << plusOrder << "\tIsRigidBodySpaces:\t"<< IsRigidBodySpaces << "\tshouldCondense:\t" << shouldCondense <<"\n";
    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
    TPZLinearAnalysis an(cmesh,RenumType::ENone);
    TPZFStructMatrix<> strmat(cmesh);
    an.SetStructuralMatrix(strmat);
    an.Run();

    {
        std::ofstream out("cmesh.txt");
        cmesh->Print(out);
        std::ofstream out2("solution.txt");
        an.Solution().Print("Solution",out2);
    }
    CheckIntegralOverDomain(cmesh,probType,h1Fam);

#ifdef USE_MAIN
    {
        TPZMultiphysicsCompMesh *mmesh = dynamic_cast<TPZMultiphysicsCompMesh *> (cmesh);
        if (mmesh)
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mmesh->MeshVector(), cmesh);
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{0};

        TPZVec<std::string> fields = {
                "Flux",
                "Pressure"};
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        vtk.Do();
    }
    std::string txt2 = "cmeshAfterSolve.txt";
    std::ofstream myfile2(txt);
    if(mcmesh) {
        mcmesh->Print(myfile2);
        for(int i =0 ; i <mcmesh->MeshVector().size(); i++) {
            std::stringstream ss; ss << "atomicmeshAfterSolve" << i << ".txt";
            std::ofstream ofs(ss.str());
            mcmesh->MeshVector()[i]->Print(ofs);
        }
    }
    else
        cmesh->Print(myfile2);
#endif

}
void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, H1Family hdivfam){

    TPZVec<std::string> fields(2);
    switch (probType){
        case ProblemType::EDarcy:
            fields[0] = "Flux";
            fields[1] = "Pressure";
            break;
        case ProblemType::EElastic:
            fields[0] = "Stress";
            fields[1] = "Displacement";
            break;

        default:
            DebugStop();
            break;
    }

    std::set<int> matids = {EDomain};
    cmesh->Reference()->ResetReference();
    TPZVec<STATE> vecint = cmesh->Integrate(fields[0], matids);

    std::cout << "\n--------------- Integral of Flux --------------" <<  std::endl;
    std::cout << "Number of components = " << vecint.size() <<  std::endl;
    for (int i = 0; i < vecint.size(); i++)
    {
        std::cout << "Integral(" << i << ") = "  << vecint[i] << std::endl;
#ifndef USE_MAIN
        REQUIRE(fabs(vecint[i]) < 1.e-4);
#endif
    }
    std::cout << std::endl;

    TPZVec<STATE> vecintp = cmesh->Integrate(fields[1], matids);
    std::cout << "\n--------------- Integral of Primal --------------" <<  std::endl;
    int ncomp = vecintp.size();
    if(ncomp > cmesh->Dimension()) ncomp = 2;
    std::cout << "Number of components = " << ncomp <<  std::endl;
    for (int i = 0; i < ncomp; i++)
    {
        std::cout << "Integral(" << i << ") = "  << vecintp[i] << std::endl;
    }
    for (int i = 0; i < ncomp; i++)
    {
#ifndef USE_MAIN
        REQUIRE(fabs(vecintp[i]) == Catch::Approx( 4.0 ).margin(1.e-4));
#endif
    }
    std::cout << std::endl;
}

void Refinement(TPZGeoMesh *gmesh){
    // children[0]->Divide(children);
    const int nel = gmesh->NElements();
    // for (int i = 0; i < nel; i++){
    int dim = gmesh->Dimension();
    TPZManVector<TPZGeoEl*,10> children;
    if (gmesh->ElementVec()[0]->Dimension() == dim) gmesh->ElementVec()[0]->Divide(children);
    for(int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(gel->Dimension() != dim-1) continue;
        if(gel->HasSubElement()) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.Neighbour();
        if(neighbour.HasSubElement()) {
            TPZManVector<TPZGeoEl*,10> children;
            gel->Divide(children);
        }
    }
}

void CheckNEqCondensedProb(TPZCompMesh* cmesh, TPZH1ApproxCreator& h1creator, MMeshType& elType) {
    
    std::cout << "\n--------------- Checking number of equations --------------" <<  std::endl;
    TPZMultiphysicsCompMesh* mpcmesh = dynamic_cast<TPZMultiphysicsCompMesh*>(cmesh);
    if (!mpcmesh) {
        std::cout << "Hybrid H1 must have a multiphysics mesh.\n";
        return;
    }
    TPZCompMesh* pmesh = mpcmesh->MeshVector()[1];
    const int lagmatid = h1creator.HybridData().fLagrangeMatId;
    const int dim = mpcmesh->Dimension();
    
    int nlag = 0;
    for(auto cel : pmesh->ElementVec()) {
        if(!cel) continue;
        TPZGeoEl* gel = cel->Reference();
        if(!gel) DebugStop();
        const int gelmatid = gel->MaterialId();
        if (gelmatid == lagmatid) {
            ++nlag;
        }
    }
    
    const int nEquations = mpcmesh->NEquations();
    const int pOrder = mpcmesh->GetDefaultOrder();
    
    int nstate = 0;
    if (h1creator.ProbType()==ProblemType::EDarcy){
        nstate = 1;
    } else if (h1creator.ProbType()==ProblemType::EElastic) {
        nstate = dim;
    } else {
        DebugStop();
    }
    
    int expNEquations = nlag * (pOrder+1) * nstate;
    if(elType == MMeshType::ETetrahedral){
        expNEquations = nstate * nlag * (pOrder+1) * (pOrder+2) / 2;
    }
    if(elType == MMeshType::EHexahedral) {
        expNEquations = nstate * nlag * (pOrder+1) * (pOrder+1);
    }
    int extrap = h1creator.GetExtraInternalOrder();
    std::cout << "Expected equations: " << expNEquations << std::endl;
    std::cout << "Mesh equations: " << nEquations << std::endl;
    bool condition = (nEquations == expNEquations) || (extrap != 0);
    if(!condition) {
        std::cout << "I should stop\n";
    }
#ifndef USE_MAIN
    REQUIRE(condition);
#endif
}

void SolveSystem(TPZCompMesh* cmesh, const bool isTestKnownSol) {
#ifdef USE_MAIN
    constexpr int nThreads{0};
#else
    constexpr int nThreads{0};
#endif
#ifdef PZ_USING_MKL
    TPZSSpStructMatrix<STATE, TPZStructMatrixOR<STATE>> matsp(cmesh);
#else
    TPZFStructMatrix<STATE> matsp(cmesh);
#endif
    matsp.SetNumThreads(nThreads);
    matsp.SetApplyConstraintInternal(true);
    
    std::cout << "\n=====> Number of equations = " << cmesh->NEquations() << std::endl << std::endl;
#ifdef USE_MAIN
    TPZLinearAnalysis an(cmesh,RenumType::ENone);
#else
    TPZLinearAnalysis an(cmesh, RenumType::ENone);
#endif
    
    an.SetStructuralMatrix(matsp);
    TPZStepSolver<STATE> step;
    step.SetDirect(ELDLt);
    an.SetSolver(step);
    
    if(isTestKnownSol){
        TestKnownSol(an,1.,cmesh);
    }
    else{
        an.Run();
    }
}

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZCompMesh* cmesh) {
    
    // Assemble matrix and rhs
    an.Assemble();
    
    // Fill the connects related with displacement with solution cteSol
    TPZMultiphysicsCompMesh* mpcmesh = dynamic_cast<TPZMultiphysicsCompMesh*>(an.Mesh());
    if (mpcmesh) {
        int64_t nc = mpcmesh->NConnects();
        TPZFMatrix<STATE> &sol = mpcmesh->Solution();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = mpcmesh->ConnectVec()[ic];
            int64_t seqnum = c.SequenceNumber();
            if(seqnum < 0) continue;
            unsigned char lagrange = c.LagrangeMultiplier();
            STATE fill = 0.;
            if(lagrange == 1)
            {
                fill = cteSol;
            }
            if(c.Order() > 1) continue;
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof ; idf++) {
                int64_t index = mpcmesh->Block().Index(seqnum, idf);
                sol(index,0) = fill;
            }
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mpcmesh->MeshVector(), mpcmesh);
        
        // Get Matrix
        TPZMatrix<STATE>* mat = an.MatrixSolver<STATE>().Matrix().operator->();
        
        // Multiply matrix by known solution vector
        const int neq = mpcmesh->NEquations();
        TPZFMatrix<STATE> res(neq,1,0.);
        mat->Multiply(mpcmesh->Solution(), res);
        res = res - an.Rhs();
        std::ofstream out("problematicElsGlob.txt");
        an.PrintVectorByElement(out, res, 1.e-6);
        
        for(int i = 0 ; i < mpcmesh->MeshVector().size() ; i++){
            std::string name = "Mesh_" + std::to_string(i) + ".txt";
            std::ofstream outmesh(name);
            mpcmesh->MeshVector()[i]->Print(outmesh);
        }

    } else {
        int64_t nc = cmesh->NConnects();
        TPZFMatrix<STATE> &sol = cmesh->Solution();
        for (int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = cmesh->ConnectVec()[ic];
            int64_t seqnum = c.SequenceNumber();
            if(seqnum < 0) continue;
            unsigned char lagrange = c.LagrangeMultiplier();
            STATE fill = 0.;
            if(lagrange == 1)
            {
                fill = cteSol;
            }
            if(c.Order() > 1) continue;
            int ndof = c.NShape() * c.NState();
            for (int idf = 0; idf < ndof ; idf++) {
                int64_t index = cmesh->Block().Index(seqnum, idf);
                sol(index,0) = fill;
            }
        }
        
        // Get Matrix
        TPZMatrix<STATE>* mat = an.MatrixSolver<STATE>().Matrix().operator->();
        
        // Multiply matrix by known solution vector
        const int neq = cmesh->NEquations();
        TPZFMatrix<STATE> res(neq,1,0.);
        mat->Multiply(cmesh->Solution(), res);
        res = res - an.Rhs();
        std::ofstream out("problematicElsGlob.txt");
        an.PrintVectorByElement(out, res, 1.e-6);
        
        std::ofstream outmesh("Mesh.txt");
        cmesh->Print(outmesh);
    }
}

void PostProcessVTK(TPZCompMesh* cmesh, ProblemType probType) {
    TPZMultiphysicsCompMesh* mpcmesh = dynamic_cast<TPZMultiphysicsCompMesh*>(cmesh);
    if (mpcmesh) {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(mpcmesh->MeshVector(), mpcmesh);
    }
    const std::string plotfile = "PostProcess"; //sem o .vtk no final
    constexpr int vtkRes{0};
    
    TPZManVector<std::string,2> fields = {"Flux","Pressure"};
    if(probType == ProblemType::EElastic){
        fields[0] = "SigmaX";
        fields[1] = "Displacement";
    }
    if (mpcmesh) {
        auto vtk = TPZVTKGenerator(mpcmesh, fields, plotfile, vtkRes);
        vtk.Do();
    } else {
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
        vtk.Do();
    }
}

void CheckError(TPZCompMesh *cmesh, TPZVec<REAL> &error, ProblemType pType){
    
    cmesh->LoadReferences();
    cmesh->EvaluateError(false,error);
    const int dim = cmesh->Dimension();
    
    std::cout << "Error = " << error << std::endl;
    if(pType == ProblemType::EDarcy){
        if(error[1] >= 1.e-5) {
            std::cout << "Deu errado\n";
        }
        REQUIRE(error[1] < 1.e-6);
    }
    if(pType == ProblemType::EElastic) {
        if(error[0] >= 1.e-6) {
            std::cout << "Deu errado\n";
        }
        REQUIRE(error[0] < 1.e-6);
    }
    /*
     for (int i = 0 ; i < error.size() ; i++) {
     #ifndef USE_MAIN
     if( pType == ProblemType::EElastic && i == 6 ) {
     // In Elastic mat i == 6 means energy error of the exact solution
     if(dim == 2){
     REQUIRE(error[i] == Catch::Approx(1.));
     }
     else if(dim == 3) {
     REQUIRE(error[i] == Catch::Approx(M_SQRT2));
     }
     continue;
     }
     
     // if elastic we dont check error in displacement for HDivConstant spaces
     if(pType == ProblemType::EElastic && i == 3 && isHDivConst){
     continue;
     }
     REQUIRE(fabs(error[i]) < 1.e-10);
     #endif
     }
     std::cout << std::endl;
     */
}

void TestH1ApproxSpaceCreator(H1Family h1Fam, ProblemType probType, int pOrder, bool isRigidBodySpaces, MMeshType mType, int extrapOrder, bool isCondensed, HybridizationType hType, bool isRef) {

    gDarcy.fExact = TLaplaceExample1::EXpY;
    gDarcy.fExact = TLaplaceExample1::EY;
    gElast2d.fProblemType = TElasticity2DAnalytic::EHomogeneous;
    gElast3d.fProblemType = TElasticity3DAnalytic::EHomogeneous;
    
    TPZSimpleTimer totaltime;
    // ==========> Initial headers <==========
    // =======================================
    static int globcount = 0;
    std::cout << "\n------------------ Starting test " << globcount++ << " ------------------" << std::endl;
    std::cout <<
    "H1Fam = " << H1FamilyToChar(h1Fam) <<
    "\nProblemType = " << ProblemTypeToChar(probType) <<
    "\nMeshType = " << mType <<
    "\npOrder = " << pOrder <<
    "\nExtra POrder = " << extrapOrder <<
    "\nHybridization type = " << HybridizationTypeToChar(hType) <<
    "\nisRBSpaces = " << std::boolalpha << isRigidBodySpaces <<
    "\nisCondensed = " << std::boolalpha << isCondensed <<
    "\nisRef = " << std::boolalpha << isRef << std::endl;
    
    // ==========> Creating GeoMesh <==========
    // ========================================
    TPZGeoMesh *gmesh;
    if (mType == MMeshType::EQuadrilateral || mType == MMeshType::ETriangular){
        gmesh = Create2DGeoMesh(probType,mType);
    } else if (mType == MMeshType::EHexahedral || mType == MMeshType::ETetrahedral){
        gmesh = Create3DGeoMesh(probType,mType);
    } else {
        DebugStop();
    }
    
    if(isRef) Refinement(gmesh);
    
    std::ofstream out("GeoMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    
    // ==========> Creating Multiphysics mesh <==========
    // ==================================================
    TPZH1ApproxCreator* h1Creator = nullptr;
    h1Creator = new TPZH1ApproxCreator(gmesh);
    h1Creator->SetProbType(probType);
    h1Creator->IsRigidBodySpaces() = isRigidBodySpaces;
    h1Creator->SetDefaultOrder(pOrder);
    int dim = gmesh->Dimension();
    h1Creator->SetHybridType(hType);
    if(hType != HybridizationType::ENone) {
        h1Creator->SetExtraInternalOrder(dim+extrapOrder);
        h1Creator->SetHybridizeBoundary();
    } else {
        h1Creator->SetExtraInternalOrder(extrapOrder);

    }
    h1Creator->SetShouldCondense(isCondensed);
    InsertMaterials(*h1Creator);
    TPZCompMesh *cmesh = nullptr;
    if (h1Creator->HybridType() == HybridizationType::ENone) {
        cmesh = h1Creator->CreateClassicH1ApproximationSpace();
    } else {
        cmesh = h1Creator->CreateApproximationSpace();
    }
    //    std::ofstream outtxt("geomeshmodified.txt");
    //    gmesh->Print(outtxt);
    
    // ==========> Check number of equations for condensed problems <==========
    // ========================================================================
    if(isCondensed && hType != HybridizationType::ENone) {
        CheckNEqCondensedProb(cmesh,*h1Creator,mType);
    }
    
    
    // ==========> Solving problem <==========
    // =======================================
    const bool isTestKnownSol = false;
    SolveSystem(cmesh,isTestKnownSol);
    if(isTestKnownSol) return; // do not perform checks if just testing known sol
    
    // ==========> Post processing <==========
    // =======================================
    //#ifdef USE_MAIN
//    PostProcessVTK(cmesh,probType);
    //#endif
    if(0)
    {
        std::ofstream out("cmesh.txt");
        int64_t neq = cmesh->Solution().Rows();
        std::map<int64_t,std::array<REAL, 3>> connectco;
        TPZFMatrix<REAL> eqco(neq,3,0.);
        int64_t nel = cmesh->NElements();
        int dim = cmesh->Dimension();
        for(int64_t el =0; el<nel; el++) {
            TPZCompEl *cel = cmesh->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() != dim) continue;
            int ncorner = gel->NCornerNodes();
            for(int i=0; i<ncorner; i++) {
                int64_t seq = cel->Connect(i).SequenceNumber();
                for(int c=0; c<3; c++) connectco[seq][c] = gel->Node(i).Coord(c);
            }
        }
        eqco.Redim(connectco.size(),3);
        int count = 0;
        for(auto it : connectco) {
            for(int c=0; c<3; c++) eqco(count,c) = it.second[c];
            count++;
        }
        cmesh->Print(out);
        eqco.Print("co = ",out, EMathematicaInput);
    }
    
    // ==========> Unit test checks <==========
    // ========================================
    // Checks if the integral over the domain is a known value (most of the cases a constant value so just the volume of the domain)
    //    CheckIntegralOverDomain(cmesh,probType,hdivFam);
    // Checks if error with respect to exact solution is close to 0
    TPZMaterial *mat = cmesh->FindMaterial(EDomain);
    TPZMatError<STATE>* materror = dynamic_cast<TPZMatError<STATE>*>(mat);
    if (!materror) DebugStop();
    int nerror = materror->NEvalErrors();
    TPZManVector<REAL,5> error(nerror);
    CheckError(cmesh,error,probType);
    
    std::cout << "\n------------------ Test ended without crashing ------------------" << std::endl;
    std::cout << "==> Total time: " << totaltime.ReturnTimeDouble()/1000. << " seconds" << std::endl << std::endl;
}
