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
#ifdef PZ_USING_MKL
#include "TPZSSpStructMatrix.h"
#endif
#include "TPZEnumApproxFamily.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZMixedElasticityND.h>
#include "TPZHDivApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzstepsolver.h"
#include "pzcheckgeom.h"
#include "pzcheckmesh.h"
#include "TPZAnalyticSolution.h"
#include "TPZMHMHDivApproxCreator.h"

#include <pzlog.h>
#include "TPZSimpleTimer.h"

// ----- Unit test includes -----
//#define USE_MAIN

#ifndef USE_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#endif

using namespace std;

/// Creates a simple mesh used for testing
TPZGeoMesh *Create2DGeoMesh(ProblemType& pType, MMeshType &mType);
TPZGeoMesh *Create3DGeoMesh(ProblemType& pType, MMeshType &mType);

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype);

void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                bool isCondensed, HybridizationType hType, bool isRef, bool isMHM);

void TestHdivApproxSpaceCreator2(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                 bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                 bool isCondensed, HybridizationType hType, bool isRef, bool isMHM);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam);

void CheckError(TPZMultiphysicsCompMesh *cmesh, TPZVec<REAL> &error, ProblemType pType);

void Refinement(TPZGeoMesh *gmesh);

void SolveSystem(TPZMultiphysicsCompMesh* cmesh, const bool isTestKnownSol);

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZMultiphysicsCompMesh* mpcmesh);

void PostProcessVTK(TPZMultiphysicsCompMesh* cmesh, ProblemType probType);

void CheckNEqCondensedProb(TPZMultiphysicsCompMesh* mpcmesh,
                           TPZHDivApproxCreator& hdivcreator,
                           MMeshType& elType);

void SetConstantPrimalSol(TPZMultiphysicsCompMesh* cmesh, ProblemType probType, bool isRigidBodySpaces);

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann,EBCDisplacementLeft,EBCDisplacementRight};

constexpr const char* HDivFamilyToChar(HDivFamily hdivfam) {
    switch (hdivfam){
        case HDivFamily::EHDivStandard: return "EHDivStandard";
        case HDivFamily::EHDivConstant: return "EHDivConstant";
        case HDivFamily::EHDivKernel: return "EHDivKernel";
        case HDivFamily::EHDivOptimized: return "EHDivOptimized";
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

#define TESTALL
#ifndef USE_MAIN
TEST_CASE("HDiv Approx Space Creator", "[hdiv_linear_solution_representation]")
{
    {
#ifdef TESTALL
             HDivFamily sType = GENERATE(HDivFamily::EHDivConstant, HDivFamily::EHDivStandard, HDivFamily::EHDivOptimized);
#else
        HDivFamily sType = GENERATE(HDivFamily::EHDivOptimized);
#endif
        SECTION("HDiv family type: " + std::string(HDivFamilyToChar(sType)))
        {
#ifdef TESTALL
            ProblemType pType = GENERATE(ProblemType::EDarcy, ProblemType::EElastic);
#else
                    ProblemType pType = GENERATE(ProblemType::EDarcy);
#endif
            SECTION("Problem type: " + std::string(ProblemTypeToChar(pType)))
            {
#ifdef TESTALL
                MMeshType mType = GENERATE(MMeshType::EQuadrilateral, MMeshType::ETriangular, MMeshType::ETetrahedral, MMeshType::EHexahedral);
#else
                MMeshType mType = GENERATE(MMeshType::EHexahedral);
#endif
                SECTION("Mesh type: " + std::string(MeshTypeToChar(mType)))
                {
#ifdef TESTALL
                    int pOrder = GENERATE(1,2);
#else
                    int pOrder = GENERATE(1);
#endif
                    SECTION("pOrder=" + std::to_string(pOrder))
                    {
#ifdef TESTALL
                        int extraporder = GENERATE(0,1);
#else
                        int extraporder = GENERATE(1);
#endif
                        SECTION("extraporder=" + std::to_string(extraporder))
                        {
#ifdef TESTALL
                            HybridizationType hType = GENERATE(HybridizationType::ENone, HybridizationType::EStandard, HybridizationType::ESemi);
#else
                            HybridizationType hType = GENERATE(HybridizationType::ENone);
#endif
                            SECTION("Hybridization type: " + std::string(HybridizationTypeToChar(hType)))
                            {
                                // bool isRBSpaces = GENERATE(true, false);
                                bool isRBSpaces = GENERATE(false,true);
                                SECTION("isRigidBodySpaces=" + std::to_string(isRBSpaces))
                                {
#ifdef TESTALL
                                    bool isCondensed = GENERATE(false, true);
#else
                                    bool isCondensed = GENERATE(false, true);
#endif
                                    //                            bool isCondensed = false;
                                    SECTION("isCondensed=" + std::to_string(isCondensed))
                                    {
#ifdef TESTALL
                                        bool isRef = GENERATE(false, true);
#else
                                        bool isRef = GENERATE(false);
#endif
                                        SECTION("isRef=" + std::to_string(isRef))
                                        {
                                            // skip extra p order for EHDivOptimized
                                            bool skipPorder = sType == HDivFamily::EHDivOptimized && extraporder > 0;
                                            skipPorder = false;
                                            // skipping 3d elements with hdivconstant and no hybridization
                                            bool skipadapt = isRef && (sType == HDivFamily::EHDivConstant || sType == HDivFamily::EHDivOptimized) && (mType == MMeshType::EPrismatic || mType == MMeshType::EHexahedral || mType == MMeshType::ETetrahedral) && hType == HybridizationType::ENone;
                                            skipadapt = false;
                                            bool isMHM = false;
                                            bool skipsemi1 = (isRef && hType == HybridizationType::ESemi);
                                            bool skipsemi2 = sType != HDivFamily::EHDivConstant && sType != HDivFamily::EHDivOptimized && hType == HybridizationType::ESemi;
                                            bool skipcondensed = (isCondensed && sType == HDivFamily::EHDivConstant);
                                            bool shouldnotskip = !skipPorder && !skipadapt && !skipsemi1 && !skipcondensed && !skipsemi2;
                                            if (shouldnotskip)
                                            {
                                                TestHdivApproxSpaceCreator(sType, pType, pOrder, isRBSpaces, mType, extraporder, isCondensed, hType, isRef, isMHM);
                                            }
                                            else {
                                                cout << "----------------- Skipping ---------------------\n" <<
                                                "HdivFam = " << HDivFamilyToChar(sType) <<
                                                "\nProblemType = " << ProblemTypeToChar(pType) <<
                                                "\nMeshType = " << mType <<
                                                "\npOrder = " << pOrder <<
                                                "\nExtra POrder = " << extraporder <<
                                                "\nHybridization type = " << HybridizationTypeToChar(hType) <<
                                                "\nisRBSpaces = " << std::boolalpha << isRBSpaces <<
                                                "\nisCondensed = " << std::boolalpha << isCondensed <<
                                                "\nisRef = " << std::boolalpha << isRef << endl;
                                                std::cout << std::boolalpha << "\nskip extra p order " << skipPorder << "\nskip adaptivity " << skipadapt <<
                                                "\nskip semi 1 (semi and refined) " << skipsemi1 << "\nskip condensed " << skipcondensed << "\nskip semi 2 (semi and hdiv family) " << skipsemi2 << std::endl
                                                << std::endl;
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
    }
}
/*
 TEST_CASE("Rigid Body", "[hdiv_space_creator_residual_test]")
 {
 HDivFamily sType = GENERATE(HDivFamily::EHDivConstant, HDivFamily::EHDivStandard, HDivFamily::EHDivOptimized);
 //    HDivFamily sType = GENERATE(HDivFamily::EHDivConstant);
 SECTION("HDiv family type: " + std::string(HDivFamilyToChar(sType)))
 {
 ProblemType pType = GENERATE(ProblemType::EDarcy, ProblemType::EElastic);
 //        ProblemType pType = GENERATE(ProblemType::EElastic);
 SECTION("Problem type: " + std::string(ProblemTypeToChar(pType)))
 {
 int pOrder = GENERATE(2,1);
 SECTION("pOrder=" + std::to_string(pOrder))
 {
 bool isRBSpaces = GENERATE(true, false);
 //                bool isRBSpaces = GENERATE(true);
 SECTION("isRigidBodySpaces=" + std::to_string(isRBSpaces))
 {
 MMeshType mType = GENERATE(MMeshType::EQuadrilateral, MMeshType::ETriangular, MMeshType::EHexahedral, MMeshType::ETetrahedral);
 //                    MMeshType mType = GENERATE(MMeshType::ETetrahedral);
 SECTION("Mesh type: " + std::string(MeshTypeToChar(mType)))
 {
 int extraporder = GENERATE(1,0);
 SECTION("extraporder=" + std::to_string(extraporder))
 {
 bool isCondensed = GENERATE(false, true);
 //                            bool isCondensed = GENERATE(true);
 SECTION("isCondensed=" + std::to_string(isCondensed))
 {
 HybridizationType hType = GENERATE(HybridizationType::ENone, HybridizationType::EStandard, HybridizationType::ESemi);
 //                                HybridizationType hType = GENERATE(HybridizationType::EStandard);
 SECTION("Hybridization type: " + std::string(HybridizationTypeToChar(hType)))
 {
 bool isRef = GENERATE(true, false);
 //                                    bool isRef = GENERATE(true);
 SECTION("isRef=" + std::to_string(isRef))
 {
 bool isMHM = false;
 if (!(isRef && hType == HybridizationType::ESemi) && !(isCondensed && sType == HDivFamily::EHDivConstant) &&
 !(sType != HDivFamily::EHDivConstant && sType != HDivFamily::EHDivOptimized && hType == HybridizationType::ESemi))
 {
 TestHdivApproxSpaceCreator2(sType, pType, pOrder, isRBSpaces, mType, extraporder, isCondensed, hType, isRef, isMHM);
 } else {
 std::cout << "Skipping test\n";
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
 }
 */
#else
int main(){
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    //    HDivFamily sType = HDivFamily::EHDivStandard;
    //    HDivFamily sType = HDivFamily::EHDivKernel;
    HDivFamily sType = HDivFamily::EHDivConstant;
    
    ProblemType pType = ProblemType::EElastic;
    // ProblemType pType = ProblemType::EDarcy;
    
    const int pord = 1;
    const bool isRBSpaces = true;
    
    MMeshType mType = MMeshType::EQuadrilateral;
    // MMeshType mType = MMeshType::ETriangular;
    //    MMeshType mType = MMeshType::EHexahedral;
    // MMeshType mType = MMeshType::ETetrahedral;
    
    int extraporder = 0;
    bool isCondensed = true;
    // bool isCondensed = false;
    //    HybridizationType hType = HybridizationType::EStandard;
    // HybridizationType hType = HybridizationType::ESemi;
    HybridizationType hType = HybridizationType::ENone;
    
    // this will create a mesh with hanging nodes
    bool isRef = true;
    // bool isRef = false;
    
    bool isMHM = false;
    
    TestHdivApproxSpaceCreator(sType,pType,pord,isRBSpaces,mType,extraporder,isCondensed,hType,isRef,isMHM);
    
    return 0;
}
#endif


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

auto exactSol = [](const TPZVec<REAL> &loc,
                   TPZVec<STATE>&u,
                   TPZFMatrix<STATE>&gradU){
    const auto &x=loc[0];
    const auto &y=loc[1];
    const auto &z=loc[2];
    
    // REAL a1 = 1./4;
    // REAL alpha = M_PI/2;
    // u[0] = x*a1*cos(x*alpha)*cosh(y*alpha) + y*a1*sin(x*alpha)*sinh(y*alpha) + x*x - y*y;
    // gradU(0,0) = -a1*(cosh(alpha*y)*(cos(alpha*x) - alpha*x*sin(alpha*x)) + alpha*y*cos(alpha*x)*sinh(alpha*y));
    // gradU(1,0) = -a1*(alpha*y*cosh(alpha*y)*sin(alpha*x) + (alpha*x*cos(alpha*x) + sin(alpha*x))*sinh(alpha*y));
    
    u[0] = x*x*x*y - y*y*y*x;
    gradU(0,0) = (3.*x*x*y - y*y*y);
    gradU(1,0) = (x*x*x - 3.*y*y*x);
    
};

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype){
    
    if (!approxCreator.GeoMesh()) {
        cout << "\nError! Please set the geomesh before inserting materials" << endl;
        DebugStop();
    }
    const int dim = approxCreator.GeoMesh()->Dimension();
    
    approxCreator.SetProbType(ptype);
    
    TPZMaterial *mat = nullptr;
    TPZMixedDarcyFlow* matdarcy = nullptr;
    TPZMixedElasticityND* matelas = nullptr;
    if (ptype == ProblemType::EDarcy) {
        matdarcy = new TPZMixedDarcyFlow(EDomain,dim);
        matdarcy->SetConstantPermeability(1.);
        matdarcy->SetExactSol(gDarcy.ExactSolution(), 1);
        gDarcy.fTensorPerm.Redim(3, 3);
        gDarcy.fInvPerm.Redim(3, 3);
        gDarcy.fTensorPerm.Identity();
        gDarcy.fInvPerm.Identity();
        mat = matdarcy;
    }
    else if (ptype == ProblemType::EElastic){
        REAL E = 1., nu = 0., fx = 0., fy = 0.;
        const int plain = 0.; //* @param plainstress = 1 \f$ indicates use of plainstress
        matelas = new TPZMixedElasticityND(EDomain, E, nu, fx, fy, plain, dim);
        matelas->SetExactSol(exactSolElastic,1);
        if(dim == 2) {
            matelas->SetExactSol(gElast2d.ExactSolution(), 1);
            gElast2d.gE = 1.;
            gElast2d.gPoisson = 0.;
        } else if (dim == 3) {
            matelas->SetExactSol(gElast3d.ExactSolution(), 1);
            gElast3d.fE = 1.;
            gElast3d.fPoisson = 0.;
        }
        mat = matelas;
    }
    
    approxCreator.InsertMaterialObject(mat);
    
    // ========> Boundary Conditions
    // -----------------------------
    
    TPZBndCondT<STATE> *BCond1 = nullptr, *BCond2 = nullptr, *BCond3 = nullptr, *BCond4 = nullptr;
    const int dirType = 0, neuType = 1;
    
    if (ptype == ProblemType::EDarcy) {
        
        // matdarcy->SetExactSol(exactSol,4);
        
        TPZFMatrix<STATE> val1(1,1,0.);
        TPZManVector<STATE> val2(1,1.);
        BCond1 = matdarcy->CreateBC(matdarcy, EBCDirichlet, dirType, val1, val2);
        BCond1->SetForcingFunctionBC(gDarcy.ExactSolution(),1);
        val2[0] = 0.;
        BCond2 = matdarcy->CreateBC(matdarcy, EBCNeumann, dirType, val1, val2);
        BCond2->SetForcingFunctionBC(gDarcy.ExactSolution(),1);
    }
    else if (ptype == ProblemType::EElastic){
        TPZFMatrix<STATE> val1(dim,dim,0.);
        TPZManVector<STATE> val2(dim,0.);
        
        val2[0] = 1.;
        val2[1] = 1.;
        if(dim == 3) val2[2] = 1.;
        BCond1 = matelas->CreateBC(matelas, EBCNeumann, dirType, val1, val2);
        BCond2 = matelas->CreateBC(matelas, EBCDisplacementRight, dirType, val1, val2);
        BCond3 = matelas->CreateBC(matelas, EBCDisplacementLeft, dirType, val1, val2);
        BCond4 = matelas->CreateBC(matelas, EBCDirichlet, dirType, val1, val2);
        if(dim == 2) {
            BCond1->SetForcingFunctionBC(gElast2d.ExactSolution(), 1);
            BCond2->SetForcingFunctionBC(gElast2d.ExactSolution(), 1);
            BCond3->SetForcingFunctionBC(gElast2d.ExactSolution(), 1);
            BCond4->SetForcingFunctionBC(gElast2d.ExactSolution(), 1);
        } else if (dim == 3) {
            BCond1->SetForcingFunctionBC(gElast3d.ExactSolution(), 1);
            BCond2->SetForcingFunctionBC(gElast3d.ExactSolution(), 1);
            BCond3->SetForcingFunctionBC(gElast3d.ExactSolution(), 1);
            BCond4->SetForcingFunctionBC(gElast3d.ExactSolution(), 1);
        }
    }
    else{
        DebugStop(); // yet not supported material
    }
    
    if(BCond1) approxCreator.InsertMaterialObject(BCond1);
    if(BCond2) approxCreator.InsertMaterialObject(BCond2);
    if(BCond3) approxCreator.InsertMaterialObject(BCond3);
    if(BCond4) approxCreator.InsertMaterialObject(BCond4);
}


void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                bool isCondensed, HybridizationType hType, bool isRef, bool isMHM){
    
    gDarcy.fExact = TLaplaceExample1::EXpY;
    gDarcy.fExact = TLaplaceExample1::EY;
    gElast2d.fProblemType = TElasticity2DAnalytic::EHomogeneous;
    gElast3d.fProblemType = TElasticity3DAnalytic::EHomogeneous;
    
    TPZSimpleTimer totaltime;
    // ==========> Initial headers <==========
    // =======================================
    static int globcount = 0;
    cout << "\n------------------ Starting test " << globcount++ << " ------------------" << endl;
    cout <<
    "HdivFam = " << HDivFamilyToChar(hdivFam) <<
    "\nProblemType = " << ProblemTypeToChar(probType) <<
    "\nMeshType = " << mType <<
    "\npOrder = " << pOrder <<
    "\nExtra POrder = " << extrapOrder <<
    "\nHybridization type = " << HybridizationTypeToChar(hType) <<
    "\nisRBSpaces = " << std::boolalpha << isRigidBodySpaces <<
    "\nisCondensed = " << std::boolalpha << isCondensed <<
    "\nisRef = " << std::boolalpha << isRef << endl;
    
    // TODO: WARNING!!!! Things to be fixed and for now we are skipping
    if(isCondensed && hdivFam == HDivFamily::EHDivConstant){
        std::cout << "\n\t======> WARNING! Condensing hdivconst spaces leads to singular K00. This happens because of the rotation space that has linear and constant functions."
        " One option is to separate the linear and constant functions is two different spaces. Please implement before using." << std::endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if(isRef && hType == HybridizationType::ESemi){
        cout << "\n\t======> This test is not working because the global matrix has wrong size!!\n" << endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if (hdivFam == HDivFamily::EHDivKernel && isRigidBodySpaces) {
        std::cout << "ERROR! Hdiv kernel currently does not support rigid body spaces \n";
        DebugStop();
        return;
    }
    if (hdivFam != HDivFamily::EHDivConstant && hdivFam != HDivFamily::EHDivOptimized && hType == HybridizationType::ESemi) {
        std::cout << "ERROR! The only HDiv space with available Semi hybridization is HDivConstant \n";
        DebugStop();
        return;
    }
    
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
    TPZHDivApproxCreator* hdivCreator = nullptr;
    if(isMHM){
        hdivCreator = new TPZMHMHDivApproxCreator(gmesh);
    }
    else{
        hdivCreator = new TPZHDivApproxCreator(gmesh);
    }
    hdivCreator->HdivFamily() = hdivFam;
    hdivCreator->SetProbType(probType);
    hdivCreator->IsRigidBodySpaces() = isRigidBodySpaces;
    hdivCreator->SetDefaultOrder(pOrder);
    hdivCreator->SetExtraInternalOrder(extrapOrder);
    hdivCreator->SetShouldCondense(isCondensed);
    hdivCreator->SetHybridType(hType);
    //    hdivCreator->SetHybridizeBoundary();
    InsertMaterials(*hdivCreator,probType);
    TPZMultiphysicsCompMesh *cmesh = hdivCreator->CreateApproximationSpace();
    //    std::ofstream outtxt("geomeshmodified.txt");
    //    gmesh->Print(outtxt);
    
    // ==========> Check number of equations for condensed problems <==========
    // ========================================================================
    if(isCondensed && hType != HybridizationType::ENone) {
        CheckNEqCondensedProb(cmesh,*hdivCreator,mType);
    }
    
    
    // ==========> Solving problem <==========
    // =======================================
    const bool isTestKnownSol = false;
    SolveSystem(cmesh,isTestKnownSol);
    if(isTestKnownSol) return; // do not perform checks if just testing known sol
    
    // ==========> Post processing <==========
    // =======================================
    //#ifdef USE_MAIN
    PostProcessVTK(cmesh,probType);
    //#endif
#ifdef PZDEBUG
    //    hdivCreator->PrintMeshElementsConnectInfo(cmesh);
    hdivCreator->PrintAllMeshes(cmesh);
#endif
    
    // ==========> Unit test checks <==========
    // ========================================
    // Checks if the integral over the domain is a known value (most of the cases a constant value so just the volume of the domain)
    //    CheckIntegralOverDomain(cmesh,probType,hdivFam);
    // Checks if error with respect to exact solution is close to 0
    TPZMaterial *mat = cmesh->FindMaterial(EDomain);
    TPZMatErrorCombinedSpaces<STATE>* materror = dynamic_cast<TPZMatErrorCombinedSpaces<STATE>*>(mat);
    if (!materror) DebugStop();
    int nerror = materror->NEvalErrors();
    TPZManVector<REAL,5> error(nerror);
    CheckError(cmesh,error,probType);
    
    cout << "\n------------------ Test ended without crashing ------------------" << endl;
    std::cout << "==> Total time: " << totaltime.ReturnTimeDouble()/1000. << " seconds" << std::endl << endl;
}

void TestHdivApproxSpaceCreator2(HDivFamily hdivFam, ProblemType probType, int pOrder,
                                 bool isRigidBodySpaces, MMeshType mType, int extrapOrder,
                                 bool isCondensed, HybridizationType hType, bool isRef, bool isMHM){
    
    gDarcy.fExact = TLaplaceExample1::EConst;
    gElast2d.fProblemType = TElasticity2DAnalytic::EDispxy;
    gElast3d.fProblemType = TElasticity3DAnalytic::EDispxyz;
    TPZSimpleTimer totaltime;
    // ==========> Initial headers <==========
    // =======================================
    static int globcount = 0;
    cout << "\n------------------ Starting test " << globcount++ << " ------------------" << endl;
    cout << "HdivFam = " << HDivFamilyToChar(hdivFam) <<
    "\nProblemType = " << ProblemTypeToChar(probType) <<
    "\npOrder = " << pOrder << "\nisRBSpaces = " << std::boolalpha << isRigidBodySpaces <<
    "\nMeshType = " << mType << "\nExtra POrder = " << extrapOrder <<
    "\nisCondensed = " << std::boolalpha << isCondensed <<
    "\nHybridization type = " << HybridizationTypeToChar(hType) <<
    "\nisRef = " << std::boolalpha << isRef << endl << endl;
    
    // TODO: WARNING!!!! Things to be fixed and for now we are skipping
    if(isCondensed && hdivFam == HDivFamily::EHDivConstant){
        std::cout << "\n\t======> WARNING! Condensing hdivconst spaces leads to singular K00. This happens because of the rotation space that has linear and constant functions."
        " One option is to separate the linear and constant functions is two different spaces. Please implement before using." << std::endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    if(isRef && hType == HybridizationType::ESemi){
        cout << "\n\t======> This test is not working because the global matrix has wrong size!!\n" << endl;
        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
        DebugStop();
        return;
    }
    //    if(hdivFam == HDivFamily::EHDivConstant && probType == ProblemType::EElastic && mType == MMeshType::ETetrahedral){
    //        cout << "\n\t======> We don't know why this configuration does not work!!\n" << endl;
    //        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
    //        DebugStop();
    //        return;
    //    }
    //    if(hdivFam == HDivFamily::EHDivConstant && extrapOrder > 0){
    //        cout << "\n\t======> WARNING! SKIPPING TEST!!\n" << endl;
    //        DebugStop();
    //        return;
    //    }
    if (hdivFam == HDivFamily::EHDivKernel && isRigidBodySpaces) {
        std::cout << "ERROR! Hdiv kernel currently does not support rigid body spaces \n";
        DebugStop();
        return;
    }
    if (hdivFam != HDivFamily::EHDivConstant && hdivFam != HDivFamily::EHDivOptimized && hType == HybridizationType::ESemi) {
        std::cout << "ERROR! The only HDiv space with available Semi hybridization is HDivConstant \n";
        DebugStop();
        return;
    }
    
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
    if(0) {
        std::ofstream out("GeoMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out);
    }
    // ==========> Creating Multiphysics mesh <==========
    // ==================================================
    TPZHDivApproxCreator* hdivCreator = nullptr;
    if(isMHM){
        hdivCreator = new TPZMHMHDivApproxCreator(gmesh);
    }
    else{
        hdivCreator = new TPZHDivApproxCreator(gmesh);
    }
    hdivCreator->HdivFamily() = hdivFam;
    hdivCreator->SetProbType(probType);
    hdivCreator->IsRigidBodySpaces() = isRigidBodySpaces;
    hdivCreator->SetDefaultOrder(pOrder);
    hdivCreator->SetExtraInternalOrder(extrapOrder);
    hdivCreator->SetShouldCondense(isCondensed);
    hdivCreator->SetHybridType(hType);
    //    hdivCreator->SetHybridizeBoundary();
    InsertMaterials(*hdivCreator,probType);
    TPZMultiphysicsCompMesh *cmesh = hdivCreator->CreateApproximationSpace();
    //    std::ofstream outtxt("geomeshmodified.txt");
    //    gmesh->Print(outtxt);
    
    // ==========> Check number of equations for condensed problems <==========
    // ========================================================================
    if(isCondensed && hType != HybridizationType::ENone) {
        CheckNEqCondensedProb(cmesh,*hdivCreator,mType);
    }
    SetConstantPrimalSol(cmesh,probType,isRigidBodySpaces);
#ifdef PZDEBUG
    //    hdivCreator->PrintMeshElementsConnectInfo(cmesh);
    //    hdivCreator->PrintAllMeshes(cmesh);
#endif
    
    TPZLinearAnalysis an(cmesh,RenumType::ENone);
    an.Assemble();
    
    // Transfering multiphysics solution to the analysis solution
    int cmesh_neq = cmesh->NEquations();
    TPZFMatrix<STATE> &cmesh_sol = cmesh->Solution();
    TPZFMatrix<STATE> &sol = an.Solution();
    for (int i = 0; i < cmesh_neq; i++)
    {
        sol.PutVal(i, 0, cmesh_sol.GetVal(i, 0));
    }
    
    TPZFMatrix<STATE> rhs = an.Rhs();
    auto mat = an.MatrixSolver<STATE>().Matrix();
    TPZFMatrix<STATE> res(rhs);
    mat->MultAdd(sol, rhs, res, -1., 1.);
    REAL normres = Norm(res);
    std::cout << "Norm res " << normres << std::endl;
    if(normres >= 1.e-10) {
        {
            std::ofstream out("ResByElement.txt");
            an.PrintVectorByElement(out,res,1.e-6);
            sol.Print("sol = ",out, EMathematicaInput);
            res.Print("res = ",out, EMathematicaInput);
        }
        std::cout << "Vai dar merda\n";
        
    }
    REQUIRE(normres < 1.e-10);
    cout << "\n------------------ Test ended without crashing ------------------" << endl;
    std::cout << "==> Total time: " << totaltime.ReturnTimeDouble()/1000. << " seconds" << std::endl << endl;
}

void SolveSystem(TPZMultiphysicsCompMesh* cmesh, const bool isTestKnownSol) {
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

void PostProcessVTK(TPZMultiphysicsCompMesh* cmesh, ProblemType probType) {
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
    const std::string plotfile = "PostProcess"; //sem o .vtk no final
    constexpr int vtkRes{0};
    
    TPZManVector<std::string,2> fields = {"Flux","Pressure"};
    if(probType == ProblemType::EElastic){
        fields[0] = "SigmaX";
        fields[1] = "Displacement";
    }
    auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);
    
    vtk.Do();
    
}

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam){
    
    int dim = cmesh -> Dimension();
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
        if (probType == ProblemType::EDarcy){
            REQUIRE(fabs(vecint[i]) < 1.e-10);
        }
#endif
    }
    std::cout << std::endl;
#ifndef USE_MAIN
    if (probType == ProblemType::EElastic){
        if (dim == 2) REQUIRE(fabs(vecint[0]) == Catch::Approx( 2.0 ));
        if (dim == 3) REQUIRE(fabs(vecint[0]) == Catch::Approx( 4.0 ));
        for (int i = 1; i < vecint.size(); i++) REQUIRE(fabs(vecint[i]) < 1.e-10);
    }
#endif
    if (hdivfam != HDivFamily::EHDivKernel){
        TPZVec<STATE> vecintp = cmesh->Integrate(fields[1], matids);
        std::cout << "\n--------------- Integral of State Var --------------" <<  std::endl;
        std::cout << "Number of components = " << vecintp.size() <<  std::endl;
        for (int i = 0; i < vecintp.size(); i++)
        {
            std::cout << "Integral(" << i << ") = "  << vecintp[i] << std::endl;
#ifndef USE_MAIN
            if (probType == ProblemType::EDarcy){
                if (dim == 2) REQUIRE(fabs(vecintp[i]) == Catch::Approx( 4.0 ));
                if (dim == 3) REQUIRE(fabs(vecintp[i]) == Catch::Approx( 8.0 ));
            }
#endif
        }
        std::cout << std::endl;
#ifndef USE_MAIN
        if (probType == ProblemType::EElastic){
            if (dim == 2) REQUIRE(fabs(vecintp[0]) == Catch::Approx( 2.0 ));
            if (dim == 3) REQUIRE(fabs(vecintp[0]) == Catch::Approx( 4.0 ));
            REQUIRE(fabs(vecintp[1]) < 1.e-10);
            REQUIRE(fabs(vecintp[2]) < 1.e-10);
        }
#endif
    }
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
    // }
    
    
    
    //    children[0]->Divide(children);
}

void CheckError(TPZMultiphysicsCompMesh *cmesh, TPZVec<REAL> &error, ProblemType pType){
    
    cmesh->LoadReferences();
    cmesh->EvaluateError(false,error);
    const int dim = cmesh->Dimension();
    
    const bool isHDivConst = cmesh->MeshVector()[0]->ApproxSpace().HDivFam() == HDivFamily::EHDivConstant;
    std::cout << "Error = " << error << std::endl;
    if(pType == ProblemType::EDarcy){
        if(error[1] >= 1.e-6) {
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

void TestKnownSol(TPZLinearAnalysis& an, const REAL cteSol, TPZMultiphysicsCompMesh* mpcmesh) {
    
    // Assemble matrix and rhs
    an.Assemble();
    
    // Fill the connects related with displacement with solution cteSol
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
        std::string name = "Mesh_" + to_string(i) + ".txt";
        std::ofstream outmesh(name);
        mpcmesh->MeshVector()[i]->Print(outmesh);
    }
}

void CheckNEqCondensedProb(TPZMultiphysicsCompMesh* mpcmesh,
                           TPZHDivApproxCreator& hdivcreator,
                           MMeshType& elType) {
    
    std::cout << "\n--------------- Checking number of equations --------------" <<  std::endl;
    TPZCompMesh* pmesh = mpcmesh->MeshVector()[1];
    const int lagmatid = hdivcreator.HybridData().fLagrangeMatId;
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
    if (hdivcreator.ProbType()==ProblemType::EDarcy){
        nstate = 1;
    } else if (hdivcreator.ProbType()==ProblemType::EElastic) {
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
    int extrap = hdivcreator.GetExtraInternalOrder();
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

#include "TPZCompElDisc.h"

void SetConstantPrimalSol(TPZMultiphysicsCompMesh *cmesh, ProblemType probType, bool isRigidBodySpaces)
{
    
    cmesh->LoadReferences();
    TPZGeoMesh *gmesh = cmesh->Reference();
    const int dim = gmesh->Dimension();
    
    TPZFMatrix<STATE> &mp_cmesh_sol = cmesh->Solution();
    TPZCompMesh *L2_cmesh = cmesh->MeshVector()[1]; // pressure for Darcy and displacement for elasticity
    TPZFMatrix<STATE> &L2_cmesh_sol = L2_cmesh->Solution();
    {
        int count = 0;
        int64_t nc = L2_cmesh->NConnects();
        for(int64_t ic = 0; ic<nc; ic++) {
            TPZConnect &c = L2_cmesh->ConnectVec()[ic];
            int ndof = c.NShape()*c.NState();
            int order = c.Order();
            if(order > 1) {
                count += ndof;
                continue;
            }
            for(int i=0; i<ndof; i++) {
                L2_cmesh_sol(count+i,0) = 1.;
            }
            count += ndof;
        }
    }
    cmesh->LoadSolutionFromMeshes();
    
    TPZManVector<STATE,3> val(3,1.);
    if (isRigidBodySpaces) // in this case we have extra avg pressure (or displacement) connects that we need to fill as well
    {
        int meshpos = (probType == ProblemType::EDarcy) ? 3 : 4; // avg pressure mesh for Darcy and avg displacement mesh for elasticity position in mesh vec
        TPZCompMesh *avg_cmesh = cmesh->MeshVector()[meshpos];
        TPZFMatrix<STATE> &avg_cmesh_sol = avg_cmesh->Solution();
        for (TPZCompEl *cel : avg_cmesh->ElementVec())
        { // Average pressure elements
            TPZGeoEl *gel = cel->Reference();
            if (gel->Dimension() != dim)
                continue;
            TPZMultiphysicsElement *mp_el = dynamic_cast<TPZMultiphysicsElement *>(gel->Reference());
#ifdef PZDEBUG
            if (!mp_el)
            {
                DebugStop();
            }
#endif
#ifdef PZDEBUG
            if (cel->NConnects() != 1)
                DebugStop();
#endif
            int cindex = mp_el->NConnects() - 1; // avg pressure connect is the last one
            TPZConnect &cloc = cel->Connect(0);
            TPZConnect &c = mp_el->Connect(cindex);
            int64_t seqloc = cloc.SequenceNumber(); //
            int64_t seq = c.SequenceNumber();
            int64_t firstEqLoc = avg_cmesh->Block().Position(seqloc);
            int64_t firstEq = cmesh->Block().Position(seq);
            int blockSize = cmesh->Block().Size(seq);
            int blockSizeloc = avg_cmesh->Block().Size(seqloc);
#ifdef PZDEBUG
            if (blockSize != blockSizeloc)
            {
                DebugStop();
            }
#endif
            if(blockSizeloc > dim) blockSizeloc = dim;
            if(blockSize > dim) blockSize = dim;
            for (int64_t eqloc = firstEqLoc; eqloc < firstEqLoc + blockSizeloc; eqloc++)
            { // atomic pressure solution
                avg_cmesh_sol.PutVal(eqloc, 0, val[0]);
            }
            for (int64_t eq = firstEq; eq < firstEq + blockSize; eq++)
            { // multiphysics pressure solution
                mp_cmesh_sol.PutVal(eq, 0, val[0]);
            }
        }
    }
}
