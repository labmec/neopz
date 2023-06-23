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

#include <pzlog.h>

// ----- Unit test includes -----
//#define USE_MAIN

#ifndef USE_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/generators/catch_generators_all.hpp>
#include <catch2/matchers/catch_matchers_all.hpp>
#endif

/// Creates a simple mesh used for testing
TPZGeoMesh *Create2DGeoMesh();

void InsertMaterials(TPZH1ApproxCreator &approxCreator);

void TestH1ApproxSpaceCreator(H1Family h1Fam, HybridizationType hybType, ProblemType probType, int pOrder, int plusOrder, bool IsRigidBodySpaces,bool shouldCondense);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, H1Family h1fam);

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann};


#ifndef USE_MAIN
TEST_CASE("Approx Space Creator", "[h1_space_creator_test]") {
    std::cout << "Testing H1 Approx Space Creator \n";

    H1Family sType = GENERATE(H1Family::EH1Standard);
    HybridizationType hybtype = GENERATE(HybridizationType::ENone);
    ProblemType pType = GENERATE(ProblemType::EDarcy);
    int pOrder = GENERATE(1,2);
    int plusOrder = GENERATE(0,1);
    bool isEnhanced = GENERATE(true,false);
    bool shouldCondense = GENERATE(true,false);

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TestH1ApproxSpaceCreator(sType, hybtype,pType,pOrder, plusOrder, isEnhanced, shouldCondense);
    std::cout << "Finish test H1 Approx Space Creator \n";
}
#else
int main(){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    HybridizationType hybtype = HybridizationType::ENone;
    TestH1ApproxSpaceCreator(H1Family::EH1Standard, hybtype, ProblemType::EDarcy,2, 1, true,true);
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

void InsertMaterials(TPZH1ApproxCreator &approxCreator){

    approxCreator.ProbType() = ProblemType::EDarcy;

    HybridizationType hybType = approxCreator.HybridType();
    TPZDarcyFlow *mat;
    if(hybType == HybridizationType::ENone){
        mat = new TPZDarcyFlow(EDomain, 2);
    }
    else {
       mat = new TPZHybridDarcyFlow(EDomain, 2);
    }
        mat->SetConstantPermeability(1.);
        // mat->SetForcingFunction()
        approxCreator.InsertMaterialObject(mat);

        //Boundary Conditions
        TPZFMatrix<STATE> val1(1, 1, 1.);
        TPZManVector<STATE> val2(1, 1.);

        //Dirichlet Boundary Conditions
        TPZBndCondT<STATE> *BCond1 = mat->CreateBC(mat, EBCDirichlet, 0, val1, val2);
        // BCond->SetForcingFunctionBC(exactSol,4);
        approxCreator.InsertMaterialObject(BCond1);

        val2[0] = 0.;
        TPZBndCondT<STATE> *BCond2 = mat->CreateBC(mat, EBCNeumann, 1, val1, val2);
        // BCond->SetForcingFunctionBC(exactSol,4);
        approxCreator.InsertMaterialObject(BCond2);
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
void TestH1ApproxSpaceCreator(H1Family h1Fam, HybridizationType hybtype ,ProblemType probType, int pOrder, int plusOrder, bool IsRigidBodySpaces, bool shouldCondense){

    TPZGeoMesh *gmesh = Create2DGeoMesh();

    TPZH1ApproxCreator h1Creator(gmesh);

    h1Creator.SetShouldCondense(shouldCondense);
    h1Creator.HybridType() = hybtype;
    h1Creator.ProbType() = probType;
    h1Creator.IsRigidBodySpaces() = IsRigidBodySpaces;
    h1Creator.SetDefaultOrder(pOrder);
    h1Creator.SetExtraInternalOrder(plusOrder);
    InsertMaterials(h1Creator);

    TPZCompMesh *cmesh;
    if(h1Creator.HybridType() == HybridizationType::ENone)
        cmesh = h1Creator.CreateClassicH1ApproximationSpace();
    else
        cmesh = h1Creator.CreateApproximationSpace();

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
    TPZLinearAnalysis an(cmesh,true);
    an.Run();

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
        REQUIRE(fabs(vecint[i]) < 1.e-10);
#endif
    }
    std::cout << std::endl;

    TPZVec<STATE> vecintp = cmesh->Integrate(fields[1], matids);
    std::cout << "\n--------------- Integral of Pressure --------------" <<  std::endl;
    std::cout << "Number of components = " << vecintp.size() <<  std::endl;
    for (int i = 0; i < vecintp.size(); i++)
    {
        std::cout << "Integral(" << i << ") = "  << vecintp[i] << std::endl;
#ifndef USE_MAIN
        REQUIRE(fabs(vecintp[i]) == Catch::Approx( 4.0 ));
#endif
    }
    std::cout << std::endl;
}

