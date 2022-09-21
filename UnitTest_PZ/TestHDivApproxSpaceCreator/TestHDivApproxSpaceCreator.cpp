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

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZHDivApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"

#include <pzlog.h>

// ----- Unit test includes -----
// #define USE_MAIN

#ifndef USE_MAIN
#include<catch2/catch.hpp>
#endif

/// Creates a simple mesh used for testing
TPZGeoMesh *Create2DGeoMesh();

void InsertMaterials(TPZHDivApproxCreator &approxCreator);

void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder, bool isEnhancedSpaces);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam);


enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann};


#ifndef USE_MAIN
TEST_CASE("Approx Space Creator", "[hdiv_space_creator_test]") {
    std::cout << "Testing HDiv Approx Space Creator \n";
    
    HDivFamily sType = GENERATE(HDivFamily::EHDivKernel,HDivFamily::EHDivConstant,HDivFamily::EHDivStandard);
    ProblemType pType = GENERATE(ProblemType::EDarcy);
    int pOrder = GENERATE(1,2,3,4);
    bool isEnhanced = GENERATE(true,false);
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TestHdivApproxSpaceCreator(sType,pType,pOrder,isEnhanced);
    std::cout << "Finish test HDiv Approx Space Creator \n";
}
#else
int main(){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TestHdivApproxSpaceCreator(HDivFamily::EHDivStandard,ProblemType::EDarcy, 1, true);
    
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

void InsertMaterials(TPZHDivApproxCreator &approxCreator){

    approxCreator.ProbType() = ProblemType::EDarcy;

    TPZMixedDarcyFlow *mat = new TPZMixedDarcyFlow(EDomain,2);
    mat->SetConstantPermeability(1.);
    // mat->SetForcingFunction()
    approxCreator.InsertMaterialObject(mat);

    //Boundary Conditions
    TPZFMatrix<STATE> val1(1,1,1.);
    TPZManVector<STATE> val2(1,1.);

    //Dirichlet Boundary Conditions
    TPZBndCondT<STATE> * BCond1 = mat->CreateBC(mat, EBCDirichlet, 0, val1, val2);
    // BCond->SetForcingFunctionBC(exactSol,4);
    approxCreator.InsertMaterialObject(BCond1);

    val2[0] = 0.;
    TPZBndCondT<STATE> * BCond2 = mat->CreateBC(mat, EBCNeumann, 1, val1, val2);
    // BCond->SetForcingFunctionBC(exactSol,4);
    approxCreator.InsertMaterialObject(BCond2);
    
}


void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder, bool isEnhancedSpaces){

    if (hdivFam == HDivFamily::EHDivKernel && isEnhancedSpaces) {
        //Hdiv kernel currently does not support enhanced spaces 
        return;
    }

    TPZGeoMesh *gmesh = Create2DGeoMesh();

    TPZHDivApproxCreator hdivCreator(gmesh);

    hdivCreator.HdivFamily() = hdivFam;
    hdivCreator.ProbType() = probType;
    hdivCreator.EnhancedSpaces() = isEnhancedSpaces;
    hdivCreator.SetDefaultOrder(pOrder);
    InsertMaterials(hdivCreator);

    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace(); 

    std::string txt = "cmesh.txt";
    std::ofstream myfile(txt);
    cmesh->Print(myfile);

    std::cout << "Number of equations = " << cmesh->NEquations() << std::endl;
    TPZLinearAnalysis an(cmesh,true);
    an.Run();

    CheckIntegralOverDomain(cmesh,probType,hdivFam);

    {
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(cmesh->MeshVector(), cmesh);
        const std::string plotfile = "myfile";//sem o .vtk no final
        constexpr int vtkRes{0};    

        TPZVec<std::string> fields = {
        "Flux",
        "Pressure"};
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        vtk.Do();
    }


}

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam){

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

    if (hdivfam != HDivFamily::EHDivKernel){
        TPZVec<STATE> vecintp = cmesh->Integrate(fields[1], matids);
        std::cout << "\n--------------- Integral of Pressure --------------" <<  std::endl;
        std::cout << "Number of components = " << vecintp.size() <<  std::endl;
        for (int i = 0; i < vecintp.size(); i++)
        {
            std::cout << "Integral(" << i << ") = "  << vecintp[i] << std::endl;
#ifndef USE_MAIN
            REQUIRE(fabs(vecintp[i]) == Approx( 4.0 ));
#endif
        }
        std::cout << std::endl;
    }
}
