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
#include "TPZEnumApproxFamily.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include <Elasticity/TPZMixedElasticityND.h>
#include "TPZHDivApproxCreator.h"
#include "TPZVTKGenerator.h"
#include "pzbuildmultiphysicsmesh.h"

#include <pzlog.h>

// ----- Unit test includes -----
//#define USE_MAIN

#ifndef USE_MAIN
#include<catch2/catch.hpp>
#endif

using namespace std;

/// Creates a simple mesh used for testing
TPZGeoMesh *Create2DGeoMesh();

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype);

void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder, bool isRigidBodySpaces);

void CheckIntegralOverDomain(TPZCompMesh *cmesh, ProblemType probType, HDivFamily hdivfam);

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann};


#ifndef USE_MAIN
TEST_CASE("Approx Space Creator", "[hdiv_space_creator_test]") {
    std::cout << "Testing HDiv Approx Space Creator \n";
    
    HDivFamily sType = GENERATE(HDivFamily::EHDivKernel,HDivFamily::EHDivConstant,HDivFamily::EHDivStandard);
    ProblemType pType = GENERATE(ProblemType::EDarcy);
    int pOrder = GENERATE(1,2,3,4);
    bool isRBSpaces = GENERATE(true,false);
    
#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    TestHdivApproxSpaceCreator(sType,pType,pOrder,isRBSpaces);
    std::cout << "Finish test HDiv Approx Space Creator \n";
}
#else
int main(){

#ifdef PZ_LOG
    TPZLogger::InitializePZLOG();
#endif
    HDivFamily sType = HDivFamily::EHDivStandard;
    ProblemType pType = ProblemType::EElastic;
    const int pord = 1;
    const bool isRBSpaces = false;
    TestHdivApproxSpaceCreator(sType,pType,pord,isRBSpaces);
    
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

void InsertMaterials(TPZHDivApproxCreator &approxCreator, ProblemType &ptype){
    
    if (!approxCreator.GeoMesh()) {
        cout << "\nError! Please set the geomesh before inserting materials" << endl;
        DebugStop();
    }
    const int dim = approxCreator.GeoMesh()->Dimension();

    approxCreator.ProbType() = ptype;

    TPZMaterial *mat = nullptr;
    TPZMixedDarcyFlow* matdarcy = nullptr;
    TPZMixedElasticityND* matelas = nullptr;
    if (ptype == ProblemType::EDarcy) {
        matdarcy = new TPZMixedDarcyFlow(EDomain,2);
        matdarcy->SetConstantPermeability(1.);
        mat = matdarcy;
    }
    else if (ptype == ProblemType::EElastic){
        REAL E = 20.59, nu = 0., fx = 0., fy = 0.;
        const int plain = 1.; //* @param plainstress = 1 \f$ indicates use of plainstress
        matelas = new TPZMixedElasticityND(EDomain, E, nu, fx, fy, plain, dim);
        mat = matelas;
    }

    approxCreator.InsertMaterialObject(mat);

    // ========> Boundary Conditions
    // -----------------------------
    
    TPZBndCondT<STATE> *BCond1 = nullptr, *BCond2 = nullptr;
    const int dirType = 0, neuType = 1;
    if (ptype == ProblemType::EDarcy) {
        TPZFMatrix<STATE> val1(1,1,1.);
        TPZManVector<STATE> val2(1,1.);
        BCond1 = matdarcy->CreateBC(matdarcy, EBCDirichlet, dirType, val1, val2);
        val2[0] = 0.;
        BCond2 = matdarcy->CreateBC(matdarcy, EBCNeumann, neuType, val1, val2);
    }
    else if (ptype == ProblemType::EElastic){
        TPZFMatrix<STATE> val1(2,2,0.);
        TPZManVector<STATE> val2(2,1.);

        BCond1 = matelas->CreateBC(matelas, EBCDirichlet, dirType, val1, val2);
        val2[0] = 0.;
        BCond2 = matelas->CreateBC(matelas, EBCNeumann, neuType, val1, val2);
    }
    else{
        DebugStop(); // yet not supported material
    }
    
    approxCreator.InsertMaterialObject(BCond1);
    approxCreator.InsertMaterialObject(BCond2);
}


void TestHdivApproxSpaceCreator(HDivFamily hdivFam, ProblemType probType, int pOrder, bool isRigidBodySpaces){

    cout << "\n------------------ Starting test ------------------" << endl;
    cout << "HdivFam = " << static_cast<std::underlying_type<HDivFamily>::type>(hdivFam) <<
    "\nProblemType = " << static_cast<std::underlying_type<ProblemType>::type>(probType) <<
    "\npOrder = " << pOrder << "\nisRBSpaces = " << isRigidBodySpaces << endl << endl;
    
    if (hdivFam == HDivFamily::EHDivKernel && isRigidBodySpaces) {
        //Hdiv kernel currently does not support enhanced spaces 
        return;
    }

    TPZGeoMesh *gmesh = Create2DGeoMesh();

    TPZHDivApproxCreator hdivCreator(gmesh);

    hdivCreator.HdivFamily() = hdivFam;
    hdivCreator.ProbType() = probType;
    hdivCreator.IsRigidBodySpaces() = isRigidBodySpaces;
    hdivCreator.SetDefaultOrder(pOrder);
    InsertMaterials(hdivCreator,probType);

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

        TPZManVector<std::string,2> fields = {"Flux","Pressure"};
        if(probType == ProblemType::EElastic){
            fields[0] = "SigmaX";
            fields[1] = "Displacement";
        }
        auto vtk = TPZVTKGenerator(cmesh, fields, plotfile, vtkRes);

        vtk.Do();
    }

    cout << "\n------------------ Test ended with serious errors ------------------" << endl << endl;
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
        std::cout << "\n--------------- Integral of State Var --------------" <<  std::endl;
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
