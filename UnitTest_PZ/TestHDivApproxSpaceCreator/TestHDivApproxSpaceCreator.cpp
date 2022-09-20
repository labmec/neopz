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

#include <pzlog.h>

// ----- Unit test includes -----
#include <catch2/catch.hpp>


TPZGeoMesh *Create2DGeoMesh();

void InsertMaterials(TPZHDivApproxCreator &approxCreator);

// ----- Run tests with or without main -----
//#define RUNWITHMAIN

enum MaterialIds {EDomain,EBCDirichlet,EBCNeumann};

int main(){

    TPZGeoMesh *gmesh = Create2DGeoMesh();

    TPZHDivApproxCreator hdivCreator(gmesh);

    InsertMaterials(hdivCreator);


    TPZMultiphysicsCompMesh *cmesh = hdivCreator.CreateApproximationSpace(); 

    std::cout << "Hello, Nathan!\n";
    return 0;
}



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