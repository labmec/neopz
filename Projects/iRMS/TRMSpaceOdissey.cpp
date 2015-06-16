//
//  TRMSpaceOdissey.cpp
//  PZ
//
//  Created by omar duran on 5/05/2015.
//
//

#include "TRMSpaceOdissey.h"
#include "TRMFlowConstants.h"
#include "TPZMatLaplacian.h"
#include "pzbndcond.h"



/// Create a H1 approximation mesh
void TRMSpaceOdissey::CreateH1Mesh()
{
    if(!fGeoMesh)
    {
        DebugStop();
    }
    fH1Mesh = new TPZCompMesh(fGeoMesh);
    fH1Mesh->SetDimModel(3);
    
    TPZMatLaplacian *material = new TPZMatLaplacian(_ReservMatId,3);
    fH1Mesh->InsertMaterialObject(material);
    
    TPZFNMatrix<1> val1(1,1,0.),val2(1,1,20);
    TPZBndCond *inflow = new TPZBndCond(material,_WellToeMatId,0,val1,val2);
    val2(0,0) = 10.;
    TPZBndCond *outflow = new TPZBndCond(material,_WellHeelMatId,0,val1,val2);
    
    fH1Mesh->InsertMaterialObject(inflow);
    fH1Mesh->InsertMaterialObject(outflow);
    
    TPZCreateApproximationSpace space;
    space.SetAllCreateFunctionsContinuous();    
    fH1Mesh->ApproxSpace() = space;
    
    fH1Mesh->AutoBuild();
    
}
