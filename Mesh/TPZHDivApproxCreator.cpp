#include "TPZHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzcmesh.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZBndCondT.h"

TPZHDivApproxCreator::TPZHDivApproxCreator(TPZGeoMesh *gmesh) : fGeoMesh(gmesh)
{ 
}

TPZHDivApproxCreator::~TPZHDivApproxCreator()
{
}


/**Insert a material object in the datastructure*/
int TPZHDivApproxCreator::InsertMaterialObject(TPZMaterial * mat) {
	if(!mat) DebugStop();
	fMaterialVec.Push(mat); 
	return fMaterialVec.size();
}

/*Due to the TPZMaterialRefactor, the type TPZBndCond
  no longer inherits from TPZMaterial. But every possible instance
  will be a TPZMaterial, so this dynamic_cast is not expected
  to fail
*/
int TPZHDivApproxCreator::InsertMaterialObject(TPZBndCond * mat) {
    return InsertMaterialObject(dynamic_cast<TPZMaterial*>(mat));
}

TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateApproximationSpace(){

    if (fProbType == ProblemType::ENone) {
        std::cout << "You have to set a proper problem type!\n";
        DebugStop();
    }

    if (!fGeoMesh) {
        std::cout << "You have to set a GeoMesh!\n";
        DebugStop();
    }

    if (fMaterialVec.size() == 0){
        std::cout << "You have to set a Material!\n";
        DebugStop();
    }
 
    if (fIsEnhancedSpaces && fHDivFam == HDivFamily::EHDivKernel){
        std::cout << "Are you sure about this?\n";
        DebugStop();
    }

    bool isElastic = fProbType == ProblemType::EElastic;
    bool isDarcy = fProbType == ProblemType::EDarcy;
    if (isElastic){
        fNumMeshes = 3;
    } else if (isDarcy) {
        fNumMeshes = 2;
    } else {
        DebugStop();
    }
    if (fIsEnhancedSpaces) fNumMeshes += 2;

    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    int countMesh = 0;
    meshvec[countMesh++] = CreateHDivSpace();
    int lagLevelCounter = 1;
    meshvec[countMesh++] = CreateL2Space(fDefaultPOrder,lagLevelCounter++);
    if (isElastic){
        meshvec[countMesh++] = CreateRotationSpace();
    }
    if (fIsEnhancedSpaces){
        meshvec[countMesh++] = CreateL2Space(0,lagLevelCounter++);
        meshvec[countMesh++] = CreateL2Space(0,lagLevelCounter++);
    }

    if (countMesh != fNumMeshes) DebugStop();

    TPZMultiphysicsCompMesh *cmeshmulti = CreateMultiphysicsSpace(meshvec);
    
    return cmeshmulti;
}

TPZCompMesh * TPZHDivApproxCreator::CreateHDivSpace(){
    
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(dim);

    for (TPZMaterial* mat:fMaterialVec)
    {
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd){
            if (mat->Dimension() != dim) DebugStop();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        } else {
            // if (bnd->Dimension() != dim-1) DebugStop();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim-1,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        }
    }
   
    //Creates computational elements       
    cmesh->ApproxSpace().SetHDivFamily(fHDivFam);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    cmesh->AutoBuild();

    return cmesh;
}

TPZCompMesh * TPZHDivApproxCreator::CreateL2Space(const int pOrder, const int lagLevel){

    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    if (fHDivFam != HDivFamily::EHDivKernel){
        for (TPZMaterial* mat:fMaterialVec)
        {
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                if (mat->Dimension() != dim) DebugStop();
                TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
                cmesh->InsertMaterialObject(nullmat);
            } 
        }
    }

    //Creates computational elements
    switch (fHDivFam)
    {
    case HDivFamily::EHDivStandard:
        if (pOrder > 0){
            cmesh->SetDefaultOrder(pOrder);
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetDefaultOrder(0);
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        break;

    case HDivFamily::EHDivConstant:
        cmesh->SetDefaultOrder(0);
        cmesh->SetAllCreateFunctionsDiscontinuous();
        break;

    case HDivFamily::EHDivKernel:
        
        break;

    default:
        DebugStop();
        break;
    }
        
    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    return cmesh;
}


TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateMultiphysicsSpace(TPZManVector<TPZCompMesh *> meshvec){
    
    int dim = fGeoMesh->Dimension();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(dim);
    
    for (TPZMaterial* mat:fMaterialVec)
    {
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *> (mat);
        if (!bnd){
            cmesh->InsertMaterialObject(mat);
        } else {
            cmesh->InsertMaterialObject(bnd);
        }
    }

    TPZManVector<int> active(fNumMeshes,1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvec);
    
    return cmesh;
}


TPZCompMesh * TPZHDivApproxCreator::CreateRotationSpace(){

   

}