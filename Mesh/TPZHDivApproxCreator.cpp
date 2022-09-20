#include "TPZHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzcmesh.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"

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

    std::cout << "Hello, Jeferson!\n";
    TPZManVector<TPZCompMesh*,7> meshvec(2);

    meshvec[0] = CreateHDivSpace();
    meshvec[1] = CreateL2Space(1);
    
    return nullptr;
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
            if (mat->Dimension() != dim-1) DebugStop();
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

TPZCompMesh * TPZHDivApproxCreator::CreateL2Space(const int lagLevel){
    
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    for (TPZMaterial* mat:fMaterialVec)
    {
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd){
            if (mat->Dimension() != dim) DebugStop();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        } 
    }
   
    //Creates computational elements
    switch (fHDivFam)
    {
    case HDivFamily::EHDivStandard:
        cmesh->SetDefaultOrder(fDefaultPOrder);
        if (fDefaultPOrder > 0){
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        break;

    case HDivFamily::EHDivConstant:
        cmesh->SetDefaultOrder(0);
        cmesh->SetAllCreateFunctionsDiscontinuous();
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


TPZCompMesh * TPZHDivApproxCreator::CreateConstantSpace(){
    return nullptr;
}

