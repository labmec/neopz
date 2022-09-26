#include "TPZHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzcmesh.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZBndCondT.h"
#include "TPZCompElDiscScaled.h"

TPZHDivApproxCreator::TPZHDivApproxCreator(TPZGeoMesh *gmesh) : TPZApproxCreator(gmesh)
{ 
}

TPZHDivApproxCreator::~TPZHDivApproxCreator()
{
}

void TPZHDivApproxCreator::CheckSetupConsistency() {
    
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
 
    if (fIsRBSpaces && fHDivFam == HDivFamily::EHDivKernel){
        std::cout << "Are you sure about this?\n";
        DebugStop();
    }
}

TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateApproximationSpace(){

    CheckSetupConsistency();

    bool isElastic = fProbType == ProblemType::EElastic;
    bool isDarcy = fProbType == ProblemType::EDarcy;
    if (isElastic){
        fNumMeshes = 3;
    } else if (isDarcy) {
        fNumMeshes = 2;
    } else {
        DebugStop();
    }
    if (fIsRBSpaces) fNumMeshes += 2;

    TPZManVector<TPZCompMesh*,7> meshvec(fNumMeshes);
    int countMesh = 0;
    meshvec[countMesh++] = CreateHDivSpace();
    int lagLevelCounter = 1;
    meshvec[countMesh++] = CreateL2Space(fDefaultPOrder,lagLevelCounter++);
    if (isElastic){
        meshvec[countMesh++] = CreateRotationSpace(fDefaultPOrder,lagLevelCounter++);
    }
    if (fIsRBSpaces){
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

    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
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
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
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
    
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
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


TPZCompMesh * TPZHDivApproxCreator::CreateRotationSpace(const int pOrder, const int lagLevel){
    const int dim = fGeoMesh->Dimension();
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(pOrder); //Insere ordem polimonial de aproximação
    cmesh->SetDimModel(dim); //Insere dimensão do modelo

    cmesh->SetAllCreateFunctionsDiscontinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    // This set will hold all the domain materials to later create the TPZCompElDiscScaled
    std::set<int> materialids;
    
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd){
            if (mat->Dimension() != dim) DebugStop();
            const int matiddomain = mat->Id();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(matiddomain,dim,mat->NStateVariables());
            if(dim == 3) nullmat->SetNStateVariables(3);
            else if(dim == 2) nullmat->SetNStateVariables(1);
            else DebugStop();
            cmesh->InsertMaterialObject(nullmat);
            materialids.insert(matiddomain);
        }
    }
    
    
    fGeoMesh->ResetReference();
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if (!gel)continue;
        if(gel->HasSubElement()) continue;
        int matid = gel->MaterialId();
        if (materialids.find(matid) == materialids.end()) {
            continue;
        }
        new TPZCompElDiscScaled(*cmesh, gel);
        gel->ResetReference();
    }

    //cmesh->LoadReferences();
    //    cmesh->ApproxSpace().CreateDisconnectedElements(false);
    //    cmesh->AutoBuild();


    int ncon = cmesh->NConnects();
    for (int i = 0; i < ncon; i++) {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    int64_t nelem = cmesh->NElements();
    for (int64_t el = 0; el < nelem; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZCompElDiscScaled *disc = dynamic_cast<TPZCompElDiscScaled *> (cel);
        if (!disc) {
            continue;
        }
        const REAL elementdim = cel->Reference()->ElementRadius();
        disc->SetTotalOrderShape();
        disc->SetFalseUseQsiEta();
        disc->SetConstC(elementdim);
        disc->SetScale(1./elementdim);
    }
    cmesh->InitializeBlock();
    return cmesh;
}

void TPZHDivApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mcmesh) {
    DebugStop(); // Implement me if needed
}
