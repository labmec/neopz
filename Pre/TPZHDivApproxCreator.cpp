#include "TPZHDivApproxCreator.h"
#include "TPZMaterial.h"
#include "TPZBndCond.h"
#include "pzcmesh.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZNullMaterialCS.h"

#include <DarcyFlow/TPZMixedDarcyFlow.h>
#include "TPZBndCondT.h"
#include "TPZCompElDiscScaled.h"
#include "TPZCompMeshTools.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

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

    if (fHDivFam == HDivFamily::EHDivKernel){
        std::cout << "Pay attention to your material, the boundary conditions \n" <<
                     "need to be checked. Its implementation is not the same as \n" <<
                     "HDivStandard and HDivConstant spaces!!!! \n";
    }

    if (fExtraInternalPOrder < 0 || fExtraInternalPOrder > 2){
        std::cout << "Extra internal pOrder can only be 0, 1 or 2.\n";
        DebugStop();
    }
    
    if(fHybridType == HybridizationType::EStandardSquared){
        std::cout << "Standard Squared hybrization not implemented. Please implement if needed" << std::endl;
        DebugStop();
    }

}

TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateApproximationSpace(){

    CheckSetupConsistency();

    if (fHybridType != HybridizationType::ENone){
        ComputePeriferalMaterialIds();
        AddHybridizationGeoElements();
        std::ofstream out("GeoMeshHybrid.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh, out);
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
        meshvec[countMesh++] = CreateConstantSpace(lagLevelCounter++);
        meshvec[countMesh++] = CreateConstantSpace(lagLevelCounter++);
    }
    
//    if(fShouldCondense && fHybridType != HybridizationType::ENone){
//        ChangeLagLevel(meshvec[1],lagLevelCounter);
//    }
    
//    std::ofstream outp("pmeshAfterLagLevelChange.txt");
//    meshvec[1]->Print(outp);

    if (countMesh != fNumMeshes) DebugStop();

    TPZMultiphysicsCompMesh *cmeshmulti = CreateMultiphysicsSpace(meshvec);

    if (fShouldCondense){  
        if (isElastic && !fIsRBSpaces){ 
            // In this case, the third (corresponding to the laglevCounter - 1) is the rotation mesh, whose can be condensed.
            // Displacements should be kept uncondensed
            TPZCompMeshTools::CondenseElements(cmeshmulti,lagLevelCounter-2,false);
        } else {
            if(fHybridType == HybridizationType::ENone){
                TPZCompMeshTools::CondenseElements(cmeshmulti,lagLevelCounter-1,false);
            }
            else{
                GroupAndCondenseElements(cmeshmulti);
            }
        }
    }
    
    return cmeshmulti;
}

TPZCompMesh * TPZHDivApproxCreator::CreateHDivSpace(){
    
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    if (fProbType == ProblemType::EElastic) {
        cmesh->SetDefaultOrder(fDefaultPOrder+1);
    }
    else if (fProbType == ProblemType::EDarcy){
        cmesh->SetDefaultOrder(fDefaultPOrder);
    }
    
    cmesh->SetDimModel(dim);

    int nstate = 0;
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
        if (!bnd){
            if (mat->Dimension() != dim) DebugStop();
            nstate = mat->NStateVariables(); // here we assume that all materials have same nstatevars
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        } else {
            // if (bnd->Dimension() != dim-1) DebugStop();
            TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim-1,mat->NStateVariables());
            cmesh->InsertMaterialObject(nullmat);
        }
    }
    
    if(nstate < 1) DebugStop();
    
    if(fHybridType != HybridizationType::ENone){
        // Then we insert the hdiv wrap materials
        const int wrapmatid = fHybridizationData.fWrapMatId;
        TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(wrapmatid,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
    }
    
    //Creates computational elements       
    cmesh->ApproxSpace().SetHDivFamily(fHDivFam);
    cmesh->ApproxSpace().SetAllCreateFunctionsHDiv(dim);
    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared) {
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
        std::set<int> volmatids = GetVolumetricMatIds();
        cmesh->AutoBuild(volmatids);
        cmesh->SetDimModel(dim);
        cmesh->SetAllCreateFunctionsHDiv();
        cmesh->ApproxSpace().CreateDisconnectedElements(false);
        std::set<int> bcmatids = GetBCMatIds();
        bcmatids.insert(fHybridizationData.fWrapMatId);
        cmesh->LoadReferences(); // So the bcs can see the volumetric elements during their creation
        cmesh->AutoBuild(bcmatids);
        FixSideOrientHydridMesh(cmesh);
    }
    else{
        cmesh->AutoBuild();
    }
    
    cmesh->SetDimModel(dim);
    cmesh->InitializeBlock();
    
    if (fExtraInternalPOrder > 0){
        ChangeInternalOrder(cmesh,fDefaultPOrder+fExtraInternalPOrder);
    }
    
#ifdef PZDEBUG
    std::ofstream out("fluxmesh.txt");
    cmesh->Print(out);
//    PrintMeshElementsConnectInfo(cmesh);
#endif

    return cmesh;
}

void TPZHDivApproxCreator::FixSideOrientHydridMesh(TPZCompMesh* cmesh) {
    const int dim = cmesh->Dimension();
    int64_t nel = cmesh->NElements();
    std::set<int> bcmatids = GetBCMatIds();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        const int matid = gel->MaterialId();
        if(bcmatids.find(matid) != bcmatids.end()) continue;
//        if(gel->Dimension() != dim) continue;
        int firstside = gel->FirstSide(dim-1);
        for(int is = firstside; is< gel->NSides()-1; is++)
        {
            TPZGeoElSide gelside(gel,is);
            auto neigh = gelside.HasNeighbour(fHybridizationData.fWrapMatId);
            if(neigh){
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                intel->SetSideOrient(is, 1.);
            }
        }
    }
}

void TPZHDivApproxCreator::PrintMeshElementsConnectInfo(TPZCompMesh* cmesh) {
    std::cout << "\n===> Printing elements connect information" << std::endl;
    for (int i = 0; i < cmesh->NElements(); i++) {
        TPZCompEl *cel = cmesh->Element(i);
        if (!cel) continue;
        if (!cel->Reference()) continue;
        cel->LoadElementReference();
        int matid = cel->Reference()->MaterialId();
        auto nconnects = cel->NConnects();
        std::cout << "Element = " << i << ", dim = " << cel->Dimension() << ", mat = " << cel->Reference()->MaterialId() << ", nstate = " << cel->Material()->NStateVariables() << ", nconnects = " << nconnects << ": ";
        for (int j = 0; j < nconnects; j++)
        {
            std::cout << cel->ConnectIndex(j) << ", ";
        }
        std::cout << std::endl;
    }
}

TPZCompMesh * TPZHDivApproxCreator::CreateConstantSpace(const int lagLevel) {
    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    int nstate = 0;
    if (fHDivFam != HDivFamily::EHDivKernel){
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                if (mat->Dimension() != dim) DebugStop();
                nstate = mat->NStateVariables();
                TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
                cmesh->InsertMaterialObject(nullmat);
                if(fProbType == ProblemType::EElastic && fIsRBSpaces && lagLevel > 2){
                    if(dim == 2) nullmat->SetNStateVariables(3);
                    else if(dim == 3) nullmat->SetNStateVariables(6);
                }
            }
        }
    }
    
    if(nstate < 1) DebugStop();

    cmesh->SetDefaultOrder(0);
    cmesh->SetAllCreateFunctionsDiscontinuous();

    cmesh->AutoBuild();

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i];
        newnod.SetLagrangeMultiplier(lagLevel);
    }

    return cmesh;
}

TPZCompMesh * TPZHDivApproxCreator::CreateL2Space(const int pOrder, const int lagLevel){

    fGeoMesh->ResetReference();
    int dim = fGeoMesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(fGeoMesh);
    cmesh->SetDimModel(dim);

    int nstate = 0;
    if (fHDivFam != HDivFamily::EHDivKernel){
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                if (mat->Dimension() != dim) DebugStop();
                nstate = mat->NStateVariables();
                TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(mat->Id(),dim,mat->NStateVariables());
                cmesh->InsertMaterialObject(nullmat);
                if(fProbType == ProblemType::EElastic && fIsRBSpaces && lagLevel > 2){
                    if(dim == 2) nullmat->SetNStateVariables(3);
                    else if(dim == 3) nullmat->SetNStateVariables(6);
                }
            } 
        }
    }
    
    if(nstate < 1) DebugStop();

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
    
    if(fHybridType == HybridizationType::EStandard){
        // Then we insert the lagrange materials into this l2 space cmesh
        const int lagmatid = fHybridizationData.fLagrangeMatId;
        TPZNullMaterial<> *nullmat = new TPZNullMaterial<>(lagmatid,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmat);
        
        const int lagpord = pOrder - 1;
        if (lagpord > 0){
            cmesh->SetDefaultOrder(lagpord);
            cmesh->SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        } else {
            cmesh->SetDefaultOrder(0);
            cmesh->SetAllCreateFunctionsDiscontinuous();
        }
        
        // Now we only autobuild the langrange materials
        std::set<int> matidtobuild;
        matidtobuild.insert(lagmatid);
        cmesh->SetDimModel(dim-1);
        cmesh->AutoBuild(matidtobuild);
    }

    int ncon = cmesh->NConnects();
    for(int i=0; i<ncon; i++)
    {
        TPZConnect &newnod = cmesh->ConnectVec()[i]; 
        newnod.SetLagrangeMultiplier(lagLevel);
    }
    
//    PrintMeshElementsConnectInfo(cmesh);

    return cmesh;
}


TPZMultiphysicsCompMesh * TPZHDivApproxCreator::CreateMultiphysicsSpace(TPZManVector<TPZCompMesh *> meshvec){
    
    int dim = fGeoMesh->Dimension();
    auto cmesh = new TPZMultiphysicsCompMesh(fGeoMesh);
    cmesh->SetDefaultOrder(fDefaultPOrder);
    cmesh->SetDimModel(dim);
    
    int nstate = 0;
    for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
        TPZMaterial* mat = matpair.second;
        TPZBndCondT<STATE> *bnd = dynamic_cast<TPZBndCondT<STATE> *> (mat);
        if (!bnd){
            nstate = mat->NStateVariables();
            cmesh->InsertMaterialObject(mat);
        } else {
            cmesh->InsertMaterialObject(bnd);
        }
    }
    
    if(nstate < 1) DebugStop();
    
    if(fHybridType != HybridizationType::ENone){
        auto * nullmatWrap = new TPZNullMaterialCS<>(fHybridizationData.fWrapMatId,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmatWrap);

        auto * nullmatLag = new TPZNullMaterialCS<>(fHybridizationData.fLagrangeMatId,dim-1,nstate);
        cmesh->InsertMaterialObject(nullmatLag);
    }

    TPZManVector<int> active(fNumMeshes,1);
    cmesh->ApproxSpace().Style() = TPZCreateApproximationSpace::EMultiphysics;
    cmesh->BuildMultiphysicsSpace(active, meshvec);
    
    if(fHybridType != HybridizationType::ENone){
        InsertInterfaceMaterialObjects(cmesh);
        AddInterfaceComputationalElements(cmesh);
    }
    
    return cmesh;
}


TPZCompMesh * TPZHDivApproxCreator::CreateRotationSpace(const int pOrder, const int lagLevel){
    const int dim = fGeoMesh->Dimension();
    //Criando malha computacional:
    TPZCompMesh * cmesh = new TPZCompMesh(fGeoMesh);
    //Insere ordem polimonial de aproximação
    if (fHDivFam == HDivFamily::EHDivConstant){
        cmesh->SetDefaultOrder(0);
    } else if (fHDivFam == HDivFamily::EHDivStandard){
        cmesh->SetDefaultOrder(pOrder);
    }else {
        DebugStop();
    }
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
    
//    PrintMeshElementsConnectInfo(cmesh);
    
    return cmesh;
}

void TPZHDivApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *cmesh) {
    
    const int dim = cmesh->Dimension();
    
    int64_t nel = cmesh->NElements();
    cmesh->LoadReferences();
    
    TPZCompEl *cel;
    TPZGeoEl *gel;
    TPZGeoElSide *gelside;
    TPZGeoElSide neigh;
    TPZCompEl *celTarget;
    
    std::set<int64_t> targetMatIds;
    targetMatIds.insert(fHybridizationData.fWrapMatId);
    targetMatIds.insert(fHybridizationData.fInterfaceMatId);
    std::set<int> bcmatids = GetBCMatIds();
    for(auto it : bcmatids){
        targetMatIds.insert(it);
    }
    
//    std::ofstream ofs1("checkGelmesh.vtk");
//    TPZVTKGeoMesh::PrintGMeshVTK(fGeoMesh,ofs1);
    
    //    std::cout << "Groups of connects " << groupindex << std::endl;
    for (int64_t el = 0; el<nel; el++) {
        cel = cmesh->Element(el);
        if(!cel){
            continue;
        }
        if (cmesh->Element(el)->Dimension() != dim) {
            continue;
        }
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
        
        gel = cel->Reference();
        
        int firstSide;
        if (gel->Dimension() >= 0) {
            firstSide = gel->FirstSide(gel->Dimension() - 1);
        } else {
            DebugStop();
        }
        
        elgr->AddElement(cel);
        for (int iside = firstSide; iside < gel->NSides() - 1; iside++) {
            gelside = new TPZGeoElSide(gel, iside);
            neigh = *gelside;
            if(gelside->HasNeighbour(bcmatids)){
                neigh = gelside->Neighbour();
                if(bcmatids.find(neigh.Element()->MaterialId()) != bcmatids.end()){
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                }
            } else{
                neigh = gelside->Neighbour();
                std::vector<int> targetMatIds= {fHybridizationData.fWrapMatId,fHybridizationData.fInterfaceMatId};
                for (int ineigh = 0; ineigh < 2; ineigh++) {
                    if (neigh.Element()->MaterialId() != targetMatIds[ineigh]) {
                        DebugStop();
                    }
                    celTarget = neigh.Element()->Reference();
                    elgr->AddElement(celTarget);
                    neigh = neigh.Neighbour();
                }
            }
            gelside = NULL;
        }
    }
    cmesh->ComputeNodElCon();
    /*if(fSpaceType == EHybridizedMixed)
     {
     int64_t nconnects = cmesh->NConnects();
     for (int64_t ic = 0; ic<nconnects; ic++) {
     TPZConnect &c = cmesh->ConnectVec()[ic];
     if(c.LagrangeMultiplier() == fConfigHybridizedMixed.lagavg) c.IncrementElConnected();
     }
     }*/
    int neoNel = cmesh->NElements();
    for (int64_t el = nel; el < neoNel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (!elgr) {
            DebugStop();
        }
        TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
        cond->SetKeepMatrix(false);
    }
    
    cmesh->CleanUpUnconnectedNodes();
    
}

void TPZHDivApproxCreator::AddInterfaceComputationalElements(TPZMultiphysicsCompMesh *mphys) {
    int numEl = mphys->NElements();

    auto fGeoMesh = mphys->Reference();
    fGeoMesh->ResetReference();
    mphys->LoadReferences();

    for(int iel = 0; iel < numEl; iel++){
        TPZCompEl *cel = mphys->Element(iel);
        if(!cel || !cel->Reference()){
            DebugStop();
        }
        TPZGeoEl *gel = cel->Reference();
        int matid = gel->MaterialId();
        if(matid != fHybridizationData.fWrapMatId){
            continue;
        }
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighLag = gelside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        if(!neighLag){
            DebugStop();
        }
        TPZCompElSide celside = gelside.Reference();
        TPZCompElSide cneigh = neighLag.Reference();

        TPZGeoElSide ginterface = gelside.Neighbour();
        if(ginterface.Element()->MaterialId() != fHybridizationData.fInterfaceMatId){
            DebugStop();
        }
        
#ifdef PZDEBUG
//        std::cout << "===> Connecting elements with matids: " << celside.Element()->Reference()->MaterialId() << "\t" << cneigh.Element()->Reference()->MaterialId() << std::endl;
//        std::cout << "===> And indexes: " << celside.Element()->Reference()->Index() << "\t" << cneigh.Element()->Reference()->Index() << std::endl;
#endif
        
        TPZMultiphysicsInterfaceElement *interface = new TPZMultiphysicsInterfaceElement(*mphys,ginterface.Element(),celside,cneigh);
    }
}

void TPZHDivApproxCreator::ChangeLagLevel(TPZCompMesh* cmesh, const int newLagLevel) {
    if(cmesh->ApproxSpace().Style() != TPZCreateApproximationSpace::MApproximationStyle::EDiscontinuous &&
       cmesh->ApproxSpace().Style() != TPZCreateApproximationSpace::MApproximationStyle::EContinuous){
        DebugStop(); // Lagrange elements can only exist in L2 meshes
    }

    // Loop over all compels and pick the ones that have matid of lagrange multiplier
    for(auto cel : cmesh->ElementVec()){
        if(!cel) continue;
        const int matid = cel->Reference()->MaterialId();
        if(matid != fHybridizationData.fLagrangeMatId) continue;
        // loop over connects changing their laglevels
        const int ncon = cel->NConnects();
        for(int ic = 0 ; ic < ncon ; ic++){
            cel->Connect(ic).SetLagrangeMultiplier(newLagLevel);
        }
    }
}
