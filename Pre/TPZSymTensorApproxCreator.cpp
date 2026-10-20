#include "TPZSymTensorApproxCreator.h"
#include "Elasticity/TPZMixedSymElasticityND.h"
#include "Elasticity/TPZInterfaceSymTensor.h"
#include "TPZRefPattern.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "TPZNullMaterial.h"
#include "TPZNullMaterialCS.h"
#include "TPZLagrangeMultiplierCS.h"
#include "TPZBndCond.h"
#include "pzgeoelbc.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include <sstream>
#include <stdexcept>

TPZSymTensorApproxCreator::TPZSymTensorApproxCreator(TPZGeoMesh *gmesh)
    : TPZApproxCreator(gmesh) {
    fProbType = ProblemType::EElastic;
    fNumMeshes = 2;
    fNState = 2;
    SetMaterialIds(fIds);
}

void TPZSymTensorApproxCreator::SetProbType(ProblemType prob) {
    if(prob != ProblemType::EElastic)
        throw std::invalid_argument("Symmetric tensor approximation requires elasticity");
    fProbType = prob;
}

void TPZSymTensorApproxCreator::SetHybridType(HybridizationType hybrid) {
    if(hybrid != HybridizationType::ENone)
        throw std::invalid_argument("Use SetLowOrder to select symmetric tensor hybridization");
    fHybridType = hybrid;
}

void TPZSymTensorApproxCreator::CheckSetupConsistency() {
    if(fCreated) throw std::logic_error("Symmetric tensor mesh has already been created");
    if(!fGeoMesh || fGeoMesh->Dimension() != 2)
        throw std::invalid_argument("Symmetric tensor approximation requires a 2D geometry");
    if(fDefaultPOrder < 1 || fLowOrder < -1 || fLowOrder > fDefaultPOrder)
        throw std::invalid_argument("Require p >= 1 and -1 <= low order <= p");
    if(fExtraInternalPOrder != 0 || fIsRBSpaces)
        throw std::invalid_argument("Extra internal order and rigid body spaces are not supported");
    auto it = fMaterialVec.find(fIds.fDomain);
    if(it == fMaterialVec.end() ||
       !dynamic_cast<TPZMixedSymElasticityND *>(it->second) || it->second->Dimension() != 2)
        throw std::invalid_argument("Insert a 2D TPZMixedSymElasticityND for the domain id");
    const std::set<int> auxiliary = {fIds.fMacroBoundary, fHybridizationData.fWrapMatId,
        fIds.fInterfacePositive, fIds.fInterfaceNegative, fIds.fDisplacement,
        fIds.fInterfacePositiveLow, fIds.fInterfaceNegativeLow,
        fIds.fDisplacementLow, fIds.fTractionLow};
    if(auxiliary.size() != 9 || auxiliary.count(fIds.fDomain))
        throw std::invalid_argument("Symmetric tensor material ids must be distinct");
    for(auto &entry : fMaterialVec) {
        if(auxiliary.count(entry.first))
            throw std::invalid_argument("Auxiliary ids collide with registered materials");
        if(entry.first != fIds.fDomain && !dynamic_cast<TPZBndCond *>(entry.second))
            throw std::invalid_argument("Only one domain material is supported");
    }
    bool hasDomain = false;
    for(auto gel : fGeoMesh->ElementVec()) {
        if(!gel || gel->HasSubElement()) continue;
        const int id = gel->MaterialId();
        if(auxiliary.count(id) && !(fGeometryIsPrepared && id == fIds.fMacroBoundary))
            throw std::invalid_argument("Geometry already contains auxiliary elements");
        if(gel->Dimension() != 2) continue;
        if(id != fIds.fDomain || (gel->Type() != ETriangle &&
           (!fGeometryIsPrepared && gel->Type() != EQuadrilateral)) ||
           (fGeometryIsPrepared && gel->Type() != ETriangle))
            throw std::invalid_argument("Expected domain triangles/quadrilaterals or prepared triangle fans");
        if(!fGeometryIsPrepared && !gel->IsRefPatternEl())
            throw std::invalid_argument("Johnson-Mercier subdivision requires refinement-pattern geometric elements");
        hasDomain = true;
    }
    if(!hasDomain) throw std::invalid_argument("Geometry contains no domain elements");
}

void TPZSymTensorApproxCreator::PrepareGeometry() {
    if(!fGeometryIsPrepared) {
        const auto nel = fGeoMesh->NElements();
        for(int64_t i = 0; i < nel; i++) {
            auto gel = fGeoMesh->Element(i);
            if(!gel || gel->HasSubElement() || gel->MaterialId() != fIds.fDomain) continue;
            for(int side = gel->FirstSide(1); side < gel->NSides()-1; side++) {
                TPZGeoElSide gs(gel,side);
                if(!gs.HasNeighbour(fIds.fMacroBoundary)) TPZGeoElBC(gs,fIds.fMacroBoundary);
            }
        }
        CreateJohnsonMercier(fGeoMesh);
    }
    // Each triangle must have exactly one center outside the macro boundary.
    for(auto gel : fGeoMesh->ElementVec()) {
        if(!gel || gel->HasSubElement() || gel->Dimension() != 2) continue;
        int centers = 0;
        for(int side = 0; side < gel->NCornerNodes(); side++)
            if(!TPZGeoElSide(gel,side).HasNeighbour(fIds.fMacroBoundary)) centers++;
        if(centers != 1) throw std::invalid_argument("Each triangle requires exactly one macroelement center");
    }
    AddWrapElements(fGeoMesh);
    if(fLowOrder >= 0) AddWrapElementsTripleHybrid(fGeoMesh);
}

TPZMultiphysicsCompMesh *TPZSymTensorApproxCreator::CreateApproximationSpace() {
    CheckSetupConsistency();
    PrepareGeometry();
    auto mesh = GenerateMultiphysicsMesh(fGeoMesh,fDefaultPOrder,fLowOrder);
    if(fShouldCondense) {
        GroupAndCondenseElements(mesh);
    } else {
        // Center trace multipliers are redundant because stress is continuous there.
        // Exclude these inactive equations even when retaining all physical unknowns.
        for(auto index : InactiveInterfaceConnects(mesh))
            mesh->ConnectVec()[index].SetCondensed(true);
        mesh->CleanUpUnconnectedNodes();
    }
    fCreated = true;
    return mesh;
}

void TPZSymTensorApproxCreator::GroupAndCondenseElements(TPZMultiphysicsCompMesh *mesh) {
    GroupElements(mesh,fLowOrder);
}

static TPZRefPattern *SymTensorTriangleRef() {
    char buf[] =
    "4     4  "
    "100       TensTri    "
    "0.     0.     0. "
    "1.     0.     0. "
    "0.     1.     0. "
    "0.3333333333333333    0.3333333333333333     0. "
    "2     3     0     1     2 "
    "2     3     0     1     3 "
    "2     3     1     2     3 "
    "2     3     2     0     3 ";
    std::istringstream str(buf);
    TPZRefPattern *refpat = new TPZRefPattern(str);
    return refpat;
}

static TPZRefPattern *SymTensorQuadRef() {
    char buf[] =
    "5     5  "
    "110       TensQuad    "
    "-1.    -1.     0. "
    " 1.    -1.     0. "
    " 1.     1.     0. "
    "-1.     1.     0. "
    " 0.     0.     0. "
    " 3    4     0     1     2     3 "
    "2     3     0     1     4 "
    "2     3     1     2     4 "
    "2     3     2     3     4 "
    "2     3     3     0     4 ";
    std::istringstream str(buf);
    TPZRefPattern *refpat = new TPZRefPattern(str);
    return refpat;
}

void TPZSymTensorApproxCreator::CreateJohnsonMercier(TPZGeoMesh *gmesh) {

    TPZAutoPointer<TPZRefPattern> Trirefpat = SymTensorTriangleRef();
    TPZAutoPointer<TPZRefPattern> Quadrefpat = SymTensorQuadRef();
    {
        int64_t nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel || gel->HasSubElement()) continue;
            if(gel->MaterialId() != fIds.fDomain) continue;
            TPZManVector<TPZGeoEl *> subels;
            if(gel->Type() == EQuadrilateral) {
                gel->SetRefPattern(Quadrefpat);
                gel->Divide(subels);
            } else if (gel->Type() == ETriangle) {
                gel->SetRefPattern(Trirefpat);
                gel->Divide(subels);
            }
        }
    }
}

TPZCompMesh * TPZSymTensorApproxCreator::CreateTensorSpace(TPZGeoMesh *gmesh, int porder, int porderlow) {
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    int dim = gmesh->Dimension();
    cmesh->SetDimModel(dim);
    cmesh->SetDefaultOrder(porder);
    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);

    TPZNullMaterial<STATE> *null = new TPZNullMaterial<STATE>(fIds.fDomain,gmesh->Dimension(),3);
    cmesh->InsertMaterialObject(null);
    TPZNullMaterial<STATE> *materialwrap = new TPZNullMaterial<STATE>(fHybridizationData.fWrapMatId,gmesh->Dimension()-1,3);

    cmesh->AutoBuild();

    cmesh->InsertMaterialObject(materialwrap);

    std::map<int64_t,int64_t> topconnectindex;
    std::set<int> matids = {fIds.fMacroBoundary};
    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != gmesh->Dimension()) continue;
        int ncorner = gel->NCornerNodes();
        int ntop = 0;
        for(int side = 0; side<ncorner; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(matids)) continue;
            ntop++;
            int64_t cindex = cel->ConnectIndex(side);
            int64_t nodeindex = gel->NodeIndex(side);
            if(topconnectindex.find(nodeindex) == topconnectindex.end()) {
                topconnectindex[nodeindex] = cindex;
            } else {
                cindex = topconnectindex[nodeindex];
            }
            cel->SetConnectIndex(side, cindex);

        }
        if(ntop != 1) DebugStop();
    }
    cmesh->ComputeNodElCon();
    cmesh->CleanUpUnconnectedNodes();

    cmesh->Reference()->ResetReference();
    nel = cmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == gmesh->Dimension()) {
            cmesh->SetDefaultOrder(porder);
            if(gel->HasSubElement()) DebugStop();
            int firstside = gel->FirstSide(dim-1);
            int nsides = gel->NSides();
            for (int side = firstside; side < nsides-1 ; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour(gelside.Neighbour());
                if(neighbour.Element()->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
                cel->LoadElementReference();
                TPZCompEl *wrap =cmesh->ApproxSpace().CreateCompEl(neighbour.Element(), *cmesh);
                gel->ResetReference();
                wrap->Reference()->ResetReference();
            }
        }
    }
    if(porderlow >= 0) {
        TPZNullMaterial<STATE> *materialtens = new TPZNullMaterial<STATE>(fIds.fTractionLow,gmesh->Dimension()-1,2);
        cmesh->InsertMaterialObject(materialtens);

        cmesh->SetDefaultOrder(porderlow);
        if(porderlow == 0) {
            cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
            cmesh->SetDimModel(dim-1);
        }
        nel = gmesh->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh->Element(el);
            if(!gel) continue;
            if(gel->MaterialId() == fIds.fTractionLow) {
                cmesh->ApproxSpace().CreateCompEl(gel, *cmesh);
                gel->ResetReference();
            }
        }
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->SetDimModel(dim);
    }
    cmesh->ExpandSolution();
    return cmesh;
}

void TPZSymTensorApproxCreator::AddWrapElements(TPZGeoMesh *gmesh) {
    int64_t nel = gmesh->NElements();
    int dim = gmesh->Dimension();
    if(dim != 2) DebugStop();
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        if(gel->Dimension() != dim || gel->MaterialId() != fIds.fDomain) continue;
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides();
        for (int side = firstside; side < nsides-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElBC gbc(gelside,fHybridizationData.fWrapMatId);
            TPZGeoElSide wrapside = gbc.CreatedElement();
            TPZGeoElSide dispsside = wrapside.HasNeighbour(fIds.fDisplacement);
            if(dispsside) {
                TPZGeoElBC gbc(wrapside,fIds.fInterfaceNegative);
                TPZGeoElSide gelsideintneg = gbc.CreatedElement();
                TPZHybridGeometry hybrid(gbc.CreatedElement()->Index(), wrapside.Element()->Index(), dispsside.Element()->Index());
                AddHybridGeometry(gelsideintneg.Element()->Index(), hybrid);
            } else { // there is no displacement boundary element
                TPZGeoElBC gbc(wrapside,fIds.fInterfacePositive);
                TPZGeoElSide gelsideintpos = gbc.CreatedElement();
                TPZGeoElBC gbc2(gelsideintpos,fIds.fDisplacement);
                TPZGeoElSide gelsidedisp = gbc2.CreatedElement();
                TPZHybridGeometry hybrid(gelsideintpos.Element()->Index(), wrapside.Element()->Index(), gelsidedisp.Element()->Index());
                AddHybridGeometry(gelsideintpos.Element()->Index(), hybrid);
            }
        }
    }
    // check the datastructure consistency
    std::set<int> matids = GetBCMatIds();
    matids.insert(fIds.fDisplacement);
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        if(gel->Dimension() != dim || gel->MaterialId() != fIds.fDomain) continue;
        int firstside = gel->FirstSide(dim-1);
        int nsides = gel->NSides();
        for (int side = firstside; side < nsides-1; side++) {
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSide neigh = gelside.Neighbour();
            if(neigh.Element()->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
            neigh++;
            int neighmat = neigh.Element()->MaterialId();
            if(neighmat != fIds.fInterfaceNegative && neighmat != fIds.fInterfacePositive) DebugStop();
            if(neighmat == fIds.fInterfacePositive) {
                if(!neigh.HasNeighbour(fIds.fInterfacePositive)) DebugStop();
            } else if (neighmat == fIds.fInterfaceNegative) {
                if(!neigh.HasNeighbour(fIds.fInterfacePositive)) DebugStop();
            }
            if(!neigh.HasNeighbour(matids)) DebugStop();
        }
    }
}

void TPZSymTensorApproxCreator::AddWrapElementsTripleHybrid(TPZGeoMesh *gmesh)
{
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->MaterialId() != fHybridizationData.fWrapMatId) continue;
        TPZGeoElSide gelwrap(gel);
        // we only add the triple hybridization elements to the wrap elements that have a neighbour on the macro boundary
        if(!gelwrap.HasNeighbour(fIds.fMacroBoundary)) continue;
        TPZGeoElSide intface = gelwrap.Neighbour();
        int matintface = intface.Element()->MaterialId();
        if(matintface == fIds.fInterfacePositive) {
          // we expect the neighbour of the interface element to be a displacement boundary element
            TPZGeoElSide lagrHigh = intface.Neighbour();
            if(lagrHigh.Element()->MaterialId() != fIds.fDisplacement) DebugStop();
            // add an interface element, then a traction boundary, then an interface element, then a displacement boundary element
            TPZGeoElBC gbcIntface(lagrHigh,fIds.fInterfaceNegativeLow);
            TPZGeoElBC gbctenslow(gbcIntface,fIds.fTractionLow);
            TPZHybridGeometry hybrid(gbcIntface.CreatedElement()->Index(), lagrHigh.Element()->Index(), gbctenslow.CreatedElement()->Index());
            AddHybridGeometry(hybrid);
            TPZGeoElBC gbcIntPos(gbctenslow,fIds.fInterfacePositiveLow);
            TPZGeoElBC gbcDispLow(gbcIntPos,fIds.fDisplacementLow);
            TPZHybridGeometry hybrid2(gbcIntPos.CreatedElement()->Index(), gbctenslow.CreatedElement()->Index(), gbcDispLow.CreatedElement()->Index());
            AddHybridGeometry(hybrid2);
        } else {
          // here we have to switch the hybridgeometry object
            TPZGeoElBC gbcDisp(intface,fIds.fDisplacement);
            TPZHybridGeometry hybrid(intface.Element()->Index(), gelwrap.Element()->Index(), gbcDisp.CreatedElement()->Index());
            AddHybridGeometry(hybrid,false);
            TPZGeoElBC gbcIntface(gbcDisp,fIds.fInterfacePositiveLow);
            TPZGeoElBC gbctenslow(gbcIntface,fIds.fTractionLow);
            TPZGeoElBC gbcIntNeg(gbctenslow,fIds.fInterfaceNegativeLow);
            TPZGeoElSide displow = gelwrap.HasNeighbour(fIds.fDisplacementLow);
            if(!displow) DebugStop();
            TPZHybridGeometry hybrid2(gbcIntNeg.CreatedElement()->Index(), gbctenslow.CreatedElement()->Index(), displow.Element()->Index());
            AddHybridGeometry(hybrid2);
        }
    }

}

TPZCompMesh * TPZSymTensorApproxCreator::CreateDisplacementSpace(TPZGeoMesh *gmesh, int porder, int porderlow) {
    int dim = gmesh->Dimension();
    TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(porder-1);
    if(porder == 1) {
        cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
    } else {
        cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
        cmesh->ApproxSpace().CreateDisconnectedElements(true);
    }

    TPZNullMaterial<STATE> *nullmat = new TPZNullMaterial<STATE>(fIds.fDomain,dim,2);
    cmesh->InsertMaterialObject(nullmat);

    std::set<int> matids= {fIds.fDomain};
    cmesh->AutoBuild(matids);

    cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh->ApproxSpace().CreateDisconnectedElements(true);
    cmesh->SetDefaultOrder(porder);
    gmesh->ResetReference();

    matids.clear();

    for (int id : GetBCMatIds()) {
        cmesh->InsertMaterialObject(new TPZNullMaterial<STATE>(id,dim-1,2));
    }
    TPZNullMaterial<> *nullLagr = new TPZNullMaterial<>(fIds.fDisplacement,dim-1,2);
    cmesh->InsertMaterialObject(nullLagr);
    matids = {fIds.fDisplacement};
    cmesh->AutoBuild(matids);
    if(porderlow >= 0) {
        TPZNullMaterial<> *nullLagr = new TPZNullMaterial<>(fIds.fDisplacementLow,dim-1,2);
        cmesh->InsertMaterialObject(nullLagr);
        if(porderlow == 0) {
            cmesh->ApproxSpace().SetAllCreateFunctionsDiscontinuous();
            cmesh->SetDimModel(dim-1); // otherwise dim-1 elements will not be created
        } else {
            cmesh->ApproxSpace().SetAllCreateFunctionsContinuous();
            cmesh->ApproxSpace().CreateDisconnectedElements(true);
        }
        cmesh->SetDefaultOrder(porderlow);
        matids = {fIds.fDisplacementLow};
        cmesh->AutoBuild(matids);

    }
    matids.clear();
    matids = GetBCMatIds();
    int64_t nel = cmesh->NElements();
    for(int64_t el = 0; el<nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(!gel || (gel->MaterialId() != fIds.fDisplacement && porderlow < 0)) continue;
        if(!gel || (gel->MaterialId() != fIds.fDisplacementLow && porderlow >= 0)) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neigh = gelside.HasNeighbour(matids);
        if(neigh) {
            cel->LoadElementReference();
            TPZCompEl *celbc = cmesh->ApproxSpace().CreateCompEl(neigh.Element(), *cmesh);
            if(porderlow == 0) {
                // Overlapping discontinuous elements may initially have no connect.
                // The boundary condition must use the existing trace unknown.
                celbc->SetConnectIndex(0, cel->ConnectIndex(0));
            }
            gel->ResetReference();
            neigh.Element()->ResetReference();
        }
    }
    if(porderlow == 0) {
        cmesh->ComputeNodElCon();
        cmesh->CleanUpUnconnectedNodes();
    }
    cmesh->SetDimModel(dim);
    cmesh->ExpandSolution();
    return cmesh;
}

TPZMultiphysicsCompMesh * TPZSymTensorApproxCreator::GenerateMultiphysicsMesh(TPZGeoMesh *gmesh, int porder, int porderlow)
{
    int dim = gmesh->Dimension();
    TPZMultiphysicsCompMesh *mfmesh = new TPZMultiphysicsCompMesh(gmesh);
    TPZManVector<TPZCompMesh *> meshvec(2,0);
    meshvec[0] = CreateTensorSpace(gmesh, porder, porderlow);
    meshvec[1] = CreateDisplacementSpace(gmesh, porder, porderlow);

    TPZNullMaterialCS<STATE> *dispinterface = new TPZNullMaterialCS<STATE> (fIds.fDisplacement,dim-1, dim);
    mfmesh->InsertMaterialObject(dispinterface);
    for (auto &entry : fMaterialVec) {
        mfmesh->InsertMaterialObject(entry.second);
    }
    TPZNullMaterialCS<STATE> *materialwrap = new TPZNullMaterialCS<STATE>(fHybridizationData.fWrapMatId,dim-1,dim);
    mfmesh->InsertMaterialObject(materialwrap);
    if(porderlow >= 0) {
        TPZNullMaterialCS<STATE> *materialtenslow = new TPZNullMaterialCS<STATE>(fIds.fTractionLow,dim-1,dim);
        mfmesh->InsertMaterialObject(materialtenslow);
        TPZNullMaterialCS<STATE> *materialdisplow = new TPZNullMaterialCS<STATE>(fIds.fDisplacementLow,dim-1,dim);
        mfmesh->InsertMaterialObject(materialdisplow);
    }
    mfmesh->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    mfmesh->BuildMultiphysicsSpace(meshvec);
    AddInterfaceElements(mfmesh,porderlow);
    return mfmesh;
}

void TPZSymTensorApproxCreator::AddInterfaceElements(TPZMultiphysicsCompMesh *mfmesh, int porderlow) {
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int dim = gmesh->Dimension();
    TPZInterfaceSymTensor *sym1 = new TPZInterfaceSymTensor(fIds.fInterfacePositive, dim-1);
    sym1->SetMultiplier(1);
    mfmesh->InsertMaterialObject(sym1);
    TPZInterfaceSymTensor *sym2 = new TPZInterfaceSymTensor(fIds.fInterfaceNegative, dim-1);
    sym2->SetMultiplier(1);
    mfmesh->InsertMaterialObject(sym2);

    if(porderlow >= 0) {
        int nstate = 2;
        TPZLagrangeMultiplierCS<STATE> *matintpos = new TPZLagrangeMultiplierCS<STATE>(fIds.fInterfacePositiveLow, dim-1, nstate);
        matintpos->SetMultiplier(1);
        mfmesh->InsertMaterialObject(matintpos);
        TPZLagrangeMultiplierCS<STATE> *matintneg = new TPZLagrangeMultiplierCS<STATE>(fIds.fInterfaceNegativeLow, dim-1, nstate);
        matintneg->SetMultiplier(-1);
        mfmesh->InsertMaterialObject(matintneg);
    }

    mfmesh->LoadReferences();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel ; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel) continue;
        if(gel->MaterialId() == fHybridizationData.fWrapMatId) {
            TPZGeoElSide gelside(gel);
            TPZCompElSide celsidewrap = gelside.Reference();
            TPZGeoElSide neighint = gelside.Neighbour();
            TPZGeoEl *neighgel = neighint.Element();
            if(neighgel->MaterialId() != fIds.fInterfaceNegative && neighgel->MaterialId() != fIds.fInterfacePositive) DebugStop();
            TPZGeoElSide neighdisp = gelside.HasNeighbour(fIds.fDisplacement);
            if(!neighdisp) DebugStop();
            auto intface = new TPZMultiphysicsInterfaceElement(*mfmesh,neighgel,celsidewrap,neighdisp.Reference());
            if(porderlow >= 0 && gelside.HasNeighbour(fIds.fMacroBoundary)) {

                TPZGeoElSide gelintfaceloworder1 = neighdisp.Neighbour();
                int intfacematid = gelintfaceloworder1.Element()->MaterialId();
                if(intfacematid != fIds.fInterfaceNegativeLow && intfacematid != fIds.fInterfacePositiveLow) DebugStop();
                TPZGeoElSide geltractionloworder = gelintfaceloworder1.Neighbour();
                if(geltractionloworder.Element()->MaterialId() != fIds.fTractionLow) DebugStop();
                new TPZMultiphysicsInterfaceElement(*mfmesh,gelintfaceloworder1.Element(),neighdisp.Reference(),geltractionloworder.Reference());
                TPZGeoElSide gelintfaceloworder2 = geltractionloworder.Neighbour();
                intfacematid = gelintfaceloworder2.Element()->MaterialId();
                if(intfacematid != fIds.fInterfaceNegativeLow && intfacematid != fIds.fInterfacePositiveLow) DebugStop();
                TPZGeoElSide geldisploworder = gelintfaceloworder2.HasNeighbour(fIds.fDisplacementLow);
                if(!geldisploworder) DebugStop();
                new TPZMultiphysicsInterfaceElement(*mfmesh,gelintfaceloworder2.Element(),geltractionloworder.Reference(),geldisploworder.Reference());
            }
        }
    }

}

std::set<int64_t> TPZSymTensorApproxCreator::InactiveInterfaceConnects(TPZMultiphysicsCompMesh *mfmesh) {
    std::set<int64_t> inactive;
    TPZGeoMesh *gmesh = mfmesh->Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        if(!gel || gel->MaterialId() != fIds.fDomain) continue;
        int nnodes = gel->NNodes();
        int firstside = 0;
        int lastside = gel->NCornerNodes();
        for (int side = firstside; side <lastside; side++) {
            TPZGeoElSide gelside(gel,side);
            if(gelside.HasNeighbour(fIds.fMacroBoundary)) continue;
            for(auto neigh = gelside.Neighbour(); neigh!=gelside; neigh++) {
                TPZGeoEl *neighgel = neigh.Element();
                if(neighgel->MaterialId() == fIds.fDisplacement) {
                    TPZCompEl *cel = neighgel->Reference();
                    int neighside = neigh.Side();
                    inactive.insert(cel->ConnectIndex(neighside));
                }
            }
        }
    }
    mfmesh->CleanUpUnconnectedNodes();
    return inactive;
}

void TPZSymTensorApproxCreator::GroupElements(TPZMultiphysicsCompMesh *mfmesh, int ploworder) {
    std::map<int64_t,std::set<int64_t> > nodetoel;
    auto inactive = InactiveInterfaceConnects(mfmesh);

    TPZGeoMesh *gmesh = mfmesh->Reference();
    int dim = gmesh->Dimension();
    int64_t nel = gmesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if(!gel || gel->HasSubElement()) continue;
        if(gel->Dimension() != dim || gel->MaterialId() != fIds.fDomain) continue;
        int nnodes = gel->NCornerNodes();
        int side = -1;
        for(side = 0 ; side<nnodes; side++) {
            TPZGeoElSide gelside(gel,side);
            auto neigh = gelside.HasNeighbour(fIds.fMacroBoundary);
            if(!neigh) break;
        }
        if(side == nnodes) DebugStop();
        int64_t nodeindex = gel->NodeIndex(side);
        if(nodetoel.find(nodeindex) != nodetoel.end()) continue;
        nodetoel[nodeindex].insert(el);
        TPZGeoElSide gelside(gel,side);
        for (TPZGeoElSide neigh = gelside.Neighbour(); neigh != gelside; neigh++) {
            TPZGeoEl *gelneigh = neigh.Element();
            if(gelneigh->HasSubElement()) continue;
            int64_t neighindex = gelneigh->Index();
            nodetoel[nodeindex].insert(neighindex);
        }
    }
    for (auto &it : nodetoel) {
        std::set<int64_t> elset = it.second;
        for (auto elindex : elset) {
            TPZGeoEl *gel = gmesh->Element(elindex);
            if(gel->Dimension() != dim || gel->MaterialId() != fIds.fDomain) continue;
            int firstside = gel->FirstSide(dim-1);
            int lastside = gel->NSides()-1;
            for(int side = firstside; side < lastside; side++) {
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neigh1 = gelside.Neighbour();
                if(neigh1.Element()->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
                TPZGeoElSide neigh2 = neigh1.Neighbour();
                int neigh2matid = neigh2.Element()->MaterialId();
                if(neigh2matid != fIds.fInterfaceNegative && neigh2matid != fIds.fInterfacePositive) DebugStop();
                it.second.insert(neigh1.Element()->Index());
                it.second.insert(neigh2.Element()->Index());
                if(ploworder >= 0 && gelside.HasNeighbour(fIds.fMacroBoundary)) {
                    TPZGeoElSide neighdisp = neigh2.Neighbour();
                    int dispmat = neighdisp.Element()->MaterialId();
                    if(dispmat != fIds.fDisplacement) DebugStop();
                    it.second.insert(neighdisp.Element()->Index());
                    TPZGeoElSide neighIntLow1 = neighdisp.Neighbour();
                    int matint1 = neighIntLow1.Element()->MaterialId();
                    if(matint1 != fIds.fInterfaceNegativeLow && matint1 != fIds.fInterfacePositiveLow) DebugStop();
                    it.second.insert(neighIntLow1.Element()->Index());
                    TPZGeoElSide neighTensionLow = neighIntLow1.Neighbour();
                    int matstress = neighTensionLow.Element()->MaterialId();
                    if(matstress != fIds.fTractionLow) DebugStop();
                    it.second.insert(neighTensionLow.Element()->Index());
                    TPZGeoElSide neighIntLow2 = neighTensionLow.Neighbour();
                    int matint2 = neighIntLow2.Element()->MaterialId();
                    if(matint2 != fIds.fInterfaceNegativeLow && matint2 != fIds.fInterfacePositiveLow) DebugStop();
                    it.second.insert(neighIntLow2.Element()->Index());
                }
            }
        }
    }
    std::set<TPZElementGroup *> groupset;
    for (auto it : nodetoel) {
        std::set<int64_t> &elset = it.second;
        TPZElementGroup *elgroup = new TPZElementGroup(*mfmesh);
        groupset.insert(elgroup);
        for(auto gelindex : elset) {
            auto cel = gmesh->Element(gelindex)->Reference();
            if(!cel) DebugStop();
            elgroup->AddElement(cel);
        }
    }
    mfmesh->ComputeNodElCon();
    for(auto it : inactive) {
        TPZConnect &c = mfmesh->ConnectVec()[it];
        int64_t ind = mfmesh->AllocateNewConnect(0,2,0);
        TPZFMatrix<STATE> rest(2,0);
        c.AddDependency(it,ind,rest,0,0,2,0);
        // c.IncrementElConnected();
    }
    for(auto it : groupset) {
        new TPZCondensedCompElT<STATE>(it);
    }
    for(auto it : inactive) {
        TPZConnect &c = mfmesh->ConnectVec()[it];
        c.SetCondensed(true);
    }
    mfmesh->CleanUpUnconnectedNodes();
}

