//
//  TPZHybridizeHDiv.cpp
//  ErrorEstimate
//
//  Created by Philippe Devloo on 16/05/18.
//

#include "TPZHybridizeHDiv.h"
#include "pzconnect.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZMaterial.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"
#include "pzintel.h"
#include "pzgeoelbc.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzgeoelside.h"
#include "pztrnsform.h"
#include "tpzintpoints.h"
#include "TPZLagrangeMultiplier.h"
#include <algorithm>

TPZHybridizeHDiv::TPZHybridizeHDiv(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    ComputePeriferalMaterialIds(meshvec_Hybrid);
    ComputeNState(meshvec_Hybrid);
}

void TPZHybridizeHDiv::ComputeNState(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    if (meshvec_Hybrid.size() > 1 && meshvec_Hybrid[1]->NMaterials() > 0) {
        NState = meshvec_Hybrid[1]->MaterialVec().begin()->second->NStateVariables();
    }
}

void TPZHybridizeHDiv::ComputePeriferalMaterialIds(TPZVec<TPZCompMesh*>& meshvec_Hybrid) {
    int maxMatId = std::numeric_limits<int>::min();
    for (auto &mesh : meshvec_Hybrid) {
        for (auto &mat : mesh->MaterialVec()) {
            maxMatId = std::max(maxMatId, mat.first);
        }
    }
    if (maxMatId == std::numeric_limits<int>::min()) {
        maxMatId = 0;
    }
    HDivWrapMatid = maxMatId + 1;
    LagrangeInterface = maxMatId + 2;
    InterfaceMatid = maxMatId + 3;
}

/// split the connect between two neighbouring elements

std::tuple<int64_t, int> TPZHybridizeHDiv::SplitConnects(const TPZCompElSide &left, const TPZCompElSide &right, TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (HDivWrapMatid == 0 || LagrangeInterface == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    //TPZCompMesh *pressuremesh = meshvec[1];
    //TPZGeoMesh *gmesh = fluxmesh->Reference();
    TPZGeoElSide gleft(left.Reference());
    TPZGeoElSide gright(right.Reference());
    TPZInterpolatedElement *intelleft = dynamic_cast<TPZInterpolatedElement *> (left.Element());
    TPZInterpolatedElement *intelright = dynamic_cast<TPZInterpolatedElement *> (right.Element());
    intelleft->SetSideOrient(left.Side(), 1);
    intelright->SetSideOrient(right.Side(), 1);
    TPZStack<TPZCompElSide> equalright;
    TPZConnect &cleft = intelleft->SideConnect(0, left.Side());

    if (cleft.HasDependency()) {
        // check whether the wrap of the large element was already created
        gright.EqualLevelCompElementList(equalright,1,0);
        // only one wrap element should exist
#ifdef PZDEBUG
        if(equalright.size() > 1)
        {
            DebugStop();
        }
        if(equalright.size()==1)
        {
            TPZGeoEl *equalgel = equalright[0].Element()->Reference();
            if(equalgel->Dimension() != fluxmesh->Dimension()-1)
            {
                DebugStop();
            }
        }
#endif
        // reset the reference of the wrap element
        if(equalright.size()) equalright[0].Element()->Reference()->ResetReference();
        cleft.RemoveDepend();
    }
    else
    {
        int64_t index = fluxmesh->AllocateNewConnect(cleft);
        TPZConnect &newcon = fluxmesh->ConnectVec()[index];
        cleft.DecrementElConnected();
        newcon.ResetElConnected();
        newcon.IncrementElConnected();
        newcon.SetSequenceNumber(fluxmesh->NConnects() - 1);

        int rightlocindex = intelright->SideConnectLocId(0, right.Side());
        intelright->SetConnectIndex(rightlocindex, index);
    }
    int sideorder = cleft.Order();
    fluxmesh->SetDefaultOrder(sideorder);
    // create HDivBound on the sides of the elements
    TPZCompEl *wrap1, *wrap2;
    {
        intelright->Reference()->ResetReference();
        intelleft->LoadElementReference();
        intelleft->SetPreferredOrder(sideorder);
        TPZGeoElBC gbc(gleft, HDivWrapMatid);
        int64_t index;
        wrap1 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        if(cleft.Order() != sideorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap1);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelleft->Reference()->ResetReference();
        wrap1->Reference()->ResetReference();
    }
    // if the wrap of the large element was not created...
    if(equalright.size() == 0)
    {
        intelleft->Reference()->ResetReference();
        intelright->LoadElementReference();
        TPZConnect &cright = intelright->SideConnect(0,right.Side());
        int rightprevorder = cright.Order();
        intelright->SetPreferredOrder(cright.Order());
        TPZGeoElBC gbc(gright, HDivWrapMatid);
        int64_t index;
        wrap2 = fluxmesh->ApproxSpace().CreateCompEl(gbc.CreatedElement(), *fluxmesh, index);
        if(cright.Order() != rightprevorder)
        {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (wrap2);
        int wrapside = gbc.CreatedElement()->NSides() - 1;
        intel->SetSideOrient(wrapside, 1);
        intelright->Reference()->ResetReference();
        wrap2->Reference()->ResetReference();
    }
    else
    {
        wrap2 = equalright[0].Element();
    }
    wrap1->LoadElementReference();
    wrap2->LoadElementReference();
    int64_t pressureindex;
    int pressureorder;
    {
        TPZGeoElBC gbc(gleft, LagrangeInterface);
        pressureindex = gbc.CreatedElement()->Index();
        pressureorder = sideorder;
    }
    intelleft->LoadElementReference();
    intelright->LoadElementReference();
    return std::make_tuple(pressureindex, pressureorder);
}

/// find an element which shares the connect and has the same dimension
// else return 0

TPZCompElSide TPZHybridizeHDiv::RightElement(TPZInterpolatedElement *intel, int side) {
    bool isrestrained = false;
    {
        TPZConnect &c = intel->SideConnect(0, side);
        if (c.HasDependency()) {
            isrestrained = true;
        }
    }
    TPZGeoEl *gel = intel->Reference();
    TPZGeoElSide gelside(gel, side);

    if (isrestrained == true) {
        /// if the side is restrained we will hybridize between the element and the larger element
        TPZCompElSide celside = gelside.LowerLevelCompElementList2(1);
        if (!celside) DebugStop();
        TPZGeoEl *neigh = celside.Element()->Reference();
        /// we assume that a larger element should not be a boundary element
        if (neigh->Dimension() != gel->Dimension()) {
            DebugStop();
        }
        return celside;
    } else {
        // if the connect is not restrained
        // - there should be only one neighbour
        //   if there is more than one neighbour the element side has already been hybridized
        // - the neighbour should be of the same dimension
        //   if the neighbour is of lower dimension it is a boundary element
        TPZStack<TPZCompElSide> celstack;
        gelside.EqualLevelCompElementList(celstack, 1, 0);
        if (celstack.size() == 1) {
            TPZGeoEl *neigh = celstack[0].Element()->Reference();
            if (neigh->Dimension() == gel->Dimension()) {
                return celstack[0];
            }
        }
    }
    return TPZCompElSide();
}


/// split the connects between flux elements and create a dim-1 pressure element

void TPZHybridizeHDiv::HybridizeInternalSides(TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    TPZGeoMesh *gmesh = fluxmesh->Reference();
    int dim = gmesh->Dimension();
    gmesh->ResetReference();
    fluxmesh->LoadReferences();
    int64_t nel = fluxmesh->NElements();
    std::list<std::tuple<int64_t, int> > pressures;
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel || intel->Reference()->Dimension() != dim) {
            continue;
        }
        // loop over the side of dimension dim-1
        TPZGeoEl *gel = intel->Reference();
        for (int side = gel->NCornerNodes(); side < gel->NSides() - 1; side++) {
            if (gel->SideDimension(side) != dim - 1) {
                continue;
            }
            TPZCompElSide celside(intel, side);
            TPZCompElSide neighcomp = RightElement(intel, side);
            if (neighcomp) {
                pressures.push_back(SplitConnects(celside, neighcomp, meshvec_Hybrid));
            }
        }
    }
    fluxmesh->InitializeBlock();
    fluxmesh->ComputeNodElCon();
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    gmesh->ResetReference();
    pressuremesh->SetDimModel(gmesh->Dimension()-1);
    for (auto pindex : pressures) {
        int64_t elindex;
        int order;
        std::tie(elindex, order) = pindex;
        TPZGeoEl *gel = gmesh->Element(elindex);
        int64_t celindex;
        TPZCompEl *cel = pressuremesh->ApproxSpace().CreateCompEl(gel, *pressuremesh, celindex);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        TPZCompElDisc *intelDisc = dynamic_cast<TPZCompElDisc *> (cel);
        if (intel){
            intel->PRefine(order);
//            intel->SetSideOrder(gel->NSides() - 1, order);
        } else if (intelDisc) {
            intelDisc->SetDegree(order);
            intelDisc->SetTrueUseQsiEta();
        } else {
            DebugStop();
        }
        int n_connects = cel->NConnects();
        for (int i = 0; i < n_connects; ++i) {
            cel->Connect(i).SetLagrangeMultiplier(2);
        }
        gel->ResetReference();
    }
    pressuremesh->InitializeBlock();
    pressuremesh->SetDimModel(gmesh->Dimension());
}

void TPZHybridizeHDiv::CreateInterfaceElements(TPZCompMesh *cmesh_Hybrid, TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (InterfaceMatid == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    int dim = pressuremesh->Dimension();
    TPZVec<TPZCompEl *> celpressure(pressuremesh->NElements(), 0);
    for (int64_t el = 0; el < pressuremesh->NElements(); el++) {
        TPZCompEl *cel = pressuremesh->Element(el);
        if (cel && cel->Reference() && cel->Reference()->Dimension() == dim - 1) {
            celpressure[el] = cel;
        }
    }
    cmesh_Hybrid->Reference()->ResetReference();
    cmesh_Hybrid->LoadReferences();
    for (auto cel : celpressure) {
        if (!cel) continue;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoEl *gel = cel->Reference();
        TPZCompEl *mphysics = gel->Reference();
        TPZGeoElSide gelside(gel, gel->NSides() - 1);
        TPZCompElSide celside = gelside.Reference();
        gelside.EqualLevelCompElementList(celstack, 0, 0);
        int count = 0;
        for (auto &celstackside : celstack) {
            if (celstackside.Reference().Element()->Dimension() == dim - 1) {
                TPZGeoElBC gbc(gelside, InterfaceMatid);
                // check if the right side has a dependency
                TPZCompEl *celneigh = celstackside.Element();
                if (celneigh->NConnects() != 1) {
                    DebugStop();
                }
                int64_t index;
                TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), index, celside, celstackside);
                count++;
            }
        }
        if (count == 1)
        {
            TPZCompElSide clarge = gelside.LowerLevelCompElementList2(false);
            if(!clarge) DebugStop();
            TPZGeoElSide glarge = clarge.Reference();
            if (glarge.Element()->Dimension() == dim) {
                TPZGeoElSide neighbour = glarge.Neighbour();
                while(neighbour != glarge)
                {
                    if (neighbour.Element()->Dimension() < dim) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if(neighbour == glarge) DebugStop();
                glarge = neighbour;
            }
            clarge = glarge.Reference();
            if(!clarge) DebugStop();
            TPZGeoElBC gbc(gelside, InterfaceMatid);

            int64_t index;
            TPZMultiphysicsInterfaceElement *intface = new TPZMultiphysicsInterfaceElement(*cmesh_Hybrid, gbc.CreatedElement(), index, celside, clarge);
            count++;
        }
        if (count != 2 && count != 0) {
            DebugStop();
        }
    }
    pressuremesh->InitializeBlock();
}

TPZCompMesh * TPZHybridizeHDiv::CreateMultiphysicsMesh(TPZCompMesh *cmesh_HDiv, TPZVec<TPZCompMesh *> &meshvector_Hybrid, double Lagrange_term_multiplier /* = 1. */) {
    TPZGeoMesh *gmesh = cmesh_HDiv->Reference();
    TPZCompMesh *cmesh_Hybrid = new TPZCompMesh(gmesh);
    cmesh_HDiv->CopyMaterials(*cmesh_Hybrid);
    InsertPeriferalMaterialObjects(cmesh_Hybrid, Lagrange_term_multiplier);
    cmesh_Hybrid->ApproxSpace().SetAllCreateFunctionsMultiphysicElem();
    cmesh_Hybrid->AutoBuild();

    TPZBuildMultiphysicsMesh::AddElements(meshvector_Hybrid, cmesh_Hybrid);
    TPZBuildMultiphysicsMesh::AddConnects(meshvector_Hybrid, cmesh_Hybrid);
    cmesh_Hybrid->LoadReferences();
    return cmesh_Hybrid;
}

/// group and condense the elements

void TPZHybridizeHDiv::GroupElements(TPZCompMesh *cmesh) {
    int64_t nel = cmesh->NElements();
    int64_t nconnects = cmesh->NConnects();
    TPZVec<TPZElementGroup *> groupindex(nconnects, 0);
    int dim = cmesh->Dimension();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference() || cel->Reference()->Dimension() != dim) {
            continue;
        }
        int64_t index;
        TPZElementGroup *elgr = new TPZElementGroup(*cmesh, index);
//        std::cout << "Created group " << index << std::endl;
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
        for (auto cindex : connectlist) {
#ifdef PZDEBUG
            if (groupindex[cindex] != 0) {
                DebugStop();
            }
#endif
            groupindex[cindex] = elgr;
        }
        elgr->AddElement(cel);
    }
//    std::cout << "Groups of connects " << groupindex << std::endl;
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        if (!cel || !cel->Reference()) {
            continue;
        }
        TPZStack<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
//        std::cout << "Analysing element " << cel->Index();
        for (auto cindex : connectlist) {
            if (groupindex[cindex] != 0) {
                groupindex[cindex]->AddElement(cel);
//                std::cout << " added to " << groupindex[cindex]->Index() << " with size " << groupindex[cindex]->GetElGroup().size();
                break;
            }
        }
//        std::cout << std::endl;
    }
    cmesh->ComputeNodElCon();
    nel = cmesh->NElements();
    for (TPZCompEl *cel : cmesh->ElementVec()) {
        TPZElementGroup *elgr = dynamic_cast<TPZElementGroup *> (cel);
        if (elgr) {
            TPZCondensedCompEl *cond = new TPZCondensedCompEl(elgr);
            cond->SetKeepMatrix(false);
        }
    }
}

/// insert the material objects for HDivWrap and LagrangeInterface in the atomic meshes

void TPZHybridizeHDiv::InsertPeriferalMaterialObjects(TPZVec<TPZCompMesh *> &meshvec_Hybrid) {
    if (LagrangeInterface == 0 || HDivWrapMatid == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }

    TPZFNMatrix<1, STATE> xk(NState, NState, 0.), xb(NState, NState, 0.), xc(NState, NState, 0.), xf(NState, 1, 0.);

    TPZCompMesh *pressuremesh = meshvec_Hybrid[1];
    int dim = pressuremesh->Dimension();
    
    if (!pressuremesh->FindMaterial(LagrangeInterface)) {
        if(dim == 2)
        {
            auto matPerif = new TPZMat1dLin(LagrangeInterface);
            matPerif->SetMaterial(xk, xc, xb, xf);
            pressuremesh->InsertMaterialObject(matPerif);
        }
        else if(dim == 3)
        {
            auto matPerif = new TPZMat2dLin(LagrangeInterface);
            matPerif->SetMaterial(xk, xc, xf);
            pressuremesh->InsertMaterialObject(matPerif);
        }
        else
        {
            DebugStop();
        }
    }
    TPZCompMesh *fluxmesh = meshvec_Hybrid[0];
    if (!fluxmesh->FindMaterial(HDivWrapMatid)) {
        if(dim == 2)
        {
            auto matPerif = new TPZMat1dLin(HDivWrapMatid);
            matPerif->SetMaterial(xk, xc, xb, xf);
            fluxmesh->InsertMaterialObject(matPerif);
        }
        else if(dim == 3)
        {
            auto matPerif = new TPZMat2dLin(HDivWrapMatid);
            matPerif->SetMaterial(xk, xc, xf);
            fluxmesh->InsertMaterialObject(matPerif);

        }
        else
        {
            DebugStop();
        }
    }
}

/// insert the material objects for HDivWrap, LagrangeInterface and InterfaceMatid in the multiphysics mesh

void TPZHybridizeHDiv::InsertPeriferalMaterialObjects(TPZCompMesh *cmesh_Hybrid, double Lagrange_term_multiplier /* = 1. */) {
    if (LagrangeInterface == 0 || HDivWrapMatid == 0 || InterfaceMatid == 0) {
        std::cerr << "Using uninitialized TPZHybridizeHDiv object. You need to call ComputePeriferalMaterialIds function first!" << std::endl;
        DebugStop();
    }
    TPZFNMatrix<1, STATE> xk(NState, NState, 0.), xb(NState, NState, 0.), xc(NState, NState, 0.), xf(NState, 1, 0.);
    int dim = cmesh_Hybrid->Dimension();

    if (!cmesh_Hybrid->FindMaterial(LagrangeInterface)) {
        if(dim == 2)
        {
            auto matPerif = new TPZMat1dLin(LagrangeInterface);
            matPerif->SetMaterial(xk, xc, xb, xf);
            cmesh_Hybrid->InsertMaterialObject(matPerif);
        }
        else if(dim == 3)
        {
            auto matPerif = new TPZMat2dLin(LagrangeInterface);
            matPerif->SetMaterial(xk, xc, xf);
            cmesh_Hybrid->InsertMaterialObject(matPerif);
        }
        else
        {
            DebugStop();
        }
    }
    if (!cmesh_Hybrid->FindMaterial(HDivWrapMatid)) {
        if(dim == 2)
        {
            auto matPerif = new TPZMat1dLin(HDivWrapMatid);
            matPerif->SetMaterial(xk, xc, xb, xf);
            cmesh_Hybrid->InsertMaterialObject(matPerif);
        }
        else if(dim == 3)
        {
            auto matPerif = new TPZMat2dLin(HDivWrapMatid);
            matPerif->SetMaterial(xk, xc, xf);
            cmesh_Hybrid->InsertMaterialObject(matPerif);
        }
        else
        {
            DebugStop();
        }
    }
    if (!cmesh_Hybrid->FindMaterial(InterfaceMatid)) {
        TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(InterfaceMatid, dim - 1, NState);
        matleft->SetMultiplier(Lagrange_term_multiplier);
        cmesh_Hybrid->InsertMaterialObject(matleft);
    }
}

std::tuple<TPZCompMesh*, TPZVec<TPZCompMesh*> > TPZHybridizeHDiv::Hybridize(TPZCompMesh* cmesh_HDiv, TPZVec<TPZCompMesh*>& meshvec_HDiv, bool group_elements, double Lagrange_term_multiplier /* = 1. */) {
    TPZManVector<TPZCompMesh *, 2> meshvec_Hybrid(meshvec_HDiv.size(), 0);
    for (int i = 0; i < meshvec_HDiv.size(); i++) {
        meshvec_Hybrid[i] = meshvec_HDiv[i]->Clone();
    }
    ComputePeriferalMaterialIds(meshvec_Hybrid);
    ComputeNState(meshvec_Hybrid);
    /// insert the material objects for HDivWrap and LagrangeInterface
    InsertPeriferalMaterialObjects(meshvec_Hybrid);
    HybridizeInternalSides(meshvec_Hybrid);
    TPZCompMesh *cmesh_Hybrid = CreateMultiphysicsMesh(cmesh_HDiv, meshvec_Hybrid, Lagrange_term_multiplier);
    CreateInterfaceElements(cmesh_Hybrid, meshvec_Hybrid);
    if (group_elements){
        GroupElements(cmesh_Hybrid);
    }
    return std::make_tuple(cmesh_Hybrid, meshvec_Hybrid);
}

static void CompareFluxes(TPZCompElSide &left, TPZCompElSide &right, std::ostream &out)
{
    TPZGeoElSide gleft = left.Reference();
    TPZGeoElSide gright = right.Reference();
    TPZTransform<REAL> trleft = gleft.Element()->SideToSideTransform(gleft.Side(),gleft.Element()->NSides()-1);
    TPZTransform<REAL> trright = gright.Element()->SideToSideTransform(gright.Side(),gright.Element()->NSides()-1);
    int sidedim = gleft.Dimension();
    int meshdim = gleft.Element()->Mesh()->Dimension();
    TPZTransform<REAL> tr(sidedim); 
    gleft.SideTransform3(gright,tr);
    TPZIntPoints *intpts = gleft.Element()->CreateSideIntegrationRule(gleft.Side(),5);
    TPZManVector<REAL,3> ptleft(sidedim,0.),ptright(sidedim,0.),ptleftvol(meshdim,0),ptrightvol(meshdim,0);
    REAL weight;
    int npoints = intpts->NPoints();
    for(int ip = 0; ip<npoints; ip++)
    {
        intpts->Point(ip,ptleft,weight);
        tr.Apply(ptleft,ptright);
        trleft.Apply(ptleft,ptleftvol);
        trright.Apply(ptright,ptrightvol);
        TPZManVector<STATE,3> fluxleft(3,0.),fluxright(3,0.);
        left.Element()->Solution(ptleftvol,0,fluxleft);
        right.Element()->Solution(ptrightvol,0,fluxright);
        TPZManVector<REAL,3> normal(3,0.);
        {
            TPZFNMatrix<9,REAL> axes(sidedim,3),jacobian(sidedim,sidedim),jacinv(sidedim,sidedim);
            REAL detjac;
            gleft.Jacobian(ptleft,jacobian,axes,detjac,jacinv);
            normal[0] = axes(0,1)*axes(1,2)-axes(0,2)*axes(1,1);
            normal[1] = -axes(0,0)*axes(1,2)+axes(0,2)*axes(1,0);
            normal[2] = axes(0,0)*axes(1,1)-axes(0,1)*axes(1,0);
        }
        STATE fluxnormalleft = 0.,fluxnormalright = 0.;
        for(int i=0; i<3; i++)
        {
            fluxnormalleft += fluxleft[i]*normal[i];
            fluxnormalright += fluxright[i]*normal[i];
        }
        REAL diff = fabs(fluxnormalleft-fluxnormalright);
        if(diff > 1.e-6)
        {
            TPZManVector<REAL,3> leftx(3),rightx(3);
            gleft.X(ptleft,leftx);
            gright.X(ptright,rightx);
            out << "Left geo element index " << gleft.Element()->Index();
            out << " Right geo element index " << gright.Element()->Index();
            out << " xleft " << leftx << " xright " << rightx << std::endl;
            out << " fluxleft " << fluxleft << std::endl;
            out << " fluxright" << fluxright << std::endl;
            out << "normal " << normal << std::endl;
            out << "fluxnormalleft " << fluxnormalleft << std::endl;
            out << "fluxnormalright " << fluxnormalright << std::endl;
        }
    }
    delete intpts;
}

/// verify the consistency of the solution of the flux mesh
void TPZHybridizeHDiv::VerifySolutionConsistency(TPZCompMesh *fluxmesh, std::ostream &out)
{
    int64_t nel = fluxmesh->NElements();
    fluxmesh->Reference()->ResetReference();
    fluxmesh->LoadReferences();
    int meshdim = fluxmesh->Dimension();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZCompEl *cel = fluxmesh->Element(el);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!intel) continue;
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() != meshdim) continue;
        int nsides = gel->NSides();
        for(int side=0; side<nsides; side++)
        {
            if(gel->SideDimension(side) != meshdim-1) continue;
            TPZStack<TPZCompElSide> celsides;
            TPZGeoElSide gelside(gel,side);
            gelside.EqualLevelCompElementList(celsides,1,0);

            {
                TPZCompElSide celneigh = gelside.LowerLevelCompElementList2(1);
                if(celneigh) 
                {
                    celsides.push_back(celneigh);
                    celneigh.Reference().EqualLevelCompElementList(celsides,1,0);
                }
            }
            TPZCompElSide celneighselect;
            int nummatch = 0;
            for(int i=0; i<celsides.size(); i++)
            {
                if(celsides[i].Reference().Element()->Dimension() == meshdim) 
                {
                    celneighselect = celsides[i];
                    nummatch++;
                }
            }
            if(nummatch > 1) DebugStop();
            TPZCompElSide celside(gelside.Reference());
            if(celneighselect)
            {
                CompareFluxes(celside,celneighselect,out);
            }
        }
    }
}
