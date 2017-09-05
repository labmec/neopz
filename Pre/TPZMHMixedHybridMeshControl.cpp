//
//  TPZMHMixedHybridMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZVecL2.h"
#include "pzbndcond.h"
#include "TPZMatLaplacian.h"
#include "TPZLagrangeMultiplier.h"


#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>
#include <algorithm>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzelementgroup.h"

#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"


// tototo
TPZMHMixedHybridMeshControl::TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<long> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices),
    fHDivWrapperMatId(105)
{
}

TPZMHMixedHybridMeshControl::TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<long> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices),
    fHDivWrapperMatId(105)
{
}

/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedHybridMeshControl::InsertPeriferalMaterialObjects()
{
    TPZMHMixedMeshControl::InsertPeriferalMaterialObjects();
    
    TPZMaterial *mat = fCMesh->MaterialVec()[1];
    if (!mat) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.), val2Flux(1,1,0.);
    int typePressure = 0;
    TPZBndCond * bcPressure = mat->CreateBC(mat, fHDivWrapperMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcPressure);
    
}



// create the elements domain per domain with approximation spaces disconnected from each other
void TPZMHMixedHybridMeshControl::CreateInternalElements()
{
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    gmesh->ResetReference();
    cmeshHDiv->LoadReferences();
    cmeshHDiv->SetDimModel(meshdim);
    cmeshHDiv->ApproxSpace().SetAllCreateFunctionsHDiv(meshdim);
    cmeshHDiv->SetDefaultOrder(fpOrderInternal);
    
    //Criar elementos computacionais malha MHM
    
    TPZGeoEl *gel = NULL;
    TPZGeoEl *gsubel = NULL;
    fConnectToSubDomainIdentifier.Expand(10000);
    
    // create the elements contained in the geometric domains associated with the coarse mesh
    for (std::map<long,long>::iterator it = fCoarseIndices.begin(); it != fCoarseIndices.end(); it++)
    {
        std::set<long> elset;
        bool LagrangeCreated = false;
        long iel = it->first;
        elset.insert(iel);
        //        std::map<long, std::pair<long,long> >::iterator it2;
        //        for (it2 = fInterfaces.begin(); it2 != fInterfaces.end(); it2++) {
        //            if (it2->first == it2->second.second && it2->second.first == *it) {
        //                elset.insert(it2->first);
        //            }
        //        }
        
        while (elset.size()) {
            std::set<long>::iterator itel = elset.begin();
            long elfirst = *itel;
            elset.erase(itel);
            gel = fGMesh->ElementVec()[elfirst];
            if(!gel) DebugStop();
            // we work with only the leaf elements
            if (! gel->HasSubElement())
            {
                long index;
                // create the flux element
                
                TPZCompEl *cel = cmeshHDiv->CreateCompEl(gel, index);
                TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
                int nsides = gel->NSides();
                for (int side=0; side<nsides; side++) {
                    if (gel->SideDimension(side) == meshdim-1) {
                        intel->SetSideOrient(side, 1);
                    }
                }
                // we create disconnected elements
                cel->Reference()->ResetReference();
                /// associate the connects with the subdomain
                SetSubdomain(cel, it->first);
            }
            else
            {
                int nsubels = gel->NSubElements();
                for (int is=0; is<nsubels; is++)
                {
                    gsubel = gel->SubElement(is);
                    elset.insert(gsubel->Index());
                }
            }
        }
        
        fGMesh->ResetReference();
        
    }
    
    fFluxMesh->ExpandSolution();
    
    CreateHDivWrappers();
}

// create the dim-1 pressure elements between the hdiv elements
void TPZMHMixedHybridMeshControl::CreatePressureInterfaces()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        long nel = fFluxMesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    long nel = fGMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fInternalWrapMatId || gel->HasSubElement()) {
            continue;
        }
        // now we have a geometric element of material id fInternalWrapMatId that doesn't have children
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) DebugStop();
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel) {
                TPZConnect &midsideconnect = intel->MidSideConnect(side);
                int order = midsideconnect.Order();
                fPressureFineMesh->SetDefaultOrder(order);
                TPZGeoElBC gelbc(gelside, fPressureSkeletonMatId);
                long index;
                TPZCompEl *sidecel = fPressureFineMesh->CreateCompEl(gelbc.CreatedElement(), index);
                int nc = sidecel->NConnects();
                for (int ic=0; ic<nc; ic++) {
                    sidecel->Connect(ic).SetLagrangeMultiplier(3);
                }
                gelbc.CreatedElement()->ResetReference();
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        if (neighbour == gelside) {
            DebugStop();
        }
    }
    fPressureFineMesh->ExpandSolution();
}

// create the approximation space associated with the skeleton and restrain the connects
void TPZMHMixedHybridMeshControl::CreateSkeleton()
{
    CreatePressureInterfaces();
//    fFluxMesh->LoadReferences();
    TPZMHMixedMeshControl::CreateSkeleton();
}

// create a boundary element around each flux element to optimize calculations
void TPZMHMixedHybridMeshControl::CreateHDivWrappers()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        long nel = fFluxMesh->NElements();
        for (long el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    long nel = fGMesh->NElements();
    for (long el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fInternalWrapMatId) {
            continue;
        }
        // now we have a geometric element of material id fInternalWrapMatId that doesn't have children
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) DebugStop();
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel) {
                intel->LoadElementReference();
                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
                int order = midsideconnect.Order();
                fFluxMesh->SetDefaultOrder(order);
                TPZGeoElBC gelbc(gelside, fHDivWrapperMatId);
                long index;
                TPZCompEl *sidecel = fFluxMesh->CreateCompEl(gelbc.CreatedElement(), index);
                // all the elements have normal pointing outwards
                intel->SetSideOrient(neighbour.Side(), 1);
                TPZInterpolatedElement *bound = dynamic_cast<TPZInterpolatedElement *>(sidecel);
                bound->SetSideOrient(side, 1);
                gelbc.CreatedElement()->ResetReference();
                neighbour.Element()->ResetReference();
            }
            neighbour = neighbour.Neighbour();
        }
    }
}

// Find the two connected flux elements to a pressure element
void TPZMHMixedHybridMeshControl::FindConnectedElements(TPZGeoElSide pressureside, TPZVec<TPZCompElSide> &fluxconnected)
{
    int meshdim = fGMesh->Dimension();
    TPZGeoElSide gelpside = pressureside;
    // find the inner skin element
    TPZGeoElSide skinside;
    TPZGeoElSide neighbour = gelpside.Neighbour();
    TPZStack<TPZCompElSide> fluxcon;
    while (neighbour != gelpside) {
        if (neighbour.Element()->MaterialId() == fHDivWrapperMatId) {
            TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(neighbour.Reference().Element());
            if(!intel) DebugStop();
            fluxcon.Push(neighbour.Reference());
        }
        else if(neighbour.Element()->MaterialId() == fInternalWrapMatId)
        {
            skinside = neighbour;
        }
        else if(neighbour.Element()->MaterialId() == fSkeletonWrapMatId)
        {
            return;
        }
        neighbour = neighbour.Neighbour();
    }
    if(!skinside) DebugStop();
    skinside = skinside.Father2();
    while (skinside) {
        neighbour = skinside.Neighbour();
        while (neighbour != skinside) {
            if (neighbour.Element()->MaterialId() == fHDivWrapperMatId) {
                TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(neighbour.Reference().Element());
                if(!intel) DebugStop();
                fluxcon.Push(neighbour.Reference());
            }
            neighbour = neighbour.Neighbour();
        }
        skinside = skinside.Father2();
    }
    if (fluxcon.size() != 2) {
        DebugStop();
    }
    fluxconnected = fluxcon;
}


/// Create the interfaces between the pressure elements of dimension dim
void TPZMHMixedHybridMeshControl::CreateMultiPhysicsInterfaceElements(int dim)
{
    int meshdim = fGMesh->Dimension();
    // this will create the interface elements for the skeleton elements
    if(dim == meshdim-1)
    {
        TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(dim);
    }
    fGMesh->ResetReference();
    fCMesh->LoadReferences();
    long nel = fCMesh->NElements();
    
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mphys) {
            continue;
        }
        TPZGeoEl *gel = mphys->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        if (mphys->Reference()->MaterialId() == fPressureSkeletonMatId) {
            TPZGeoElSide gelpressureside(gel,gel->NSides()-1);
            TPZManVector<TPZCompElSide,4> connected;
            FindConnectedElements(gelpressureside, connected);
            int ncon = connected.size();
            for (int icon=0; icon<ncon; icon++) {
                long index;
                /** @brief Constructor */
                //                TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, long &index, TPZCompElSide left, TPZCompElSide right);

                TPZGeoElBC gbc(gelpressureside,fLagrangeMatIdRight);
                new TPZMultiphysicsInterfaceElement(fCMesh,gbc.CreatedElement(),index,gelpressureside.Reference(),connected[icon]);
            }
        }
    }
}

/// group and condense the elements
void TPZMHMixedHybridMeshControl::GroupandCondenseElements()
{
    for (std::map<long,long>::iterator it=fCoarseIndices.begin(); it != fCoarseIndices.end(); it++) {
        TPZCompEl *cel = fCMesh->Element(it->second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        subcmesh->ComputeNodElCon();
        
        GroupElements(subcmesh);
        subcmesh->ComputeNodElCon();
        
        bool keeplagrange = false;
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, keeplagrange);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 0;
        int preconditioned = 0;
        TPZAutoPointer<TPZGuiInterface> guiInterface;
        
        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
    }
    
}

/// group element H(div) elements with surrounding interface elements
void TPZMHMixedHybridMeshControl::GroupElements(TPZCompMesh *cmesh)
{
    
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int dim = cmesh->Dimension();
    long nel = cmesh->NElements();
    // all elements that have been put in a group
    std::set<long> grouped;
    for (long el=0; el<nel; el++) {
        if (grouped.find(el) != grouped.end()) {
            continue;
        }
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZMultiphysicsElement *mphel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mphel) {
            DebugStop();
        }
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(mphel->Element(0));
        if(!intel)
        {
            DebugStop();
        }
        std::set<long> elgroup;
        
        std::set<long> connectlist;
        cel->BuildConnectList(connectlist);
        //        std::cout << "cel " << cel->Index() << " connects ";
        //        for (std::set<long>::iterator it=connectlist.begin(); it != connectlist.end(); it++) {
        //            std::cout << *it << " ";
        //        }
        //        std::cout << std::endl;
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            //            std::cout << "side " << is << std::endl;
            TPZGeoElSide gelside(gel,is);
            if (gelside.Dimension() != dim-1) {
                continue;
            }
            // we assume that we are working on H(div) elements
            if (intel->NSideConnects(is) != 1) {
                DebugStop();
            }
            int clocindex = intel->SideConnectLocId(0, is);
            
            long cindex = mphel->ConnectIndex(clocindex);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZMultiphysicsInterfaceElement *mphys = dynamic_cast<TPZMultiphysicsInterfaceElement *>(neighbour.Element()->Reference());
                if (mphys) {
                    TPZCompEl *right = mphys->RightElement();
                    TPZCompEl *left = mphys->LeftElement();
                    if (right->NConnects() == 1 && right->ConnectIndex(0)== cindex) {
                        elgroup.insert(mphys->Index());
                    }
                    if (left->NConnects() == 1 && left->ConnectIndex(0)== cindex) {
                        elgroup.insert(mphys->Index());
                    }
                }
                neighbour = neighbour.Neighbour();
            }
        }
        for (int is=0; is<ns; is++) {
            //            std::cout << "side " << is << std::endl;
            TPZGeoElSide gelside(gel,is);
            TPZStack<TPZCompElSide> celstack;
            gelside.ConnectedCompElementList(celstack, 0, 0);
            long nelstack = celstack.size();
            for (long elst=0; elst<nelstack; elst++) {
                TPZCompElSide celst=celstack[elst];
                //                TPZGeoElSide gelst =celst.Reference();
                TPZCompEl *celsidelement = celst.Element();
                if (grouped.find(celsidelement->Index()) != grouped.end()) {
                    continue;
                }
                std::set<long> smallset;
                celsidelement->BuildConnectList(smallset);
                //                std::cout << "neigh " << celsidelement->Index() << " connects ";
                //                for (std::set<long>::iterator it=smallset.begin(); it != smallset.end(); it++) {
                //                    std::cout << *it << " ";
                //                }
                //                std::cout << std::endl;
                if (std::includes(connectlist.begin(), connectlist.end(), smallset.begin(), smallset.end()))
                {
                    //                    std::cout << "Is included\n";
                    elgroup.insert(celsidelement->Index());
                }
            }
        }
        if (elgroup.size()) {
            elgroup.insert(el);
            long grindex;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,grindex);
            for (std::set<long>::iterator it = elgroup.begin(); it != elgroup.end(); it++) {
                elgr->AddElement(cmesh->Element(*it));
                grouped.insert(*it);
            }
            grouped.insert(grindex);
        }
    }
    cmesh->ComputeNodElCon();

}

