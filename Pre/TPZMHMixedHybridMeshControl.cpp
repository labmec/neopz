//
//  TPZMHMixedHybridMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMixedHybridMeshControl.h"
#include "TPZBndCond.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"


#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzelementgroup.h"

#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"

#include "TPZVTKGeoMesh.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mhmixedhybridmeshcontrol");
#endif

using namespace std;
#include <algorithm>

/*
TPZMHMixedHybridMeshControl::TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, std::set<int64_t> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices)
{
}
*/

TPZMHMixedHybridMeshControl::TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh, TPZVec<int64_t> &coarseindices) : TPZMHMixedMeshControl(gmesh,coarseindices)
{
}

TPZMHMixedHybridMeshControl::TPZMHMixedHybridMeshControl(TPZAutoPointer<TPZGeoMesh> gmesh) : TPZMHMixedMeshControl(gmesh)
{
}

TPZMHMixedHybridMeshControl::~TPZMHMixedHybridMeshControl()
{
    
}




// create the elements domain per domain with approximation spaces disconnected from each other
void TPZMHMixedHybridMeshControl::CreateInternalFluxElements()
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
    fConnectToSubDomainIdentifier[cmeshHDiv].Expand(10000);
    int64_t nel = fGMesh->NElements();
    
    for (auto it= this->fMHMtoSubCMesh.begin(); it != this->fMHMtoSubCMesh.end(); it++)
    {
        fGMesh->ResetReference();
        int64_t MHMIndex = it->first;
        for (int64_t el = 0; el< nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel) continue;
            int geldim = gel->Dimension();
            int64_t gelMHM = fGeoToMHMDomain[el];
            if(!gel || gel->Dimension() != fGMesh->Dimension() || gel->HasSubElement() || gelMHM != MHMIndex)
            {
                continue;
            }
            // create the flux element
            
            TPZCompEl *cel = cmeshHDiv->CreateCompEl(gel);
            /// associate the connects with the subdomain
            SetSubdomain(cel, it->first);
        }
    }
    
    fGMesh->ResetReference();    
    fFluxMesh->ExpandSolution();
    
    SplitFluxElementsAroundFractures();

    CreateHDivWrappers();
}

/// Create all data structures for the computational mesh
void TPZMHMixedHybridMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    CreateHDivMHMMesh();
    
    // hybridize the fluxes between macro elements
    InsertPeriferalPressureMaterialObjects();
    /// hybridize the flow objects on the skeleton in all cases
    for (auto it = fSkeletonWithFlowMatId.begin(); it != fSkeletonWithFlowMatId.end(); it++) {
        HybridizeSkeleton(*it, fSkeletonWithFlowPressureMatId);
    }
    CreateSkeletonAxialFluxes();
    if(fHybridize)
    {
        /// hybridize the normal fluxes
        TPZMHMixedMeshControl::HybridizeSkeleton(fSkeletonMatId,fPressureSkeletonMatId);
    }
    
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    // create the discontinuous pressure mesh
    CreatePressureMHMMesh();
    
    CreateInternalAxialFluxes();
#ifdef PZDEBUG
    if(1)
    {
        ofstream out("cmeshflux.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fFluxMesh.operator->(), out);
    }
#endif
#ifdef PZDEBUG
    if(1)
    {
        ofstream outp("cmeshpres.vtk");
        TPZVTKGeoMesh::PrintCMeshVTK(fPressureFineMesh.operator->(), outp);
    }
#endif

    AdjustBoundaryConditionsOfFractures();
    fPressureFineMesh->ExpandSolution();
    fFluxMesh->ExpandSolution();
    

    InsertPeriferalMaterialObjects();
    
    /// create the multiphysics mesh
    CreateHDivPressureMHMMesh();


    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
    fCMesh->ComputeNodElCon();
#ifdef PZDEBUG
    if(1)
    {
        fFluxMesh->ComputeNodElCon();
        std::ofstream outfl("FluxMesh.txt");
        fFluxMesh->Print(outfl);
        fPressureFineMesh->ComputeNodElCon();
        std::ofstream outpr("PressureMesh.txt");
        fPressureFineMesh->Print(outpr);
        std::ofstream outmf("MFMesh.txt");
        fCMesh->Print(outmf);
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif
    
    if (usersubstructure) {
        HideTheElements();
    }
    fNumeq = fCMesh->NEquations();

}



/// hybridize the flux elements - each flux element becomes 5 elements
void TPZMHMixedHybridMeshControl::HybridizeSkeleton(int skeletonmatid, int pressurematid)
{
    TPZMHMixedMeshControl::HybridizeSkeleton(skeletonmatid,pressurematid);
    fFluxMesh->ExpandSolution();
    fPressureFineMesh->ExpandSolution();
}


// create the dim-1 pressure elements between the hdiv elements
void TPZMHMixedHybridMeshControl::CreatePressureInterfaces()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fFluxMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if(!cel) continue;
            TPZGeoEl *gel = cel->Reference();
            if (gel) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    int64_t nel = fGMesh->NElements();
#ifdef PZDEBUG
    // look for neighbouring fracture elements
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int matid = gel->MaterialId();
        bool found = fFractureFlowDim1MatId.find(matid) != fFractureFlowDim1MatId.end();
        if(!found) continue;
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while (neighbour != gelside) {
            int neighmatid = neighbour.Element()->MaterialId();
            if (fFractureFlowDim1MatId.find(neighmatid) != fFractureFlowDim1MatId.end()) {
                DebugStop();
            }
            neighbour = neighbour.Neighbour();
        }
    }
#endif
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        int matid = gel->MaterialId();
        bool found = fFractureFlowDim1MatId.find(matid) != fFractureFlowDim1MatId.end();
        if (!gel || !found || gel->HasSubElement()) {
            continue;
        }
        
        int subdomain = this->fGeoToMHMDomain[el];
        if (subdomain == -1) {
            DebugStop();
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            int side = neighbour.Side();
            int subdomain2 = WhichSubdomain(intel);
#ifdef PZDEBUG
            if (subdomain != subdomain2) {
                DebugStop();
            }
#endif
            if (intel) {
                TPZGeoElBC pressure_gelbc(gelside,fPressureDim1MatId);
                TPZConnect &midsideconnect = intel->MidSideConnect(side);
                int order = midsideconnect.Order();
                fPressureFineMesh->SetDefaultOrder(order);
                TPZCompEl *sidecel = fPressureFineMesh->CreateCompEl(pressure_gelbc.CreatedElement());
#ifdef PZ_LOG
                if(logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Created a pressure interface element " << sidecel->Index() << " with gel index " << pressure_gelbc.CreatedElement()->Index()
                     << " and material id " << pressure_gelbc.CreatedElement()->MaterialId();
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
                SetSubdomain(sidecel, subdomain);
                int nc = sidecel->NConnects();
                for (int ic=0; ic<nc; ic++) {
                    sidecel->Connect(ic).SetLagrangeMultiplier(3);
                }
                pressure_gelbc.CreatedElement()->ResetReference();
                // only one element will be created on a side
                break;
            }
            neighbour = neighbour.Neighbour();
        }
        // if we did not find a flux element stop the code!
        if (neighbour == gelside) {
            DebugStop();
        }
    }
    fPressureFineMesh->ExpandSolution();
}

/// Create the pressure mesh which is dual to the flux mesh
void TPZMHMixedHybridMeshControl::CreatePressureMHMMesh()
{
    TPZMHMixedMeshControl::CreatePressureMHMMesh();
    // will create a pressure element associated with internal interfaces only
    // the fluxes associated with the skeleton interfaces will be restrained
    CreatePressureInterfaces();
    CreateLowerDimensionPressureElements();
    fPressureFineMesh->ExpandSolution();
}



// create the approximation space associated with the skeleton and restrain the connects
void TPZMHMixedHybridMeshControl::CreateSkeleton()
{
    TPZMHMixedMeshControl::CreateSkeleton();
}

/// verifies if a HDiv wrapper needs to be created for a given element/side
// the method will look if there is a neighbouring geometric element with material id fPressureDim1MatId
bool TPZMHMixedHybridMeshControl::NeedsHDivWrapper(TPZGeoElSide gelside)
{
    TPZGeoElSide wrap = gelside.Neighbour();
    while(wrap != gelside && wrap.Element()->MaterialId() != fInternalWrapMatId)
    {
        wrap = wrap.Neighbour();
    }
    if (wrap.Element()->MaterialId() != fInternalWrapMatId) {
        DebugStop();
    }
    // we go all the way up
    while (wrap.Father2()) {
        wrap = wrap.Father2();
    }
    // now we look for a neighbour one level at a time
    TPZStack<TPZGeoElSide> elstack;
    elstack.Push(wrap);
    while (elstack.size()) {
        TPZGeoElSide search = elstack.Pop();
        TPZGeoElSide neighbour = search.Neighbour();
        while (neighbour != search) {
            if (fFractureFlowDim1MatId.find(neighbour.Element()->MaterialId()) != fFractureFlowDim1MatId.end()) {
                return true;
            }
            neighbour = neighbour.Neighbour();
        }
        if (search.HasSubElement()) {
            int nsub = search.Element()->NSubElements();
            for (int isub = 0; isub < nsub; isub++)
            {
                TPZGeoEl *subel = search.Element()->SubElement(isub);
                elstack.Push(TPZGeoElSide(subel,subel->NSides()-1));
            }
        }
    }
    return false;
}

/// split the fluxes between the flux elements adjacent to a fracture
void TPZMHMixedHybridMeshControl::SplitFluxElementsAroundFractures()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fFluxMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fInternalWrapMatId) {
            continue;
        }
        // now we have a geometric element of material id fInternalWrapMatId
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) continue;
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        if (!NeedsHDivWrapper(gelside)) {
            continue;
        }
        // create the HDiv wrapper element for all connected flux elements
        TPZGeoElSide neighbour = gelside.Neighbour();
        int first = 0;
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel && first == 0) {
                intel->LoadElementReference();
#ifdef PZDEBUG
                int midsideconnectlocid = intel->MidSideConnectLocId(neighbour.Side());
                int cindex = intel->ConnectIndex(midsideconnectlocid);
                // check whether a wrapper with this connect already exists
                // if affirmative - DebugStop
#endif
                int locindex = intel->MidSideConnectLocId(neighbour.Side());
                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
                if(midsideconnect.NElConnected() != 2)
                {
                    DebugStop();
                }
                // create a new connect
                // associate with the element
                int64_t index = fFluxMesh->AllocateNewConnect(midsideconnect.NShape(), midsideconnect.NState(), midsideconnect.Order());
                int domain = WhichSubdomain(intel);
                SetSubdomain(fFluxMesh.operator->(), index, domain);
                intel->SetConnectIndex(locindex, index);
                midsideconnect.DecrementElConnected();
                fFluxMesh->ConnectVec()[index].IncrementElConnected();
                intel->SetSideOrient(neighbour.Side(), 1);
                first++;
            }
            else if(intel && first == 1)
            {
                intel->SetSideOrient(neighbour.Side(), 1);
                first++;
            }
            else if(intel)
            {
                DebugStop();
            }
            neighbour = neighbour.Neighbour();
        }
    }
}



// create a boundary element around each flux element of the mesh dimension to optimize calculations
void TPZMHMixedHybridMeshControl::CreateHDivWrappers()
{
    fGMesh->ResetReference();
    int meshdim = fGMesh->Dimension();
    // build a vector which for each geoel index points to the flux element
    TPZVec<TPZInterpolatedElement *> fluxelements(fGMesh->NElements(),0);
    {
        int64_t nel = fFluxMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fFluxMesh->Element(el);
            if (cel && cel->Reference()) {
                fluxelements[cel->Reference()->Index()] = dynamic_cast<TPZInterpolatedElement *>(cel);
            }
        }
    }
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel || gel->MaterialId() != fInternalWrapMatId) {
            continue;
        }
        // now we have a geometric element of material id fInternalWrapMatId
        int nsides = gel->NSides();
        int side = nsides-1;
        if(gel->SideDimension(side) != meshdim-1) continue;
        // look for a neighbouring flux element, determine its side order
        TPZGeoElSide gelside(gel,side);
        if (!NeedsHDivWrapper(gelside)) {
            continue;
        }
        // create the HDiv wrapper element for all connected flux elements
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            TPZInterpolatedElement *intel = fluxelements[neighbour.Element()->Index()];
            if (intel) {
                intel->LoadElementReference();
#ifdef PZDEBUG
                int midsideconnectlocid = intel->MidSideConnectLocId(neighbour.Side());
                int cindex = intel->ConnectIndex(midsideconnectlocid);
                // check whether a wrapper with this connect already exists
                // if affirmative - DebugStop
#endif
                TPZConnect &midsideconnect = intel->MidSideConnect(neighbour.Side());
                int order = midsideconnect.Order();
                fFluxMesh->SetDefaultOrder(order);
                TPZGeoElBC gelbc(gelside, fHDivWrapperMatId);
                TPZCompEl *sidecel = fFluxMesh->CreateCompEl(gelbc.CreatedElement());
                // all the elements have normal pointing outwards
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Creating cap on element " << intel->Index() << " side " << neighbour.Side() << " element created " << sidecel->Index();
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
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

// Find the connected flux elements to a pressure element
void TPZMHMixedHybridMeshControl::FindConnectedElements(TPZGeoElSide &pressureside, int domain, TPZVec<TPZCompElSide> &fluxconnected)
{
    // the multiphysics mesh has to be loaded
    int meshdim = fGMesh->Dimension();
    int pressuredim = pressureside.Dimension();
    TPZGeoElSide gelpside = pressureside;
    // find the inner skin element
    TPZGeoElSide skinside;
    TPZGeoElSide neighbour = gelpside.Neighbour();
    TPZStack<TPZCompElSide> fluxcon;
    while (neighbour != gelpside) {
        if (neighbour.Element()->MaterialId() == fHDivWrapperMatId && neighbour.Element()->Dimension() == pressuredim) {
            TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(neighbour.Reference().Element());
            if(!intel) DebugStop();
            TPZCompEl *fluxel = intel->Element(0);
#ifdef PZDEBUG
            if (!fluxel) {
                DebugStop();
            }
#endif
            int fluxdomain = WhichSubdomain(fluxel);
            if (fluxdomain == domain || domain == -1) {
                fluxcon.Push(neighbour.Reference());
            }
        }
        else if(neighbour.Element()->MaterialId() == fInternalWrapMatId)
        {
            skinside = neighbour;
        }
        else if(neighbour.Element()->MaterialId() == fSkeletonWrapMatId && pressuredim == meshdim-1)
        {
            return;
        }
        else if(neighbour.Element()->MaterialId() == fSkeletonWrapMatId)
        {
            skinside = neighbour;
        }
        neighbour = neighbour.Neighbour();
    }
    if(!skinside) DebugStop();
    skinside = skinside.Father2();
    while (skinside) {
        neighbour = skinside.Neighbour();
        while (neighbour != skinside) {
            if (neighbour.Element()->MaterialId() == fHDivWrapperMatId && neighbour.Element()->Dimension()== pressuredim) {
                TPZMultiphysicsElement *intel = dynamic_cast<TPZMultiphysicsElement *>(neighbour.Reference().Element());
                if(!intel) DebugStop();
                TPZCompEl *fluxel = intel->Element(0);
#ifdef PZDEBUG
                if (!fluxel) {
                    DebugStop();
                }
#endif
                int fluxdomain = WhichSubdomain(fluxel);
                if (domain == -1 || fluxdomain == domain) {
                    fluxcon.Push(neighbour.Reference());
                }
            }
            neighbour = neighbour.Neighbour();
        }
        skinside = skinside.Father2();
    }
    if (fluxcon.size() != 2 && pressuredim == meshdim-1) {
        DebugStop();
    }
    if (fluxcon.size() == 0) {
        TPZManVector<REAL,3> co(3);
        pressureside.Element()->NodePtr(0)->GetCoordinates(co);
        std::cout << "Element dimension " << pressureside.Element()->Dimension() << " node 0 co " << co << std::endl;
        std::cout << "Neighbouring material ids ";
        TPZGeoElSide neighbour = pressureside.Neighbour();
        while (neighbour != pressureside) {
            std::cout << neighbour.Element()->MaterialId() << " ";
            neighbour = neighbour.Neighbour();
        }
        std::cout << std::endl;
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
        std::pair<int,int> skelmatid(fSkeletonMatId,fSecondSkeletonMatId);
        TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(dim,fPressureSkeletonMatId, skelmatid);
        std::pair<int,int> skelflowmatid(fSkeletonMatId,fSecondSkeletonMatId);
        TPZMHMixedMeshControl::CreateMultiPhysicsInterfaceElements(dim,fSkeletonWithFlowPressureMatId, skelflowmatid);
    }
    fGMesh->ResetReference();
//    fCMesh->LoadReferences();
    
    // Reference only elements of dimension dim
    // otherwise the function that returns connected element list will return elements with wrong dimension
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (gel && gel->Dimension() == dim) {
            cel->LoadElementReference();
        }
    }
    
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        if (!mphys) {
            continue;
        }
        TPZGeoEl *gel = mphys->Reference();
        if (gel->Dimension() != dim) {
            continue;
        }
        // work only with the pressure elements
        if ((mphys->Reference()->MaterialId() == fPressureDim1MatId && dim==meshdim-1) || (mphys->Reference()->MaterialId() == fPressureDim2MatId && dim == meshdim-2)) {
            TPZGeoElSide gelpressureside(gel,gel->NSides()-1);
            // there cannot be a flux element associated with this multiphysics element
            if (mphys->Element(0)) {
                DebugStop();
            }
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Investigating pressure element " << mphys->Index() << " dim " << dim;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            TPZCompEl *presel = mphys->Element(1);
            int domain = WhichSubdomain(presel);
            // look for all connected flux elements
            TPZManVector<TPZCompElSide,4> connected;
            FindConnectedElements(gelpressureside, domain, connected);
            int ncon = connected.size();
            for (int icon=0; icon<ncon; icon++) {
                /** @brief Constructor */
                //                TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);

                TPZCompElSide HDivWrap = connected[icon];
                TPZMultiphysicsElement *mphyscon = dynamic_cast<TPZMultiphysicsElement *>(HDivWrap.Element());
                if (!mphyscon) {
                    DebugStop();
                }
                // there can be no pressure element associate with the right element
                if (mphyscon->Element(1)) {
                    DebugStop();
                }
                TPZGeoElBC gbc(gelpressureside,fLagrangeMatIdRight);
                // what if the pressure skeleton element has a flux associated with it?
                
                TPZCompEl *intface = new TPZMultiphysicsInterfaceElement(fCMesh,gbc.CreatedElement(),gelpressureside.Reference(),connected[icon]);
#ifdef PZ_LOG
                if (logger.isDebugEnabled())
                {
                    std::stringstream sout;
                    sout << "Created an interface element index " << intface->Index() << " with cap mfys index " << mphyscon->Index() << " flux cap element " << mphyscon->Element(0)->Index();
                    LOGPZ_DEBUG(logger, sout.str())
                }
#endif
            }
        }
    }
}

/// group and condense the elements
void TPZMHMixedHybridMeshControl::GroupandCondenseElements()
{
    for (auto it:fMHMtoSubCMesh) {
        TPZCompEl *cel = fCMesh->Element(it.second);
        TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *>(cel);
        if (!subcmesh) {
            DebugStop();
        }
        subcmesh->ComputeNodElCon();
        GroupElements(subcmesh);
        subcmesh->InitializeBlock();
#ifdef PZ_LOG
        if(logger.isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        bool keeplagrange = true;
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, keeplagrange);
        subcmesh->CleanUpUnconnectedNodes();
        int numthreads = 0;
        int preconditioned = 0;
        TPZAutoPointer<TPZGuiInterface> guiInterface;

        subcmesh->SetAnalysisSkyline(numthreads, preconditioned, guiInterface);
    }
    fCMesh->ComputeNodElCon();
    fCMesh->CleanUpUnconnectedNodes();
}

/// group element H(div) elements with surrounding interface elements
void TPZMHMixedHybridMeshControl::GroupElements(TPZCompMesh *cmesh)
{
    
    cmesh->Reference()->ResetReference();
    cmesh->LoadReferences();
    int dim = cmesh->Dimension();
    int64_t nel = cmesh->NElements();
    // all the elements that have been put in a group
    std::set<int64_t> grouped;
    for (int64_t el=0; el<nel; el++) {
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
        std::set<int64_t> elgroup;
        
        std::set<int64_t> connectlist;
        cel->BuildConnectList(connectlist);
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
        {
            std::stringstream sout;
            sout << "cel " << cel->Index() << " connects ";
            for (std::set<int64_t>::iterator it=connectlist.begin(); it != connectlist.end(); it++) {
                sout << *it << " ";
            }
            sout << std::endl;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        int ns = gel->NSides();
        for (int is=gel->NCornerNodes(); is<gel->NSides(); is++) {
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
            
            int64_t cindex = mphel->ConnectIndex(clocindex);
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                TPZMultiphysicsInterfaceElement *mphys = dynamic_cast<TPZMultiphysicsInterfaceElement *>(neighbour.Element()->Reference());
                if (mphys) {
                    TPZCompEl *right = mphys->RightElement();
                    TPZCompEl *left = mphys->LeftElement();
                    if (right->NConnects() == 1 && right->ConnectIndex(0)== cindex) {
#ifdef PZ_LOG
                        if (logger.isDebugEnabled())
                        {
                            std::stringstream sout;
                            sout << "Adding interface " << mphys->Index() << " right connect index " << right->ConnectIndex(0);
                            LOGPZ_DEBUG(logger, sout.str())
                        }
#endif
                        elgroup.insert(mphys->Index());
                    }
                    if (left->NConnects() == 1 && left->ConnectIndex(0)== cindex) {
#ifdef PZ_LOG
                        if (logger.isDebugEnabled())
                        {
                            std::stringstream sout;
                            sout << "Adding interface " << mphys->Index() << " left connect index " << left->ConnectIndex(0);
                            LOGPZ_DEBUG(logger, sout.str())
                        }
#endif
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
            int64_t nelstack = celstack.size();
            for (int64_t elst=0; elst<nelstack; elst++) {
                TPZCompElSide celst=celstack[elst];
                //                TPZGeoElSide gelst =celst.Reference();
                TPZCompEl *celsidelement = celst.Element();
                if (grouped.find(celsidelement->Index()) != grouped.end()) {
                    continue;
                }
                std::set<int64_t> smallset;
                celsidelement->BuildConnectList(smallset);
                //                std::cout << "neigh " << celsidelement->Index() << " connects ";
                //                for (std::set<int64_t>::iterator it=smallset.begin(); it != smallset.end(); it++) {
                //                    std::cout << *it << " ";
                //                }
                //                std::cout << std::endl;
                if (std::includes(connectlist.begin(), connectlist.end(), smallset.begin(), smallset.end()))
                {
                    //                    std::cout << "Is included\n";
#ifdef PZ_LOG
                    if (logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Adding connected element " << celsidelement->Index() << " connects are ";
                        for (auto it=smallset.begin(); it != smallset.end(); it++) {
                            sout << *it << ' ';
                        }
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    elgroup.insert(celsidelement->Index());
                }
            }
        }
#ifdef PZDEBUG
        {
            for (auto it = elgroup.begin(); it != elgroup.end(); it++) {
                if (cmesh->Element(*it) == 0) {
                    std::cout << "Element " << *it << " does not exist?\n";
                    DebugStop();
                }
            }
        }
#endif
        if (elgroup.size()) {
            elgroup.insert(el);
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh);
            const int64_t grindex = elgr->Index();
            for (std::set<int64_t>::iterator it = elgroup.begin(); it != elgroup.end(); it++) {
                elgr->AddElement(cmesh->Element(*it));
                grouped.insert(*it);
            }
            grouped.insert(grindex);
        }
    }
    cmesh->ComputeNodElCon();

}

#ifdef MACOSX
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wreturn-type"
#endif
static TPZInterpolatedElement *FindPressureSkeleton(TPZGeoEl *InternalWrap, int pressmatid, TPZVec<TPZInterpolatedElement *> pressureelements)
{
    TPZGeoElSide gelside(InternalWrap,InternalWrap->NSides()-1);
    TPZGeoElSide neighbour(gelside.Neighbour());
    while (neighbour != gelside) {
        if (neighbour.Element()->MaterialId() == pressmatid) {
            TPZInterpolatedElement *intel = pressureelements[neighbour.Element()->Index()];
            if (!intel) {
                DebugStop();
            }
            return intel;
        }
        neighbour = neighbour.Neighbour();
    }
    TPZGeoElSide father = gelside.Father2();
    if (father) {
        return FindPressureSkeleton(father.Element(), pressmatid, pressureelements);
    }
    else
    {
        DebugStop();
        //silencing warning
        return nullptr;
    }
}
#ifdef MACOSX
#pragma clang diagnostic pop
#endif

/// Create lower dimension pressure interfaces (dim-2)
/// The (dim-2) geometric elements have already been created.
// when dim==2 create a point pressure element for each skeleton element
void TPZMHMixedHybridMeshControl::CreateLowerDimensionPressureElements()
{
    fGMesh->ResetReference();
    // create a vector that points a geometric element index to a pressure computational element
    TPZManVector<TPZInterpolatedElement *> pressureelements(fGMesh->NElements(),0);
    {
        int64_t nel = fPressureFineMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZCompEl *cel = fPressureFineMesh->Element(el);
            if (!cel) {
                continue;
            }
            TPZInterpolatedElement *press = dynamic_cast<TPZInterpolatedElement *>(cel);
            if(!press) DebugStop();
#ifdef PZDEBUG
            if (!cel->Reference()) {
                DebugStop();
            }
#endif
            TPZGeoEl *gel = cel->Reference();
            pressureelements[gel->Index()] = press;
        }
    }
    int64_t nel = fGMesh->NElements();
    int meshdim = fGMesh->Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        // only interested in elements of dimension dim-2
        if (!gel || gel->Dimension() != meshdim-2) {
            continue;
        }
        // initially we will only create fracture flow aint64_t the skeleton elements
        if (gel->MaterialId() == fSkeletonWrapMatId || gel->MaterialId() == fInternalWrapMatId) {
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmatid = neighbour.Element()->MaterialId();
            while(neighbour != gelside)
            {
                TPZCompEl *cel = pressureelements[neighbour.Element()->Index()];
                if (neighmatid == fSkeletonWithFlowPressureMatId && cel) {
                    break;
                }
                neighbour = neighbour.Neighbour();
                neighmatid = neighbour.Element()->MaterialId();
            }
            if (neighbour == gelside) {
                neighbour = neighbour.Neighbour();
                neighmatid = neighbour.Element()->MaterialId();
                while(neighbour != gelside)
                {
                    TPZCompEl *cel = pressureelements[neighbour.Element()->Index()];
                    if (neighmatid == fPressureDim1MatId && cel) {
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                    neighmatid = neighbour.Element()->MaterialId();
                }
            }
            if (neighbour == gelside) {
                continue;
            }
            TPZCompEl *pressel = pressureelements[neighbour.Element()->Index()];
            if(!pressel)
            {
                DebugStop();
            }
            TPZGeoElBC gelbc(gel,gel->NSides()-1,fPressureDim2MatId);
            TPZGeoEl *newgel = gelbc.CreatedElement();
            int64_t newgelindex = newgel->Index();
            TPZCompEl *cel = fPressureFineMesh->CreateCompEl(newgel);
            TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
            newgel->ResetReference();
#ifdef PZDEBUG
            if (!intel) {
                DebugStop();
            }
#endif
            int nc = cel->NConnects();
            for (int ic=0; ic<nc; ic++) {
                cel->Connect(ic).SetLagrangeMultiplier(3);
            }
            SetSubdomain(cel, -1);
            if (gel->MaterialId() != fSkeletonWrapMatId && neighmatid == fPressureDim1MatId) {
                int64_t neighindex = neighbour.Element()->Index();
                int subdomain = fGeoToMHMDomain[neighindex];
                SetSubdomain(cel, subdomain);
            }
            // do not apply the restraints
            if (0 && gel->MaterialId() == fInternalWrapMatId)
            {
                // if the element is of type internal, it may be neighbouring a skeleton element
                // in this case we need to constrain the lower dimensional pressure element to the skeleton pressure space
                TPZGeoElSide gelside(newgel,newgel->NSides()-1);
                TPZGeoElSide neighbour = gelside.Neighbour();
                while (neighbour != gelside) {
                    if (neighbour.Element()->MaterialId() == fSkeletonWrapMatId) {
                        // find the pressure skeleton element and restrain the pressure element
#ifdef PZDEBUG
                        if (neighbour.Element()->Dimension() != meshdim-1) {
                            DebugStop();
                        }
#endif
                        TPZInterpolatedElement *pressure_skeleton = FindPressureSkeleton(neighbour.Element(), fPressureSkeletonMatId, pressureelements);
                        // find out left and right subdomain
                        int64_t skelindex = pressure_skeleton->Reference()->Index();
#ifdef PZDEBUG
                        if (fInterfaces.find(skelindex) == fInterfaces.end()) {
                            DebugStop();
                        }
#endif
                        int64_t left = fInterfaces[skelindex].first;
                        int64_t right = fInterfaces[skelindex].second;
                        TPZGeoElBC gelbc(gel,gel->NSides()-1,fPressureDim2MatId);
                        TPZGeoEl *newgel = gelbc.CreatedElement();
                        TPZCompEl* mycel = fPressureFineMesh->CreateCompEl(newgel);
                        const int64_t index = mycel->Index();
                        newgel->ResetReference();
                        TPZInterpolatedElement *intel2 = dynamic_cast<TPZInterpolatedElement *>(fPressureFineMesh->Element(index));
                        SetSubdomain(intel, left);
                        SetSubdomain(intel2, right);
                        intel->RestrainSide(newgel->NSides()-1, pressure_skeleton, pressure_skeleton->Reference()->NSides()-1);
                        intel2->RestrainSide(newgel->NSides()-1, pressure_skeleton, pressure_skeleton->Reference()->NSides()-1);
                        break;
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (neighbour == gelside) {
                    // we DIDNT find a neighbouring skeleton element
                    // figure out in which subdomain we are
                    neighbour = gelside.Neighbour();
                    while(neighbour != gelside)
                    {
                        TPZGeoEl *gelneigh = neighbour.Element();
                        int matneigh = gelneigh->MaterialId();
                        if (matneigh == fPressureDim1MatId) {
                            TPZInterpolatedElement *pressel = pressureelements[gelneigh->Index()];
#ifdef PZDEBUG
                            if (!pressel) {
                                DebugStop();
                            }
#endif
                            int domain = WhichSubdomain(pressel);
                            if(domain == -1) DebugStop();
                            SetSubdomain(intel, domain);
                            break;
                        }
                        neighbour = neighbour.Neighbour();
                    }
                    // I don't understand
                    if (neighbour==gelside) {
                        DebugStop();
                    }
                }
            }
        }
    }
}


static bool HasNeighbour(TPZGeoElSide gelside, int materialid)
{
    if (gelside.Element()->MaterialId() == materialid) {
        return true;
    }
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        if (neighbour.Element()->MaterialId() == materialid) {
            return true;
        }
        neighbour = neighbour.Neighbour();
    }
    return false;
}
static int64_t HasNeighbour(TPZGeoElSide gelside, std::set<int> &materialid)
{
    if(materialid.size() == 0) return -1;
    if (materialid.find(gelside.Element()->MaterialId()) != materialid.end()) {
        return gelside.Element()->Index();
    }
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        if (materialid.find(neighbour.Element()->MaterialId()) != materialid.end()) {
            return neighbour.Element()->Index();
        }
        neighbour = neighbour.Neighbour();
    }
    return -1;
}
// create HDiv approximation spaces associated with the pressure skeleton elements
// initially we will create only HDiv elements associated with the coarse skeleton elements
void TPZMHMixedHybridMeshControl::CreateSkeletonAxialFluxes()
{
    fGMesh->ResetReference();
    // create a vector that points a geometric element index to a pressure computational element
    int64_t nel = fPressureFineMesh->NElements();
    TPZManVector<TPZInterpolatedElement *> pressureelements(fGMesh->NElements(),0);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        if (!cel) {
            continue;
        }
        TPZInterpolatedElement *press = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!press) DebugStop();
#ifdef PZDEBUG
        if (!cel->Reference()) {
            DebugStop();
        }
#endif
        TPZGeoEl *gel = cel->Reference();
        pressureelements[gel->Index()] = press;
    }
    /// indices of the skeleton elements and their left/right geometric elements of the skeleton mesh
//    std::map<int64_t, std::pair<int64_t,int64_t> > fInterfaces;
    for (auto it = fInterfaces.begin(); it != fInterfaces.end(); it++) {
        // skip the boundary elements
        if (it->first == it->second.second) {
            continue;
        }
        TPZInterpolatedElement *press = pressureelements[it->first];
        if (!press) {
            continue;
        }
        int pressmatid = press->Reference()->MaterialId();
        TPZGeoEl *gel = fGMesh->Element(it->first);
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        int64_t elindex = -1;
        elindex = HasNeighbour(gelside, fSkeletonWithFlowMatId);
        if (pressmatid == fSkeletonWithFlowPressureMatId && elindex == -1) {
            DebugStop();
        }
        if (elindex == -1) {
            continue;
        }
        TPZGeoEl *gelflux = fGMesh->Element(elindex);
        int gelfluxmatid = gelflux->MaterialId();
#ifdef PZDEBUG
        {
            // if two fracture flow elements would be created, an error in the datastructure occurred
            TPZGeoEl *gel = fGMesh->Element(it->first);
            TPZGeoElSide gelside(gel,gel->NSides()-1);
            if (HasNeighbour(gelside, fFractureFlowDim1MatId) != -1) {
                DebugStop();
            }
        }
#endif
        gelflux->SetMaterialId(fSkeletonMatId);
        CreateAxialFluxElement(press, gelfluxmatid);
    }
}

// create HDiv approximation spaces associated with the pressure elements on internal side
// initially we will create only HDiv elements associated with the coarse skeleton elements
void TPZMHMixedHybridMeshControl::CreateInternalAxialFluxes()
{
    if (fFractureFlowDim1MatId.size() == 0)
    {
        // we dont have any fractures in the mesh
        return;
    }
    fGMesh->ResetReference();
    int64_t nel = fPressureFineMesh->NElements();
    int meshdim = fGMesh->Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(fPressureFineMesh->Element(el));
        if (!intel) {
            continue;
        }
        TPZGeoEl *gel = intel->Reference();
        if (!gel || gel->Dimension() != meshdim-1) {
            continue;
        }
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        // we dont handle skeleton elementos
        if (gel->MaterialId() == fSkeletonWithFlowPressureMatId) {
            continue;
        }
        if (HasNeighbour(gelside, fSkeletonWithFlowMatId) != -1) {
            continue;
        }
        int64_t gelfracindex = HasNeighbour(gelside, fFractureFlowDim1MatId);
        if(gelfracindex == -1)
        {
            DebugStop();
        }
        int gelfracmatid = fGMesh->Element(gelfracindex)->MaterialId();
        CreateAxialFluxElement(intel,gelfracmatid);
    }
}


/// Add axial flux to a pressure element
// the material id of the pressure element can be either fPressureSkeletonMatId or fPressureDim1MatId
// An H(div) and a pressure element will be created with material Id fFractureFlowDim1MatId
// HDivWrapper boundary elements will also be created
void TPZMHMixedHybridMeshControl::CreateAxialFluxElement(TPZInterpolatedElement *PressureElement, int gelfluxmatid)
{
    TPZGeoEl *gel = PressureElement->Reference();
    if (!gel || gel->Dimension() != fGMesh->Dimension()-1) {
        DebugStop();
    }
    if (gel->MaterialId() != fSkeletonWithFlowPressureMatId && gel->MaterialId() != fPressureDim1MatId &&
        gel->MaterialId() != fPressureSkeletonMatId) {
        DebugStop();
    }
    int subdomain = WhichSubdomain(PressureElement);
    int meshdim = fGMesh->Dimension();
    // duplicate the pressure element
    PressureElement->LoadElementReference();
    fPressureFineMesh->SetAllCreateFunctions(*PressureElement);
    int order = PressureElement->PreferredSideOrder(gel->NSides()-1);
    fPressureFineMesh->SetDefaultOrder(order);
    TPZGeoElBC gelfracbc(gel,gel->NSides()-1,gelfluxmatid);
    int64_t gelfracindex = gelfracbc.CreatedElement()->Index();
    TPZGeoEl *gelfrac = fGMesh->Element(gelfracindex);
    TPZCompEl* loccel = fPressureFineMesh->CreateCompEl(gelfrac);
    const int64_t indexnewpressure = loccel->Index();
    TPZCompEl *presclone = fPressureFineMesh->Element(indexnewpressure);
    // a pressure element is created using the same connects
#ifdef PZDEBUG
    int nc = PressureElement->NConnects();
    for (int ic=0; ic<nc; ic++) {
        if(PressureElement->ConnectIndex(ic) != presclone->ConnectIndex(ic))
        {
            DebugStop();
        }
    }
#endif
    gel->ResetReference();
    presclone->Reference()->ResetReference();
    
    // Create a Flux Element and its Hdiv Wrapper Elements
    fFluxMesh->SetDefaultOrder(order);
    fFluxMesh->SetDimModel(meshdim-1);
    fFluxMesh->SetAllCreateFunctionsHDiv();
    fFluxMesh->SetDefaultOrder(order);
    TPZCompEl* mycel = fFluxMesh->CreateCompEl(gelfrac);
    TPZInterpolatedElement *fluxcel = dynamic_cast<TPZInterpolatedElement *>(mycel);
    SetSubdomain(fluxcel, subdomain);
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        std::stringstream sout;
        sout << "Created a " << meshdim-1 << " dimensional flux element on interface index " << gel->Index() << " and material id " << gelfrac->MaterialId();
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    int nsides = gel->NSides();
    for (int is = 0; is<nsides; is++) {
        if (gel->SideDimension(is) != meshdim-2) {
            continue;
        }
        fluxcel->SetSideOrient(is, 1);
        TPZGeoElBC gelcapbc(gel,is,fHDivWrapperMatId);
        TPZGeoEl *gelcap = gelcapbc.CreatedElement();
        TPZCompEl *cel = fFluxMesh->CreateCompEl(gelcap);
        TPZInterpolatedElement *celcap = dynamic_cast<TPZInterpolatedElement *>(cel);
        celcap->SetSideOrient(gelcap->NSides()-1, 1);
        if (!celcap) {
            DebugStop();
        }
        if (fluxcel->SideConnectIndex(0, is) != celcap->ConnectIndex(0)) {
            DebugStop();
        }
#ifdef PZ_LOG
        if (logger.isDebugEnabled()) {
            std::stringstream sout;
            sout << "Created a " << meshdim-2 << "-dimensional cap with fluxmesh index " << celcap->Index() << " connect index " << celcap->ConnectIndex(0);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        celcap->Reference()->ResetReference();
    }
    gelfrac->ResetReference();
    fFluxMesh->SetDimModel(meshdim);
}



/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedHybridMeshControl::InsertPeriferalMaterialObjects()
{
    TPZMHMixedMeshControl::InsertPeriferalMaterialObjects();
    
    int matid = *fMaterialIds.begin();
    auto *mat =
        dynamic_cast<TPZMaterialT<STATE> *>(fCMesh->FindMaterial(matid));
    if (!mat) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.);
    TPZManVector<STATE,1> val2Flux(1,0.);
    int typeFlux = 1;
    int typePressure = 0;
    {
        TPZBndCond * bcPressure = mat->CreateBC(mat, fPressureSkeletonMatId, typePressure, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcPressure);
    }
    {
        TPZBndCond * bcSecondFlux = mat->CreateBC(mat, fSecondSkeletonMatId, typePressure, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcSecondFlux);
    }
    TPZBndCond * bcPressure = mat->CreateBC(mat, fHDivWrapperMatId, typePressure, val1, val2Flux);
    //    bcN->SetForcingFunction(0,force);
    fCMesh->InsertMaterialObject(bcPressure);


    if (fFractureFlowDim1MatId.size() != 0)
    {
        bcPressure = mat->CreateBC(mat, fPressureDim1MatId, typePressure, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcPressure);
        
        bcPressure = mat->CreateBC(mat, fSkeletonWithFlowPressureMatId, typePressure, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcPressure);
        
        bcPressure = mat->CreateBC(mat, fPressureDim2MatId, typePressure, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcPressure);

        auto *matflux =
        dynamic_cast<TPZMaterialT<STATE> *>(fCMesh->FindMaterial(*fFractureFlowDim1MatId.begin()));
        if (!matflux) {
            DebugStop();
        }
        bcPressure = matflux->CreateBC(matflux, fHomogeneousNeumannBcMatId, typeFlux, val1, val2Flux);
        //    bcN->SetForcingFunction(0,force);
        fCMesh->InsertMaterialObject(bcPressure);
    }

    int nstate = 1;
    int dim = fGMesh->Dimension();
    auto *matleft = new TPZLagrangeMultiplier<STATE>(fLagrangeMatIdLeft,dim,nstate);
    auto *matright = new TPZLagrangeMultiplier<STATE>(fLagrangeMatIdRight,dim,nstate);
    
    fCMesh->InsertMaterialObject(matleft);
    fCMesh->InsertMaterialObject(matright);

}


/// Insert the necessary H(div) material objects to create the flux mesh
void TPZMHMixedHybridMeshControl::InsertPeriferalHdivMaterialObjects()
{
    TPZMHMixedMeshControl::InsertPeriferalHdivMaterialObjects();
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int meshdim = gmesh->Dimension();
    TPZCompMesh * cmeshHDiv = fFluxMesh.operator->();
    TPZNullMaterial<STATE> *matl2 = 0;
    for(auto iter:cmeshHDiv->MaterialVec())
    {
        TPZMaterial *mat = iter.second;
        matl2 = dynamic_cast<TPZNullMaterial<STATE> *>(mat);
        if(matl2) break;
    }
    if (!matl2) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(1,1,0.);
    TPZManVector<STATE,1> val2(1,0.);
    TPZBndCond *bc;
    // totototo
    bc = matl2->CreateBC(matl2, fHDivWrapperMatId, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    

    int matid = *fFractureFlowDim1MatId.begin();
    auto *fracdim1 =
        dynamic_cast<TPZMaterialT<STATE> *>(fFluxMesh->FindMaterial(matid));
    if (fracdim1)
    {
        bc = fracdim1->CreateBC(fracdim1, fHomogeneousNeumannBcMatId, 1, val1, val2);
        cmeshHDiv->InsertMaterialObject(bc);
    }
    
}

/// Insert the necessary pressure material objects to create the pressure mesh
void TPZMHMixedHybridMeshControl::InsertPeriferalPressureMaterialObjects()
{
    auto *matdim1 = new TPZNullMaterial<STATE>(fPressureDim1MatId);
    matdim1->SetDimension(fGMesh->Dimension()-1);
    fPressureFineMesh->InsertMaterialObject(matdim1);
    
    auto *matdim2 = new TPZNullMaterial<STATE>(fPressureDim2MatId);
    matdim2->SetDimension(fGMesh->Dimension()-2);
    fPressureFineMesh->InsertMaterialObject(matdim2);
    
    
    auto *matpres = new TPZNullMaterial<STATE>(fSkeletonWithFlowPressureMatId);
    matpres->SetDimension(fGMesh->Dimension()-1);
    fPressureFineMesh->InsertMaterialObject(matpres);
    
    for(auto it = fFractureFlowDim1MatId.begin(); it != fFractureFlowDim1MatId.end(); it++)
    {
        int matid = *it;
        auto *mat1d = new TPZNullMaterial<STATE>(matid);
        fPressureFineMesh->InsertMaterialObject(mat1d);
    }
    
    for(auto it = fSkeletonWithFlowMatId.begin(); it != fSkeletonWithFlowMatId.end(); it++)
    {
        int matid = *it;
        auto *mat1d = new TPZNullMaterial<STATE>(matid);
        fPressureFineMesh->InsertMaterialObject(mat1d);
    }
    
    TPZMHMixedMeshControl::InsertPeriferalPressureMaterialObjects();
}

// find a neighbouring element whose connect is equal to my side connect
static TPZGeoElSide FindCapElement(TPZGeoElSide fracel)
{
    TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *>(fracel.Element()->Reference());
    if(!fracel) DebugStop();
    if (intel->NSideConnects(fracel.Side()) != 1) {
        DebugStop();
    }
    int64_t connectindex = intel->SideConnectIndex(0, fracel.Side());
    TPZGeoElSide neighbour = fracel.Neighbour();
    while (neighbour != fracel) {
        TPZGeoEl *gelneigh = neighbour.Element();
        TPZInterpolatedElement *intelneigh = dynamic_cast<TPZInterpolatedElement *>(gelneigh->Reference());
        if (intelneigh && intelneigh->NConnects() == 1 && intelneigh->ConnectIndex(0) == connectindex) {
            return neighbour;
        }
        neighbour = neighbour.Neighbour();
    }
    return TPZGeoElSide();
}
/// for a given boundary wrap element
/// determine the appropriate boundary condition
/// determine if it has a neighbouring fracture element
/// change the material id of the geometric element to apply the boundary condition
void TPZMHMixedHybridMeshControl::ApplyNeighbourBoundaryCondition(TPZGeoEl *gel)
{
    // we assume the Flux mesh has been loaded
    if (fGMesh->Reference() != fFluxMesh.operator->()) {
        DebugStop();
    }
    int meshdim = fGMesh->Dimension();
#ifdef PZDEBUG
    if (!gel || gel->MaterialId() != fBoundaryWrapMatId || gel->Dimension() != meshdim-2) {
        DebugStop();
    }
#endif
    // only change the boundary condition if there is a neighbouring fracture element
    TPZGeoElSide gelside(gel,gel->NSides()-1);
    TPZGeoElSide neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        int neighmatid = neighbour.Element()->MaterialId();
        if (neighmatid == fPressureDim1MatId || neighmatid == fSkeletonWithFlowPressureMatId) {
            break;
        }
        neighbour = neighbour.Neighbour();
    }
    if (neighbour == gelside) {
        return;
    }
    {
        TPZManVector<REAL,3> co(3);
        gel->Node(0).GetCoordinates(co);
        std::cout << "ApplyNeighbourBoundaryCondition co = " << co << std::endl;
    }
    neighbour = gelside.Neighbour();
    // find a boundary wrap of dimension dim-1
    while (neighbour != gelside) {
        if (neighbour.Element()->Dimension() == meshdim-1 && neighbour.Element()->MaterialId() == fBoundaryWrapMatId) {
            break;
        }
        neighbour = neighbour.Neighbour();
    }
    if (neighbour == gelside) {
        DebugStop();
    }
    TPZGeoElSide boundside(neighbour.Element(),neighbour.Element()->NSides()-1);
    // we take the largest ancestor
    boundside = boundside.LowestFatherSide();
    // find an element associated with a boundary element
    TPZBndCond *bcfound = 0;
    while(boundside)
    {
        // look for a neighbour that qualifies
        neighbour = boundside.Neighbour();
        while (neighbour != boundside) {
            if (neighbour.Element()->Dimension() == meshdim-1 && neighbour.Element()->Reference()) {
                // a candidate
                int neighmatid = neighbour.Element()->MaterialId();
                if (fFractureFlowDim1MatId.find(neighmatid) != fFractureFlowDim1MatId.end()) {
                    continue;
                }
                if (fSkeletonWithFlowMatId.find(neighmatid) != fSkeletonWithFlowMatId.end()) {
                    continue;
                }
                if (neighmatid == fSkeletonMatId || neighmatid == fSecondSkeletonMatId || neighmatid == fHDivWrapperMatId) {
                    neighbour = neighbour.Neighbour();
                    continue;
                }
                TPZBndCond *bc = dynamic_cast<TPZBndCond *>(neighbour.Element()->Reference()->Material());
                if (bc) {
                    bcfound = bc;
                    break;
                }
            }
            neighbour = neighbour.Neighbour();
        }
        if (bcfound) {
            break;
        }
        // if not found, look for a neighbour of one of the subelements
        if (boundside.Element()->HasSubElement()) {
            TPZGeoEl *sub = boundside.Element()->SubElement(0);
            boundside = TPZGeoElSide(sub,sub->NSides()-1);
        }
        else
        {
            boundside = TPZGeoElSide();
        }
    }
    if (!bcfound) {
        DebugStop();
    }
    int bctype = bcfound->Type();
    
    // find next fracture flux geometric element
    neighbour = gelside.Neighbour();
    while (neighbour != gelside) {
        int neighmatid = neighbour.Element()->MaterialId();
        bool isfracture = fFractureFlowDim1MatId.find(neighmatid) != fFractureFlowDim1MatId.end();
        if(!isfracture) isfracture = fSkeletonWithFlowMatId.find(neighmatid) != fSkeletonWithFlowMatId.end();
        if (neighbour.Element()->Dimension() == meshdim-1 && neighbour.Element()->Reference() && isfracture) {
            // we found a fracture element
            TPZGeoElSide cap = FindCapElement(neighbour);
            if (!cap) {
                DebugStop();
            }
            if (bctype == 0) {
                cap.Element()->SetMaterialId(bcfound->Id());
            }
            else
            {
                cap.Element()->SetMaterialId(fHomogeneousNeumannBcMatId);
            }
        }
        neighbour = neighbour.Neighbour();
    }
}

/// Change the material id of the boundary elements associated with fracture flow
// if the element is neighbour of a pressure boundary condition, maintain the pressure bc
// if the element is neighbour of a flux boundary condition, apply a homogeneous flux condition
void TPZMHMixedHybridMeshControl::AdjustBoundaryConditionsOfFractures()
{
    fGMesh->ResetReference();
    fFluxMesh->LoadReferences();
    int meshdim = fGMesh->Dimension();
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(gel && gel->Dimension() == meshdim-2 && gel->MaterialId() == fBoundaryWrapMatId)
        {
            ApplyNeighbourBoundaryCondition(gel);
        }
    }
}

/// verify the consistency of the datastructure
//  only implemented in the hybrid Hdiv version
void TPZMHMixedHybridMeshControl::CheckMeshConsistency()
{
    // all the elements with material id fHDivWrapperMatId should have 3 elements connected
    // -> an Hdiv element - itself - HdivBound element - interface element
    int meshdim = fGMesh->Dimension();
    int64_t nel = fCMesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        if(!cel) continue;
        TPZMultiphysicsElement *mphys = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZMultiphysicsInterfaceElement *mphinterface = dynamic_cast<TPZMultiphysicsInterfaceElement *>(cel);
        if (!mphys && !mphinterface) {
            std::cout << "I dont understand\n";
            DebugStop();
        }
        if (!mphys) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        // check the consistency of the material id of the geometric element and the referred elements are equal
        int matid = mphys->Reference()->MaterialId();
        int nelref = mphys->NMeshes();
        if (nelref != 2) {
            std::cout << "I dont understand\n";
            DebugStop();
        }
        int nactive = 0;
        for (int iref=0; iref<nelref; iref++) {
            TPZCompEl *referred = mphys->Element(iref);
            if (referred) {
                nactive++;
                if (referred->Reference()->MaterialId() != matid) {
                    DebugStop();
                }
            }
        }
        if (matid == fHDivWrapperMatId) {
            int ncon = mphys->NConnects();
            if (ncon != 1) {
                std::cout << "I dont understand\n";
                DebugStop();
            }
            TPZConnect &c = mphys->Connect(0);
            if (c.NElConnected() < 3) {
                int ncorner = gel->NCornerNodes();
                TPZManVector<REAL,3> co(3);
                std::cout << "gel coordinates\n";
                for (int i=0; i<ncorner; i++) {
                    gel->Node(i).GetCoordinates(co);
                    std::cout << co << std::endl;
                }
                std::cout << "I dont understand nelconnected = " << c.NElConnected() <<  "\n";
                DebugStop();
            }
        }
        if (matid == fPressureDim1MatId) {
            int ncon = mphys->NConnects();
            if (meshdim == 2 && ncon != 3) {
                std::cout << "I dont understand\n";
                DebugStop();
            }
            if (gel->Dimension() != meshdim-1)
            {
                DebugStop();
            }
            for(int ic=0; ic<ncon; ic++)
            {
                TPZConnect &c = mphys->Connect(ic);
                if (c.NElConnected() != 4) {
                    std::cout << "I dont understand\n";
                    DebugStop();
                }
            }
        }
        if (matid == fPressureDim2MatId) {
            int ncon = mphys->NConnects();
            if (meshdim == 2 && ncon != 1) {
                std::cout << "point pressure should have only one connect\n";
                DebugStop();
            }
            TPZConnect &c = mphys->Connect(0);
            if (c.NElConnected() < 2) {
                std::cout << "point pressure element isolated\n";
                DebugStop();
            }
        }
        if (gel->Dimension() == meshdim && nactive == 2) {
            // the number of nelconnected should be 1
            int nc = mphys->NConnects();
            int nfluxconnects = mphys->Element(0)->NConnects();
            for (int ic=nfluxconnects; ic<nc; ic++) {
                TPZConnect &c = mphys->Connect(ic);
                if (c.NElConnected() != 1) {
                    std::cout << "I dont understand\n";
                    DebugStop();
                }
            }
        }
        
    }
    
    // the connects of the pressure elements should have NelConnected > 1 -> pressure + at least one interface element
    
    
}

/// print the elements in a readable format
void TPZMHMixedHybridMeshControl::PrintFriendly(std::ostream &out)
{
    int64_t nel = fCMesh->NElements();
    int meshdim = fGMesh->Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = fCMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != meshdim) {
            continue;
        }
        int ncorners = gel->NCornerNodes();
        TPZFNMatrix<12,REAL> coords(gel->NCornerNodes(),meshdim);
        for (int i=0; i<ncorners; i++) {
            TPZManVector<REAL,3> co(3);
            gel->Node(i).GetCoordinates(co);
            for (int j=0; j<meshdim; j++) {
                coords(i,j) = co[j];
            }
        }
        TPZMultiphysicsElement *mfel = dynamic_cast<TPZMultiphysicsElement *>(cel);
        TPZManVector<TPZCompEl *> els(mfel->NMeshes(),0);
        for(int i=0; i<mfel->NMeshes(); i++) els[i] = mfel->Element(i);
        out << "MF-El index " << el << " referred el indices";
        for (int i=0; i<mfel->NMeshes(); i++) {
            if (els[i]) {
                out << els[i]->Index() << " ";
            }
        }
        out << " connect indexes/numelcon ";
        for (int i=0; i<mfel->NConnects(); i++) {
            out << mfel->ConnectIndex(i) << "/" << mfel->Connect(i).NElConnected() << " ";
        }
        out << std::endl;
        coords.Print(out);
    }
    
}

