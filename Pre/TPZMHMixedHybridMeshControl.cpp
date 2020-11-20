//
//  TPZMHMixedHybridMeshControl.cpp
//  PZ
//
//  Created by Philippe Devloo on 09/10/16.
//
//

#include "TPZMHMixedHybridMeshControl.h"
#include "pzbndcond.h"
#include "TPZNullMaterial.h"
#include "TPZLagrangeMultiplier.h"


#include <iostream>
#include <sstream>
#include <iterator>
#include <numeric>

#include "TPZGeoElSideAncestors.h"
#include "TPZGeoElSidePartition.h"
#include "pzsubcmesh.h"

#include "pzbuildmultiphysicsmesh.h"
#include "TPZCompMeshTools.h"
#include "pzelementgroup.h"
#include "pzcondensedcompel.h"

#include "pzinterpolationspace.h"
#include "pzintel.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmultiphysicselement.h"

#include "pzmat1dlin.h"

#include "TPZVTKGeoMesh.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mhmixedhybridmeshcontrol"));
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




// create discontinuous elements and increase the internal flux orders
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
    cmeshHDiv->ApproxSpace().CreateDisconnectedElements(true);
    //Criar elementos computacionais malha MHM
    cmeshHDiv->AutoBuild(fMaterialIds);
    {
        int64_t nel = cmeshHDiv->NElements();
        for (int64_t i = 0; i < nel; i++) {
            TPZCompEl *cel = cmeshHDiv->Element(i);
            if(!cel) continue;
            // Create the flux element
            TPZGeoEl *gel = cel->Reference();
            int MHMDomain = fGeoToMHMDomain[gel->Index()];
            // Associate the connects with the subdomain for the flux mesh
            SetSubdomain(cel, MHMDomain);
        }
    }
    // create the boundary elements with the same connect as the volumetric element
    // we assume the boundary elements have the same level as the flux elements
    {
        int64_t nel = cmeshHDiv->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZCompEl *cel = cmeshHDiv->Element(el);
            TPZGeoEl *gel = cel->Reference();
            if(gel->Dimension() != meshdim) DebugStop();
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            for(int side = gel->FirstSide(meshdim-1); side < gel->NSides()-1; side++)
            {
                intel->SetSideOrient(side, 1);
                // create the boundary element with the same connect as the volumetric element
                TPZGeoElSide gelside(gel,side);
                TPZGeoElSide neighbour = gelside.HasNeighbour(fMaterialBCIds);
                if(neighbour)
                {
                    cel->LoadElementReference();
                    int64_t index;
                    TPZCompEl *bcel = cmeshHDiv->CreateCompEl(neighbour.Element(), index);
                    int domain = WhichSubdomain(cel);
                    SetSubdomain(bcel, domain);
                    gel->ResetReference();
                    neighbour.Element()->ResetReference();
                }
            }
        }
    }
    
    int64_t nelprev = cmeshHDiv->NElements();
    cmeshHDiv->SetDefaultOrder(fpOrderSkeleton);
    std::set<int> skelmatids({fSkeletonMatId,fSecondSkeletonMatId});
    cmeshHDiv->AutoBuild(skelmatids);
    {
        int64_t nel = cmeshHDiv->NElements();
        for (int el = nelprev; el<nel; el++) {
            TPZCompEl *cel = cmeshHDiv->Element(el);
            SetSubdomain(cel, -1);
        }
    }
    // duplicate the skeleton elements depending on whether the macro fluxes are hybridized or not
//    if(fHybridize == false)
//    {
//        int64_t nel = cmeshHDiv->NElements();
//        for (int64_t el = 0; el<nel; el++) {
//            TPZCompEl *cel = cmeshHDiv->Element(el);
//            TPZGeoEl *gel = cel->Reference();
//            if(gel->MaterialId() == fSkeletonMatId)
//            {
//                if(this->fInterfaces.find(gel->Index()) == fInterfaces.end()) DebugStop();
//                auto domains = this->fInterfaces[gel->Index()];
//                TPZGeoElSide gelside(gel);
//                TPZGeoElBC gbc1(gelside,fSkeletonMatId);
//                TPZGeoElBC gbc2(gelside,fSkeletonMatId);
//                int64_t index;
//                TPZCompEl *cel1 = cmeshHDiv->CreateCompEl(gbc1.CreatedElement(), index);
//                TPZCompEl *cel2 = cmeshHDiv->CreateCompEl(gbc2.CreatedElement(), index);
//                SetSubdomain(cel1, domains.first);
//                SetSubdomain(cel2, domains.second);
//            }
//        }
//    }
//    else
//    {
//        // the two flux elements already exist
//        // each one should stay in the super element
//        DebugStop();
//    }
    cmeshHDiv->LoadReferences();
    int64_t nel = cmeshHDiv->NElements();
    
    if (fGeoToMHMDomain.size() != gmesh->NElements()) DebugStop();

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    fFluxMesh->ExpandSolution();
    
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    {
        TPZVec<int64_t> eldata(fGeoToMHMDomain);
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel || gel->HasSubElement()) eldata[el] = -999;
        }
        std::ofstream out("MHMDomainNoWrap.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh.operator->(), out, eldata);
    }


    
    /// create dim-1 elements, either of material id fHDivWrapperMatId or fWrapMatID
    CreateHDivWrappers();
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// Create geometric Skeleton elements for flux and pressure in case of hybridization of skeleton fluxes
void TPZMHMixedHybridMeshControl::CreateSecondSkeletonGeometricElements()
{
    /// loop over the skeleton elements and create two geometric neighbouring elements
    int64_t nel = fGMesh->NElements();
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel || gel->MaterialId() != fSkeletonMatId || gel->HasSubElement())
        {
            continue;
        }
#ifdef PZDEBUG
        {
            // the element has to correspond a an interface element
            if(fInterfaces.find(el) == fInterfaces.end())
            {
                DebugStop();
            }
            // it cannot be a boundary interface (they keep the matid of the boundary)
            if(fInterfaces[el].second == el) DebugStop();
        }
#endif
        TPZGeoElSide gelside(gel);
        TPZGeoElBC gbc(gelside,fPressureSkeletonMatId);
        TPZGeoElSide presside(gbc.CreatedElement());
        TPZGeoElBC(presside,fSecondSkeletonMatId);
    }
}

/// Create geometric elements for the pressure interface elements with mat id fPressureDim1MatId
// this method relies upon the existence of the HDivWrapMatId elements
void TPZMHMixedHybridMeshControl::CreatePressureInterfaceGeometricElements()
{
    // for each HDivWrapElement, decide if a neighbouring pressure element
    // needs to be created
    // always create :
    // - when there is a skeleton element (larger or equal)
    // when NOT TO CREATE :
    //  - when there is a neighbouring pressure element of fPressureDim1MatId
    //  - when there is a smaller HDivWrapMatiD element
    int64_t nel = fGMesh->NElements();
    fGMesh->ResetReference();
    int dim = fGMesh->Dimension();
    std::set<int> SkelMatids(fMaterialBCIds);
    SkelMatids.insert(fSkeletonMatId);
    for(int64_t el=0; el<nel; el++)
    {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel || gel->Dimension() != dim-1 || gel->MaterialId() != fHDivWrapperMatId) continue;
        TPZGeoElSide gelside(gel);
        int domain = fGeoToMHMDomain[el];
        TPZGeoElSideAncestors ancest(gelside);
        TPZGeoElSidePartition part(gelside);
        bool has_smaller = part.HasHigherLevelNeighbour(fHDivWrapperMatId);
        bool has_skel_neighbour = ancest.HasLargerorEqual(SkelMatids);
        bool has_large = ancest.HasLarger(fPressureDim1MatId);
        if(has_large) DebugStop();
        bool has_equal = gelside.HasNeighbour(fPressureDim1MatId);
        if((!(has_equal) && !has_smaller) || has_skel_neighbour)
        {
            TPZGeoElBC gbc(gelside,fPressureDim1MatId);
            SetSubdomain(gbc.CreatedElement(), domain);
//            int64_t index;
//            fPressureFineMesh->CreateCompEl(gbc.CreatedElement(), index);
//            gbc.CreatedElement()->ResetReference();
        }
    }
    /*
    /// adjust the polynomial orders
    TPZVec<TPZInterpolatedElement *> intel(nel,0);
    /// create a vector pointing to the flux elements
    int64_t nelflux = fFluxMesh->NElements();
    for (int64_t el = 0; el<nelflux; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        TPZInterpolatedElement *inter = dynamic_cast<TPZInterpolatedElement *>(cel);
        if(!inter) DebugStop();
        intel[gel->Index()] = inter;
    }
    /// for each pressure interface, set the order of the pressure to the flux neighbour
    // if there are two neighbours, it is assumed they have the same polynomial order
    int64_t nelpres = fPressureFineMesh->NElements();
    for (int64_t el = 0; el<nelpres; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        if(!cel) DebugStop();
        TPZGeoEl *gel = cel->Reference();
        if(gel->MaterialId() != fPressureDim1MatId) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.HasNeighbour(fMaterialIds);
        if(!neighbour) DebugStop();
        TPZInterpolatedElement *cintel = intel[neighbour.Element()->Index()];
        if(!cintel) DebugStop();
        TPZConnect &c = cintel->MidSideConnect(neighbour.Side());
        int order = c.Order();
        cintel->PRefine(order);
    }
    fPressureFineMesh->ExpandSolution();
     */
}



/// Create all data structures for the computational mesh
void TPZMHMixedHybridMeshControl::BuildComputationalMesh(bool usersubstructure)
{
    if (fpOrderInternal == 0 || fpOrderSkeleton == 0) {
        DebugStop();
    }
    // when the approximation space is HDiv, it makes no sense to have "large boundary skeleton elements
    DivideBoundarySkeletonElements();
    if(fHybridizeSkeleton)
    {
        /// hybridize the normal fluxes
        CreateSecondSkeletonGeometricElements();
    }
    InsertPeriferalMaterialObjects();

#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        sout << "Subdomain indices " << fGeoToMHMDomain << std::endl;
        sout << "Recognized domains ";
        for (auto it = fMHMtoSubCMesh.begin(); it != fMHMtoSubCMesh.end(); it++) {
            sout << it->first << ' ';
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    // this method will create the flux elements, include HDivWratMatId elements
    CreateHDivMHMMesh();
    
    {
        TPZVec<int64_t> eldata(fGeoToMHMDomain);
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel || gel->HasSubElement()) eldata[el] = -999;
        }
        std::ofstream out("MHMDomainNoPress.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh.operator->(), out, eldata);
    }

    // create a pressure dim1 geometric element next to each fHDivWrapperMatId element
    CreatePressureInterfaceGeometricElements();
    
    {
        TPZVec<int64_t> eldata(fGeoToMHMDomain);
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel || gel->HasSubElement()) eldata[el] = -999;
        }
        std::ofstream out("MHMDomainWithPress.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh.operator->(), out);
    }
    // hybridize the fluxes between macro elements
    
#ifdef PZDEBUG
    if (fFluxMesh->Dimension() != fGMesh->Dimension()) {
        DebugStop();
    }
#endif
    // create the discontinuous pressure mesh
    // this method creates the internal pressure elements (super class)
    // this method will create the dim-1 skeleton elements as well
    CreatePressureMHMMesh();
#ifdef PZDEBUG
    {
        int numc = fPressureFineMesh->NConnects();
        int numc2 = fConnectToSubDomainIdentifier[fPressureFineMesh.operator->()].size();
        if(numc != numc2) DebugStop();
    }
#endif
    if(fNState > 1)
    {
        CreateRotationMesh();
    }
    // create meshes with distributed fluxes and average solutions for each element
    CreateAverageSolutionMeshes();
    
    // create the meshes with one rigid body mode per subdomain
    CreateLagrangeMultiplierMesh();
    
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

    {
        TPZVec<int64_t> eldata(fGeoToMHMDomain);
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel || gel->HasSubElement()) eldata[el] = -999;
        }
        std::ofstream out("MHMDomainAfterPress.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh.operator->(), out, eldata);
    }

    
    // Set Lagrange multiplier levels of the connects of the atomic meshes
    SetLagrangeMultiplierLevels();
    
    /// create the multiphysics mesh
    CreateMultiphysicsMHMMesh();


    std::cout << "Total number of equations " << fCMesh->Solution().Rows() << std::endl;
    fGlobalSystemWithLocalCondensationSize = fCMesh->NEquations();
    fGlobalSystemSize = fCMesh->Solution().Rows();
    
    fCMesh->ComputeNodElCon();
    /*
    {
        int neq = fPressureFineMesh->NEquations();
        for(int eq=0; eq<neq; eq++) fPressureFineMesh->Solution()(eq,0) = 1.;
    }
    {
        int neq = fAverageSolution->NEquations();
        for(int eq=0; eq<neq; eq++) fAverageSolution->Solution()(eq,0) = 1.;
    }
    {
        int neq = fCMeshConstantPressure->NEquations();
        for(int eq=0; eq<neq; eq++) fCMeshConstantPressure->Solution()(eq,0) = 1.;
    }
    fCMesh->LoadSolutionFromMeshes();
     */
#ifdef PZDEBUG
    if(1)
    {
        fFluxMesh->ComputeNodElCon();
        std::ofstream outfl("FluxMesh.txt");
        fFluxMesh->Print(outfl);
        fPressureFineMesh->ComputeNodElCon();
        std::ofstream outpr("PressureMesh.txt");
        fPressureFineMesh->Print(outpr);
        if(fRotationMesh)
        {
            std::ofstream outrot("RotationMesh.txt");
            fRotationMesh->Print(outrot);
        }
        std::ofstream outdistflux("DistrFluxMesh.txt");
        fDistributedFlux->Print(outdistflux);
        std::ofstream outavSol("AveragePressureMesh.txt");
        fAverageSolution->Print(outavSol);
        std::ofstream outdistFluxDomain("DistrFluxByDomain.txt");
        fCMeshLagrange->Print(outdistFluxDomain);
        std::ofstream outAvPresDomain("AveragePressureByDomainMesh.txt");
        fCMeshConstantPressure->Print(outAvPresDomain);

        std::ofstream outmf("MFMesh.txt");
        fCMesh->Print(outmf);
        {
            int64_t nel = fCMesh->NElements();
            for (int64_t el = 0; el < nel; el++) {
                TPZCompEl *cel = fCMesh->Element(el);
                TPZGeoEl *gel = cel->Reference();
                int domain = WhichSubdomain(cel);
                outmf << "El index " << el << " matid " << gel->MaterialId() << " domain " << domain << std::endl;
            }
        }
        std::ofstream out("Friendly.txt");
        PrintFriendly(out);
    }
    CheckMeshConsistency();
#endif
    
    if (usersubstructure) {
        char average_displacement_subdomain = 10;
        HideTheElements(average_displacement_subdomain);
    }
    fNumeq = fCMesh->NEquations();

}

/// create the multiphysics mesh
void TPZMHMixedHybridMeshControl::CreateMultiphysicsMHMMesh()
{
    TPZManVector<TPZCompMesh *,9> meshvec;
    GetMeshVec(meshvec);
    // build the multiphysics mesh
    TPZVec<int64_t> gelindices;
    this->GetGeometricElementPartition(gelindices);
    
    fCMesh->BuildMultiphysicsSpace(meshvec,gelindices);
    // identify connects with submeshes
    JoinSubdomains(meshvec, fCMesh.operator->());
    
    // create the interface elements
    int dim = fCMesh->Dimension();
    CreateMultiPhysicsInterfaceElements(dim-1);

    {
        TPZVec<int64_t> eldata(fGeoToMHMDomain);
        int64_t nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if(!gel || gel->HasSubElement()) eldata[el] = -999;
        }
        std::ofstream out("MHMDomain.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fGMesh.operator->(), out, eldata);
    }

    // put elements in substructures
    // -> the macro fluxes stay in the coarse domain
    // group and condense elements in substructures (2 levels)
    // -> internal fluxes, pressure, distributed flux onto boundary flux, average
    // -> boundary flux, average onto skeleton pressure
    // at the level of macro domain
    // condense the fluxes, average onto the pressure skeleton
//    DebugStop();
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
    std::set<int> matids = {fPressureSkeletonMatId};
    int nelbefore = fPressureFineMesh->NElements();
    fPressureFineMesh->ApproxSpace().SetAllCreateFunctionsContinuous();
    fPressureFineMesh->ApproxSpace().CreateDisconnectedElements(true);
    fPressureFineMesh->SetDefaultOrder(fpOrderSkeleton);
    fPressureFineMesh->AutoBuild(matids);
    int Nelafter = fPressureFineMesh->NElements();
    for (int64_t el = nelbefore; el< Nelafter; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        SetSubdomain(cel, -1);
    }
    
    matids.clear();
    matids.insert(fPressureDim1MatId);
    fPressureFineMesh->SetDefaultOrder(fpOrderInternal);
    nelbefore = Nelafter;
    fPressureFineMesh->AutoBuild(matids);
    Nelafter = fPressureFineMesh->NElements();
    for (int64_t el = nelbefore; el< Nelafter; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        int domain = fGeoToMHMDomain[gel->Index()];
        SetSubdomain(cel, domain);
    }

    fPressureFineMesh->ExpandSolution();
}

/// Create the pressure mesh which is dual to the flux mesh
void TPZMHMixedHybridMeshControl::CreatePressureMHMMesh()
{
    // this will create all pressure elements of the dimension of the geometric
    // mesh
    TPZMHMixedMeshControl::CreatePressureMHMMesh();
    // will create a pressure element associated with internal interfaces only
    {
        int numc = fPressureFineMesh->NConnects();
        int numc2 = fConnectToSubDomainIdentifier[fPressureFineMesh.operator->()].size();
        if(numc != numc2) DebugStop();
    }
    

    CreatePressureInterfaces();
}



// create the approximation space associated with the skeleton and restrain the connects
void TPZMHMixedHybridMeshControl::CreateSkeleton()
{
    // do nothing, the skeleton elements were generated with the regular flux
    // elements
    
//    TPZMHMixedMeshControl::CreateSkeleton();
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
        if(fluxelements[el] == 0) continue;
        TPZGeoEl *gel = fGMesh->Element(el);
        if(gel->Dimension() != meshdim) continue;
        for(int side = gel->NCornerNodes(); side < gel->NSides()-1; side++)
        {
            if(gel->SideDimension(side) != meshdim-1) continue;
            TPZGeoElSide gelside(gel,side);
            TPZGeoElSideAncestors ancestor(gelside);
            // create an element on the side of the H(div) element
            // verify if there a BC neighbouring
            TPZGeoElSide neighbour;
            neighbour = ancestor.HasLargerorEqual(fMaterialBCIds);
            // if the element is on the domain boundary, no wrapper is created
            if(neighbour) continue;
            // create a neighbouring element HDivWrap or BC matid
            TPZGeoElBC gbc(gelside,fHDivWrapperMatId);
            neighbour = TPZGeoElSide(gbc.CreatedElement());
            int64_t neighindex = neighbour.Element()->Index();
            int neighmatid = neighbour.Element()->MaterialId();
            fluxelements[el]->LoadElementReference();
            TPZConnect &c = fluxelements[el]->SideConnect(0, side);
            fFluxMesh->SetDefaultOrder(c.Order());
            int64_t index;
            TPZCompEl *cel = fFluxMesh->ApproxSpace().CreateCompEl(neighbour.Element(), fFluxMesh, index);
            int64_t mhm_domain = fGeoToMHMDomain[gel->Index()];
            // if the geometric element already existed, keep its mhm domain index
            SetSubdomain(cel, mhm_domain);
            gel->ResetReference();
            cel->Reference()->ResetReference();
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
        else if(neighbour.Element()->MaterialId() == fWrapMatId && pressuredim == meshdim-1)
        {
            return;
        }
        else if(neighbour.Element()->MaterialId() == fWrapMatId)
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
    int problemdim = fGMesh->Dimension();
    if(dim != problemdim-1)
    {
        // untested
        DebugStop();
    }
    if(fGMesh->Reference()) fGMesh->ResetReference();
    fCMesh->LoadReferences();
    int64_t nel = fGMesh->NElements();
    std::set<int> allSkeleton(fMaterialBCIds);
    allSkeleton.insert(fSkeletonMatId);
    // loop over the hybridized pressure elements
    for (int64_t el = 0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if(!gel || gel->MaterialId() != fPressureDim1MatId) continue;
        TPZGeoElSide gelside(gel);
        int domain = fGeoToMHMDomain[gel->Index()];
        TPZGeoElSideAncestors ancestor(gelside);
        TPZGeoElSide skelside = ancestor.HasLargerorEqual(allSkeleton);
        TPZStack<TPZGeoElSide> connected;
        // by default the sign of the lagrange interface corresponds to left matid
        int interfacematid = fLagrangeMatIdLeft;
        if(skelside && !fHybridizeSkeleton)
        {
            connected.Push(skelside);
            int64_t gelindex = skelside.Element()->Index();
            if(fInterfaces.find(gelindex) == fInterfaces.end()) DebugStop();
            auto twodomains = fInterfaces[gelindex];
            if(twodomains.first == domain)
            {
                interfacematid = fLagrangeMatIdLeft;
            }
            else if (twodomains.second == domain)
            {
                interfacematid = fLagrangeMatIdRight;
            }
            else DebugStop();
        }
        else if(fHybridizeSkeleton)
            DebugStop();
        TPZGeoElSide neighbour = gelside.Neighbour();
        while(neighbour != gelside)
        {
            int neighdomain = WhichSubdomain(neighbour.Element());
            bool samedomain = neighdomain == domain;
            if(neighbour.Element()->MaterialId() == fHDivWrapperMatId && samedomain) connected.Push(neighbour);
            neighbour = neighbour.Neighbour();
        }
        if(connected.size() < 2)
        {
            TPZGeoElSide large = ancestor.HasLarger(fHDivWrapperMatId);
            if(large) connected.Push(large);
        }
        if(connected.size() != 2)
        {
            std::cout << "the pressure element didn't find two flux elements to connect\n";
            gel->Print();
            DebugStop();
        }
        // generate the interface element
        TPZCompElSide celpressuredim1 = gelside.Reference();
        TPZCompElSide celconnect1 = connected[0].Reference();
        TPZCompElSide celconnect2 = connected[1].Reference();
        if(!celpressuredim1 || !celconnect1 || !celconnect2) DebugStop();
        TPZGeoElBC gbc1(gelside,interfacematid);
        TPZGeoElBC gbc2(gelside,fLagrangeMatIdLeft);
        int64_t index;
        TPZCompEl *interf1 = new TPZMultiphysicsInterfaceElement(fCMesh,gbc1.CreatedElement(),index,celpressuredim1,celconnect1);
        TPZCompEl *interf2 = new TPZMultiphysicsInterfaceElement(fCMesh,gbc2.CreatedElement(),index,celpressuredim1,celconnect2);
        SetSubdomain(interf2, domain);
        if(!skelside) SetSubdomain(interf1, domain);
        else SetSubdomain(gbc1.CreatedElement(), domain);
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
        TPZCompMeshTools::GroupElements(subcmesh);
        TPZStack<int64_t> grouped_elements;
        {
            int64_t nel = subcmesh->NElements();
            for (int64_t el = 0; el<nel; el++) {
                TPZCompEl *cel = subcmesh->Element(el);
                TPZElementGroup *group = dynamic_cast<TPZElementGroup *>(cel);
                if(group)
                {
                    grouped_elements.Push(group->Index());
                }
            }
        }
//        int lagrangeLevelPressureInterface = 6;
//        TPZCompMeshTools::CondenseElements(subcmesh, lagrangeLevelPressureInterface);
//        TPZStack<int64_t> condensed_elements;
//        {
//            int64_t nel = subcmesh->NElements();
//            for (int64_t el = 0; el<nel; el++) {
//                TPZCompEl *cel = subcmesh->Element(el);
//                TPZCondensedCompEl *condense = dynamic_cast<TPZCondensedCompEl *>(cel);
//                if(condense) condensed_elements.Push(condense->Index());
//            }
//        }
        unsigned char boundaryflux_lagrangelevel = 5;
        TPZCompMeshTools::AgglomerateElements(subcmesh, grouped_elements,boundaryflux_lagrangelevel);
        subcmesh->ComputeNodElCon();
        bool keeplagrange = false;
        bool keepmatrix = false;
        TPZCompMeshTools::CreatedCondensedElements(subcmesh, keeplagrange, keepmatrix);
        subcmesh->InitializeBlock();
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            std::stringstream sout;
            subcmesh->Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
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
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
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
#ifdef LOG4CXX
                        if (logger->isDebugEnabled())
                        {
                            std::stringstream sout;
                            sout << "Adding interface " << mphys->Index() << " right connect index " << right->ConnectIndex(0);
                            LOGPZ_DEBUG(logger, sout.str())
                        }
#endif
                        elgroup.insert(mphys->Index());
                    }
                    if (left->NConnects() == 1 && left->ConnectIndex(0)== cindex) {
#ifdef LOG4CXX
                        if (logger->isDebugEnabled())
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
#ifdef LOG4CXX
                    if (logger->isDebugEnabled())
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
            int64_t grindex;
            TPZElementGroup *elgr = new TPZElementGroup(*cmesh,grindex);
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
    }
}
#ifdef MACOSX
#pragma clang diagnostic pop
#endif


/// Insert Boundary condition objects that do not perform any actual computation
void TPZMHMixedHybridMeshControl::InsertPeriferalMaterialObjects()
{
    TPZMHMixedMeshControl::InsertPeriferalMaterialObjects();
    int dim = fGMesh->Dimension();
    {
        TPZNullMaterial *mat = new TPZNullMaterial(fHDivWrapperMatId);
        mat->SetNStateVariables(fNState);
        mat->SetDimension(dim-1);
        fCMesh->InsertMaterialObject(mat);
    }
    {
        TPZNullMaterial *mat = new TPZNullMaterial(fPressureDim1MatId);
        mat->SetNStateVariables(fNState);
        mat->SetDimension(dim-1);
        fCMesh->InsertMaterialObject(mat);
    }


    TPZLagrangeMultiplier *matleft = new TPZLagrangeMultiplier(fLagrangeMatIdLeft,dim,fNState);
    TPZLagrangeMultiplier *matright = new TPZLagrangeMultiplier(fLagrangeMatIdRight,dim,fNState);
    matright->SetMultiplier(-1.);
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
    TPZNullMaterial *matl2 = 0;
    for(auto iter:cmeshHDiv->MaterialVec())
    {
        TPZMaterial *mat = iter.second;
        matl2 = dynamic_cast<TPZNullMaterial *>(mat);
        if(matl2) break;
    }
    if (!matl2) {
        DebugStop();
    }
    TPZFNMatrix<1,STATE> val1(fNState,fNState,0.),val2(fNState,1,0.);
    TPZBndCond *bc;
    bc = matl2->CreateBC(matl2, fHDivWrapperMatId, 0, val1, val2);
    cmeshHDiv->InsertMaterialObject(bc);
    

}

/// Insert the necessary pressure material objects to create the pressure mesh
void TPZMHMixedHybridMeshControl::InsertPeriferalPressureMaterialObjects()
{
    TPZMHMixedMeshControl::InsertPeriferalPressureMaterialObjects();
    
    if(fPressureDim1MatId == 0) DebugStop();
    if(fPressureFineMesh->FindMaterial(fPressureDim1MatId) != 0) DebugStop();
    TPZNullMaterial *nullmat = new TPZNullMaterial(fPressureDim1MatId);
    nullmat->SetDimension(fPressureFineMesh->Dimension());
    nullmat->SetNStateVariables(fNState);
    fPressureFineMesh->InsertMaterialObject(nullmat);
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
/// verify the consistency of the datastructure
//  only implemented in the hybrid Hdiv version
void TPZMHMixedHybridMeshControl::CheckMeshConsistency()
{
    // all the elements with material id fHDivWrapperMatId should have 3 elements connected
    // -> an Hdiv element - itself - HdivBound element - interface element
    int meshdim = fGMesh->Dimension();
    int64_t nel = fCMesh->NElements();
    fCMesh->ComputeNodElCon();
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
        // in the near future this will be seven...
        if ((fNState > 1 && nelref != 7) || (fNState == 1 && nelref != 6)) {
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
                std::cout << "gel index " << gel->Index() << std::endl;
                std::cout << "material id " << gel->MaterialId() << std::endl;
                std::cout << "Connect index " << mphys->ConnectIndex(0) << std::endl;
                gel->Print();
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
                if (c.NElConnected() != 3) {
                    std::cout << "I dont understand\n";
                    DebugStop();
                }
            }
        }
        if (gel->Dimension() == meshdim && nactive == 7) {
            // the number of nelconnected of all connects except the fluxconnects should be 1
            int nc = mphys->NConnects();
            int nfluxconnects = mphys->Element(0)->NConnects();
            // the last two connects receive contributions of all elements in a macro domain
            for (int ic=nfluxconnects; ic<nc-2; ic++) {
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

/// Control the numbering of the equations based on the lagrange multiplier level of the connects
void TPZMHMixedHybridMeshControl::SetLagrangeMultiplierLevels()
{
//    flux internal 0 - DONE
//    1, 2 or 3 displacement connects 1 DONE
//    distributed flux 2 DONE
//    other displacement 3 DONE
//    rotation connect 4 DONE
//    element boundary flux 5 DONE
//    element average displacement 6 DONE
//    all hybrid displacement - 1 7 DONE
//    dist flux subdomain 8 DONE
//    last hybrid displacement 9 DONE
//    skeleton fluxes 10 DONE
//    av displ subdomain 11 DONE
//    hybrid displacement skeleton 12 DONE
    
    int dim = fCMesh->Dimension();
    int64_t nelflux = fFluxMesh->NElements();
    for (int64_t el=0; el<nelflux; el++) {
        TPZCompEl *cel = fFluxMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == dim)
        {
            // internal (last) connect level 0
            // other connects 5
            int nc = cel->NConnects();
            for(int ic=0; ic<nc-1; ic++) cel->Connect(ic).SetLagrangeMultiplier(5);
            cel->Connect(nc-1).SetLagrangeMultiplier(0);
        }
        int matid = gel->MaterialId();
        if(matid == fSkeletonMatId || matid == fSecondSkeletonMatId)
        {
            if(cel->NConnects() != 1) DebugStop();
            cel->Connect(0).SetLagrangeMultiplier(10);
        }
    }
    int nfirst_disp = 1;
    if(fNState == 2) nfirst_disp = 2;
    if(fNState == 3) nfirst_disp = 3;
    int64_t neldisp = fPressureFineMesh->NElements();
    std::map<int64_t,int64_t> last_hybrid_dispel;
    for (int64_t el = 0; el<neldisp; el++) {
        TPZCompEl *cel = fPressureFineMesh->Element(el);
        TPZGeoEl *gel = cel->Reference();
        if(gel->Dimension() == dim)
        {
            int nc = cel->NConnects();
            if(nc <= nfirst_disp) DebugStop();
            for(int ic=0; ic<nfirst_disp; ic++) cel->Connect(ic).SetLagrangeMultiplier(1);
            for(int ic=nfirst_disp; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(3);
        }
        int matid = gel->MaterialId();
        if(matid == fPressureDim1MatId)
        {
            int domain = fGeoToMHMDomain[gel->Index()];
            if(domain == -1) DebugStop();
            last_hybrid_dispel[domain] = el;
            int nc = cel->NConnects();
            for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(7);
        }
        if(matid == fPressureSkeletonMatId)
        {
            int nc = cel->NConnects();
            for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(12);
        }
    }
    for(auto it : last_hybrid_dispel)
    {
        TPZCompEl *cel = fPressureFineMesh->Element(it.second);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(9);
    }
    if(fRotationMesh)
    {
        int64_t nelrot = fRotationMesh->NElements();
        for (int64_t el = 0; el<nelrot; el++) {
            TPZCompEl *cel = fRotationMesh->Element(el);
            int nc = cel->NConnects();
            for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(4);
        }
    }
    int neldistflux = fDistributedFlux->NElements();
    for(int64_t el = 0; el<neldistflux; el++)
    {
        TPZCompEl *cel = fDistributedFlux->Element(el);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(2);

    }
    int nelavdisp = fAverageSolution->NElements();
    for(int64_t el = 0; el<nelavdisp; el++)
    {
        TPZCompEl *cel = fAverageSolution->Element(el);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(6);

    }
    int neldistfluxDomain = fCMeshLagrange->NElements();
    for(int64_t el = 0; el<neldistfluxDomain; el++)
    {
        TPZCompEl *cel = fCMeshLagrange->Element(el);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(8);

    }
    int nelavdispDomain = fCMeshConstantPressure->NElements();
    for(int64_t el = 0; el<nelavdispDomain; el++)
    {
        TPZCompEl *cel = fCMeshConstantPressure->Element(el);
        int nc = cel->NConnects();
        for(int ic=0; ic<nc; ic++) cel->Connect(ic).SetLagrangeMultiplier(11);

    }
}
