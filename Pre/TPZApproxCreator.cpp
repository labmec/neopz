//
// Created by victor on 21/09/2022.
//

#include "TPZApproxCreator.h"
#include "TPZBndCond.h"
#include "pzgeoelbc.h"
#include "TPZGeoElSidePartition.h"
#include "TPZNullMaterialCS.h"
#include "TPZMultiphysicsCompMesh.h"
#include "pzintel.h"
#include "TPZLagrangeMultiplierCS.h"
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.CreateMultiphysicsSpace");
#endif

TPZApproxCreator::TPZApproxCreator(TPZGeoMesh *gmesh) : fGeoMesh(gmesh)
{
}

/**Insert a material object in the datastructure*/
int TPZApproxCreator::InsertMaterialObject(TPZMaterial * mat) {
    if(!mat) DebugStop();
    const int matid = mat->Id();
    int nstate = mat->NStateVariables();
    SetNState(nstate);
    if(fMaterialVec.find(matid) != fMaterialVec.end()) DebugStop();
    fMaterialVec.insert(std::make_pair(matid, mat));
    return (int)fMaterialVec.size();
}

int TPZApproxCreator::InsertMaterialObject(TPZBndCond *mat)  {
    return InsertMaterialObject(dynamic_cast<TPZMaterial*>(mat));
}

void TPZApproxCreator::ComputePeriferalMaterialIds(int base) {
    
    // TODO: Victor. Please explain the logic here. For instance, why is base by default 10 and if base < 2 then base = 2? Also, why aren't the matids in order of +1?
    if(base < 2) base = 2;
    int max_matid = 0;
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) continue;
        max_matid = std::max(max_matid,gel->MaterialId());
    }
    int remain = max_matid % base;
    int matid_base = max_matid - remain + base;
    fHybridizationData.fWrapMatId = matid_base;
    fHybridizationData.fLeftInterfaceMatId = matid_base+base;
    fHybridizationData.fRightInterfaceMatId = matid_base+2*base;
    fHybridizationData.fLagrangeMatId = matid_base + 3*base;
    fHybridizationData.fSecondLeftInterfaceMatId = matid_base + 4*base;
    fHybridizationData.fSecondRightInterfaceMatId = matid_base + 5*base;
    fHybridizationData.fSecondLagrangeMatId = matid_base + 6*base;
}

void TPZApproxCreator::AddHybridizationGeoElements(){
    /// For fHybridType == HybridizationType::EStandard, wrap, interface and lagrange geometric elements are created;
    ///     if fHybridizeBCLevel == 0, the boundary is not hybridized, if fHybridizeBCLevel == 1, it is.
    /// For fHybridType == HybridizationType::EStandardSquared, wrap, interface, lagrange, second interface and second lagrange
    /// geometric elements are created;
    ///     BC can be hybridized 0, 1 or 2 times in this set up;
    /// fHybridType == HybridizationType::Semi has yet to be implemented.

    if(fHybridType == HybridizationType::EStandardSquared) {
        std::cout << "Hybridization type squared should be implemented here\n";
        DebugStop();
    }
    if(fHybridizationData.fWrapMatId == -123456){
        std::cout << "\nERROR! Please call TPZApproxCreator::ComputePeriferalMaterialIds() before TPZApproxCreator::AddHybridizationGeoElements()" << std::endl;
        DebugStop();
    }

#ifdef PZ_LOG
    std::map<int,int> numcreated;
#endif
    // create the wrap and interface geometric elements
    // the creation of wrap and interface needs to happen in separate loops when there are hanging nodes
    int64_t nel = fGeoMesh->NElements();
    int dim = fGeoMesh->Dimension();
    std::set<int> bcMatIds = GetBCMatIds();
    // Wrap and interface creation for standard hyb.
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        int side = gel->FirstSide(dim-1);
        // loop over the sides of dimension dim-1
        for(; side < nsides-1; side++)
        {
            TPZGeoElSide gelside(gel,side);
            // we want to create side elements of type
            // first fMatWrapId
            TPZGeoElSide neighbour = gelside.Neighbour();
#ifdef PZDEBUG
            {
                int neighMatId = neighbour.Element()->MaterialId();
                if(neighMatId == fHybridizationData.fWrapMatId)
                {
                    std::cout << __PRETTY_FUNCTION__ << " should be called only once!\n";
                    DebugStop();
                }
            }
#endif
            // if there is a neighbour that is a boundary condition
            // do not create the wrap layers
            bool hasBCNeighbour = gelside.HasNeighbour(bcMatIds);
            if(hasBCNeighbour && fHybridizationData.fHybridizeBCLevel == 0)
            {
                // no interface will be created between the element and a flux space
                continue;
            }
            // create the wrap material id
            TPZGeoElBC(gelside, fHybridizationData.fWrapMatId);
#ifdef PZ_LOG
            numcreated[fHybridizationData.fWrapMatId]++;
#endif
#ifdef PZDEBUG
            neighbour = gelside.Neighbour();
            if(neighbour.Element()->MaterialId() != fHybridizationData.fWrapMatId)
            {
                DebugStop();
            }
#endif
        }
    }
    
    // Creates the lagrange geometric elements
    // we create the interface geometric elements here
    nel = fGeoMesh->NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim-1) continue;
        int gelmatid = gel->MaterialId();
        if(gelmatid != fHybridizationData.fWrapMatId) continue;
        TPZGeoElSide gelside(gel);
        TPZGeoElSide neighbour = gelside.Neighbour();
        // if the neighbour is a boundary condition and no hybridization is applied
        // do not create the wrap layers
        int neighmat = neighbour.Element()->MaterialId();
        int leftinterface = fHybridizationData.fLeftInterfaceMatId;
        int rightinterface = fHybridizationData.fRightInterfaceMatId;
        
        TPZGeoElSide hasLagrangeElement = gelside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        
        TPZGeoElSide hasLargeElementNeighbour = gelside.HasLowerLevelNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide wrapNeighbour = neighbour.HasNeighbour(fHybridizationData.fWrapMatId);
        // the haswrapneighbour is false
        // if on a boundary and if bc hybridization is 0
        // if there are smaller elements connected
        bool haswrapNeighbour = (wrapNeighbour != gelside);
        
        TPZGeoElSide hasbcNeighbour = gelside.HasNeighbour(bcMatIds);
        
        if(hasLargeElementNeighbour && haswrapNeighbour) DebugStop();
        
        if(haswrapNeighbour) { // neighbouring elements of the same level
            if(hasLagrangeElement) {
                // we create the interface geometric element (necessarily)
                TPZGeoElBC gbc(gelside,rightinterface);
                int64_t lagr_id = hasLagrangeElement.Element()->Index();
                int64_t wrap_id = wrapNeighbour.Element()->Index();
                int64_t int_id = gbc.CreatedElement()->Index();
                TPZHybrid hybridconfig(int_id,wrap_id,lagr_id);
                fHybridizationData.fInterfaces.Push(hybridconfig);
#ifdef PZ_LOG
                numcreated[rightinterface]++;
#endif
            } else {
                TPZGeoElBC gbc(gelside,leftinterface);
                // we create the lagrange geometric element (necessarily)
                TPZGeoElBC gbcLagrange(gbc.CreatedElement(),fHybridizationData.fLagrangeMatId);
                int64_t int_id = gbc.CreatedElement()->Index();
                int64_t wrap_id = gel->Index();
                int64_t lagr_id = gbcLagrange.CreatedElement()->Index();
                TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                fHybridizationData.fInterfaces.Push(hybrid);
#ifdef PZ_LOG
                numcreated[leftinterface]++;
                numcreated[fHybridizationData.fLagrangeMatId]++;
#endif
            }

        } else { // it has no wrap neighbour
                // handling hanging sides
            if(hasLargeElementNeighbour) {
                // we create the lagrange geometric element (necessarily)
                // from small to large -> leftinterface
                TPZGeoElBC gbc(gelside,leftinterface);
                int64_t lagr_id = hasLargeElementNeighbour.Element()->Index();
                int64_t int_id = gbc.CreatedElement()->Index();
                int64_t wrap_id = gel->Index();
                TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                fHybridizationData.fInterfaces.Push(hybrid);
#ifdef PZ_LOG
                numcreated[leftinterface]++;
#endif
            } else { // there is no wrap neighbour AND there is no large element -> I am the large element
                // large element -> rightinterface
                TPZGeoElBC gbc(gelside,rightinterface);
                int64_t int_id = gbc.CreatedElement()->Index();
                int64_t wrap_id = gel->Index();
                
                // we create the lagrange geometric element if there is no boundary element as neighbour
                if(!hasbcNeighbour) { // if there is a BC element, it will act as lagrange element
                    TPZGeoElBC gbcLagrange(gbc.CreatedElement(),fHybridizationData.fLagrangeMatId);
                    int64_t lagr_id = gbc.CreatedElement()->Index();
                    TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                    fHybridizationData.fInterfaces.Push(hybrid);
                } else {
                    int64_t lagr_id = hasbcNeighbour.Element()->Index();
                    TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                    fHybridizationData.fInterfaces.Push(hybrid);
                }
#ifdef PZ_LOG
                numcreated[leftinterface]++;
                if(!hasbcNeighbour) {
                    numcreated[fHybridizationData.fLagrangeMatId]++;
                }
#endif
            }
        }
    }
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << __PRETTY_FUNCTION__ << std::endl;
        for(auto it : numcreated)
        {
            sout << "Material id " << it.first << " number of elements created " << it.second << std::endl;
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef PZDEBUG
    {
//        std::ofstream out("gmesh.txt");
//        fGeoMesh->Print(out);
    }
    // Check if neighborhood is correct
    CheckNeighborhoodHybridization();
#endif
}

void TPZApproxCreator::CheckNeighborhoodHybridization() const {
    const int dim = fGeoMesh->Dimension();
    for(auto gel : fGeoMesh->ElementVec()) {
        if(gel->Dimension() != dim || gel->HasSubElement()) continue;
        
        const int firstFace = gel->FirstSide(dim-1);
        for (int iside = firstFace; iside < gel->NSides()-1; iside++) {
            // verify if the first neighbour is matwrap
            TPZGeoElSide gelside(gel,iside);
            if (!gelside.HasNeighbour(fHybridizationData.fWrapMatId)) continue;
            
            TPZGeoElSide gelsidewrap = gelside.Neighbour();
            TPZGeoEl* gelWrap = gelsidewrap.Element();
            if (gelWrap->MaterialId() != fHybridizationData.fWrapMatId) {
                DebugStop(); // Wrap has to be the first neighbor
            }
            
            TPZGeoElSide neighwrap = gelsidewrap.HasNeighbour(fHybridizationData.fWrapMatId);
            if(neighwrap == gelsidewrap) continue;
            
            // must have an interface element of same size
            TPZGeoEl* gelInter = gelside.Neighbour().Neighbour().Element();
            if (gelInter->MaterialId() != fHybridizationData.fLeftInterfaceMatId &&
                gelInter->MaterialId() != fHybridizationData.fRightInterfaceMatId) {
                DebugStop(); // Wrap has to be the first neighbor
            }
            
            // Now check if the wrap is the correct first neighbor by checking if the wrap nodes are contained in the element itself
            std::set<int64_t> wrapNodes,geoNodes;
            const int gelWrapNSideNodes = gelWrap->NSideNodes(gelWrap->NSides()-1),
                gelSideNSideNodes = gel->NSideNodes(iside);
            
            if(gelWrapNSideNodes != gelSideNSideNodes) DebugStop();
            
            for(int in = 0 ; in < gelWrap->NSideNodes(gelWrap->NSides()-1) ; in++) {
                wrapNodes.insert(gelWrap->SideNodeIndex(gelWrap->NSides()-1, in));
                geoNodes.insert(gel->SideNodeIndex(iside, in));
            }
            std::set<int> diffset;
            std::set_difference(wrapNodes.begin(), wrapNodes.end(), geoNodes.begin(), geoNodes.end(),inserter(diffset,diffset.begin()));
            
            if(diffset.size()) DebugStop();
            
        }
    }
//    std::cout << "\nElement index " << neighgel->Reference()->Index() << " matid " << neighgel->MaterialId();
}

void TPZApproxCreator::InsertWrapAndLagrangeMaterialObjects(TPZMultiphysicsCompMesh *mphys){
    int periphelDim = mphys->Dimension()-1;
    int nstate = NState();
    auto Fwd = [](int matid, int dim, int nState,  TPZMultiphysicsCompMesh *cmesh){
        TPZMaterial *mat = cmesh->FindMaterial(matid);
        if(mat) return;
        auto nullmat = new TPZNullMaterialCS(matid);
        nullmat->SetDimension(dim);
        nullmat->SetNStateVariables(nState);
        cmesh->InsertMaterialObject(nullmat);
    };

    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi)
    {
        Fwd(fHybridizationData.fLagrangeMatId,periphelDim,nstate,mphys);
        Fwd(fHybridizationData.fWrapMatId,periphelDim,nstate,mphys);
    }
    else if (fHybridType == HybridizationType::EStandardSquared) {
        Fwd(fHybridizationData.fLagrangeMatId,periphelDim,nstate,mphys);
        Fwd(fHybridizationData.fWrapMatId,periphelDim,nstate,mphys);
        Fwd(fHybridizationData.fSecondLagrangeMatId,periphelDim,nstate,mphys);
    }
}

void TPZApproxCreator::InsertInterfaceMaterialObjects(TPZMultiphysicsCompMesh *mphys){
    if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::EStandardSquared  || fHybridType == HybridizationType::ESemi) {

        int nstate = NState();
                
        //TODO Should we support more general interface classes? Support CSTATE?
        TPZLagrangeMultiplierCS<STATE> *interface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fLeftInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
        if(fProbType == ProblemType::EElastic){
            interface->SetMultiplier(-1.);
        }
        mphys->InsertMaterialObject(interface);

        if(fHybridizationData.fLeftInterfaceMatId != fHybridizationData.fRightInterfaceMatId) {
            TPZLagrangeMultiplierCS<STATE> *interface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fRightInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            if(fProbType == ProblemType::EElastic){
                interface->SetMultiplier(-1.);
            }
            mphys->InsertMaterialObject(interface);
        }
        
        if (fHybridType == HybridizationType::EStandardSquared) {
            TPZLagrangeMultiplierCS<STATE> *secondLeftInterface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fSecondLeftInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            secondLeftInterface->SetMultiplier(-1.);
            mphys->InsertMaterialObject(secondLeftInterface);

            TPZLagrangeMultiplierCS<STATE> *secondRightInterface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fSecondRightInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            secondRightInterface->SetMultiplier(-1.);
            mphys->InsertMaterialObject(secondRightInterface);

        }
    }
    else{
        std::cout << __PRETTY_FUNCTION__ << "Should this be called\n?";
        DebugStop();
    }
}

std::set<int> TPZApproxCreator::GetBCMatIds(){
    std::set<int> bcMatIdsVec;
    for(std::pair<int,TPZMaterial*> mat: fMaterialVec){
        TPZBndCond *bcmat = dynamic_cast<TPZBndCond*>(mat.second);
        if(bcmat){
            bcMatIdsVec.insert(mat.first);
        }
    }
    return bcMatIdsVec;
}

std::set<int> TPZApproxCreator::GetVolumetricMatIds(){
    const int domaindim = fGeoMesh->Dimension();
    std::set<int> volumeMatIds;
    for(std::pair<int,TPZMaterial*> mat: fMaterialVec){
        const int matdim = mat.second->Dimension();
        TPZBndCond *bcmat = dynamic_cast<TPZBndCond*>(mat.second);
        if(!bcmat && domaindim == matdim){
            volumeMatIds.insert(mat.first);
        }
    }
    return volumeMatIds;
}

void TPZApproxCreator::Print(std::ostream &ofs){
    std::stringstream ss;
    ss << "\n";
    int numStars = 70;
    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    std::string hybridization = [](HybridizationType hybtype){
        switch (hybtype){
            case HybridizationType::ENone:
                return  "ENone";
            case HybridizationType::EStandard:
                return "EStandard";
            case HybridizationType::EStandardSquared:
                return "EStandardSquared";
            case HybridizationType::ESemi:
                return "ESemi";
            default:
                DebugStop();
                return "garbage";//avoids compiler warning on gcc
        }
    }(this->fHybridType);

    std::string problem = [](ProblemType probtype){
        switch (probtype){
            case ProblemType::ENone:
                return  "ENone";
            case ProblemType::EElastic:
                return "EElastic";
            case ProblemType::EDarcy:
                return "EDarcy";
            case ProblemType::EStokes:
                return "EStokes";
            default:
                DebugStop();
                return "garbage";//avoids compiler warning on gcc
        }
    }(this->fProbType);

    ss << __PRETTY_FUNCTION__ << "\nApproximation space creation set up:\n\n";
    ss << "fHybridType:\t" << hybridization << "fProbType:\t" << problem;
    ss <<"\nExtra Internal Pressure Order: " << fExtraInternalPOrder;
    ss << "\nDefault Pressure Order:\t" << fDefaultPOrder;
    if(this->fHybridType != HybridizationType::ENone) {
        ss << "\nfHybridizeBCLevel: " << fHybridizationData.fHybridizeBCLevel;
    }
    ss << "\nfShouldCondense:\t" << fShouldCondense;
    ss << "\nfIsRBSpaces" << fIsRBSpaces;
    ss << "\n\nmatids:\n\t{";
    for(auto it : fMaterialVec){
        ss << it.first <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    if(this->fHybridType == HybridizationType::EStandard || this->fHybridType == HybridizationType::EStandardSquared) {
        ss << "peripheralMatids:\n";
        ss << "\tfWrapMatId = " << fHybridizationData.fWrapMatId;
        ss << "\n\tfLeftInterfaceMatId = " << fHybridizationData.fLeftInterfaceMatId;
        ss << "\n\tfRightInterfaceMatId = " << fHybridizationData.fRightInterfaceMatId;
        ss << "\n\tfLagrangeMatId = " << fHybridizationData.fLagrangeMatId;
        if( this->fHybridType == HybridizationType::EStandardSquared) {
            ss << "\n\tfSecondLeftInterfaceMatId = " << fHybridizationData.fSecondLeftInterfaceMatId;
            ss << "\n\tfSecondRightInterfaceMatId = " << fHybridizationData.fSecondRightInterfaceMatId;
            ss << "\n\tfSecondLagrangeMatId = " << fHybridizationData.fSecondLagrangeMatId;
        }
        ss << "\n";
    }

    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    ofs << ss.str();
}

void TPZApproxCreator::ChangeInternalOrder(TPZCompMesh *cmesh, int pOrder) const {

    int64_t nel = cmesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZCompEl *cel = cmesh->Element(el);
        if (!cel) continue;

        TPZGeoEl *gel = cel->Reference();
        if (gel->Dimension() != cmesh->Dimension()) {//Only elements with the same dimension of the mesh
            continue;
        }
        int nc = cel->NConnects();
        //Gets the volumetric connect
        int64_t conIndex = cel->ConnectIndex(nc-1);
        TPZInterpolatedElement *intel = dynamic_cast<TPZInterpolatedElement *> (cel);
        if (!intel) DebugStop();

        TPZConnect &c = cmesh->ConnectVec()[conIndex];
        int64_t index;
        //Sets the new connect order
        c.SetOrder(pOrder,index);

        //Gets connect information to update block size (stiffness matrix data structure)
        int64_t seqnum = c.SequenceNumber();
        int nvar = 1;
        TPZMaterial * mat = cel->Material();
        if (mat) nvar = mat->NStateVariables();
        int nshape = intel->NConnectShapeF(nc-1,pOrder);
        c.SetNShape(nshape);
        // c.SetNState(nvar);
        cmesh->Block().Set(seqnum, nvar * nshape);

        cel->SetIntegrationRule(2*pOrder);
    }
    cmesh->InitializeBlock();
}

void TPZApproxCreator::SetMeshElementType(){

    MElementType firstElement;
    const int dim = fGeoMesh->Dimension();
    bool isfirst = true;
    bool isDiffElMesh = false;

    for (int iEl = 0; iEl < fGeoMesh->NElements(); iEl++){
        
        TPZGeoEl *gel = fGeoMesh->ElementVec()[iEl];
        if (!gel) continue;
        if (gel->Dimension() != dim) continue;

        if (isfirst){
            firstElement = gel->Type();
            isfirst = false;
        } else {
            if (firstElement != gel->Type()){
                isDiffElMesh = true;
                if (fExtraInternalPOrder != 0){
                    std::cout << "Error! Please implemet extra internal order for different element types" << std::endl;
                    DebugStop();
                };
            }
        }
    }

    if (!isDiffElMesh){
        fElementType = firstElement;
    }   

}
