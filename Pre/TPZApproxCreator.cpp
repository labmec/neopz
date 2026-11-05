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
static TPZLogger logger("TPZApproxCreator");
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

    if(fHybridType == HybridizationType::ENone) {
        std::cout << "Hybridization type is none, this method should not be called\n";
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
                int64_t wrap_id = el;
                int64_t int_id = gbc.CreatedElement()->Index();
                TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                AddHybridGeometry(int_id,hybrid);
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
                AddHybridGeometry(int_id,hybrid);
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
                if(hasLagrangeElement) DebugStop();
                TPZGeoElSide lagrange = hasLargeElementNeighbour.HasNeighbour(fHybridizationData.fLagrangeMatId);
                if(!lagrange) {
                    TPZGeoElBC gbcLagrange(hasLargeElementNeighbour,fHybridizationData.fLagrangeMatId);
                    lagrange = TPZGeoElSide(gbcLagrange.CreatedElement());
                }
                TPZGeoElBC gbc(gelside,leftinterface);
                int64_t lagr_id = lagrange.Element()->Index();
                int64_t int_id = gbc.CreatedElement()->Index();
                int64_t wrap_id = gel->Index();
                TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                AddHybridGeometry(int_id,hybrid);
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
                    int64_t lagr_id = -1;
                    if(hasLagrangeElement) {
                        lagr_id = hasLagrangeElement.Element()->Index();
                    } else {
                        TPZGeoElBC gbcLagrange(gbc.CreatedElement(),fHybridizationData.fLagrangeMatId);
                        lagr_id = gbcLagrange.CreatedElement()->Index();
                    }
                    TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                    AddHybridGeometry(int_id,hybrid);
                } else { // connected to a boundary element
                    int64_t lagr_id = hasbcNeighbour.Element()->Index();
                    TPZHybrid hybrid(int_id,wrap_id,lagr_id);
                    AddHybridGeometry(int_id,hybrid);
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
  if(fHybridType == HybridizationType::EStandardSquared) {
    this->AddHybridSquareGeoElements();
  }
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
        interface->SetMultiplier(fHybridizationData.fMultipliers[0]);
        mphys->InsertMaterialObject(interface);

        if(fHybridizationData.fLeftInterfaceMatId != fHybridizationData.fRightInterfaceMatId) {
            TPZLagrangeMultiplierCS<STATE> *interface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fRightInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            interface->SetMultiplier(fHybridizationData.fMultipliers[1]);
            mphys->InsertMaterialObject(interface);
        }
        
        if (fHybridType == HybridizationType::EStandardSquared) {
            TPZLagrangeMultiplierCS<STATE> *secondLeftInterface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fSecondLeftInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            secondLeftInterface->SetMultiplier(fHybridizationData.fMultipliers[2]);
            mphys->InsertMaterialObject(secondLeftInterface);

            TPZLagrangeMultiplierCS<STATE> *secondRightInterface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fSecondRightInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            secondRightInterface->SetMultiplier(fHybridizationData.fMultipliers[3]);
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

void TPZApproxCreator::HybridizationData::Print(std::ostream &out) const {
    out << "peripheralMatids:\n";
    out << "\tfWrapMatId = " << fWrapMatId << "\n";
    out << "\tfLeftInterfaceMatId = " << fLeftInterfaceMatId << "\n";
    out << "\tfRightInterfaceMatId = " << fRightInterfaceMatId << "\n";
    out << "\tfLagrangeMatId = " << fLagrangeMatId << "\n";
    out << "\tfSecondLeftInterfaceMatId = " << fSecondLeftInterfaceMatId << "\n";
    out << "\tfSecondRightInterfaceMatId = " << fSecondRightInterfaceMatId << "\n";
    out << "\tfSecondLagrangeMatId = " << fSecondLagrangeMatId << "\n";
    out << "\tfHybridizeBCLevel = " << fHybridizeBCLevel << "\n";
    out << "\tfMultipliers = {" << fMultipliers[0] << ", " << fMultipliers[1] << ", " << fMultipliers[2] << ", " << fMultipliers[3] << "}\n";
    out << "\tfInterfaces:\n";
    for (auto &it : fInterfaces) {
        out << "\t\tinterface " << it.first << ": (lagrange=" << it.second.fLagrange << ", left=" << it.second.fLeft << ", right=" << it.second.fRight << ")\n";
    }
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
        fHybridizationData.Print(ss);
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


/// Add a hybrid geometric configuration
void TPZApproxCreator::AddHybridGeometry(int64_t interface_id, TPZHybrid &hybrid) {
//    std::cout << "AddHybridG adding interface " << interface_id << " left " << hybrid.fLeft << " right " << hybrid.fRight << std::endl;
#ifdef PZDEBUG
    if(fHybridizationData.fInterfaces.find(interface_id) != fHybridizationData.fInterfaces.end()) {
        DebugStop();
    }
#endif
    fHybridizationData.fInterfaces[interface_id] = hybrid;
}

/// Add a hybrid geometric configuration
static void AddHybrid(std::map<int64_t,TPZHybrid> &tree, int64_t interface_id, TPZHybrid &hybrid) {
//    std::cout << "Hybrid adding interface " << interface_id << " left " << hybrid.fLeft << " right " << hybrid.fRight << std::endl;
#ifdef PZDEBUG
    if(tree.find(interface_id) != tree.end()) {
        DebugStop();
    }
#endif
    tree[interface_id] = hybrid;
}



// the right interface is created next to a large element
/// Add geometric elements to represent squared hybridization
void TPZApproxCreator::AddHybridSquareGeoElements() {
  if(fHybridType != HybridizationType::EStandardSquared) DebugStop();
    if(0) {
        std::cout << "Hybridization structure before hybrid square\n";
        for(auto &it : fHybridizationData.fInterfaces) {
            TPZHybrid &hybridlarge = it.second;
            std::cout << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
        }
    }
#ifdef PZDEBUG
    std::set<int64_t> largeflux;
#endif


    std::map<int64_t,TPZHybrid> hybridsq;
    int bchybridlevel = fHybridizationData.fHybridizeBCLevel;
    for(auto &it : fHybridizationData.fInterfaces) {
        TPZGeoEl *gelinterface = fGeoMesh->Element(it.first);
        TPZGeoEl *gelwrap = fGeoMesh->Element(it.second.fLeft);
        if(gelwrap->MaterialId() != fHybridizationData.fWrapMatId) DebugStop();
        TPZGeoElSide wrapside(gelwrap);
        
        bool hasLargeElementNeighbour = wrapside.HasLowerLevelNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide wrapneighbour = wrapside.Neighbour().HasNeighbour(fHybridizationData.fWrapMatId);
        TPZGeoElSide lagrangeneighbour = wrapside.HasNeighbour(fHybridizationData.fLagrangeMatId);
        bool haswrapNeighbour = (wrapneighbour != wrapside);
        bool hasbcNeighbour = wrapside.HasNeighbour(GetBCMatIds());
        bool isLargeElement = false;
        if(!haswrapNeighbour && !hasLargeElementNeighbour && !hasbcNeighbour) {
            // we are in the case where there is no wrap neighbour, no large element neighbour and no bc neighbour, so we are the large element and we need to create the second interface and second lagrange elements
            isLargeElement = true;
        }
        // find the element configuration :
        // - it is an equal level element
        // - it is a small element linked to a large element
        // - it is a large element
        if(isLargeElement) {
          // copy the TPZHybrid object
          TPZHybrid hybridlarge =(it.second);
          AddHybrid(hybridsq,it.first,hybridlarge);
#ifdef PZDEBUG
            largeflux.insert(it.second.fRight);
#endif
          continue;
        }
        if(hasbcNeighbour) {
          // we do not doubly hybridize boundary conditions
          // copy the TPZHybrid object
          TPZHybrid bchybrid(it.second);
          AddHybrid(hybridsq,it.first,bchybrid);
          continue;
        }
        if(haswrapNeighbour) {
          TPZManVector<int64_t> geoelem(5,-1);
          geoelem[0] = gelwrap->Index();
          geoelem[1] = gelinterface->Index();
          int interfacematid = gelinterface->MaterialId();
          int secondinterfacematid = 0;
          if(interfacematid == fHybridizationData.fLeftInterfaceMatId) {
            secondinterfacematid = fHybridizationData.fSecondLeftInterfaceMatId;
          } else if(interfacematid == fHybridizationData.fRightInterfaceMatId) {
            secondinterfacematid = fHybridizationData.fSecondRightInterfaceMatId;
          } else {
            DebugStop();
          }
          if(interfacematid == fHybridizationData.fLeftInterfaceMatId) {
            int nsides = gelinterface->NSides();
            auto newflux = gelinterface->CreateBCGeoEl(nsides-1,fHybridizationData.fLagrangeMatId);
            geoelem[2] = newflux->Index();
            lagrangeneighbour = TPZGeoElSide(newflux);      
          } else {
              int64_t lagrangeindex = it.second.fRight;
              auto gel = fGeoMesh->Element(lagrangeindex);
              if(gel->MaterialId() != fHybridizationData.fLagrangeMatId) DebugStop();
              lagrangeneighbour = TPZGeoElSide(gel);
              if(!lagrangeneighbour) DebugStop();
              geoelem[2] = lagrangeneighbour.Element()->Index();
          }
          {
            TPZGeoElBC gbc(lagrangeneighbour,secondinterfacematid);
            geoelem[3] = gbc.CreatedElement()->Index();
          }
          // we need to create a new interface element
          // create a secondlagrange if it doesnt exist
          TPZGeoElSide secondlagrange = wrapside.HasNeighbour(fHybridizationData.fSecondLagrangeMatId);
          if(!secondlagrange) {
            TPZGeoEl *gelinterface = fGeoMesh->Element(geoelem[3]);
            TPZGeoElSide gelint(gelinterface);
            int nsides = gelinterface->NSides();
            secondlagrange = TPZGeoElSide(gelinterface->CreateBCGeoEl(nsides-1,fHybridizationData.fSecondLagrangeMatId));
            geoelem[4] = secondlagrange.Element()->Index();
          } else {
            geoelem[4] = secondlagrange.Element()->Index();
          }
          TPZHybrid hybr1(geoelem[1],geoelem[0],geoelem[2]);
          TPZHybrid hybr2(geoelem[3],geoelem[2],geoelem[4]);
          AddHybrid(hybridsq,geoelem[1],hybr1);
          AddHybrid(hybridsq,geoelem[3],hybr2);
          // create a new lagrange element if a left interface
          // make a copy of the iterator
          // create a second TPZHybrid object
          // create a left/right second interface
          // the left element is the lagrange
          // the right element is the secondlagrange
          continue;
        }
        if(hasLargeElementNeighbour) {
          TPZManVector<int64_t,7> geoelem(7,-1);
          geoelem[0] = it.second.fLeft;
          geoelem[1] = it.first;
          geoelem[6] = it.second.fRight;
          TPZGeoEl *gelinterface = fGeoMesh->Element(it.first);
          if(gelinterface->MaterialId() != fHybridizationData.fLeftInterfaceMatId) {
            DebugStop();
          }
          int nsides = gelinterface->NSides();
          TPZGeoElBC gbclagrange(gelinterface,nsides-1,fHybridizationData.fLagrangeMatId);
          geoelem[2] = gbclagrange.CreatedElement()->Index();
          TPZGeoElBC gbcsecinterface(gbclagrange,fHybridizationData.fSecondLeftInterfaceMatId);
          geoelem[3] = gbcsecinterface.CreatedElement()->Index();
          TPZGeoElBC gbcsecondlagrange(gbcsecinterface,fHybridizationData.fSecondLagrangeMatId);
          geoelem[4] = gbcsecondlagrange.CreatedElement()->Index();
          TPZGeoElBC gbcsecondinterfacer(gbcsecondlagrange,fHybridizationData.fSecondRightInterfaceMatId);
          geoelem[5] = gbcsecondinterfacer.CreatedElement()->Index();
//            std::cout << "geoelem " << geoelem << std::endl;
          TPZHybrid h1(geoelem[1],geoelem[0],geoelem[2]);
          TPZHybrid h2(geoelem[3],geoelem[2],geoelem[4]);
          TPZHybrid h3(geoelem[5],geoelem[6],geoelem[4]);
          AddHybrid(hybridsq,geoelem[1],h1);
          AddHybrid(hybridsq,geoelem[3],h2);
          AddHybrid(hybridsq,geoelem[5],h3);
          continue;
        }

    }
#ifdef PZDEBUG
    {
        // number of times an element appears in a position (only once)
        std::set<int64_t> asleft, asright, aslagrange;
        // number of times an element is referenced (max 2)
        std::map<int64_t,int> numref;
        int error = 0;
        for(auto &it : hybridsq) {
            int64_t lag = it.first;
            int64_t left = it.second.fLeft;
            int64_t right = it.second.fRight;
            bool leftlarge = largeflux.find(left) != largeflux.end();
            bool rightlarge = largeflux.find(right) != largeflux.end();
            if(!leftlarge && asleft.find(left) != asleft.end()) {
                std::cout << "left " << left << " appears more than once\n";
                error++;
            }
            if(rightlarge && asright.find(right) != asright.end()) {
                TPZGeoEl *gel = fGeoMesh->Element(right);
                int matid = gel->MaterialId();
                if(matid != fHybridizationData.fSecondLagrangeMatId) {
                    std::cout << "right " << right << " appears more than once\n";
                    error++;
                }
            }
            if(aslagrange.find(lag) != aslagrange.end()) {
                std::cout << "lagrange " << lag << " appears more than once\n";
                error++;
            }
            asleft.insert(left);
            asright.insert(right);
            aslagrange.insert(lag);
            numref[left]++;
            numref[right]++;
            numref[lag]++;
        }
        for(auto &it : numref) {
            bool islarge = largeflux.find(it.first) != largeflux.end();
            if(!islarge && it.second > 2) {
                std::cout << "element " << it.first << " appears more than twice\n";
                error++;
            }
        }
        if(error > 0) {
            for(auto &it : hybridsq) {
                TPZHybrid &hybridlarge = it.second;
                std::cout << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
            }
            DebugStop();
        }
    }
#endif
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                std::stringstream out;
                for(auto &it : hybridsq) {
                    TPZHybrid &hybridlarge = it.second;
                    out << "lag " << it.first << " left " << hybridlarge.fLeft << " right " << hybridlarge.fRight << std::endl;
                }
                LOGPZ_DEBUG(logger, out.str())
            }
#endif
    fHybridizationData.fInterfaces = hybridsq;
}

/// compute the lagrange multiplier coeficients as a function of the problem type and hybridization
void TPZApproxCreator::HybridizationData::SetProblemHybridH1(ProblemType prob, HybridizationType hybrid) {
    if(prob == ProblemType::ENone || hybrid == HybridizationType::ENone) return;
    if(prob == ProblemType::EDarcy) {
        fMultipliers[0] = 1.;
        fMultipliers[1] = 1.;
        fMultipliers[2] = 1.;
        fMultipliers[3] = -1.;
    } else if(prob == ProblemType::EElastic) {
        fMultipliers[0] = -1.;
        fMultipliers[1] = -1.;
        fMultipliers[2] = -1.;
        fMultipliers[3] = 1.;
    } else {
        DebugStop();
    }
}

/// compute the lagrange multiplier coeficients as a function of the problem type and hybridization
void TPZApproxCreator::HybridizationData::SetProblemHybridHDiv(ProblemType prob, HybridizationType hybrid) {
    if(prob == ProblemType::ENone || hybrid == HybridizationType::ENone) return;
    if(prob == ProblemType::EDarcy){
        if(hybrid == HybridizationType::EStandard) {
            fMultipliers[0] = 1.;
            fMultipliers[1] = 1.;
            fMultipliers[2] = 1.;
            fMultipliers[3] = -1.;
        } else if (hybrid == HybridizationType::ESemi) {
            fMultipliers[0] = 1.;
            fMultipliers[1] = -1.;
            fMultipliers[2] = 1.;
            fMultipliers[3] = -1.;
        } else {
            DebugStop();
        }
    } else if(prob == ProblemType::EElastic) {
        if(hybrid == HybridizationType::EStandard) {
            fMultipliers[0] = -1.;
            fMultipliers[1] = -1.;
            fMultipliers[2] = -1.;
            fMultipliers[3] = 1.;
        } else if(hybrid == HybridizationType::ESemi) {
            fMultipliers[0] = -1.;
            fMultipliers[1] = 1.;
            fMultipliers[2] = -1.;
            fMultipliers[3] = 1.;
        }

    } else {
        DebugStop();
    }
}
