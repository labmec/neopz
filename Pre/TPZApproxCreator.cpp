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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("CreateMultiphysicsSpace"));
#endif

TPZApproxCreator::TPZApproxCreator(TPZGeoMesh *gmesh) : fGeoMesh(gmesh)
{
}

/**Insert a material object in the datastructure*/
int TPZApproxCreator::InsertMaterialObject(TPZMaterial * mat) {
    if(!mat) DebugStop();
    const int matid = mat->Id();
    if(fMaterialVec.find(matid) != fMaterialVec.end()) DebugStop();
    fMaterialVec.insert(std::make_pair(matid, mat));
    return fMaterialVec.size();
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
    fHybridizationData.fInterfaceMatId = matid_base+base;
    fHybridizationData.fLagrangeMatId = matid_base + 2*base;
    fHybridizationData.fSecondInterfaceMatId = matid_base + 3*base;
    fHybridizationData.fSecondLagrangeMatId = matid_base + 4*base;
}

void TPZApproxCreator::AddHybridizationGeoElements(){
    /// For fHybridType == HybridizationType::EStandard, wrap, interface and lagrange geometric elements are created;
    ///     if fHybridizeBCLevel == 0, the boundary is not hybridized, if fHybridizeBCLevel == 1, it is.
    /// For fHybridType == HybridizationType::EStandardSquared, wrap, interface, lagrange, second interface and second lagrange
    /// geometric elements are created;
    ///     BC can be hybridized 0, 1 or 2 times in this set up;
    /// fHybridType == HybridizationType::Semi has yet to be implemented.

    if(fHybridizationData.fWrapMatId == -123456){
        std::cout << "\nERROR! Please call TPZApproxCreator::ComputePeriferalMaterialIds() before TPZApproxCreator::AddHybridizationGeoElements()" << std::endl;
        DebugStop();
    }

#ifdef LOG4CXX
    std::map<int,int> numcreated;
#endif
    int64_t nel = fGeoMesh->NElements();
    int dim = fGeoMesh->Dimension();
    std::set<int> bcMatIds = GetBCMatIds();
    // Wrap and interface creation for standard hyb., wrap, interface and lag. creation for double hyb.
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim) continue;
        int nsides = gel->NSides();
        int side = 0;
        for(int d=0; d<dim-1; d++) side += gel->NSides(d);
        // loop over the sides of dimension dim-1
        for(; side < nsides-1; side++)
        {
            TPZGeoElSide gelside(gel,side);
            // we want to create side elements of type
            // first fMatWrapId
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighMatId = neighbour.Element()->MaterialId();
#ifdef PZDEBUG
            {
                if(neighMatId == fHybridizationData.fWrapMatId)
                {
                    std::cout << __PRETTY_FUNCTION__ << " should be called only once!\n";
                    DebugStop();
                }
            }
#endif
            // if the neighbour is a boundary condition and no hybridization is applied
            // do not create the wrap layers
            bool HasBCNeighbour = (bcMatIds.find(neighMatId) != bcMatIds.end());
            if(fHybridizationData.fHybridizeBCLevel == 0 && HasBCNeighbour)
            {
                // no interface will be created between the element and a flux space
                continue;
            }
            // first fH1Hybrid.fMatWrapId
            TPZGeoElBC(gelside, fHybridizationData.fWrapMatId);
            neighbour = gelside.Neighbour();
#ifdef LOG4CXX
            numcreated[fH1Hybrid.fMatWrapId]++;
            if(neighbour.Element()->MaterialId() != fH1Hybrid.fMatWrapId)
            {
                DebugStop();
            }
#endif
            if(fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi)
            {
                // then, depending on orientation fH1Hybrid.fLagrangeMatid.first or second
                int dir = gel->NormalOrientation(side);
                TPZGeoElBC(neighbour,fHybridizationData.fInterfaceMatId);
#ifdef LOG4CXX
                numcreated[fHybridizationData.fInterfaceMatId]++;
#endif
            } else if(fHybridType == HybridizationType::EStandardSquared) {
                TPZGeoElBC(neighbour,fHybridizationData.fInterfaceMatId);
#ifdef LOG4CXX
                numcreated[fH1Hybrid.fLagrangeMatid.first]++;
#endif
                neighbour = neighbour.Neighbour();
                if(!HasBCNeighbour || fHybridizationData.fHybridizeBCLevel == 2)
                {
                    TPZGeoElBC(neighbour,fHybridizationData.fLagrangeMatId);
#ifdef LOG4CXX
                    numcreated[fH1Hybrid.fFluxMatId]++;
#endif
                }
            } else {
                // only two cases handled so far
                DebugStop();
            }
        }
    }
    nel = fGeoMesh->NElements();
    for(int64_t el = 0; el<nel; el++)
    {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel || gel->HasSubElement() || gel->Dimension() != dim-1) continue;
        int matid = gel->MaterialId();
        TPZGeoElSide gelside(gel,gel->NSides()-1);
        TPZGeoElSide neighbour(gelside.Neighbour());
        // if the neighbour is a boundary condition and no hybridization is applied
        // do not create the wrap layers
        int neighmat = neighbour.Element()->MaterialId();
        bool HasBCNeighbour = (bcMatIds.find(neighmat) != bcMatIds.end());

        int interface = fHybridizationData.fInterfaceMatId;
        bool isinterface = (matid == interface);
        bool islag = (matid == fHybridizationData.fLagrangeMatId);
        if((fHybridType == HybridizationType::EStandard || fHybridType == HybridizationType::ESemi) && isinterface)
        {
            bool ShouldCreateLag =  !HasBCNeighbour;
            if(gelside.HasLowerLevelNeighbour(fHybridizationData.fLagrangeMatId))
            {
                // create the flux element
                //                TPZGeoElBC(gelside,fH1Hybrid.fFluxMatId);
#ifdef LOG4CXX
                //                numcreated[fH1Hybrid.fFluxMatId]++;
#endif
            }
            else if(ShouldCreateLag && !gelside.HasNeighbour(fHybridizationData.fLagrangeMatId))
            {
                // create the flux element
                TPZGeoElBC(gelside,fHybridizationData.fLagrangeMatId);
#ifdef LOG4CXX
                numcreated[fH1Hybrid.fFluxMatId]++;
#endif
            }
        }
        if(fHybridType == HybridizationType::EStandardSquared && islag)
        {
            // we need to include three layers of geometric elements
            // if there is fH1Hybrid.fInterfacePressure element as neighbour dont create
            if(gelside.HasNeighbour(fHybridizationData.fSecondLagrangeMatId)) continue;
            // if there is a higher level connected element, dont create
            TPZGeoElSidePartition partition(gelside);
            if(partition.HasHigherLevelNeighbour(fHybridizationData.fLagrangeMatId)) continue;
            // now we can create the elements
            TPZGeoElSide neighbour = gelside.Neighbour();
            int neighmat = neighbour.Element()->MaterialId();
            bool HasBCNeighbour = (bcMatIds.find(neighmat) != bcMatIds.end());

            TPZGeoElBC gbc1(gelside, fHybridizationData.fSecondInterfaceMatId);
            // if the flux is neighbour, the interface will be between lag and the
            // boundary condition
#ifdef LOG4CXX
            numcreated[fH1Hybrid.fSecondLagrangeMatid]++;
#endif
            if(HasBCNeighbour) continue;
            neighbour = gelside.Neighbour();
            // pressure element
            TPZGeoElBC gbc2(neighbour,fHybridizationData.fSecondLagrangeMatId);
            neighbour = neighbour.Neighbour();
            // lagrange element
            TPZGeoElBC gbc3(neighbour,fHybridizationData.fSecondInterfaceMatId);
#ifdef LOG4CXX
            numcreated[fH1Hybrid.fInterfacePressure]++;
            numcreated[fH1Hybrid.fSecondLagrangeMatid]++;
#endif
        }
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
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

}

void TPZApproxCreator::InsertWrapAndLagrangeMaterialObjects(TPZMultiphysicsCompMesh *mphys){
    int periphelDim = mphys->Dimension()-1;
    int nstate = 1;
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

        int nstate = 0;
        for (std::pair<int,TPZMaterial*> matpair : fMaterialVec) {
            TPZMaterial* mat = matpair.second;
            TPZBndCond *bnd = dynamic_cast<TPZBndCond *> (mat);
            if (!bnd){
                nstate = mat->NStateVariables(); // here we assume that all materials have same nstatevars
                break;
            }
        }
        if(nstate < 1) DebugStop();
                
        //TODO Should we support more general interface classes? Support CSTATE?
        TPZLagrangeMultiplierCS<STATE> *interface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
        if(fProbType == ProblemType::EElastic){
            interface->SetMultiplier(-1.);
        }
        mphys->InsertMaterialObject(interface);

        if (fHybridType == HybridizationType::EStandardSquared) {

            //TODO Should we support more general interface classes? Support CSTATE?
            TPZLagrangeMultiplierCS<STATE> *secondInterface = new TPZLagrangeMultiplierCS<STATE>(fHybridizationData.fSecondInterfaceMatId, fGeoMesh->Dimension()-1, nstate);
            secondInterface->SetMultiplier(-1.);
            mphys->InsertMaterialObject(secondInterface);

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
        ss << "\n\tfInterfaceMatId = " << fHybridizationData.fInterfaceMatId;
        ss << "\n\tfLagrangeMatId = " << fHybridizationData.fLagrangeMatId;
        if( this->fHybridType == HybridizationType::EStandardSquared) {
            ss << "\n\tfSecondInterfaceMatId = " << fHybridizationData.fSecondInterfaceMatId;
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
