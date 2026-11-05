//
//  TPZBuildSBFem.cpp
//  PZ
//
//  Created by Philippe Devloo on 06/01/17.
//
//

#include "TPZBuildSBFem.h"
#include "TPZSBFemVolume.h"
#include "TPZSBFemElementGroup.h"
#include "pzcompel.h"
#include "tpzgeoblend.h"
#include "TPZGeoLinear.h"

#include "tpzgeoelrefpattern.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.tpzbuildsbfem");
#endif

TPZBuildSBFem::TPZBuildSBFem(const TPZBuildSBFem &copy) : fGMesh(copy.fGMesh), fMatIdTranslation(copy.fMatIdTranslation),
    fBoundaryMatIds(copy.fBoundaryMatIds), fSkeletonMatId(copy.fSkeletonMatId),
    fElementPartition(copy.fElementPartition), fPartitionCenterNode(copy.fPartitionCenterNode),
    fPOrderBubbleFunctions(copy.fPOrderBubbleFunctions), fSkeletonPOrder(copy.fSkeletonPOrder),
    fNRefSkeleton(copy.fNRefSkeleton)
{

}

/// standard configuration means each element is a partition and a center node is created
void TPZBuildSBFem::StandardConfiguration()
{
    TPZVec<int64_t> elindices(fGMesh->NElements());
    int64_t nel = elindices.size();
    for (int64_t el=0; el<nel; el++) {
        elindices[el] = el;
    }
    StandardConfiguration(elindices);
}

/// standard configuration means each element is a partition and a center node is created for the indicated elements
void TPZBuildSBFem::StandardConfiguration(TPZVec<int64_t> &elementindices)
{
    CreateElementCenterNodes(elementindices);
    AddSkeletonElements();
}

/// build element groups according to the id of the scaling centers
void TPZBuildSBFem::Configure(TPZVec<int64_t> &scalingcenters)
{
    std::map<int64_t,int64_t> nodetogroup;
    //    int maxpartition = 4;
    int64_t count = 0;
    for (int64_t el=0; el<scalingcenters.size(); el++) {
        if (scalingcenters[el] == -1) {
            continue;
        }
        if (nodetogroup.find(scalingcenters[el]) == nodetogroup.end()) {
            nodetogroup[scalingcenters[el]] = count++;
            //            if (count >= maxpartition) {
            //                break;
            //            }
        }
    }
    if(fPartitionCenterNode.size())
    {
        DebugStop();
    }
    int dim = fGMesh->Dimension();
    fPartitionCenterNode.resize(nodetogroup.size());
    for (std::map<int64_t,int64_t>::iterator it = nodetogroup.begin(); it != nodetogroup.end(); it++) {
        fPartitionCenterNode[it->second] = it->first;
    }
    int64_t nel = fGMesh->NElements();
    for (int64_t el=0; el < nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (gel->Dimension() != dim && scalingcenters[el] != -1) {
            DebugStop();
        }
        else if(gel->Dimension() != dim)
        {
            continue;
        }
        if (scalingcenters[el] == -1) {
            continue;
        }
        if (nodetogroup.find(scalingcenters[el]) == nodetogroup.end()) {
            DebugStop();
        }
        int64_t partition = nodetogroup[scalingcenters[el]];
        //        if(partition < maxpartition)
        {
            fElementPartition[el] = partition;
        }
    }
    AddSkeletonElements();
    
}


/// add a partition manually
int64_t TPZBuildSBFem::AddPartition(TPZVec<int64_t> &elindices, int64_t centernodeindex)
{
    {
        int64_t nel = fGMesh->NElements();
        if(fElementPartition.size() < nel) {
            fElementPartition.Resize(nel, -1);
        }
    }
    int64_t npart = fPartitionCenterNode.size();
    fPartitionCenterNode.resize(npart+1);
    if (fGMesh->NodeVec().NElements() <= centernodeindex) {
        DebugStop();
    }
    fPartitionCenterNode[npart] = centernodeindex;
    int64_t nel = elindices.size();
    for (int64_t el=0; el<nel; el++) {
        int64_t elindex = elindices[el];
        if (fElementPartition[elindex] != -1 || fElementPartition[elindex] >= npart) {
            DebugStop();
        }
        fElementPartition[elindex] = npart;
    }
    return npart;
}

/// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFem::BuildComputationMesh(TPZCompMesh &cmesh)
{
    // create the lower dimensional mesh
    std::set<int> matids;
    std::set<int> sbfemmatids = GetMaterialIds();
    int dim = cmesh.Dimension();
    TPZGeoMesh *gmesh = cmesh.Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        int gelmatid = gel->MaterialId();
        if(sbfemmatids.find(gelmatid) != sbfemmatids.end()) continue;
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    // create the boundary elements
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.SetDefaultOrder(fSkeletonPOrder);
    cmesh.AutoBuild(matids);
    // at this point all elements of lower dimension have been created
    CreateVolumetricElements(cmesh);
    CreateElementGroups(cmesh);
    
}

/// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFem::BuildComputationalMeshFromSkeleton(TPZCompMesh &cmesh)
{
    // create the lower dimensional mesh
    std::set<int> matids;
    int dim = cmesh.Dimension();
    TPZGeoMesh *gmesh = cmesh.Reference();
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() < dim) {
            matids.insert(gel->MaterialId());
        }
    }
    // create the boundary elements
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.AutoBuild(matids);
    CreateVolumetricElementsFromSkeleton(cmesh);
    CreateElementGroups(cmesh);
    
}

/// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFem::BuildComputationMesh(TPZCompMesh &cmesh, const std::set<int> &volmatids, const std::set<int> &boundmatids)
{
    // create the boundary elements
    cmesh.ApproxSpace().SetAllCreateFunctionsContinuous();
    cmesh.AutoBuild(boundmatids);
    CreateVolumetricElements(cmesh,volmatids);
    CreateElementGroups(cmesh);
    
}

static void DivideToFinest(TPZGeoEl *gelskel) {
    if(gelskel->HasSubElement()) return;
    int matid = gelskel->MaterialId();
    TPZGeoElSide gelside(gelskel);
    TPZGeoElSide neighbour = gelside.Neighbour();
    while(neighbour != gelside) {
        if(neighbour.HasSubElement()) {
            TPZStack<TPZGeoEl *> subels;
            gelskel->Divide(subels);
            for(int i=0; i<subels.NElements(); i++) {
                DivideToFinest(subels[i]);
            }
            break;
        }
        neighbour = neighbour.Neighbour();
    }
}


/// create the geometric skeleton elements
void TPZBuildSBFem::AddSkeletonElements()
{
    // create a lower dimension element on each boundary
    int dim = fGMesh->Dimension();
    
    int64_t nel = fGMesh->NElements();
    // first create skeleton elements at the level of the root mesh
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh->Element(el);
        if (!gel) {
            continue;
        }
        if (gel->Dimension() != dim) {
            continue;
        }
        // the element doesnt belong to any partition, do not create a skeleton element
        if (fElementPartition[el] == -1) {
            continue;
        }
        if(gel->HasSubElement()) {
            continue;
        }
        int64_t elpartition = fElementPartition[el];
        int nsides = gel->NSides();
        for (int is=gel->FirstSide(dim-1); is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            TPZGeoElSide thisside(gel,is);
            // we do not create skeleton elements on the boundary of the domain without boundary condition
            TPZGeoElSide neighbour = thisside.Neighbour();
            if (neighbour == thisside) {
                continue;
            }
            int64_t neighbourelpartition = -1;
            // look for a neighbour with mesh dimension
            while(neighbour != thisside && neighbour.Element()->Dimension() != dim)
            {
                neighbour = neighbour.Neighbour();
            }
            // we found one!
            if(neighbour != thisside)
            {
                neighbourelpartition = fElementPartition[neighbour.Element()->Index()];
            } else {
                // if there is no neighbour of same dimension
                // -> my side is fine (the neighbour will create the skeleton and divide it
                continue;
            }
            if(thisside.HasNeighbour(this->fBoundaryMatIds))
            {
                continue;
            }
            if(thisside.HasNeighbour(fSkeletonMatId))
            {
                continue;
            }

            if (elpartition != neighbourelpartition) {
                gel->CreateBCGeoEl(is,fSkeletonMatId);
            }
        }
        /// divide the skeleton elements until they have no divided neighbour
        nel = fGMesh->NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = fGMesh->Element(el);
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            if(gel->MaterialId() != fSkeletonMatId) continue;
            DivideToFinest(gel);
        }

    }
}

/// create a geometric node at the center of each partition
void TPZBuildSBFem::CreateElementCenterNodes(TPZVec<int64_t> &elindices)
{
    if(fPartitionCenterNode.size())
    {
        DebugStop();
    }
    int dim = fGMesh->Dimension();
    int64_t nel = elindices.size();
    fPartitionCenterNode.resize(nel);
    int64_t count = 0;
    for (int64_t el=0; el<elindices.size(); el++) {
        TPZGeoEl *gel = fGMesh->Element(elindices[el]);
        if (!gel || gel->Dimension() != dim || gel->HasSubElement()) {
            continue;
        }
        int nsides = gel->NSides();
        TPZManVector<REAL,3> xicenter(dim),xcenter(3);
        gel->CenterPoint(nsides-1,xicenter);
        gel->X(xicenter,xcenter);
        int64_t middlenode = fGMesh->NodeVec().AllocateNewElement();
        fGMesh->NodeVec()[middlenode].Initialize(xcenter,fGMesh);
        fPartitionCenterNode[count] = middlenode;
        fElementPartition[elindices[el]] = count;
        count++;
    }
    fPartitionCenterNode.resize(count);
}

/// create geometric volumetric elements
// the lower dimensional elements already exist (e.g. all connects have been created
void TPZBuildSBFem::GenerateCollapsedGeometricElements(const std::set<int> &matids)
{
    TPZGeoMesh *gmesh = fGMesh.operator->();
    int dim = gmesh->Dimension();
    std::set<int> generatingskel = fBoundaryMatIds;
    generatingskel.insert(fSkeletonMatId);

    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement() || gel->Reference()) {
            continue;
        }
        // we create SBFemVolume elements by partitioning the volume elements
        if (gel->Dimension() != dim) {
            continue;
        }
        if (fElementPartition[el] == -1) {
            continue;
        }
        if (matids.find(gel->MaterialId()) == matids.end()) {
            continue;
        }
        int nsides = gel->NSides();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != dim-1) {
                continue;
            }
            
            TPZGeoElSide gelside(gel,is);
            TPZGeoElSide gelskel = gelside.HasNeighbour(generatingskel);
            if(!gelskel) {
                continue;
            }
            int matid = fMatIdTranslation[gel->MaterialId()];
            TPZStack<TPZGeoElSide> gelstack;
            gelskel.YoungestChildren(gelstack);
            // we create a volume element based on all smaller elements linked to this side
            int64_t ngelstack = gelstack.NElements();
            for (int igel=0; igel<ngelstack; igel++) {
                TPZGeoElSide subgelside = gelstack[igel];
                // we are only interested in faces
                if (subgelside.Dimension() != dim-1) {
                    continue;
                }
                // verify is there is a neighbour with the partition id and material id of the collapsed element
                // if so debugstop()??
                TPZGeoElSide neighbour = subgelside.Neighbour();
                while(neighbour != subgelside)
                {
                    if (neighbour.Element()->Dimension() == dim) {
                        int64_t neighbourelpartition = fElementPartition[neighbour.Element()->Index()];
                        if (neighbourelpartition == fElementPartition[el] && neighbour.Element()->MaterialId() == matid) {
                            std::cout << __PRETTY_FUNCTION__ << "This is very strange, investigate please\n";
                            break;
                        }
                    }
                    neighbour = neighbour.Neighbour();
                }
                if (neighbour != subgelside)
                {
                    continue;
                }

                int nnodes = subgelside.NSideNodes();
                TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
                int64_t index;
                for (int in=0; in<nnodes; in++) {
                    Nodes[in] = subgelside.SideNodeIndex(in);
                }
                int elpartition = fElementPartition[el];
                for (int in=nnodes; in < 2*nnodes; in++) {
                    Nodes[in] = fPartitionCenterNode[elpartition];
                }
                if (subgelside.IsLinearMapping())
                {
                    switch(nnodes)
                    {
                        case 2:
                            gmesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index, gel->IsRefPatternEl());
                            break;
                        case 4:
                            gmesh->CreateGeoElement(ECube, Nodes, matid, index, gel->IsRefPatternEl());
                            break;
                        case 3:
                            gmesh->CreateGeoElement(EPrisma, Nodes, matid, index, gel->IsRefPatternEl());
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                    
                }
                else
                {
                    int64_t elementid = gmesh->NElements()+1;
                    switch(nnodes)
                    {
                        case 2:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);
                            break;
                        case 4:
                            new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmesh,index);
                            break;
                        case 3:
                            gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                            break;
                        default:
                            std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                            DebugStop();
                    }
                }
                if (index >= fElementPartition.size()) {
                    fElementPartition.Resize(index+1,-1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    gmesh->BuildConnectivity();
}

/// create geometric volumetric elements
// the lower dimensional elements already exist (e.g. all connects have been created
void TPZBuildSBFem::CreateVolumetricElements(TPZCompMesh &cmesh) {
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
        // all computational elements have been loaded
    std::set<int> matids, matidstarget;
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmesh.FindMaterial(mat)) {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }
    if(matidstarget.size() == 0) {
        std::cout << __PRETTY_FUNCTION__ << " matidstarget has zero size\n";
        DebugStop();
    }
    cmesh.LoadReferences();
    GenerateCollapsedGeometricElements(matids);
    cmesh.ApproxSpace().SetAllCreateFunctionsSBFem(dim);
    cmesh.AutoBuild(matidstarget);
}

/// create geometric volumetric elements
void TPZBuildSBFem::CreateVolumetricElementsFromSkeleton(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    cmesh.LoadReferences();
    std::set<int> matids, matidstarget;
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if (cmesh.FindMaterial(mat)) {
            matids.insert(it->first);
            matidstarget.insert(it->second);
        }
    }
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement() ) {
            continue;
        }
        if (!gel->Reference()) {
            continue;
        }
        if (fElementPartition[el] == -1) {
            continue;
        }
        if (gel->Dimension() > dim - 1) {
            DebugStop();
        }
        if (matids.find(gel->MaterialId()) == matids.end()) {
            continue;
        }
        int nsides = gel->NSides();
        int is = nsides - 1;
        TPZStack<TPZCompElSide> celstack;
        TPZGeoElSide gelside(gel,is);
        int nnodes = gelside.NSideNodes();
        TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
        int matid = fMatIdTranslation[gel->MaterialId()];
        int64_t index;
        for (int in=0; in<nnodes; in++) {
            Nodes[in] = gelside.SideNodeIndex(in);
        }
        int elpartition = fElementPartition[el];
        for (int in=nnodes; in < 2*nnodes; in++) {
            Nodes[in] = fPartitionCenterNode[elpartition];
        }
        if (gelside.IsLinearMapping())
        {
            switch(nnodes)
            {
                case 1:
                    gmesh->CreateGeoElement(EOned, Nodes, matid, index);
                    break;
                case 2:
                    gmesh->CreateGeoElement(EQuadrilateral, Nodes, matid, index);
                    break;
                case 4:
                    gmesh->CreateGeoElement(ECube, Nodes, matid, index);
                    break;
                case 3:
                    gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                    break;
                default:
                    std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                    DebugStop();
            }
            
        }
        else
        {
            int64_t elementid = gmesh->NElements()+1;
            switch(nnodes)
            {
                case 2:
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);
                    break;
                case 4:
                    new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmesh,index);
                    break;
                case 3:
                    gmesh->CreateGeoElement(EPrisma, Nodes, matid, index);
                    break;
                default:
                    std::cout << "Don't understand the number of nodes per side : nnodes " << nnodes << std::endl;
                    DebugStop();
            }
        }
        if (index >= fElementPartition.size()) {
            fElementPartition.Resize(index+1,-1);
        }
        fElementPartition[index] = elpartition;
        
        
    }
    gmesh->BuildConnectivity();
    cmesh.ApproxSpace().SetAllCreateFunctionsSBFem(dim);
    cmesh.AutoBuild(matidstarget);
}

/// create geometric volumetric elements for all elements with the matid
void TPZBuildSBFem::CreateVolumetricElements(TPZCompMesh &cmesh, const std::set<int> &matids)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    cmesh.LoadReferences();
    std::set<int> matidsorig,matidstarget;
    for (std::map<int,int>::iterator it = fMatIdTranslation.begin(); it!= fMatIdTranslation.end(); it++) {
        int64_t mat = it->second;
        if(matids.find(mat) != matids.end())
        {
            matidsorig.insert(it->first);
        }
    }
    int64_t nel = gmesh->NElements();
    for (int64_t el=0; el<nel; el++) {
        TPZGeoEl *gel = gmesh->Element(el);
        if (!gel || gel->HasSubElement() || gel->Reference()) {
            continue;
        }
        if (fElementPartition[el] == -1) {
            continue;
        }
        if (matidsorig.find(gel->MaterialId()) == matidsorig.end()) {
            continue;
        }
        int nsides = gel->NSides();
        int geldim = gel->Dimension();
        for (int is = 0; is<nsides; is++) {
            if (gel->SideDimension(is) != geldim-1) {
                continue;
            }
            // we only create an SBFem volume element if it is connected to a skeleton element?
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel,is);
            int onlyinterpolated = true;
            int removeduplicates = true;
            gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);
            int ncelstack = celstack.NElements();
            for (int icel=0; icel<ncelstack; icel++) {
                TPZGeoElSide subgelside = celstack[icel].Reference();
                // we are only interested in faces
                if (subgelside.Dimension() != geldim-1) {
                    continue;
                }
                int nnodes = subgelside.NSideNodes();
//                if (nnodes != 2) {
//                    std::cout << "Please extend the code to higher dimensions\n";
//                }
                TPZManVector<int64_t,4> Nodes(nnodes*2,-1);
                int matid = fMatIdTranslation[gel->MaterialId()];
                int64_t index;
                for (int in=0; in<nnodes; in++) {
                    Nodes[in] = subgelside.SideNodeIndex(in);
                }
                int elpartition = fElementPartition[el];
                for (int in=nnodes; in < nnodes*2; in++) {
                    Nodes[in] = fPartitionCenterNode[elpartition];
                }
                MElementType targettype = ENoType;
                switch (nnodes) {
                    case 2:
                        targettype = EQuadrilateral;
                        break;
                    case 1:
                        targettype = EOned;
                        break;
                    case 3:
                        targettype = EPrisma;
                        break;
                    case 4:
                        targettype = ECube;
                        break;
                        
                    default:
                        DebugStop();
                        break;
                }
                if (subgelside.IsLinearMapping())
                {
                    gmesh->CreateGeoElement(targettype, Nodes, matid, index);
#ifdef PZ_LOG
                    if(logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Created element of type " << targettype << " with nodes " << Nodes << " index " << index;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                }
                else
                {
                    int64_t elementid = gmesh->NElements()+1;
                    TPZGeoEl *gblend = 0;
                    switch (targettype) {
                        case EOned:
                            gblend = gmesh->CreateGeoElement(targettype, Nodes, matid, index);
                            break;
                        case EQuadrilateral:
                            gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoQuad> > (Nodes, matid, *gmesh,index);
                            break;
                        case EPrisma:
                            gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoPrism> > (Nodes, matid, *gmesh,index);
                            break;
                        case ECube:
                            gblend = new TPZGeoElRefPattern< pzgeom::TPZGeoBlend<pzgeom::TPZGeoCube> > (Nodes, matid, *gmesh,index);
                            break;
                        default:
                            DebugStop();
                            break;
                    }
#ifdef PZ_LOG
                    if(logger.isDebugEnabled())
                    {
                        std::stringstream sout;
                        sout << "Created element of type " << targettype << " with nodes " << Nodes << " index " << index;
                        LOGPZ_DEBUG(logger, sout.str())
                    }
#endif
                    
                    
                }
                if (index >= fElementPartition.size()) {
                    fElementPartition.Resize(index+1,-1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    gmesh->BuildConnectivity();
    cmesh.ApproxSpace().SetAllCreateFunctionsSBFem(dim);
    cmesh.AutoBuild(matids);
    
}


/// put the sbfem volumetric elements in element groups
void TPZBuildSBFem::CreateElementGroups(TPZCompMesh &cmesh)
{
    cmesh.LoadReferences();
    int64_t numgroups = fPartitionCenterNode.size();
    int64_t groupelementindices(numgroups);
    
    TPZVec<int64_t> elementgroupindices(numgroups);
    TPZSBFemElementGroup::SetDefaultPolynomialOrder(this->fPOrderBubbleFunctions);
    for (int64_t el=0; el<numgroups; el++) {
        TPZCompEl* cel = new TPZSBFemElementGroup(cmesh);
        elementgroupindices[el] = cel->Index();
    }
    
    
    int64_t nel = cmesh.NElements();
    int dim = cmesh.Dimension();
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = cmesh.Element(el);
        if (!cel) {
            continue;
        }
        TPZSBFemVolume *sbfem = dynamic_cast<TPZSBFemVolume *>(cel);
        if (sbfem) {
            TPZGeoEl *gel = sbfem->Reference();
            int64_t gelindex = gel->Index();
            int side = -1;
            if (gel->Type() == EQuadrilateral) {
                side = 4;
            }
            else if(gel->Type() == EOned)
            {
                side = 0;
            }
            else if(gel->Type() == ECube)
            {
                side = 20;
            }
            else if(gel->Type() == EPrisma)
            {
                side = 15;
            }
            else
            {
                DebugStop();
            }
            TPZGeoElSide gelside(gel,side);
            int geldim = gel->Dimension();
            int nsidenodes = gel->NSideNodes(side);
            TPZManVector<int64_t,8> cornernodes(nsidenodes);
            for (int node = 0; node<nsidenodes; node++) {
                cornernodes[node] = gel->SideNodeIndex(side, node);
            }
            
            TPZGeoElSide neighbour = gelside.Neighbour();
            while (neighbour != gelside) {
                if(neighbour.Element()->Dimension() == geldim-1 && neighbour.Element()->Reference())
                {
                    int nsidenodesneighbour = neighbour.Element()->NCornerNodes();
                    if (nsidenodesneighbour == nsidenodes)
                    {
                        TPZManVector<int64_t,8> neighnodes(nsidenodesneighbour);
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            neighnodes[i] = neighbour.Element()->NodeIndex(i);
                        }
                        int numequal = 0;
                        for (int i=0; i<nsidenodesneighbour; i++) {
                            if (neighnodes[i] == cornernodes[i]) {
                                numequal++;
                            }
                        }
                        if (numequal == nsidenodesneighbour) {
                            break;
                        }
                    }
                }
                neighbour = neighbour.Neighbour();
            }
            if (neighbour == gelside) {
                // we are not handling open sides (yet)
                DebugStop();
            }
            int64_t skelindex = neighbour.Element()->Reference()->Index();
            sbfem->SetSkeleton(skelindex);
            
            int64_t gelgroup = fElementPartition[gelindex];
            if (gelgroup == -1) {
                DebugStop();
            }
            int64_t celgroupindex = elementgroupindices[gelgroup];
            TPZCompEl *celgr = cmesh.Element(celgroupindex);
            TPZSBFemElementGroup *sbfemgr = dynamic_cast<TPZSBFemElementGroup *>(celgr);
            if (!sbfemgr) {
                DebugStop();
            }
            sbfemgr->AddElement(sbfem);
            //            sbfem->SetElementGroupIndex(celgroupindex);
        }
    }
    
    for (int64_t el=0; el<numgroups; el++) {
        int64_t index;
        
        index = elementgroupindices[el];
        TPZCompEl *cel = cmesh.Element(index);
        TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
        if (!sbfemgroup) {
            DebugStop();
        }
        const TPZVec<TPZCompEl *> &subgr = sbfemgroup->GetElGroup();
        int64_t nsub = subgr.NElements();
        for (int64_t is=0; is<nsub; is++) {
            TPZCompEl *cel = subgr[is];
            TPZSBFemVolume *femvol = dynamic_cast<TPZSBFemVolume *>(cel);
            if (!femvol) {
                DebugStop();
            }
            femvol->SetElementGroupIndex(index);
        }
        if (nsub == 0)
        {
            delete sbfemgroup;
        }
        
    }

    if (this->fPOrderBubbleFunctions != 0)
    {
        for (int64_t el=0; el<numgroups; el++)
        {
            int64_t index;
            
            index = elementgroupindices[el];
            TPZCompEl *cel = cmesh.Element(index);
            TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
            if (!sbfemgroup || sbfemgroup->NConnects() == 0) {
                continue;
            }
            sbfemgroup->InitializeInternalConnect();
        }
    }
    // for (int64_t el=0; el<numgroups; el++) {
    //     int64_t index;
        
    //     index = elementgroupindices[el];
    //     TPZCompEl *cel = cmesh.Element(index);
    //     TPZSBFemElementGroup *sbfemgroup = dynamic_cast<TPZSBFemElementGroup *>(cel);
    //     if (!sbfemgroup || sbfemgroup->NConnects() == 0) {
    //         continue;
    //     }
        

        // sbfemgroup->ComputeEigenvalues();
    // }
}

/// Divide de skeleton elements
void TPZBuildSBFem::DivideSkeleton(int nref)
{
    int dim = fGMesh->Dimension();
    for (int ir=0; ir<nref; ir++)
    {
        TPZAdmChunkVector<TPZGeoEl *> elvec = fGMesh->ElementVec();
        int64_t nel = elvec.NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = elvec[el];
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            if (gel->Dimension() != dim-1) {
                continue;
            }
            TPZManVector<TPZGeoEl *,10> subel;
            gel->Divide(subel);
            int64_t partition = fElementPartition[el];
            int nsub = subel.size();
            for (int isub=0; isub<nsub; isub++) {
                while (fElementPartition.size() <= subel[isub]->Index()) {
                    fElementPartition.Resize(fElementPartition.size()*2, -1);
                }
                fElementPartition[subel[isub]->Index()] = partition;
            }
        }
    }
}

/// Divide de skeleton elements
void TPZBuildSBFem::DivideSkeleton(int nref, const std::set<int> &volmatids)
{
    std::set<int> exclude;
    for (auto it = fMatIdTranslation.begin(); it != fMatIdTranslation.end(); it++) {
        if(volmatids.find(it->second) != volmatids.end())
        {
            exclude.insert(it->first);
        }
    }
    int dim = fGMesh->Dimension();
    for (int ir=0; ir<nref; ir++)
    {
        TPZAdmChunkVector<TPZGeoEl *> elvec = fGMesh->ElementVec();
        int64_t nel = elvec.NElements();
        for (int64_t el=0; el<nel; el++) {
            TPZGeoEl *gel = elvec[el];
            if (!gel || gel->HasSubElement()) {
                continue;
            }
            int matid = gel->MaterialId();
            // skip the elements with matid volmatids
            if (exclude.find(matid) != exclude.end()) {
                continue;
            }
            /// skip the elements which do not have material translation
            if (fMatIdTranslation.find(matid) == fMatIdTranslation.end()) {
                continue;
            }
            if (gel->Dimension() != dim-1) {
                continue;
            }
            TPZManVector<TPZGeoEl *,10> subel;
            gel->Divide(subel);
            int64_t partition = fElementPartition[el];
            int nsub = subel.size();
            for (int isub=0; isub<nsub; isub++) {
                while (fElementPartition.size() <= subel[isub]->Index()) {
                    fElementPartition.Resize(fElementPartition.size()*2, -1);
                }
                fElementPartition[subel[isub]->Index()] = partition;
            }
            
        }
    }
}

/// @brief print the SBFem configuration
/// @param out output stream
void TPZBuildSBFem::Print(std::ostream &out) const
{
    out << "SBFem Configuration:" << std::endl;
    out << "Skeleton Material Id: " << fSkeletonMatId << std::endl;
    out << "Material Id Translation:" << std::endl;
    for (const auto &item : fMatIdTranslation) {
        out << "  From: " << item.first << " To: " << item.second << std::endl;
    }
    out << "fElement partition " << fElementPartition << std::endl;
    out << "Element Partition:" << std::endl;
    for (int64_t el=0; el<fElementPartition.size(); el++) {
        auto elpart = fElementPartition[el];
        out << "  Element: " << el << " Partition: " << elpart;
        TPZGeoEl *gel = fGMesh->Element(el);
        if(gel) {
            out << " matid " << gel->MaterialId();
        }
        if(elpart != -1)
        {
            out << " Center node: " << fPartitionCenterNode[elpart];
        }
        out << std::endl;
    }
}

/// @brief Identify the partition corresponding to the element group
/// @param elgr SBFemElementGroup that contains the elements of a given partition
/// @return index of the partition
int64_t TPZBuildSBFem::GetPartition(TPZSBFemElementGroup *elgr) {
    auto elvec = elgr->GetElGroup();
    int64_t nel = elvec.size();
    if(nel == 0) DebugStop();
    TPZCompEl *cel = elvec[0];
    TPZGeoEl *gel = cel->Reference();
    int64_t index = gel->Index();
    int64_t partition = fElementPartition[index];
    if(partition == -1) DebugStop();
    return partition;
}
