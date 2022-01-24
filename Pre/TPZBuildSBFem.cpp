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
void TPZBuildSBFem::AddPartition(TPZVec<int64_t> &elindices, int64_t centernodeindex)
{
    int64_t npart = fPartitionCenterNode.size();
    fPartitionCenterNode.resize(npart+1);
    if (fGMesh->NodeVec().NElements() <= centernodeindex) {
        DebugStop();
    }
    fPartitionCenterNode[npart] = centernodeindex;
    int64_t nel = elindices.size();
    for (int64_t el=0; el<nel; el++) {
        int64_t elindex = elindices[el];
        if (fElementPartition[elindex] != -1) {
            DebugStop();
        }
        fElementPartition[elindex] = npart;
    }
}

/// add the sbfem elements to the computational mesh, the material should exist in cmesh
void TPZBuildSBFem::BuildComputationMesh(TPZCompMesh &cmesh)
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



/// create the geometric skeleton elements
void TPZBuildSBFem::AddSkeletonElements()
{
    // create a lower dimension element on each boundary
    int dim = fGMesh->Dimension();
    
    int64_t nel = fGMesh->NElements();
    
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
        int64_t elpartition = fElementPartition[el];
        int nsides = gel->NSides();
        for (int is=0; is<nsides; is++) {
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
            }
            
            // look for a neighbour of dimension dim-1
            neighbour = thisside.Neighbour();
            while (neighbour != thisside && neighbour.Element()->Dimension() != dim-1) {
                neighbour = neighbour.Neighbour();
            }
            // if we didnt find a lower dimension element and the neighbouring element of same dimension belong to a different partition
            if (thisside == neighbour && elpartition != neighbourelpartition) {
                gel->CreateBCGeoEl(is,fSkeletonMatId);
            }
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
        if (gel->Dimension() != dim) {
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
void TPZBuildSBFem::CreateVolumetricElements(TPZCompMesh &cmesh)
{
    TPZGeoMesh *gmesh = cmesh.Reference();
    gmesh->ResetReference();
    int dim = gmesh->Dimension();
    cmesh.LoadReferences();
    // all computational elements have been loaded
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
            TPZStack<TPZCompElSide> celstack;
            TPZGeoElSide gelside(gel,is);
            int onlyinterpolated = true;
            int removeduplicates = true;
            // we identify all computational elements connected to this element side
            gelside.EqualorHigherCompElementList2(celstack, onlyinterpolated, removeduplicates);
            // we create a volume element based on all smaller elements linked to this side
            int ncelstack = celstack.NElements();
            for (int icel=0; icel<ncelstack; icel++) {
                TPZGeoElSide subgelside = celstack[icel].Reference();
                // we are only interested in faces
                if (subgelside.Dimension() != dim-1) {
                    continue;
                }
                int nnodes = subgelside.NSideNodes();
                TPZManVector<int64_t,8> Nodes(nnodes*2,-1);
                int matid = fMatIdTranslation[gel->MaterialId()];
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
                    fElementPartition.resize(index+1);
                }
                fElementPartition[index] = elpartition;
            }
        }
    }
    gmesh->BuildConnectivity();
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
            fElementPartition.resize(index+1);
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
                    fElementPartition.resize(index+1);
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
    int64_t numgroups = fPartitionCenterNode.size();
    int64_t groupelementindices(numgroups);
    
    TPZVec<int64_t> elementgroupindices(numgroups);
    
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
    }

    if (TPZSBFemElementGroup::gDefaultPolynomialOrder != 0)
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
