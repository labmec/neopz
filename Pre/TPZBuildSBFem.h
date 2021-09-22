//
//  TPZBuildSBFem.hpp
//  PZ
//
//  Created by Philippe Devloo on 06/01/17.
//
//

#ifndef TPZBuildSBFem_hpp
#define TPZBuildSBFem_hpp

#include <stdio.h>
#include "pzmanvector.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzautopointer.h"
#include <map>

class TPZBuildSBFem
{
protected:
    /// geometric mesh
    TPZAutoPointer<TPZGeoMesh> fGMesh;
    
    /// The volumetric elements with Mat Id will spawn SBFemVolume elements with MatId
    std::map<int,int> fMatIdTranslation;
    
    /// Material Id associated with the skeleton elements
    int fSkeletonMatId;
    
    /// partition to which each element belongs
    TPZManVector<int64_t> fElementPartition;
    
    /// center node id for each partition
    TPZManVector<int64_t> fPartitionCenterNode;
    
public:
    
    /// simple constructor
    TPZBuildSBFem(TPZAutoPointer<TPZGeoMesh> & gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : fGMesh(gmesh), fMatIdTranslation(matidtranslation), fSkeletonMatId(skeletonmatid)
    {
        fElementPartition.resize(fGMesh->NElements());
        fElementPartition.Fill(-1);
    }
    
    /// set the matid translation
    void SetMatIdTranslation(const std::map<int,int> &matidtranslation)
    {
        fMatIdTranslation = matidtranslation;
    }
    
    /// standard configuration means each element is a partition and a center node is created
    void StandardConfiguration();

    /// standard configuration means each element is a partition and a center node is created for the indicated elements
    void StandardConfiguration(TPZVec<int64_t> &elementindices);
    
    /// build element groups according to the id of the scaling centers
    void Configure(TPZVec<int64_t> &scalingcenters);
    
    /// add a partition manually
    void AddPartition(TPZVec<int64_t> &elids, int64_t centernodeindex);
    
    /// define the partition index of each element and the ids of the scaling centers
    void SetPartitions(TPZVec<int64_t> &gelpartitionids, TPZVec<int64_t> &partition_nodeindices)
    {
#ifdef PZDEBUG
        if(gelpartitionids.size() != fGMesh->NElements())
        {
            DebugStop();
        }
#endif
        fElementPartition = gelpartitionids;
        fPartitionCenterNode = partition_nodeindices;
    }
    
    /// add the sbfem elements to the computational mesh, the material should exist in cmesh
    void BuildComputationMesh(TPZCompMesh &cmesh);
    
    /// add the sbfem elements to the computational mesh, the material should exist in cmesh
    void BuildComputationMesh(TPZCompMesh &cmesh, const std::set<int> &volmatids, const std::set<int> &boundmatids);
    
    /// build the computational elements of the skeleton and build the volume elements directly from the skeleton elements
    void BuildComputationalMeshFromSkeleton(TPZCompMesh &cmesh);
    
    /// Divide the skeleton elements
    void DivideSkeleton(int nref);

    /// Divide the skeleton elements - elements of dimension dim-1 which are not in volmatids
    void DivideSkeleton(int nref, const std::set<int> &volmatids);

private:
    /// create the geometric skeleton elements
    void AddSkeletonElements();
    
    /// create a geometric node at the center of each partition
    void CreateElementCenterNodes(TPZVec<int64_t> &elindices);
    
    /// create geometric volumetric elements
    void CreateVolumetricElements(TPZCompMesh &cmesh);
    
    /// create geometric volumetric elements
    void CreateVolumetricElementsFromSkeleton(TPZCompMesh &cmesh);
    
    /// create geometric volumetric elements for all elements with the matid
    void CreateVolumetricElements(TPZCompMesh &cmesh, const std::set<int> &matids);
    
    /// put the sbfem volumetric elements in element groups
    void CreateElementGroups(TPZCompMesh &cmesh);
};

#endif /* TPZBuildSBFem_hpp */
