//
//  TPZBuildSBFemMultiphysics.h
//  PZ
//
//  Created by Karolinne Coelho on 18/01/21.
//
//

#pragma once

#include "TPZBuildSBFem.h"
#include "TPZMultiphysicsCompMesh.h"

using namespace std;

class TPZBuildSBFemMultiphysics : public TPZBuildSBFem
{
    int fDifPressure, fInternal, fLeftFlux, fRightFlux, fInterface;

    set<int> fCondensedMatids;

    map<int64_t,TPZCompEl *> fGeltocel;

    TPZManVector<int64_t> fElGroupIndexes;
    
public:
    
    /// simple constructor
    TPZBuildSBFemMultiphysics(TPZAutoPointer<TPZGeoMesh> & gmesh, int skeletonmatid, std::map<int,int> &matidtranslation) : TPZBuildSBFem(gmesh,skeletonmatid, matidtranslation)
    {
        // fSkeletonMatId represents the external average pressure;
        fInterface = fSkeletonMatId+1; //7 

        // fRightFlux represents the external right flux (next to fSkeletonMatId)
        fRightFlux = fInterface+1; //8

        // fInternal will gather internal flux and pressures (next to fRightFlux)
        fInternal = fRightFlux+1; //9
        
        // fLeftFlux represents the external left flux (next to fInternal)
        fLeftFlux = fInternal +1; //10

        // fDifPressure represents the differential of the pressure (next to fLeftFlux)
        fDifPressure = fLeftFlux + 1; //11
        
        // Order of the elements, from left to right:
        // fDifPressure -> fInterface -> fLeftFlux -> fInternal -> fRightFlux -> fInterface -> fSkeletonMatId
        fCondensedMatids = {fInterface, fLeftFlux, fInternal, fRightFlux};
    }

    int GetSideSkeletonEl(TPZGeoEl * gel);

    int GetSideCollapsedEl(TPZGeoEl * gel);

    void BuildMultiphysicsCompMesh(TPZMultiphysicsCompMesh &cmesh);

    void CreateExternalElements(TPZAutoPointer<TPZGeoMesh> & gmesh, set<int> & matidtarget, set <int> & matids1d);

    void CreateCollapsedGeoEls(TPZCompMesh & cmeshpressure, set<int> & matidstarget, set<int> & matids1d);

    void CreateCompElPressure(TPZCompMesh &cmeshpressure, set<int> & matids1d);

    void CreateCompElFlux(TPZCompMesh &cmeshflux, set<int> & matidtarget, set<int> & matid1d);

    void CreateSBFemMultiphysicsMesh(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

    void AddInterfaceElements(TPZMultiphysicsCompMesh & cmeshm, set<int> &matids1d);

    void CreateSBFemVolumeFlux(TPZCompMesh & cmesh, set<int> &matids1d, set<int> & matidtarget);

    void CreateSBFemVolumePressure(TPZCompMesh & cmesh, set<int> &matids1d, set<int> & matidtarget);
    
    void CreateSBFemMultiphysicsVol(TPZMultiphysicsCompMesh & cmeshm, set<int> &matid1d, set<int> & matidtarget);

    void CreateSBFemMultiphysicsElGroups(TPZMultiphysicsCompMesh & cmeshm, set<int> & matidtarget);

    void GroupandCondense(TPZMultiphysicsCompMesh & cmeshm);

    void StandardConfigurationHdiv();

    void AddInternalElements();

// NOT READY YET
    void BuildMultiphysicsCompMeshfromSkeleton(TPZCompMesh &cmesh);
};