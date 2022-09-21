//
// Created by victor on 21/09/2022.
//

#ifndef PZ_TPZAPPROXCREATOR_H
#define PZ_TPZAPPROXCREATOR_H

#include "TPZEnumApproxFamily.h"
#include "pzstack.h"
#include "pzgmesh.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

enum class HybridizationType {ENone, EStandard, EStandardSquared, ESemi};
enum class ProblemType {ENone, EElastic, EDarcy, EStokes};

class TPZApproxCreator {
protected:

    HybridizationType fHybridType = HybridizationType::ENone;

    ProblemType fProbType = ProblemType::ENone;

    TPZStack<std::pair<int, TPZMaterial * > > fMaterialVec;

    TPZStack<std::pair<int, TPZBndCond * > > fBCMaterialVec;

    TPZGeoMesh *fGeoMesh = nullptr;

    int fDefaultPOrder = 1;

    /// @brief This sets hdiv+/hdiv++ internal/bubble POrder: 0- hdiv, 1- hdiv+, 2- hdiv++
    int fExtraInternalPOrder = 0;

    int fNumMeshes = 2;

    bool fShouldCondense = 1;

    bool fIsEnhancedSpaces = false;

    struct HybridizationData{
        /// Material ids of additional elements generated due hybridization
        int fWrapMatId;

        int fInterfaceMatId;

        int fLagrangeMatId;

        /// Material ids generated due second hybridization
        int fSecondInterfaceMatId;

        int fSecondLagrangeMatId;

        /// indicates whether the boundary conditions should be hybridized and how many times it should be
        int fHybridizeBCLevel = 0;
    };

    HybridizationData fHybridizationData;

public:
    TPZApproxCreator() = default;

    TPZApproxCreator(TPZGeoMesh * gmesh);

    ~TPZApproxCreator() = default;

    int InsertMaterialObject(int matid, TPZMaterial * mat);

    int InsertBCMaterialObject(int matid, TPZBndCond * mat);

    TPZGeoMesh *GeoMesh(){return fGeoMesh;}
    const TPZGeoMesh *GeoMesh() const {return fGeoMesh;}

    HybridizationType &HybridType(){return fHybridType;}
    const HybridizationType &HybridType() const {return fHybridType;}

    ProblemType &ProbType(){return fProbType;}
    const ProblemType &ProbType() const {return fProbType;}

    void SetDefaultOrder(const int ord){
        fDefaultPOrder = ord;
    }
    int &GetDefaultOrder(){return fDefaultPOrder;}

    void SetExtraInternalOrder(const int ord){
        fExtraInternalPOrder = ord;
    }

    void SetShouldCondense( bool isCondensed){
        fShouldCondense = isCondensed;
    }

    bool GetShouldCondense( bool isCondensed){
        return fShouldCondense;
    }

    int &GetExtraInternalOrder(){return fExtraInternalPOrder;}

    bool &EnhancedSpaces(){return fIsEnhancedSpaces;}
    const bool &EnhancedSpaces() const {return fIsEnhancedSpaces;}

    virtual TPZMultiphysicsCompMesh *CreateApproximationSpace() = 0;

    void Print(std::ostream &out = std::cout);

protected:
    /// Compute Periferal Material ids for additional elements generated during hybridization
    void ComputePeriferalMaterialIds(int base = 10);

    void AddHybridizationGeoElements();

    ///This method checks if the current configuration is valid
    virtual void CheckSetupConsistency() = 0;
};


#endif //PZ_TPZAPPROXCREATOR_H
