

#ifndef TPZHDIVAPPROXCREATOR_H
#define TPZHDIVAPPROXCREATOR_H

#include "TPZEnumApproxFamily.h"
#include "pzstack.h"
#include "pzgmesh.h"

class TPZMaterial;
class TPZMultiphysicsCompMesh;

enum class HybridizationType {ENone, EStandard, ESemi};
enum class ProblemType {ENone, EElastic, EDarcy, EStokes};

class TPZHDivApproxCreator
{
protected:
    
    HDivFamily fHDivFam;
    
    HybridizationType fHybridType = HybridizationType::ENone;
    
    ProblemType fProbType = ProblemType::ENone;

    TPZStack<TPZMaterial * > fMaterialVec;

    TPZGeoMesh *fGeoMesh = nullptr;

    int fDefaultPOrder = 1;
    
    /// @brief This sets hdiv+/hdiv++ internal/bubble POrder: 0- hdiv, 1- hdiv+, 2- hdiv++
    int fExtraInternalPOrder = 0; 

public:
    TPZHDivApproxCreator() = default;

    TPZHDivApproxCreator(TPZGeoMesh * gmesh);

    ~TPZHDivApproxCreator();

    int InsertMaterialObject(TPZMaterial * mat);
    int InsertMaterialObject(TPZBndCond * mat);

    TPZGeoMesh *GeoMesh(){return fGeoMesh;}
    const TPZGeoMesh *GeoMesh() const {return fGeoMesh;}

    HDivFamily &HdivFamily(){return fHDivFam;}
    const HDivFamily &HdivFamily() const {return fHDivFam;}

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
    int &GetExtraInternalOrder(){return fExtraInternalPOrder;}

    TPZMultiphysicsCompMesh *CreateApproximationSpace();

private:

    TPZCompMesh * CreateHDivSpace();

    TPZCompMesh * CreateL2Space(const int lagLevel);

    TPZCompMesh * CreateConstantSpace();
   

};




#endif