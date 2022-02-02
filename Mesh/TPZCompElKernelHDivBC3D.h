/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCompElKernelHDivBC3D_H
#define TPZCompElKernelHDivBC3D_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"
#include "TPZCompElHCurlNoGrads.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDivBC3D : public TPZCompElHCurlNoGrads<TSHAPE>  {
   
private:
    int fSideOrient = 1;

    /// Family of the HDiv/HCurl space being used. Changing this will change the shape generating class
    // The values are set as Default HDiv or HCurl, but the class will only work for HDivKernel or HCurlNoGrads
    HDivFamily fhdivfam = DefaultFamily::fHDivDefaultValue;
    HCurlFamily fhcurlfam = DefaultFamily::fHCurlDefaultValue;
    
public:
	    
	TPZCompElKernelHDivBC3D();
    
    TPZCompElKernelHDivBC3D(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue, 
                            const HCurlFamily hcurlfam = DefaultFamily::fHCurlDefaultValue);
	
    virtual void InitMaterialData(TPZMaterialData &data) override;

    
    void ComputeRequiredDataT(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi);
    
    void SetSideOrient(int orient);

    int GetSideOrient();
    void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override{
        ComputeRequiredDataT(data,qsi);
    }

};


#endif
