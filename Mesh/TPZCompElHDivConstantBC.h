/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELHDIVCONSTANTBC_H
#define TPZCOMPELHDIVCONSTANTBC_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"
#include "pzelchdivbound2.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivConstantBC : public TPZCompElHDivBound2<TSHAPE>  {



public:
    // Type of HDiv Space
    enum MShapeType {EHDivKernel, EHDivConstant, ECurlNoGrads};
	    
	TPZCompElHDivConstantBC();
    
    TPZCompElHDivConstantBC(TPZCompMesh &mesh, TPZGeoEl *gel, int shapetype = EHDivKernel);
	
	virtual ~TPZCompElHDivConstantBC();

    virtual void InitMaterialData(TPZMaterialData &data) override;

	void ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data) override;
    
    virtual int NConnectShapeF(int connect, int order) const override;

    void AdjustConnects();
    
private:
    /// the type of space this object will generate
    int fShapeType;

};


#endif
