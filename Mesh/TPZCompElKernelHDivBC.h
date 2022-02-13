/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCOMPELKERNELHDIVBC_H
#define TPZCOMPELKERNELHDIVBC_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"
#include "TPZCompElH1.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDivBC : public TPZCompElH1<TSHAPE>  {
    int fSideOrient;

public:	    
	TPZCompElKernelHDivBC();
    
    TPZCompElKernelHDivBC(TPZCompMesh &mesh, TPZGeoEl *gel);
	
	virtual ~TPZCompElKernelHDivBC();

    virtual void InitMaterialData(TPZMaterialData &data) override;

	void ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data) override;

    virtual void SetSideOrient(int side, int sideorient) override;
    virtual int GetSideOrient(int side) override;
};


#endif
