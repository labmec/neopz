/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCompElKernelHDiv_H
#define TPZCompElKernelHDiv_H

#include "pzelctemp.h"
#include "TPZBndCond.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"
#include "TPZCompElH1.h"
#include "TPZCompElHCurl.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElKernelHDiv : public TPZCompElH1<TSHAPE>  {

protected:
  ///! Indexes of the connects associated with the elements
  TPZManVector<int64_t,TSHAPE::NSides> fConnectIndexes =
    TPZManVector<int64_t,TSHAPE::NSides>(TSHAPE::NSides,-1);

public:
	    
	TPZCompElKernelHDiv();
    
    TPZCompElKernelHDiv(TPZCompMesh &mesh, TPZGeoEl *gel);
	
    virtual void InitMaterialData(TPZMaterialData &data) override;

    // void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override;
    void ComputeShape(TPZVec<REAL> &qsi,TPZMaterialData &data) override;

	/** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
    
    int NConnectShapeF(int connect, int order) const override;

    /// Return the maximum order
    virtual int MaxOrder() override;

protected:
	 //@{
    /** @brief Compute the solution using Hdiv structure */
	void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override{
        ComputeSolutionKernelHdivT(data);
    }
    void ReallyComputeSolution(TPZMaterialDataT<CSTATE> &data) override{
        ComputeSolutionKernelHdivT(data);
    }

	template<class TVar>
    void ComputeSolutionKernelHdivT(TPZMaterialDataT<TVar> &data);
    template<class TVar>
    void ComputeRequiredDataT(TPZMaterialDataT<TVar> &data, TPZVec<REAL>&qsi);
    
};

#include "tpzquadrilateral.h"
#include "tpztriangle.h"

/** @brief Creates computational quadrilateral element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);

/** @brief Creates computational linear element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);

#endif
