/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef TPZCompElKernelHDiv3D_H
#define TPZCompElKernelHDiv3D_H

#include "pzelctemp.h"
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
class TPZCompElKernelHDiv3D : public TPZCompElHCurl<TSHAPE> {
    
    /// Family of the HDiv space being used. Only HDivFamily::EHDivKernel is supported
    HDivFamily fhdivfam = HDivFamily::EHDivKernel;

public:
	    
	TPZCompElKernelHDiv3D(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = DefaultFamily::fHDivDefaultValue);
	
	TPZCompElKernelHDiv3D(){};
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
		
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	// virtual void InitMaterialData(TPZMaterialData &data) override;
    void ComputeShape(TPZVec<REAL> &qsi, TPZMaterialData &data) override;

    /** @brief Initialize a material data and its attributes based on element dimension, number
    * of state variables and material definitions */
    void InitMaterialData(TPZMaterialData &data) override;

	/** @brief Returns the unique identifier for reading/writing objects to streams */
    int ClassId() const override;

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

};


template<class TSHAPE>
int TPZCompElKernelHDiv3D<TSHAPE>::ClassId() const{
    return Hash("TPZCompElKernelHDiv3D") ^ TPZCompElHCurl<TSHAPE>::ClassId() << 1;
}

#include "pzcmesh.h"

template<class TSHAPE>
void TPZCompElKernelHDiv3D<TSHAPE>::SetCreateFunctions(TPZCompMesh* mesh) {
    mesh->ApproxSpace().SetHDivFamily(HDivFamily::EHDivKernel);
    mesh->ApproxSpace().SetAllCreateFunctionsHDiv(3);
}

#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "tpzquadrilateral.h"

/** @brief Creates computational cube element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational prismal element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational tetrahedral element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);

/** @brief Creates computational quadrilateral element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDivKernel approximate space */
TPZCompEl *CreateHDivKernelBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh, const HDivFamily hdivfam);


#endif
