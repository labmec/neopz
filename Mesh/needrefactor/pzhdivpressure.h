/**
 * @file
 * @brief Contains declaration of TPZCompElHDivPressure class which implements a generic computational element (HDiv scope).
 */

#ifndef PZHDIVPRESSUREHTT
#define PZHDIVPRESSUREHTT

#include "pzelchdiv.h"
#include "TPZEnumApproxFamily.h"

/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivPressure : public TPZCompElHDiv<TSHAPE> {
	
	/** @brief Defines the interpolation order for pressure variable*/
	int fPressureOrder;
	
	/** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
public:
	
	TPZCompElHDivPressure(TPZCompMesh &mesh, TPZGeoEl *gel, const HDivFamily hdivfam = HDivFamily::EDefault);
	
	TPZCompElHDivPressure(TPZCompMesh &mesh, const TPZCompElHDivPressure<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDivPressure(TPZCompMesh &mesh,
						  const TPZCompElHDivPressure<TSHAPE> &copy,
						  std::map<int64_t,int64_t> & gl2lcConMap,
						  std::map<int64_t,int64_t> & gl2lcElMap);
	
	TPZCompElHDivPressure();
	
	virtual ~TPZCompElHDivPressure();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override {
		return new TPZCompElHDivPressure<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int64_t,int64_t> & gl2lcConMap,std::map<int64_t,int64_t>&gl2lcElMap) const override
	{
		return new TPZCompElHDivPressure<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
  virtual void SetCreateFunctions(TPZCompMesh *mesh) override {
		mesh->SetAllCreateFunctionsHDivPressure();
	}
	
	virtual MElementType Type() override;
	
	virtual int NConnects() const override;
	
	virtual void SetConnectIndex(int i, int64_t connectindex) override;
	
	virtual int NConnectShapeF(int connect, int order) const override;
	
	virtual int Dimension() const override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const override {
		return 0;
	}
	
	virtual int64_t ConnectIndex(int node) const override;
    
    /** @brief returns the index of the pressure connect
     * returns -1 if their is no pressure connect
     */
    virtual int PressureConnectIndex() const override
    {
        return NConnects()-1;
    }
	
	/** @brief Identifies the interpolation order for pressure variable*/
	virtual void SetPressureOrder(int ord);
	/** @brief Returns the interpolation order to dual variable */
	int DualOrder();
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Sets the preferred interpolation order along a side
	 
	 This method only updates the datastructure of the element
	 In order to change the interpolation order of an element, use the method PRefine*/
	virtual void SetPreferredOrder(int order) override;
	
	virtual int ConnectOrder(int connect) const override;
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data) override;
	
	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
	/** 
	 * @brief Compute the shape functions corresponding to the dual space
	 */
	virtual void ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
    //@{
    ///Compute the solution for a given variable
	void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override{
		SolutionT(qsi,var,sol);
	}
	void Solution( TPZVec<REAL> &qsi,int var,TPZVec<CSTATE> &sol) override{
		SolutionT(qsi,var,sol);
	}
	//@}
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension) override;
	
	/** @brief Returns the transformation which transform a point from the side to the interior of the element */
	TPZTransform<> TransformSideToElement(int side) override;
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */

	int ClassId() const override;

	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;

protected:
	/** @brief Compute the solution using Hdiv structure */
	void ReallyComputeSolution(TPZMaterialDataT<STATE> &data) override{
		ComputeSolutionHDivPressureT(data);
	}
	void ReallyComputeSolution(TPZMaterialDataT<CSTATE> &data) override{
		ComputeSolutionHDivPressureT(data);
	}
	template<class TVar>
	void ComputeSolutionHDivPressureT(TPZMaterialDataT<TVar>&data);
	template<class TVar>
	void SolutionT( TPZVec<REAL> &qsi,int var,TPZVec<TVar> &sol);
};

template<class TSHAPE>
int TPZCompElHDivPressure<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivPressure") ^ TPZCompElHDiv<TSHAPE>::ClassId() << 1;
}

/** @brief Creates computational point element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational linear element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational quadrilateral element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational triangular element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational cube element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational prismal element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational pyramidal element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressurePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);
/** @brief Creates computational tetrahedral element for HDivPressure approximate space */
TPZCompEl *CreateHDivPressureTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,const HDivFamily hdivfam);

/** @} */

#endif
