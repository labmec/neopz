/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */

#ifndef PZELCHDIVHTT
#define PZELCHDIVHTT

#include "pzelctemp.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @addtogroup CompElement
 * @{
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDiv : public TPZIntelGen<TSHAPE> {
	
    TPZManVector<int, TSHAPE::NFaces> fSideOrient;
    
	/** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
public:
	
	TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
	
	TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDiv(TPZCompMesh &mesh,
				  const TPZCompElHDiv<TSHAPE> &copy,
				  std::map<long,long> & gl2lcConMap,
				  std::map<long,long> & gl2lcElMap);
	
	TPZCompElHDiv();
	
	virtual ~TPZCompElHDiv();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElHDiv<TSHAPE> (mesh, *this);
	}
	
	/**
	 * @brief Create a copy of the given element. The clone copy have the connect indexes
	 * mapped to the local clone connects by the given map
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<long,long> & gl2lcConMap,std::map<long,long>&gl2lcElMap) const
	{
		return new TPZCompElHDiv<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsHDiv();
	}
	
    /** @brief Prints the relevant data of the element to the output stream */
	virtual void Print(std::ostream &out = std::cout) const;
	

	
	virtual MElementType Type();
	
	virtual int NConnects() const;
	
	virtual void SetConnectIndex(int i, long connectindex);
	
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const {
		return 0;
	}
	/** 
     * @brief return the number of shape for flux(just for flux)
	 **/
	virtual int NFluxShapeF() const;
	
	virtual int NSideConnects(int side) const;
    
	/** 
     * @brief return the local index for connect
	 **/
	virtual int SideConnectLocId(int node, int side) const;
    
    /** 
     * @brief return the local index for side
     **/
	virtual int ConnectSideLocId(int connect) const;
	
	virtual long ConnectIndex(int con) const;
    
	
	virtual void SetIntegrationRule(int ord);
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/*
     * @brief Sets the preferred interpolation order along a side \n
	 * This method only updates the datastructure of the element
	 * In order to change the interpolation order of an element, use the method PRefine
	 */
	virtual void SetPreferredOrder(int order);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int SideOrder(int side) const;
	
    /**
     * @brief return the interpolation order of the polynomial for connect
     **/
	virtual int ConnectOrder(int connect) const;
	/**
     * @brief return the number of continuous functions 
     **/
	int NShapeContinuous(TPZVec<int> &order);
    
    /// Fill the polynomial order needed from the continuous shape functions
    void FillOrder(TPZVec<int> &order) const ;
	
    /// Return the maximum order??
    virtual int MaxOrder();
    
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data);
    
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi);

	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<long> &shapeindex);
	
	/** 
	 * @brief Returns the vector index  of the first index shape associate to to each side 
	 * Special implementation to Hdiv
	 */
	void FirstShapeIndex(TPZVec<long> &Index) const;
    
	/**
     * @brief Returns a matrix index of the shape and vector  associate to element
     * @param[in] VectorSide Indicates the side associated with each vector
     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
	 */
	void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);
	
	/**
     * @brief Returns a matrix index of the shape and vector  associate to element
     * @param[in] VectorSide Indicates the side associated with each vector
     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
	 */
	void IndexShapeToVec(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);
    void IndexShapeToVec2(TPZVec<int> &VectorSide, TPZVec<int> &bilinear, TPZVec<int> &direction, TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);

	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
    
    /** @brief Compute the solution for a given variable */
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol);
	
private:
    virtual	void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
    
public:
    
    /** 
	 * @brief Compute shape functions based on master element in the classical FEM manne. 
	 * @param[in] qsi point in master element coordinates 
	 * @param[in] data stores all input data
	 */
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZMaterialData &data);
	
    void ComputeSolutionHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data);
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);	
    /** @brief Compute the solution using Hdiv structure */
	void ComputeSolutionHDiv(TPZMaterialData &data);
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	
	/** Jorge 09/06/2001
	 * @brief Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform TransformSideToElement(int side);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
		/** @brief Refinement along the element */
		virtual void PRefine(int order);
	
};

/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational prismal element for HDiv approximate space */
TPZCompEl *CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational pyramidal element for HDiv approximate space */
TPZCompEl *CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);

/** @} */

#endif
