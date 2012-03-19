/**
 * @file
 * @brief Contains declaration of TPZCompElHDiv class which implements a generic computational element (HDiv scope).
 */
//$Id: pzelctemp.h,v 1.14 2008-10-08 02:13:33 phil Exp $

#ifndef PZELCHDIVHTT
#define PZELCHDIVHTT

#include "pzelctemp.h"


/**
 * @brief This class implements a "generic" computational element to HDiv scope. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
/** 
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDiv : public TPZIntelGen<TSHAPE> {
	/** @brief Defines the interpolation order for pressure variable*/
	int fPressureOrder;
	
	/** @brief To append vectors */
	void Append(TPZFMatrix<REAL> &u1, TPZFMatrix<REAL> &u2, TPZFMatrix<REAL> &u12);
public:
	
	TPZCompElHDiv(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	TPZCompElHDiv(TPZCompMesh &mesh, const TPZCompElHDiv<TSHAPE> &copy);
	
	/**
	 * @brief Constructor used to generate patch mesh... generates a map of connect index from
	 * global mesh to clone mesh
	 */
	TPZCompElHDiv(TPZCompMesh &mesh,
				  const TPZCompElHDiv<TSHAPE> &copy,
				  std::map<int,int> & gl2lcConMap,
				  std::map<int,int> & gl2lcElMap);
	
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
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const
	{
		return new TPZCompElHDiv<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}
	
    /** @brief Set create function in TPZCompMesh to create elements of this type
	 */
	virtual void SetCreateFunctions(TPZCompMesh *mesh){
		mesh->SetAllCreateFunctionsHDiv();
	}

	
	virtual MElementType Type();
	
	virtual int NConnects() const;
	
	virtual void SetConnectIndex(int i, int connectindex);
	
	virtual int NConnectShapeF(int connect) const;
	
	virtual int Dimension() const {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const {
		return 0;
	}
	
	virtual int NSideConnects(int side) const;
	
	virtual int SideConnectLocId(int node, int side) const;
	
	virtual int ConnectIndex(int node) const;
    
    /** @brief returns the index of the pressure connect
     * returns -1 if their is no pressure connect
     */
    virtual int PressureConnectIndex() const
    {
        return NConnects()-1;
    }

	
	virtual void SetIntegrationRule(int ord);
	/** @brief Identifies the interpolation order for pressure variable*/
	virtual void SetPressureOrder(int ord);
	/** @brief Returns the interpolation order to dual variable */
	int DualOrder();
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/** @brief Sets the preferred interpolation order along a side
	 
	 This method only updates the datastructure of the element
	 In order to change the interpolation order of an element, use the method PRefine*/
	virtual void SetPreferredOrder(int order);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int SideOrder(int side) const;
	
	virtual int ConnectOrder(int connect) const;
	
	/** Transform a point in the parameter space of the side into a point in the space
     of the master element*/
	//  virtual void SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
	
	/** Transform a point in the parameter space of the master element into a point in the
     space of the side*/
	//  virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);
	
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data);
	
	/**
	 * @brief Compute the correspondence between the normal vectors and the shape functions
	 */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int> &shapeindex);
	
	/** 
	 * @brief Returns the vector index  of the first index shape associate to element 
	 *
	 * Special implementation to Hdiv
	 */
	void FirstShapeIndex(TPZVec<int> &Index);
	/** @brief Returns a matrix index of the shape and vector  associate to element
     * @param VectorSide input : indicates the side associated with each vector
     * @param IndexVecShape output : indicates for the vector/shape function for the approximation space
     */
	void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,int> > & IndexVecShape);
	
	/** @brief Computes the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/** 
	 * @brief Compute the shape functions corresponding to the dual space
	 */
	virtual void ShapeDual(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
	
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
    
    ///Compute the solution for a given variable
	virtual void Solution( TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);
    
    
    virtual	void ComputeSolution(TPZVec<REAL> &qsi, TPZSolVec &sol, TPZGradSolVec &dsol,TPZFMatrix<REAL> &axes);
	
    void ComputeSolutionHDiv(TPZVec<REAL> &qsi, TPZMaterialData &data);
    virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                                 const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);	
    /**   
	 * @brief Compute the solution using Hdiv structure
	 */
	void ComputeSolutionHDiv(TPZMaterialData &data);
	
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	//virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);
	
	/** Jorge 09/06/2001
	 * @brief Returns the transformation which transform a point from the side to the interior of the element
	 */
	TPZTransform TransformSideToElement(int side);
	
	/**
	 * @brief Returns the unique identifier for reading/writing objects to streams
	 */
	virtual int ClassId() const;
	/**
	 @brief Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/**
	 @brief Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
};

/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateHDivPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational prismal element for HDiv approximate space */
TPZCompEl *CreateHDivPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational pyramidal element for HDiv approximate space */
TPZCompEl *CreateHDivPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif
