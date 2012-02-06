/**
 * @file
 * @brief Contains declaration of TPZCompElHDivBound2 class which implements a generic computational element (HDiv scope variant).
 */

#ifndef PZELCHDIVBOUND_H_
#define PZELCHDIVBOUND_H_

#include "pzelctemp.h"

/** \addtogroup CompElement */
/** @{ */
/**
 @brief Implements a generic computational element to HDiv scope. \ref CompElement "Computational Element"
 @author Philippe Devloo
 @since Sep 29, 2009.
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivBound2 : public TPZIntelGen<TSHAPE> {
	
	/** @brief Method to append vectors */
	void Append(TPZFMatrix &u1, TPZFMatrix &u2, TPZFMatrix &u12);
public:
	
	TPZCompElHDivBound2(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	TPZCompElHDivBound2(TPZCompMesh &mesh, const TPZCompElHDivBound2<TSHAPE> &copy);

	/** 
	 * @brief Constructor used to generate patch mesh...\n 
	 * Generates a map of connect index from global mesh to clone mesh
	 */
	TPZCompElHDivBound2(TPZCompMesh &mesh,
						const TPZCompElHDivBound2<TSHAPE> &copy,
						std::map<int,int> & gl2lcConMap,
						std::map<int,int> & gl2lcElMap);
	
	/** @brief Default constructor */
	TPZCompElHDivBound2();
	/** @brief Default destructor */
	virtual ~TPZCompElHDivBound2();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
		return new TPZCompElHDivBound2<TSHAPE> (mesh, *this);
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
		return new TPZCompElHDivBound2<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
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
	
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord);
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
	virtual int PreferredSideOrder(int iside);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int SideOrder(int side) const;
	
	virtual int ConnectOrder(int connect) const;
	
	/* *transform a point in the parameter space of the side into a point in the space
     of the master element*/
	//  virtual void SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);
	
	/* *transform a point in the parameter space of the master element into a point in the
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
	
	/** @brief Returns the vector index  of the first index shape associate to element */
	/** Special implementation to Hdiv */
	void FirstShapeIndex(TPZVec<int> &Index);
	/* * return a matrix index of the shape and vector  associate to element*/
	//	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int> > & IndexShapeToVec);
	
	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi);
	
	/** @brief Compute the shape function at the integration point */
	void Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
	
	/** @brief Returns a matrix index of the shape and vector  associate to element*/
	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int> > & IndexVecShape);

	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

};


//TPZCompEl *CreateHDivBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
//TPZCompEl *CreateHDivBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
//TPZCompEl *CreateHDivBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
//TPZCompEl *CreateHDivBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

/** @} */

#endif /* PZELCHDIVBOUND_H_ */
