//
//  pzhdivpressurebound.h
//  PZ
//
//  Created by Agnaldo Farias on 25/06/13.
//
//

/**
 * @file
 * @brief Contains declaration of TPZCompElHDivPressureBound class which implements a generic computational element (HDiv-Pressure scope variant).
 */


#ifndef __PZ__pzhdivpressurebound__
#define __PZ__pzhdivpressurebound__

#include <iostream>

#include "pzelchdivbound2.h"

/** \addtogroup CompElement */
/** @{ */
/**
 * @brief Implements a generic computational element to HDiv-Pressure scope. \ref CompElement "Computational Element"
 * @author Agnaldo Farias
 * @since Jun 25, 2013.
 */
/**
 * By varying the classes passed as template arguments, the complete family of computational elements are implemented
 */
template<class TSHAPE>
class TPZCompElHDivPressureBound : public TPZCompElHDivBound2<TSHAPE> {

    /** @brief Defines the interpolation order for pressure variable*/
	int fPressureOrder;
    
public:
    
    TPZCompElHDivPressureBound(TPZCompMesh &mesh, TPZGeoEl *gel);
    
    
    /** @brief Default constructor */
	TPZCompElHDivPressureBound();
    
    /**
	 * @brief Constructor used to generate patch mesh...\n
	 * Generates a map of connect index from global mesh to clone mesh
	 */
	TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivPressureBound<TSHAPE> &copy, std::map<int64_t,int64_t> & gl2lcConMap, std::map<int64_t,int64_t> & gl2lcElMap);
    
    
    TPZCompElHDivPressureBound(TPZCompMesh &mesh, const TPZCompElHDivPressureBound<TSHAPE> &copy);
   
    
	/** @brief Default destructor */
	virtual ~TPZCompElHDivPressureBound();
	
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const  override {
		return new TPZCompElHDivPressureBound<TSHAPE> (mesh, *this);
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
		return new TPZCompElHDivPressureBound<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
	}

    
    
    virtual int NConnects() const override;
    
    /** @brief Identifies the interpolation order for pressure variable*/
	virtual void SetPressureOrder(int ord);
    
    virtual void SetConnectIndex(int i, int64_t connectindex) override;
    
    virtual int NConnectShapeF(int connect, int order) const override;
    
    virtual int SideConnectLocId(int node, int side) const override;
    
    /** @brief Set create function in TPZCompMesh to create elements of this type */
	virtual void SetCreateFunctions(TPZCompMesh *mesh) override;
	
    //	virtual MElementType Type();
	
	
	virtual int Dimension() const  override {
		return TSHAPE::Dimension;
	}
	
	virtual int NCornerConnects() const  override {
		return 0;
	}
	
//	virtual int NSideConnects(int side) const;
	
	/** @brief Identifies the interpolation order on the interior of the element*/
	virtual void GetInterpolationOrder(TPZVec<int> &ord) override;
	
	/** @brief Returns the preferred order of the polynomial along side iside*/
    //	virtual int PreferredSideOrder(int iside);
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order) override;
	
	/** @brief Returns the actual interpolation order of the polynomial along the side*/
	virtual int EffectiveSideOrder(int side) const override;
	
	virtual int ConnectOrder(int connect) const override;
    
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data) override;
	
	/** @brief Compute the correspondence between the normal vectors and the shape functions */
	void ComputeShapeIndex(TPZVec<int> &sides, TPZVec<int64_t> &shapeindex);
	
	/** @brief Returns the vector index  of the first index shape associate to element */
	/** Special implementation to Hdiv */
	void FirstShapeIndex(TPZVec<int64_t> &Index);
	
	/** @brief Compute the values of the shape function of the side*/
	virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) override;
	
	/** @brief Compute the shape function at the integration point */
	void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) override;
	
	/** @brief Returns a matrix index of the shape and vector  associate to element*/
	void IndexShapeToVec(TPZVec<int> &fVectorSide,TPZVec<std::pair<int,int64_t> > & IndexVecShape);
    
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
int ClassId() const override;

    
	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
    void Read(TPZStream &buf, void *context) override;
};

template<class TSHAPE>
int TPZCompElHDivPressureBound<TSHAPE>::ClassId() const{
    return Hash("TPZCompElHDivPressureBound") ^ TPZCompElHDivBound2<TSHAPE>::ClassId() << 1;
}

#endif /* defined(__PZ__pzhdivpressurebound__) */
