//
//  pzcompelhdivcurved.h
//  PZ
//
//  Created by Douglas Castro on 6/25/14.
//
//

#ifndef __PZ__pzcompelhdivcurved__
#define __PZ__pzcompelhdivcurved__

#include <iostream>
#include "pzelctemp.h"
#include "pzelchdiv.h"

template<class TSHAPE>
class TPZCompElHdivCurved : public TPZCompElHDiv<TSHAPE>  //????????????????????????
{
    
    // computematerialdata
    // initmaterialdata
    
public:
    
    TPZCompElHdivCurved(TPZCompMesh &mesh, TPZGeoEl *gel, long &index);
	
	TPZCompElHdivCurved();
	
	virtual ~TPZCompElHdivCurved();
	
    /** @brief Set create function in TPZCompMesh to create elements of this type */
//	virtual void SetCreateFunctions(TPZCompMesh *mesh){
//		mesh->SetAllCreateFunctionsHDivCurved(); // Implementar ?????? veja hdivfull
//	}
	
	virtual MElementType Type();
	
    /**
     * @brief Number of shapefunctions of the connect associated
     * @param connect connect number
     * @return number of shape functions
     */
	virtual int NConnectShapeF(int connect) const;
	
	/**
     * @brief return the number of shape for flux(just for flux)
	 **/
	virtual int NFluxShapeF() const;
	
	/** @brief Sets the interpolation order of side to order*/
	virtual void SetSideOrder(int side, int order);
	
	/**
     * @brief return the number of continuous functions
     **/
	virtual void NShapeContinuous(TPZVec<int> &order, int &nshape  );
    
	
	/** @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions */
	virtual void InitMaterialData(TPZMaterialData &data);
    
	/**
	 * @brief Returns the vector index  of the first index shape associate to to each side
	 * Special implementation to Hdiv
	 */
	virtual void FirstShapeIndex(TPZVec<long> &Index);
    
	/**
     * @brief Returns a matrix index of the shape and vector  associate to element
     * @param[in] VectorSide Indicates the side associated with each vector
     * @param[out] IndexVecShape Indicates for the vector/shape function for the approximation space
	 * @param[in] pressureorder Order of the pressure (to select shape functions?)
	 */
	virtual void IndexShapeToVec(TPZVec<int> &VectorSide,TPZVec<std::pair<int,long> > & IndexVecShape, int pressureorder);
	
	void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
    
    virtual void Shape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
    

};

/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational cube element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational prismal element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational pyramidal element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational tetrahedral element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational point element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedBoundPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational linear element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedBoundLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational quadrilateral element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedBoundQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);
/** @brief Creates computational triangular element for HDiv approximate space */
TPZCompEl *CreateHDivCurvedBoundTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,long &index);



#endif /* defined(__PZ__pzcompelhdivcurved__) */
