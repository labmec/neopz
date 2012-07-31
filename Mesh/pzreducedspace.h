//
//  pzreducedspace.h
//  PZ
//
//  Created by Philippe Devloo on 7/30/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef PZ_pzreducedspace_h
#define PZ_pzreducedspace_h

#include "pzinterpolationspace.h"

class TPZReducedSpace : public TPZInterpolationSpace
{
    /** @brief Default constructor */
	TPZReducedSpace();
	
	/** @brief Default destructor */
	virtual ~TPZReducedSpace();
	
	/** @brief Puts a copy of the element in the referred mesh */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy);
	
	/** @brief Puts a copy of the element in the patch mesh */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, std::map<int,int> &gl2lcElMap);
	
	/** @brief Copy of the element in the new mesh whit alocated index */
	TPZReducedSpace(TPZCompMesh &mesh, const TPZReducedSpace &copy, int &index);
	
	/**
	 * @brief Create a computational element within mesh
	 * @param mesh mesh wher will be created the element
	 * @param gel geometrical element to insert
	 * @param index new elemen index
	 */
	/** Inserts the element within the data structure of the mesh */
	TPZReducedSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const 
    {
        return 1;
    }
	
	/** @brief It returns the shapes number of the element */
	virtual int NShapeF() const;
	
	/** @brief Returns the number of shapefunctions associated with a connect*/
	virtual int NConnectShapeF(int inod) const;
	
	/** @brief Returns the max order of interpolation. */
	virtual int MaxOrder();
	
	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	
	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphix matrix of derivatives of shapefunctions, dimension (dim,numshape)
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	virtual void ShapeX(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphix, TPZFMatrix<REAL> &axes);

	/** 
	 * @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data);
	
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi);

    /**
     * @brief Computes solution and its derivatives in local coordinate qsi
     * @param qsi master element coordinate
     * @param phi matrix containing shape functions compute in qsi point
     * @param dphix matrix containing the derivatives of shape functions in the direction of the axes
     * @param axes [in] axes indicating the direction of the derivatives
     * @param sol finite element solution
     * @param dsol solution derivatives
     */
    void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphix,
                         const TPZFMatrix<REAL> &axes, TPZSolVec &sol, TPZGradSolVec &dsol);
    
	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/** @brief Initialize element matrix in which is computed in CalcResidual */
	void InitializeElementMatrix(TPZElementMatrix &ef);
	
	/** @brief Save the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Read the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
private:

    TPZInterpolationSpace *ReferredIntel();
};


#endif
