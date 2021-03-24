/**
 * @file
 * @brief Contains TPZShapePiram class which implements the shape functions of a pyramid element.
 * author: Nathan Shauer
 * date: 01/04/2015
 */

#ifndef SHAPEPIRAMHDIV_H
#define SHAPEPIRAMHDIV_H

#include "pzshapepiram.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {

	/** 
	 * @brief Implements the shape functions of a pyramid (3D) element for Hdiv space. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is -1,1 for the base and [0,1] for the height
	 */
	class TPZShapePiramHdiv  : public TPZShapePiram{

    public:
		
        int ClassId() const override;
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a pyramid element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions   
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
        /**
         * @brief Number of shapefunctions of the connect associated with the side, considering the order
         * of interpolation of the element
         * @param side associated side
         * @param order vector of integers indicating the interpolation order of the element
         * @return number of shape functions
         */
        static int NConnectShapeF(int side, int order);
	};
        
};

#endif
