/**
 * @file
 * @brief Contains TPZShapeQuad class which implements the shape functions of a quadrilateral element.
 */
// $Id: pzshapequad.h,v 1.9 2008-03-26 20:17:34 phil Exp $

#ifndef SHAPEQUADHPP
#define SHAPEQUADHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzquadrilateral.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/// groups all classes dedicated to the computation of shape functions
namespace pzshape{

	/** 
	 * @brief Implements the shape functions of a quadrilateral (2D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/** 
	 * The quadrilateral shape functions are also used in 3D elements \n
	 * The range of the master element is -1,1
	 */
	class TPZShapeQuad  : public pztopology::TPZQuadrilateral{

	public:
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a quadrilateral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions		 
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points \n
		 * The shapefunction computation uses the shape functions of the linear element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
							  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

		/**
		 * @brief Computes the corner shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (4) shape functions
		 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the (4) shape functions
		 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/** 
		 * @brief Compute the internal functions of the quadrilateral shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 * @param quad_transformation_index determines the transformation applied to the internal shape
		 * functions. \n This parameter is computed by the GetTransformId2dQ method
		 * @see GetTransformId2dQ
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape2dQuadInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,
								  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,int quad_transformation_index);
		
#ifdef _AUTODIFF
		/**
		 * @brief Compute the internal functions of the quadrilateral shape function at a point
		 * @param x coordinate of the point (already setup with derivatives)
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values (with derivatives)
		 * @param quad_transformation_index determines the transformation applied to the internal shape
		 * functions. This parameter is computed by the GetTransformId2dQ method
		 * @see GetTransformId2dQ
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape2dQuadInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void Shape2dQuadInternal(TPZVec<FADREAL> &x, int order,
										TPZVec<FADREAL> &phi,int quad_transformation_index);
#endif

		/**
		 * @brief Transform the derivatives of num shapefunctions in place for a quadrilateral
		 * @param transid identifier of the transformation of the quad element as obtained by the GetTransformId2dQ method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the derivatives of the shapefunctions as a row vector \n
		 * the values of the derivatives contained in this matrix are modified upon return
		 */
		static void TransformDerivative2dQ(int transid, int num, TPZFMatrix<REAL> &in);
		
#ifdef _AUTODIFF
		/**
		 * Transform the derivatives of num shapefunctions in place for a quadrilateral
		 * @param transid identifier of the transformation of the quad element as obtained by the GetTransformId2dQ method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the shapefunctions and derivatives.n return
		 */
		//  static void TransformDerivative2dQ(int transid, int num, TPZVec<FADREAL> &in);
#endif
		
		/**
		 * @brief Transform the coordinates of the point in the space of the quadrilateral
		 * master element based on the transformation id
		 * @param transid identifier of the transformation of the element as obtained by the GetTransformId2dQ method
		 * @param in coordinates of the variational parameter
		 * @param out coordinates of the transformed parameter
		 */
		static void TransformPoint2dQ(int transid,TPZVec<REAL> &in,TPZVec<REAL> &out);
		
#ifdef _AUTODIFF
		
		/**
		 * @brief Transform the coordinates of the point in the space of the quadrilateral
		 * master element based on the transformation id
		 * @param transid identifier of the transformation of the element as obtained by the GetTransformId2dQ method
		 * @param in coordinates of the variational parameter (with derivatives)
		 * @param out coordinates of the transformed parameter (with derivatives)
		 */
		static void TransformPoint2dQ(int transid,TPZVec<FADREAL> &in,TPZVec<FADREAL> &out);
		
#endif
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		
		static void ProjectPoint2dQuadToRib(int rib, TPZVec<REAL> &in, REAL &out);
		/**
		 * @brief Method which identifies the quadrilateral transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point
		 */
		static int GetTransformId2dQ(TPZVec<int> &id);

		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the two dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (2,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions
		 */
		static void TransformDerivativeFromRibToQuad(int rib,int num,TPZFMatrix<REAL> &dphi);

#ifdef _AUTODIFF

		/**
		 * Transforms the derivative of a shapefunction computed on the rib into the two dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (2,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @num number of shapefunction derivatives which need to be transformed
		 * @phi values of the shapefunctions (with derivatives)
		 */
		//  static void TransformDerivativeFromRibToQuad(int rib,int num,TPZVec<FADREAL> &phi);
#endif

		/** @brief Data structure which defines the quadrilateral transformations */
		static REAL gTrans2dQ[8][2][2];
		/** @brief Data structure which defines the quadrilateral transformations */
		static REAL gFaceTr2dQ[6][2][3];
		/** @brief Data structure which defines the quadrilateral transformations */
		static REAL gRibTrans2dQ1d[4][2];

		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NConnectShapeF(int side, int order);

		/**
		 * @brief Total number of shapefunctions, considering the order
		 * of interpolation of the element
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NShapeF(TPZVec<int> &order);
		
	};
	
};

#endif
