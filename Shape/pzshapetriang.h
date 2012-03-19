/**
 * @file
 * @brief Contains TPZShapeTriang class which implements the shape functions of a triangular element.
 */
// $ Id: $
#ifndef SHAPETRIANGHPP
#define SHAPETRIANGHPP

#include "tpztriangle.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a triangular (2D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The triangular shape functions are also used in 3D elements \n
	 * The range of the master element is 0,1
	 */
	class TPZShapeTriang : public pztopology::TPZTriangle {
		
	public:

		/**
		 * @brief Computes the values of the shape functions and their derivatives for a triangular element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (4 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions		 
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
							  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the corner shape functions for a triangular element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (3) shape functions
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
		 * @brief Compute the internal functions of the triangle shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 * @param triangle_transformation_index determines the transformation applied to the internal shape
		 * functions. \n This parameter is computed by the GetTransformId2dT method
		 * @see GetTransformId2dQ
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions. \n
		 * Shape2dTriangleInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,int triangle_transformation_index);

		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out);

		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the two dimensional derivative \n
		 * of the function with respect to the element. The parameter dphi should be dimensioned (2,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions
		 */
		static void TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Method which identifies the triangular transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point
		 */
		static int GetTransformId2dT(TPZVec<int> &id);
		
		/**
		 * @brief Transform the coordinates of the point in the space of the triangle
		 * master element based on the transformation id
		 * @param transid identifier of the transformation of the element as obtained by the GetTransformId2dT method
		 * @param in coordinates of the variational parameter
		 * @param out coordinates of the transformed parameter
		 */
		static void TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Transform the derivatives of num shapefunctions in place for a triangle
		 * @param transid identifier of the transformation of the triangle element as obtained by the GetTransformId2dT method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the derivatives of the shapefunctions as a row vector \n
		 * the values of the derivatives contained in this matrix are modified upon return
		 */
		static void TransformDerivative2dT(int transid, int num, TPZFMatrix<REAL> &in);
		
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gTrans2dT[6][2][2];
		
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gVet2dT[6][2];
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gRibTrans2dT1d[3][2];
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gVet1dT[3];
		
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
