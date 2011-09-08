/**
 * @file
 * @brief Contains TPZShapeCube class which implements the shape functions of a hexaedral element.
 */
// $Id: pzshapecube.h,v 1.8 2008-03-26 20:17:33 phil Exp $
#ifndef SHAPECUBEHPP
#define SHAPECUBEHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzcube.h"
//#include "pzgeoelside.h"


#ifdef _AUTODIFF
#include "fadType.h"
#endif
/** @brief Groups all classes dedicated to the computation of shape functions */
namespace pzshape {
	
	/**
	 * @brief Implements the shape functions of a hexahedral (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/** 
	 * The range of the master element is [-1 ,1]
	 */
	class TPZShapeCube : public pztopology::TPZCube {
		
	public:
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a hexahedral element
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
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi);
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi);
		
#ifdef _AUTODIFF
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a hexahedral element
		 * @param point (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions and derivatives
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
		 */
		static void ShapeCube(TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZVec<FADREAL> &phi);
#endif
		
		/**
		 * @brief Computes the corner shape functions for a hexahedral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (8) shape functions
		 * @param dphi (output) value of the derivatives of the (8) shape functions holding the derivatives in a column
		 */
		static void ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
		
	private:
		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input) value of the (4) shape functions
		 * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
		
#ifdef _AUTODIFF
		/**
		 * @brief Computes the corner shape functions for a hexahedral element
		 * @param pt (input) point where the shape function is computed already setup with derivatives
		 * @param phi (output) value of the (8) shape functions and derivatives
		 */
		static void ShapeCornerCube(TPZVec<FADREAL> &pt, TPZVec<FADREAL> &phi);
#endif
		
		/**
		 * @brief Compute the internal functions of the hexahedral shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation \n
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
								  TPZFMatrix &dphi);//,int quad_transformation_index
#ifdef _AUTODIFF
		/**
		 * @brief Compute the internal functions of the hexahedral shape function at a point
		 * @param x coordinate of the point (with derivatives already setup)
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values (and derivatives)
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner shape functions\n
		 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation \n
		 * determined by the transformation index
		 */
		static void Shape3dCubeInternal(TPZVec<FADREAL> &x, int order,TPZVec<FADREAL> &phi);//,int quad_transformation_index
#endif
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeToRib(int face, TPZVec<REAL> &in, REAL &outval);
		
#ifdef _AUTODIFF
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element already setup with derivatives
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeToRib(int face, TPZVec<FADREAL> &in, FADREAL &outval);
#endif
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeSide(int face, TPZVec<REAL> &in, REAL &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinates of the point on the face
		 */
		static void ProjectPoint3dCubeFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dCubeToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
#ifdef _AUTODIFF
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element (with derivatives)
		 * @param outval coordinates of the point on the face (with derivatives)
		 */
		static void ProjectPoint3dCubeToFace(int face, TPZVec<FADREAL> &in, TPZVec<FADREAL> &outval);
#endif
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToCube(int rib,int num,TPZFMatrix &dphi);
#ifdef _AUTODIFF
		/**
		 * Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param phi values of the derivatives of the shapefunctions (modified in place)
		 */
		//static void TransformDerivativeFromRibToCube(int rib,int num,TPZVec<FADREAL> &phi);
#endif
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToCube(int rib,int num,TPZFMatrix &dphi);
		
#ifdef _AUTODIFF
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param phi values of the shapefunctions (modified in place) with derivatives.
		 */
		//static void TransformDerivativeFromFaceToCube(int rib,int num,TPZVec<FADREAL> &phi);
#endif
		
		/** @brief Data structure which defines the hexahedral transformations */
		static REAL gFaceTrans3dCube2d[6][2][3];
		/** @brief Data structure which defines the hexahedral transformations */
		static REAL gRibTrans3dCube1d[12][3];
		
	public:
		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NConnectShapeF(int side, int order);
		
		/**
		 * @brief Total number of shapefunctions, considering the order of interpolation of the element
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NShapeF(TPZVec<int> &order);

	};
	
};

#endif
