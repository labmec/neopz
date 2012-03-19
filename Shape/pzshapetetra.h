/**
 * @file
 * @brief Contains TPZShapeTetra class which implements the shape functions of a tetrahedral element.
 */
// $Id: pzshapetetra.h,v 1.6 2008-03-26 20:17:34 phil Exp $
#ifndef SHAPETETRAHPP
#define SHAPETETRAHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpztetrahedron.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a tetrahedral (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is 0,1
	 */
	class TPZShapeTetra : public pztopology::TPZTetrahedron{
		
	public:

		/**
		 * @brief Computes the values of the shape functions and their derivatives for a tetrahedral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (11 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and triangular element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the corner shape functions for a tetrahedral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (4) shape functions
		 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void CornerShape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the (4) shape functions
		 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		
		/** 
		 * @brief Compute the internal functions of the tetrahedral shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dTetraInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index (also on the faces)
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dTetraToRib(int rib, TPZVec<REAL> &in, REAL &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param side side to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint3dTetrSide(int side, TPZVec<REAL> &in, REAL &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dTetraToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinates of the point on the face
		 */
		static void ProjectPoint3dTetrFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToTetra(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param face face index to which the point should be projected
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToTetra(int face,int num,TPZFMatrix<REAL> &dphi);
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gVet1dTetr[6];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceTrans3dTetr2d[4][2][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gVet2dTetr[4][2];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibTrans3dTetr1d[6][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceSum3dTetra2d[4][2];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceTrans3dTetra2d[4][2][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibSum3dTetra1d[6];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibTrans3dTetra1d[6][3];
		
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
