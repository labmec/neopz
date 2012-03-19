/**
 * @file
 * @brief Contains TPZShapePrism class which implements the shape functions of a prism element.
 */
// $ Id: $
#ifndef SHAPEPRISMHPP
#define SHAPEPRISMHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzprism.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a prism (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is \f$ [-1,1] \f$ for the \f$ z \f$ direction and \f$ [0,1] \f$ for the triangular plane
	 */
	class TPZShapePrism : public pztopology::TPZPrism {
		
	public:
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a prism element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and quadrilateral and triangular element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
							  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the corner shape functions for a prism element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (6) shape functions
		 * @param dphi (output) value of the derivatives of the (6) shape functions holding the derivatives in a column
		 */
		static void CornerShape(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input) value of the (4) shape functions
		 * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/** 
		 * @brief Compute the internal functions of the prism shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dPrismaInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
								  TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dPrismaToRib(int rib, TPZVec<REAL> &in, REAL &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dPrismaToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
		
		/**
		 * @brief Transforms a point on the face by the corresponding transformation
		 * @param transid transformation index of the face
		 * @param face face to which the coordinates refer
		 * @param in coordinate of the point on the face before transformation
		 * @param out coordinate of the point after transformation
		 */
		/**
		 * This method applies two dimensional transformations which are implemented by the classes \n
		 * associated with triangles and quadrilaterals
		 */
		static void TransformPoint3dPrismaFace(int transid, int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToPrisma(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (3,num), at least
		 * @param face rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToPrisma(int face,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transform the derivatives of the shapefunction on the shape (i.e. two dimensional derivative) to acount
		 * for the transformation on the corresponding face.
		 * @param transid id of the transformation which needs to be applied to the face
		 * @param face face to which the coordinates of the point refers
		 * @param num number of derivatives of shapefunctions which need to be transformed
		 * @param in matrix of derivatives of dimension (2,num)
		 */
		static void TransformDerivativeFace3dPrisma(int transid, int face, int num, TPZFMatrix<REAL> &in);
		
		/** @brief Data structure which defines the hexahedral transformations and topology */
		static REAL gFaceTrans3dPrisma2d[5][2][3];
		
		/** @brief Data structure which defines the hexahedral transformations and topology */
		static REAL gFaceSum3dPrisma2d[5][2];
		
		/** @brief Data structure which defines the hexahedral transformations and topology */
		static REAL gRibTrans3dPrisma1d[9][3];
		
		/** @brief Data structure which defines the hexahedral transformations and topology */
		static REAL gRibSum3dPrisma1d[9];
		
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
