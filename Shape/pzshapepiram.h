/**
 * @file
 * @brief Contains TPZShapePiram class which implements the shape functions of a pyramid element.
 */

#ifndef SHAPEPIRAMHPP
#define SHAPEPIRAMHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzpyramid.h"
#include "pzshtmat.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a pyramid (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is -1,1 for the base and [0,1] for the height
	 */
	class TPZShapePiram  : public pztopology::TPZPyramid{
		
	public:
		
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
        
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the side associated with each shapefunction
         */
        static void ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders);//, TPZVec<int64_t> &sides
        
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the internal shapefunctions of a side
         * @param sides is a vector with copy of side as much as needed, it depends on the order
         */
        static void SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders);
        
		/** 
		 * @brief Compute the internal functions of the pyramid shape function at a point\n
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
								  TPZFMatrix<REAL> &dphi);
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the (4) shape functions
		 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
        /**
         * @brief Computes the generating shape functions for a quadrilateral element
         * @param pt (input) point where the shape function is computed
         * @param phi (input/output) value of the (4) shape functions
         * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
        
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dPiramToRib(int rib, TPZVec<REAL> &in, REAL &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dPiramToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
		
		/**
		 * @brief Transforms a point on the face by the corresponding transformation
		 * @param transid transformation index of the face
		 * @param face face to which the coordinates refer
		 * @param in coordinate of the point on the face before transformation
		 * @param out coordinate of the point after transformation
		 */
		/**
		 * This method applies two dimensional transformations which are implemented by the classes
		 * associated with triangles and quadrilaterals
		 */
		static void TransformPoint3dPiramFace(int transid, int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToPiram(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param face face to which the coordinates refer
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToPiram(int face,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transform the derivatives of the shapefunction on the shape (i.e. two dimensional derivative) to acount
		 * for the transformation on the corresponding face
		 * @param transid id of the transformation which needs to be applied to the face
		 * @param face face to which the coordinates of the point refers
		 * @param num number of derivatives of shapefunctions which need to be transformed
		 * @param in matrix of derivatives of dimension (2,num)
		 */
		static void TransformDerivativeFace3dPiram(int transid, int face, int num, TPZFMatrix<REAL> &in);
		/** @brief Data structure which defines the pyramid transformations and topology */
		static REAL gFaceTrans3dPiram2d[5][2][3];
		
		/** @brief Data structure which defines the pyramid transformations and topology */
		static REAL gRibTrans3dPiram1d[8][3];
		
		/** @brief Data structure which defines the pyramid transformations and topology */
		static REAL gFaceSum3dPiram2d[5][2];
		
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
		static int NShapeF(const TPZVec<int> &order);
		
        int ClassId() const override;
        static void ShapeInternal(int side, TPZVec<REAL> &x, int order,  TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
        static void ShapeCorner(const TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
	};
	
};

#endif
