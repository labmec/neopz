/**
 * @file
 * @brief Contains TPZShapeLinear class which implements the shape functions of a linear one-dimensional element.
 */
 // $ Id: $
#ifndef SHAPELINEARHPP
#define SHAPELINEARHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzline.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/**
	 * @brief Implements the shape functions of a linear (1D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/** 
	 * The linear shape functions form the basis of all other shape function computations \n
	 * The range of the master element is -1,1 \n
	 * The orthogonal function which generates the linear shape functions can be modified
	 * by changing the function pointer fOrthogonal \n
	 * all static tables and functions concerning one-d elements will be grouped in this class
	 */
	class TPZShapeLinear : public pztopology::TPZLine{
		
	public:
		
		/**
		 * @brief Pointer to function which returns num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void (*fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		
#ifdef _AUTODIFF
		/**
		 *	@brief Pointer to function which returns num orthogonal functions at the point x
		 * @param x coordinate of the point with derivatives already setup
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values with derivatives
		 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		static void (*FADfOrthogonal)(FADREAL&,int ,TPZVec<FADREAL> &);
#endif
		/**
		 * @brief Chebyshev orthogonal function, computes num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void Chebyshev(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		/**
		 * @brief Legendre orthogonal function, computes num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void Legendre(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		/**
		 * @brief Legendre function computing several derivatives.
		 * @see Legendre(Real, int, TPZFMatrix, TPZFMatrix)
		 */
		static void Legendre(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int nderiv);
		
#ifdef _AUTODIFF
		/**
		 * @brief Chebyshev orthogonal function, computes num orthogonal functions at the point x
		 * @param x coordinate of the point (with derivative already setup)
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values and derivatives
		 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		static void Chebyshev(FADREAL & x,int num,TPZVec<FADREAL> &phi);
#endif
	public:
		
		/**
		 * @brief Functions which computes the shapefunctions of a one-d element
		 * @param pt coordinate of the point
		 * @param order order of the shape functions to be computed 0<= order
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 * @param id determines the orientation of the shape functions
		 */
		/**
		 * The orientation of the shapefunctions depend on the order of the id parameters
		 * if \f$ id[0] < id[1] \f$ the shapefunctions are unchanged
		 * if \f$ id[0] > id[1] \f$ the odd ordered shapefunctions are inverted
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Computes the values of the orthogonal shapefunctions before multiplying them by the
		 * corner shapefunctions
		 * @param x coordinate of the point
		 * @param ord order of the shape functions to be computed 0<= order
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 * @param transformation_index determines the transformation applied to the internal shape
		 * functions. \n This parameter is computed by the GetTransformId1d method
		 * @see GetTransformId1d
		 */
		/**
		 * The shape1dInternal function is extensively used by the shapefunction computation of
		 * the other elements
		 */
		static void ShapeInternal(TPZVec<REAL> &x,int ord,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,int transformation_index);
		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the  shape functions
		 * @param dphi (input/output) value of the derivatives of the shape functions holding the derivatives in a column
		 */
		static void ShapeGenerating(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
#ifdef _AUTODIFF
		/**
		 * @brief Computes the values of the orthogonal shapefunctions before multiplying them by the
		 * corner shapefunctions
		 * @param x coordinate of the point (with derivative already setup)
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values and derivatives
		 * @param transformation_index determines the transformation applied to the internal shape
		 * functions. \n This parameter is computed by the GetTransformId1d method
		 * @see GetTransformId1d
		 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		/**
		 * The shape1dInternal function is extensively used by the shapefunction computation of the other elements
		 */
		static void ShapeInternal(FADREAL & x,int num,TPZVec<FADREAL> & phi,int transformation_index);
		
#endif
		/**
		 * @brief Computes the transformation applied to the variational parameter of the one-d element
		 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
		 * @param in coordinate of the variational parameter
		 * @param out transformed parameter
		 */
		/**
		 * The transformation applied to compensate for odd/even brokerage between elements can also
		 * be viewed by the transformation of a variational parameter.
		 */
		static void TransformPoint1d(int transid,REAL in,REAL &out);
#ifdef _AUTODIFF
		/**
		 * @brief Computes the transformation applied to the variational parameter of the one-d element
		 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
		 * @param in coordinate of the variational parameter (with derivatives)
		 * @param out transformed parameter (with derivatives)
		 * @note REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		/*
		 * The transformation applied to compensate for odd/even brokerage between elements can also
		 * be viewed by the transformation of a variational parameter.
		 */
		static void TransformPoint1d(int transid,FADREAL & in,FADREAL &out);
#endif
		/**
		 * @brief Applies the transformation on the values of the derivatives of the shape functions of the
		 * internal shape functions
		 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the derivatives of the shapefunctions as a row vector \n
		 * the values of the derivatives contained in this matrix are modified upon return
		 */
		static void TransformDerivative1d(int transid,int num,TPZFMatrix<REAL> &in);
#ifdef _AUTODIFF
		/**
		 * Applies the transformation on the values of the derivatives of the shape functions of the
		 * internal shape functions
		 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
		 * @param num number of shapefunctions needed to transform
		 * @param in TPZVec of a differentiable type containing the values of the derivatives of the shapefunctions.
		 * @note REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		//static void TransformDerivative1d(int transid,int num,TPZVec<FADREAL> &in);
#endif
		/**
		 * @brief Computes the id of the transformation which will be applied on the parameter of the element\n
		 * @param id contains two distinct integer numbers which determine the orientation of the element
		 * @return index of the tranformation
		 */
		/**
		 * The return value is used in several methods of this class
		 */
		static int GetTransformId1d(TPZVec<int> &id);
		
		
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

