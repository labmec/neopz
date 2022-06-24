/**
 * @file
 * @brief Contains TPZShapeLinear class which implements the shape functions of a linear one-dimensional element.
 */

#ifndef SHAPELINEARHPP
#define SHAPELINEARHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzline.h"
#include "pzshtmat.h"

#include "fadType.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	

 
/**
 * @brief Implements the shape functions of a linear (1D) element. \ref shape "Shape"
 * @ingroup shape
 */
///
///  The linear shape functions form the basis of all other shape function computations \n
///  The range of the master element is -1,1 \n
///  The orthogonal function which generates the linear shape functions can be modified
///  by changing the function pointer fOrthogonal \n
///  all static tables and functions concerning one-d elements will be grouped in this
///
	class TPZShapeLinear : public pztopology::TPZLine{
		
	public:

		/** @{
		 * @name Orthogonal polynomials family 
		 */
		
		/**
		 * @brief Pointer to function which returns num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void (*fOrthogonal)(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		/**
		 * @brief Chebyshev orthogonal polynomial, computes num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void Chebyshev(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
        /**
         * @brief Exponential polynomial, computes num orthogonal functions at the point x
         * @param x coordinate of the point
         * @param num number of shape functions to be computed
         * @param phi shapefunction values
         * @param dphi values of the derivatives of the shape functions
         */
        static void Expo(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
        
		/**
		 * @brief Legendre orthogonal polynomial, computes num orthogonal functions at the point x
		 * @param x coordinate of the point
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		static void Legendre(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		/**
		 * @brief Jacobi orthogonal polynomials
		 */
		static void Jacobi(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
        
		/**
		 * @brief Jacobi parameters orthogonal polynomials
		 */
        static REAL JacobiA;
        static REAL JacobiB;
        
		/**
		 * @brief Hermite orthogonal polynomials
		 */
		static void Hermite(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi);
		
		/**
		 * @brief Legendre function computing several derivatives.
		 * @see Legendre(Real, int, TPZFMatrix, TPZFMatrix)
		 */
		static void Legendre(REAL x,int num,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int nderiv);

		/**
		 *	@brief Pointer to function which returns num orthogonal functions at the point x
		 * @param x coordinate of the point with derivatives already setup
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values with derivatives
		 * REMARK: The Derivative classes MUST store at least one derivative - 1d problem
		 */
        static void (*FADfOrthogonal)(const FADREAL& x,int num,TPZFMatrix<FADREAL> & phi,TPZFMatrix<FADREAL> & dphi);
		/**
		 * @brief Chebyshev orthogonal polynomial, computes num orthogonal functions at the point x
		 * @param x coordinate of the point (with derivative already setup)
		 * @param num number of shape functions to be computed
		 * @param phi shapefunction values and derivatives
		 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
		 */
		static void Chebyshev(const FADREAL & x,int num,TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi);

		/** @} */
		
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
		static void Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
        template<class T>
        static void ShapeCorner(const TPZVec<T> &pt,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
            phi(0,0) = (1-pt[0])/2.;
            phi(1,0) = (1+pt[0])/2.;
            dphi(0,0) = -0.5;
            dphi(0,1)= 0.5;
        }

		
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
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
         * @brief Computes the values of the orthogonal shapefunctions before multiplying them by the
         * corner shapefunctions
         * @param x coordinate of the point
         * @param ord order of the shape functions to be computed 0<= order
         * @param phi shapefunction values
         * @param dphi values of the derivatives of the shape functions
         * @param transformation_index = 0;
         * functions. \n This parameter is computed by the GetTransformId1d method
         * @see GetTransformId1d
         */
        /**
         * The shape1dInternal function is extensively used by the shapefunction computation of
         * the other elements
         */
        static void ShapeInternal(TPZVec<REAL> &x,int ord,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {    //versao nova-> o ponto vem transformado
            // Quadratic or higher shape functions
            int num = ord-1;
            if(num <= 0) return;
            REAL y;
            fOrthogonal(x[0],num,phi,dphi);
        }
        static void ShapeInternal(TPZVec<Fad<REAL>> &x,int ord,TPZFMatrix<Fad<REAL>> &phi,TPZFMatrix<Fad<REAL>> &dphi)
        {    //versao nova-> o ponto vem transformado
            // Quadratic or higher shape functions
            int num = ord-1;
            if(num <= 0) return;
            REAL y;
            FADfOrthogonal(x[0],num,phi,dphi);
        }

            

		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the  shape functions
		 * @param dphi (input/output) value of the derivatives of the shape functions holding the derivatives in a column
		 */
        template<class T>
		static void ShapeGenerating(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
                
                phi(2,0) = 4.*phi(0,0)*phi(1,0);
                dphi(0,2) = 4.*(dphi(0,0)*phi(1,0)+phi(0,0)*dphi(0,1));
                
        }

        /**
         * @brief Computes the generating shape functions for a quadrilateral element
         * @param pt (input) point where the shape function is computed
         * @param phi (input/output) value of the  shape functions
         * @param dphi (input/output) value of the derivatives of the shape functions holding the derivatives in a column
         */
        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            ShapeGenerating(pt, phi, dphi);
        }

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
		static void ShapeInternal(FADREAL & x,int num,TPZFMatrix<FADREAL> & phi,int transformation_index);

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
		/**
		 * @brief Applies the transformation on the values of the derivatives of the shape functions of the
		 * internal shape functions
		 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the derivatives of the shapefunctions as a row vector \n
		 * the values of the derivatives contained in this matrix are modified upon return
		 */
		static void TransformDerivative1d(int transid,int num,TPZFMatrix<REAL> &in);

		/**
		 * @brief Computes the id of the transformation which will be applied on the parameter of the element\n
		 * @param id contains two distinct integer numbers which determine the orientation of the element
		 * @return index of the tranformation
		 */
		/**
		 * The return value is used in several methods of this class
		 */
		static int GetTransformId1d(TPZVec<int64_t> &id);

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
        static void TransformPoint1d(int transid,double &val) ;
        static TPZTransform<REAL> ParametricTransform(int transid);
        
        // compute internal shape function of the side, x was already projected
        static void ShapeInternal(int side, TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            if (side != 2) {
                phi(0,0) = 1.;
                return;
            }
            ShapeInternal(x, order, phi, dphi);
        }
        static void ShapeInternal(int side, TPZVec<FADREAL> &x, int order,TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi)
        {
            if (side != 2) {
                phi(0,0) = 1.;
                return;
            }
            ShapeInternal(x, order, phi, dphi);
        }

		
	};
	
};

#endif

