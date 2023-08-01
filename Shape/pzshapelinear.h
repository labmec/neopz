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
        using TTOPOL=pztopology::TPZLine;
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
		

        template<class T>
        static void ShapeCorner(const TPZVec<T> &pt,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
            phi(0,0) = (1-pt[0])/2.;
            phi(1,0) = (1+pt[0])/2.;
            dphi(0,0) = -0.5;
            dphi(0,1)= 0.5;
        }

		
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the side associated with each shapefunction
         */
        static void ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders);//, TPZVec<int64_t> &sides
        
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the side associated with each shapefunction
         */
        static void InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders);

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
        template<class T>
        static void ShapeInternal(TPZVec<T> &x,int ord,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi)
        {    //versao nova-> o ponto vem transformado
            // Quadratic or higher shape functions
            int num = ord-1;
            if(num <= 0) return;
            REAL y;
            if constexpr (std::is_same_v<FADREAL, T>)
            {
                TPZShapeLinear::FADfOrthogonal(x[0],num,phi,dphi);
            }
            else{
                TPZShapeLinear::fOrthogonal(x[0],num,phi,dphi);
            }
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
        
	};
	
};

#endif

