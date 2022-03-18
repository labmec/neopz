/** 
 * @file
 * @brief Contains the declaration of the shape function discontinuous
 */

#ifndef SHAPEDISCHPP
#define SHAPEDISCHPP

#include "pzfmatrix.h"
#include "pzvec.h"

class TPZShapeData;

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
    
    
    /** 
	 * @brief Implements the shape functions discontinuous elements. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
     * This class computes the shape functions for n-dimensional elements
     * the shape functions can be tensor based or interpolation order based
     * The Shape2dFull method also computes higher order derivatives
     */
    class TPZShapeDisc {
        
    public:
        
        /**
		 * @brief Polynomial shape function
         * @param C:      normalizing constant, it provided by the discontinous element
         * @param x0:     interior point of the deformed element, it can be, for example, the center of mass of the element
         * @param x:      coordinate of the point
         * @param degree: degree of interpolation of the element
         * @param phi:    shapefunction values
         * @param dphi:   values of the derivatives of the shape functions
         * @param n:      number of derivatives to be computed
         */
        static void Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n = 1);
        
        static void PolynomialWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n = 1);
        
        static void IntegratedLegendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, int n = 1);
        
        static void Legendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n = 1);
        
        static void ChebyshevWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n = 1);
        
        static void LegendreWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix<REAL> & phi,TPZFMatrix<REAL> & dphi, int n = 1);
        
        /**
		 * @brief Pointer to orthogonal function to use as shape function
         * @note UseOrthoShape = 1 means it will be used Legendre polynomial as shape function. 
         * @since Mar 31, 2004
         */
        static void (*fOrthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix<REAL> & phi, TPZFMatrix<REAL> & dphi, int n);
        
    public:
        
        enum MShapeType {ETensorial, EOrdemTotal, ETensorialFull, EOrdemTotalFull};
        
        TPZShapeDisc();
        
        TPZShapeDisc(const TPZShapeDisc &copy);
        
        ~TPZShapeDisc();
        
        static void Shape(int dimension, REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type);
		
		/**
		 * @brief Computes the shape function set at the point qsi for dual space.
		 * @param dimension dimension of element
		 * @param degree degree of polynomial
		 * @param qsi point in master element coordinates
		 * @param phi vector of values of shapefunctions, dimension (numshape,1)
		 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
		 * @param type type of polinomial space
		 */
		
		static void Shape(int &dimension,int &degree,TPZVec<REAL> &qsi, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,MShapeType type);	
		
        /// initialize the TPZShapeData datastructure
        static void Initialize(int order, int dimension, MShapeType shapetype, TPZShapeData &data);
        

        
        /** @brief Number of shapefunctions dependent on the dimension and order of interpolation */
        static int NShapeF(int degree, int dimension, MShapeType type);
        
        /** @brief Discontinous polynomials of the line element */
        static void Shape0D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
        /** @brief Discontinous polynomials of the line element */
        static void Shape1D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
        /** @brief Discontinous bases of the two-dimensional elements */
        static void Shape2D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type);
        
        /** @brief Discontinous bases of the two-dimensional elements with many derivatives! */
        static void Shape2DFull(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type);
        
        /** @brief Discontinous bases of the three-dimensional elements */
        static void Shape3D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, MShapeType type);
		
    };
    
};

#endif
