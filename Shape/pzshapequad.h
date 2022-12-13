/**
 * @file
 * @brief Contains TPZShapeQuad class which implements the shape functions of a quadrilateral element.
 */

#ifndef SHAPEQUADHPP
#define SHAPEQUADHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzquadrilateral.h"
#include "pzshtmat.h"
#include "pzshapelinear.h"

#include "fadType.h"

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
		 * @brief Computes the corner shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (4) shape functions
		 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
        template<class T>
		static void ShapeCorner(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            REAL dx[2],dy[2];
            T x[2], y[2];
            x[0]  =  (1.-pt[0])/2.;
            x[1]  =  (1.+pt[0])/2.;
            dx[0] = -0.5;
            dx[1] =  0.5;
            y[0]  =  (1.-pt[1])/2.;
            y[1]  =  (1.+pt[1])/2.;
            dy[0] = -0.5;
            dy[1] =  0.5;
            phi(0,0)  = x[0]*y[0];
            phi(1,0)  = x[1]*y[0];
            phi(2,0)  = x[1]*y[1];
            phi(3,0)  = x[0]*y[1];
            dphi(0,0) = dx[0]*y[0];
            dphi(1,0) = x[0]*dy[0];
            dphi(0,1) = dx[1]*y[0];
            dphi(1,1) = x[1]*dy[0];
            dphi(0,2) = dx[1]*y[1];
            dphi(1,2) = x[1]*dy[1];
            dphi(0,3) = dx[0]*y[1];
            dphi(1,3) = x[0]*dy[1];
        }
		
		/**
		 * @brief Computes the generating shape functions for a quadrilateral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (input/output) value of the (4) shape functions
		 * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
        template<class T>
		static void ShapeGenerating(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            int is;
            for(is=4; is<8; is++)
            {
                phi(is,0) = phi(is%4,0)*phi((is+1)%4,0);
                dphi(0,is) = dphi(0,is%4)*phi((is+1)%4,0)+phi(is%4,0)*dphi(0,(is+1)%4);
                dphi(1,is) = dphi(1,is%4)*phi((is+1)%4,0)+phi(is%4,0)*dphi(1,(is+1)%4);
            }
            phi(8,0) = phi(0,0)*phi(2,0);
            dphi(0,8) = dphi(0,0)*phi(2,0)+phi(0,0)*dphi(0,2);
            dphi(1,8) = dphi(1,0)*phi(2,0)+phi(0,0)*dphi(1,2);

            // Make the generating shape functions linear and unitary
            for(is=4; is<8; is++)
            {
                phi(is,0) += phi(8,0);
                dphi(0,is) += dphi(0,8);
                dphi(1,is) += dphi(1,8);
                phi(is,0) *= 4.;
                dphi(0,is) *= 4.;
                dphi(1,is) *= 4.;
            }
            phi(8,0) *= 16.;
            dphi(0,8) *= 16.;
            dphi(1,8) *= 16.;
        }
		
        
        static void ShapeInternal(TPZVec<REAL> &x, int order,
                                  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {
            if((order-2) < 0) return;
            int ord1 = order - 1;
            int numshape = (order-1)*(order-1);
            if(numshape > phi.Rows() || phi.Cols() < 1) phi.Resize(numshape,1);
            if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
            TPZFNMatrix<20, REAL> phi0(ord1,1),phi1(ord1,1),dphi0(1,ord1),dphi1(1,ord1);
            TPZShapeLinear::fOrthogonal(x[0],ord1,phi0,dphi0);
            TPZShapeLinear::fOrthogonal(x[1],ord1,phi1,dphi1);
            for (int i=0;i<ord1;i++) {
                for (int j=0;j<ord1;j++) {
                    int index = i*ord1+j;
                    phi(index,0) =  phi0(i,0)* phi1(j,0);
                    dphi(0,index) = dphi0(0,i)* phi1(j,0);
                    dphi(1,index) =  phi0(i,0)*dphi1(0,j);
                }
            }
        }
        static void ShapeInternal(TPZVec<FADREAL> &x, int order,
                                  TPZFMatrix<FADREAL> &phi,TPZFMatrix<FADREAL> &dphi)
        {
            if((order-2) < 0) return;
            int ord1 = order - 1;
            int numshape = (order-1)*(order-1);
            if(numshape > phi.Rows() || phi.Cols() < 1) phi.Resize(numshape,1);
            if(dphi.Rows() < 2 || dphi.Cols() < numshape) dphi.Resize(2,numshape);
            TPZFNMatrix<20, FADREAL> phi0(ord1,1),phi1(ord1,1),dphi0(1,ord1),dphi1(1,ord1);
            TPZShapeLinear::FADfOrthogonal(x[0],ord1,phi0,dphi0);
            TPZShapeLinear::FADfOrthogonal(x[1],ord1,phi1,dphi1);
            for (int i=0;i<ord1;i++) {
                for (int j=0;j<ord1;j++) {
                    int index = i*ord1+j;
                    phi(index,0) =  phi0(i,0)* phi1(j,0);
                    dphi(0,index) = dphi0(0,i)* phi1(j,0);
                    dphi(1,index) =  phi0(i,0)*dphi1(0,j);
                }
            }
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
        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
#ifdef PZDEBUG
            if (side < 4 || side > 8) {
                DebugStop();
            }
#endif
            switch (side) {

                case 4:
                case 5:
                case 6:
                case 7:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 8:
                {
                    ShapeInternal(x, order, phi, dphi);
                }
                    break;
                default:
                    std::cout<< "Wrong side parameter" << std::endl;
                    DebugStop();
                    break;
            }
        }
        
		
	};
	
};

#endif
