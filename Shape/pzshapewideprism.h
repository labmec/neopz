/**
 * @file
 * @brief Contains TPZShapeWidePrism class which implements the shape functions of a prism element.
 */

#ifndef SHAPEWIDEPRISMHPP
#define SHAPEWIDEPRISMHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzprism.h"
#include "pzshtmat.h"
#include "pzshapelinear.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a prism (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is \f$ [-1,1] \f$ for the \f$ z \f$ direction and \f$ [0,1] \f$ for the triangular plane
	 */
	class TPZShapeWidePrism : public pztopology::TPZPrism {
		
	public:
		
        
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the side associated with each shapefunction
         */
        static void InternalShapeOrder(const TPZVec<int64_t> &id, int order, TPZGenMatrix<int> &shapeorders);

        
        /**
         * @brief Computes the corner shape functions for a prism element
         * @param pt (input) point where the shape function is computed
         * @param phi (output) value of the (6) shape functions
         * @param dphi (output) value of the derivatives of the (6) shape functions holding the derivatives in a column
         */
        template<class T>
        static void ShapeCorner(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            phi(0,0)  = .5*(1.-pt[0]-pt[1])*(1.-pt[2]);
            phi(1,0)  = .5*pt[0]*(1.-pt[2]);
            phi(2,0)  = .5*pt[1]*(1.-pt[2]);
            phi(3,0)  = .5*(1.-pt[0]-pt[1])*(1.+pt[2]);
            phi(4,0)  = .5*pt[0]*(1.+pt[2]);
            phi(5,0)  = .5*pt[1]*(1.+pt[2]);
            
            dphi(0,0) = -.5*(1.-pt[2]);
            dphi(1,0) = -.5*(1.-pt[2]);
            dphi(2,0) = -.5*(1.-pt[0]-pt[1]);
            
            dphi(0,1) =  .5*(1.-pt[2]);
            dphi(1,1) =  .0*pt[1];
            dphi(2,1) = -.5*pt[0];
            
            dphi(0,2) =  .0*pt[1];
            dphi(1,2) =  .5*(1.-pt[2]);
            dphi(2,2) = -.5*pt[1];
            
            dphi(0,3) = -.5*(1.+pt[2]);
            dphi(1,3) = -.5*(1.+pt[2]);
            dphi(2,3) =  .5*(1.-pt[0]-pt[1]);
            
            dphi(0,4) =  .5*(1.+pt[2]);
            dphi(1,4) =  .0*pt[1];
            dphi(2,4) =  .5*pt[0];
            
            dphi(0,5) =  .0*pt[1];
            dphi(1,5) =  .5*(1.+pt[2]);
            dphi(2,5) =  .5*pt[1];
        }

        /**
         * @brief Computes the generating shape functions for the prism element
         * @param pt (input) point where the shape function is computed
         * @param phi (input) value of the (4) shape functions
         * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        template<class T>
        static void ShapeGenerating(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            int is;
            // 9 ribs
            for(is=6; is<NSides; is++)
            {
                int nsnodes = NSideNodes(is);
                switch(nsnodes)
                {
                    case 2:
                    {
                        int is1 = SideNodeLocId(is,0);
                        int is2 = SideNodeLocId(is,1);
                        phi(is,0) = phi(is1,0)*phi(is2,0);
                        dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                        dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                        dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
                    }
                        break;
                    case 3:
                    {
                        int is0 = 0;
                        int is1 = 1;
                        int is2 = 2;
                        int is3 = 3;
                        int is4 = 4;
                        if(is == 19)
                        {
                            is2 = 5;
                        }
                        phi(is,0) = (phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*phi(is2,0);
                        int d;
                        for(d=0; d<3; d++)
                        {
                            dphi(d,is) =
                            (dphi(d,is0)+dphi(d,is3))*(phi(is1,0)+phi(is4,0))*phi(is2,0) +
                            (phi(is0,0)+phi(is3,0))*(dphi(d,is1)+dphi(d,is4))*phi(is2,0) +
                            (phi(is0,0)+phi(is3,0))*(phi(is1,0)+phi(is4,0))*dphi(d,is2);
                        }
                    }
                        break;
                    case 4:
                    {
                        int is1 = SideNodeLocId(is,0);
                        int is2 = SideNodeLocId(is,2);
                        phi(is,0) = phi(is1,0)*phi(is2,0);
                        dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                        dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                        dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
                    }
                        break;
                    case 6:
                    {
                        int is1 = 0;
                        int is2 = 4;
                        int is3 = 2;
                        int is4 = 5;
                        phi(is,0) = phi(is1,0)*phi(is2,0)*(phi(is3,0)+phi(is4,0));
                        dphi(0,is) = dphi(0,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(0,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(0,is3)+dphi(0,is4));
                        dphi(1,is) = dphi(1,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(1,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(1,is3)+dphi(1,is4));
                        dphi(2,is) = dphi(2,is1)*phi(is2,0)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*dphi(2,is2)*(phi(is3,0)+phi(is4,0))+phi(is1,0)*phi(is2,0)*(dphi(2,is3)+dphi(2,is4));
                    }
                        break;
                    default:
                        DebugStop();
                }
            }
            // Make the generating shape functions linear and unitary
            for(is=6; is<NSides; is++)
            {
                TPZStack<int> highsides;
                HigherDimensionSides(is,highsides);
                int h, nh = highsides.NElements();
                for(h=0; h<nh; h++)
                {
                    int hs = highsides[h];
                    if(NSideNodes(hs) != 4) continue;
                    phi(is,0) += phi(hs,0);
                    dphi(0,is) += dphi(0,hs);
                    dphi(1,is) += dphi(1,hs);
                    dphi(2,is) += dphi(2,hs);
                }
            }
            REAL mult[] = {1.,1.,1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,4.,4.,4.,27.,16.,16.,16.,27.,8.};
            for(is=6;is<NSides; is++)
            {
                phi(is,0) *= mult[is];
                dphi(0,is) *= mult[is];
                dphi(1,is) *= mult[is];
                dphi(2,is) *= mult[is];
            }
            
        }

        
        /**
         * @brief Compute the bubble functions of the prism shape function at a point
         * @param x coordinate of the point
         * @param order maximum order of shape functions to be computed
         * @param phi shapefunction values
         * @param dphi values of the derivatives of the shape functions
         */
        static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
                                  TPZFMatrix<REAL> &dphi)
        {
            if(order < 3) return;
            int ord1 = order-1;
            int ord2 = order-1;
            
            TPZFNMatrix<20,REAL> phi0(ord1,1),phi1(ord1,1),phi2(ord2,1),
            dphi0(1,ord1),dphi1(1,ord1),dphi2(1,ord2);
            TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
            TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
            TPZShapeLinear::fOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
            int index = 0;
            for (int i=0;i<ord1;i++) {
                for (int j=0;j<ord1;j++) {
                    for (int k=0;k<ord2;k++) {
                        if( i+j < ord1 && k < ord2) {
                            phi(index,0) =     phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) = 2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =  2.*phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) =     phi0(i,0)* phi1(j,0)*dphi2(0,k);
                            index++;
                        }
                    }
                }
            }
        }
        static void ShapeInternal(TPZVec<FADREAL> &x, int order,TPZFMatrix<FADREAL> &phi,
                                  TPZFMatrix<FADREAL> &dphi)
        {
            if(order < 3) return;
            int ord1 = order-1;
            int ord2 = order-1;
            
            TPZFNMatrix<20,FADREAL> phi0(ord1,1),phi1(ord1,1),phi2(ord2,1),
            dphi0(1,ord1),dphi1(1,ord1),dphi2(1,ord2);
            TPZShapeLinear::FADfOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
            TPZShapeLinear::FADfOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
            TPZShapeLinear::FADfOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
            int index = 0;
            for (int i=0;i<ord1;i++) {
                for (int j=0;j<ord1;j++) {
                    for (int k=0;k<ord2;k++) {
                        if( i+j < ord1 && k < ord2) {
                            phi(index,0) =     phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) = 2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =  2.*phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) =     phi0(i,0)* phi1(j,0)*dphi2(0,k);
                            index++;
                        }
                    }
                }
            }
        }

		
        /*
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
        
//       static void ShapeInternal(int side, TPZVec<REAL> &x, int order, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
	};
	
};

#endif
