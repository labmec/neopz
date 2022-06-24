/**
 * @file
 * @brief Contains TPZShapePrism class which implements the shape functions of a prism element.
 */

#ifndef SHAPEPRISMHPP
#define SHAPEPRISMHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzprism.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshtmat.h"

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
		static void Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,
							  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
        
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
         * @brief Computes the generating shape functions for the prism element
         * @param pt (input) point where the shape function is computed
         * @param phi (input) value of the (4) shape functions
         * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
        
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
            int ord1 = order-2;
            int ord2 = order-1;
            
            TPZFNMatrix<20,REAL> phi0(ord1,1),phi1(ord1,1),phi2(ord2,1),
            dphi0(1,ord1),dphi1(1,ord1),dphi2(1,ord2);
            TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
            TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
            TPZShapeLinear::fOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
            int index = 0;//x � ponto de integra��o dentro da pir�mide
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
            int ord1 = order-2;
            int ord2 = order-1;
            
            TPZFNMatrix<20,FADREAL> phi0(ord1,1),phi1(ord1,1),phi2(ord2,1),
            dphi0(1,ord1),dphi1(1,ord1),dphi2(1,ord2);
            TPZShapeLinear::FADfOrthogonal(2.*x[0]-1.,ord1,phi0,dphi0);//f e df       0<=x0<=1 -> -1<=2*x0-1<=1
            TPZShapeLinear::FADfOrthogonal(2.*x[1]-1.,ord1,phi1,dphi1);//g e dg             0<=x1<=1 -> -1<=2*x1-1<=1
            TPZShapeLinear::FADfOrthogonal(x[2],ord2,phi2,dphi2);//h e dh      -1<=x3<=1
            int index = 0;//x � ponto de integra��o dentro da pir�mide
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
		static int NShapeF(const TPZVec<int> &order);
        static void CornerShape(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
         

        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            if (side < 6 || side > 20) {
                DebugStop();
            }
            
            switch (side) {
                    
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                    
                case 15:
                case 19:
                    
                {
                
                    pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                    
                case 16:
                case 17:
                case 18:
                {
                    pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
                }
                    break;

                case 20:
                {
                    ShapeInternal(x, order, phi, dphi);
                }
                    break;
                default:
                    std::cout << "Wrong side parameter side " << side << std::endl;
                    DebugStop();
                    break;
            }
         
            
        }
	};
	
};

#endif
