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
#include "fadType.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
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

        /// Compute the corner shape functions
        static void ShapeCorner(const TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {
            if(fabs(pt[0])<1.e-10 && fabs(pt[1])<1.e-10 && pt[2]==1.) {
                //para testes com transformaçoes geometricas-->>Que é o que faz o RefPattern!!
                //(0,0,1) nunca é um ponto de integração
                phi(0,0)  = 0.;
                phi(1,0)  = 0.;
                phi(2,0)  = 0.;
                phi(3,0)  = 0.;
                phi(4,0)  = pt[2];
                dphi(0,0)  = -0.25;
                dphi(1,0)  = -0.25;
                dphi(2,0)  = -0.25;
                dphi(0,1)  = 0.25;
                dphi(1,1)  = -0.25;
                dphi(2,1)  = -0.25;
                dphi(0,2)  = 0.25;
                dphi(1,2)  = 0.25;
                dphi(2,2)  = -0.25;
                dphi(0,3)  = -0.25;
                dphi(1,3)  = 0.25;
                dphi(2,3)  = -0.25;
                dphi(0,4)  = 0;
                dphi(1,4)  = 0;
                dphi(2,4)  = 1.;
                return;
            }
            
            REAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
            REAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
            REAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
            REAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
            REAL lmez = (1.-pt[2]);
            phi(0,0)  = T0xz*T0yz*lmez;
            phi(1,0)  = T1xz*T0yz*lmez;
            phi(2,0)  = T1xz*T1yz*lmez;
            phi(3,0)  = T0xz*T1yz*lmez;
            phi(4,0)  = pt[2];
            REAL lmexmez = 1.-pt[0]-pt[2];
            REAL lmeymez = 1.-pt[1]-pt[2];
            REAL lmaxmez = 1.+pt[0]-pt[2];
            REAL lmaymez = 1.+pt[1]-pt[2];
            dphi(0,0) = -.25*lmeymez / lmez;
            dphi(1,0) = -.25*lmexmez / lmez;
            dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
            
            dphi(0,1) =  .25*lmeymez / lmez;
            dphi(1,1) = -.25*lmaxmez / lmez;
            dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
            
            dphi(0,2) =  .25*lmaymez / lmez;
            dphi(1,2) =  .25*lmaxmez / lmez;
            dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
            
            dphi(0,3) = -.25*lmaymez / lmez;
            dphi(1,3) =  .25*lmexmez / lmez;
            dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
            
            dphi(0,4) =  0.0;
            dphi(1,4) =  0.0;
            dphi(2,4) =  1.0;
        }
        static void ShapeCorner(const TPZVec<FADREAL> &pt,TPZFMatrix<FADREAL> &phi,TPZFMatrix<FADREAL> &dphi)
        {
            const int dim = pt[0].size();
            FADREAL one(dim,1.);
            if(fabs(pt[0].val())<1.e-10 && fabs(pt[1].val())<1.e-10 && pt[2].val()==1.) {
                //para testes com transformaçoes geometricas-->>Que é o que faz o RefPattern!!
                //(0,0,1) nunca é um ponto de integração
                phi(0,0)  = one*0.;
                phi(1,0)  = one*0.;
                phi(2,0)  = one*0.;
                phi(3,0)  = one*0.;
                phi(4,0)  = one*pt[2];
                dphi(0,0)  = one*-0.25;
                dphi(1,0)  = one*-0.25;
                dphi(2,0)  = one*-0.25;
                dphi(0,1)  = one*0.25;
                dphi(1,1)  = one*-0.25;
                dphi(2,1)  = one*-0.25;
                dphi(0,2)  = one*0.25;
                dphi(1,2)  = one*0.25;
                dphi(2,2)  = one*-0.25;
                dphi(0,3)  = one*-0.25;
                dphi(1,3)  = one*0.25;
                dphi(2,3)  = one*-0.25;
                dphi(0,4)  = one*0;
                dphi(1,4)  = one*0;
                dphi(2,4)  = one*1.;
                return;
            }
            
            FADREAL T0xz = .5*(1.-pt[2]-pt[0]) / (1.-pt[2]);
            FADREAL T0yz = .5*(1.-pt[2]-pt[1]) / (1.-pt[2]);
            FADREAL T1xz = .5*(1.-pt[2]+pt[0]) / (1.-pt[2]);
            FADREAL T1yz = .5*(1.-pt[2]+pt[1]) / (1.-pt[2]);
            FADREAL lmez = (1.-pt[2]);
            phi(0,0)  = T0xz*T0yz*lmez;
            phi(1,0)  = T1xz*T0yz*lmez;
            phi(2,0)  = T1xz*T1yz*lmez;
            phi(3,0)  = T0xz*T1yz*lmez;
            phi(4,0)  = pt[2];
            FADREAL lmexmez = 1.-pt[0]-pt[2];
            FADREAL lmeymez = 1.-pt[1]-pt[2];
            FADREAL lmaxmez = 1.+pt[0]-pt[2];
            FADREAL lmaymez = 1.+pt[1]-pt[2];
            dphi(0,0) = -.25*lmeymez / lmez;
            dphi(1,0) = -.25*lmexmez / lmez;
            dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
            
            dphi(0,1) =  .25*lmeymez / lmez;
            dphi(1,1) = -.25*lmaxmez / lmez;
            dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
            
            dphi(0,2) =  .25*lmaymez / lmez;
            dphi(1,2) =  .25*lmaxmez / lmez;
            dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
            
            dphi(0,3) = -.25*lmaymez / lmez;
            dphi(1,3) =  .25*lmexmez / lmez;
            dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
            
            dphi(0,4) =  0.0;
            dphi(1,4) =  0.0;
            dphi(2,4) =  1.0;
        }
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
								  TPZFMatrix<REAL> &dphi)
        {
            // calculate the values of the function and derivatives of the product of the orthogonal functions
            if(order < 3) return;
            int ord = order-2;
            TPZFNMatrix<20, REAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),
            dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
            TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);//f e df            -1<=x0<=1
            TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);//g e dg            -1<=x1<=1
            TPZShapeLinear::fOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);//h e dh       0<=x3<=1 -> -1<=2*x2-1<=1
            int index = 0;   // integration point internal to pyramid
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        if( i+j+k < ord ) {
                            //int index = ord*(ord*i+j)+k; //i,j,k is the order of the orthogonal functions
                            phi(index,0) =    phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) =   dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =    phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) = 2.*phi0(i,0)* phi1(j,0)*dphi2(0,k);
                            index++;
                        }
                    }
                }
            }
        }
        static void ShapeInternal(TPZVec<FADREAL> &x, int order,TPZFMatrix<FADREAL> &phi,
                                  TPZFMatrix<FADREAL> &dphi)
        {
            // calculate the values of the function and derivatives of the product of the orthogonal functions
            if(order < 3) return;
            int ord = order-2;
            TPZFNMatrix<20, FADREAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),
            dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
            TPZShapeLinear::FADfOrthogonal(x[0],ord,phi0,dphi0);//f e df            -1<=x0<=1
            TPZShapeLinear::FADfOrthogonal(x[1],ord,phi1,dphi1);//g e dg            -1<=x1<=1
            TPZShapeLinear::FADfOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);//h e dh       0<=x3<=1 -> -1<=2*x2-1<=1
            int index = 0;   // integration point internal to pyramid
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        if( i+j+k < ord ) {
                            //int index = ord*(ord*i+j)+k; //i,j,k is the order of the orthogonal functions
                            phi(index,0) =    phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) =   dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =    phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) = 2.*phi0(i,0)* phi1(j,0)*dphi2(0,k);
                            index++;
                        }
                    }
                }
            }
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
            // contribute the ribs
            for(is=NCornerNodes; is<NCornerNodes+8; is++)
            {
                int is1,is2;
                is1 = ContainedSideLocId(is,0);
                is2 = ContainedSideLocId(is,1);
                phi(is,0) = phi(is1,0)*phi(is2,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
            }
            // contribution of the faces
            // quadrilateral face
            is = 13;
            {
                int is1,is2;
                is1 = ShapeFaceId[0][0];// SideConnectLocId(is,0);
                is2 = ShapeFaceId[0][2];// SideConnectLocId(is,2);
                phi(is,0) = phi(is1,0)*phi(is2,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
            }
            is++;
            // triangular faces
            for(;is<18; is++)
            {
                int is1,is2,is3;
                is1 = ContainedSideLocId(is,0);
                is2 = ContainedSideLocId(is,1);
                is3 = ContainedSideLocId(is,2);
                phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
                
            }
            {
                int is1 = 0;
                int is2 = 2;
                int is3 = 4;
                phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
            }
            
            // Make the generating shape functions linear and unitary
            // contribute the ribs
            for(is=NCornerNodes; is<NCornerNodes+4; is++)
            {
                int isface = 13;
                phi(is,0) += phi(isface,0);
                dphi(0,is) += dphi(0,isface);
                dphi(1,is) += dphi(1,isface);
                dphi(2,is) += dphi(2,isface);
            }
            // scaling the shapefunctions
            {
                REAL sidescale[] = {1.,1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,4.,4.,16.,27.,27.,27.,27.,64.};
                for(is=5; is<NSides; is++)
                {
                    phi(is,0) *= sidescale[is];
                    dphi(0,is) *= sidescale[is];
                    dphi(1,is) *= sidescale[is];
                    dphi(2,is) *= sidescale[is];
                }
            }
            
            
        }
		
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
        
        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order,  TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            if (side < 5 || side > 18) {
                DebugStop();
            }
            
            switch (side) {
                    
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 13:
                {
                    pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 14:
                case 15:
                case 16:
                case 17:
                {
                    pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                
                case 18:
                {
                    ShapeInternal(x, order, phi, dphi);
                }
                    break;
                default:
                    std::cout << __PRETTY_FUNCTION__ << " Wrong side parameter side " << side << std::endl;
                    DebugStop();
                    break;
            }
           
            
        }
	};
	
};

#endif
