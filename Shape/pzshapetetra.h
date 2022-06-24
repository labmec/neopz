/**
 * @file
 * @brief Contains TPZShapeTetra class which implements the shape functions of a tetrahedral element.
 */

#ifndef SHAPETETRAHPP
#define SHAPETETRAHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpztetrahedron.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "fadType.h"
#include "pzshtmat.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a tetrahedral (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The range of the master element is 0,1
	 */
	class TPZShapeTetra : public pztopology::TPZTetrahedron{
		
	public:

		/**
		 * @brief Computes the values of the shape functions and their derivatives for a tetrahedral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (11 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and triangular element for its implementation
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
		 * @brief Computes the corner shape functions for a tetrahedral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (4) shape functions
		 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeCorner(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
            phi(1,0)  = pt[0];
            phi(2,0)  = pt[1];
            phi(3,0)  = pt[2];
            
            dphi(0,0) = -1.0;
            dphi(1,0) = -1.0;
            dphi(2,0) = -1.0;
            dphi(0,1) =  1.0;
            dphi(1,1) =  0.0;
            dphi(2,1) =  0.0;
            dphi(0,2) =  0.0;
            dphi(1,2) =  1.0;
            dphi(2,2) =  0.0;
            dphi(0,3) =  0.0;
            dphi(1,3) =  0.0;
            dphi(2,3) =  1.0;
        }
        static void ShapeCorner(const TPZVec<FADREAL> &pt, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi)
        {
            phi(0,0)  = 1-pt[0]-pt[1]-pt[2];
            phi(1,0)  = pt[0];
            phi(2,0)  = pt[1];
            phi(3,0)  = pt[2];
            const int dim = pt[0].size();
            dphi(0,0) = FADREAL(dim,-1.0);
            dphi(1,0) = FADREAL(dim,-1.0);
            dphi(2,0) = FADREAL(dim,-1.0);
            dphi(0,1) = FADREAL(dim, 1.0);
            dphi(1,1) = FADREAL(dim, 0.0);
            dphi(2,1) = FADREAL(dim, 0.0);
            dphi(0,2) = FADREAL(dim, 0.0);
            dphi(1,2) = FADREAL(dim, 1.0);
            dphi(2,2) = FADREAL(dim, 0.0);
            dphi(0,3) = FADREAL(dim, 0.0);
            dphi(1,3) = FADREAL(dim, 0.0);
            dphi(2,3) = FADREAL(dim, 1.0);
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
            // 6 ribs
            for(is=4; is<NSides; is++)
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
                        //int face = is-10;
                        int is1 = SideNodeLocId(is,0); //ShapeFaceId[face][0];
                        int is2 = SideNodeLocId(is,1); //ShapeFaceId[face][1];
                        int is3 = SideNodeLocId(is,2); //ShapeFaceId[face][2];
                        phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
                        dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
                        dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);
                        dphi(2,is) = dphi(2,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(2,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(2,is3);
                    }
                        break;
                    case 4:
                    {
                        phi(is,0) = phi(0,0)*phi(1,0)*phi(2,0)*phi(3,0);
                        for(int xj=0;xj<3;xj++) {
                            dphi(xj,is) = dphi(xj,0)* phi(1 ,0)* phi(2 ,0)* phi(3 ,0) +
                            phi(0, 0)*dphi(xj,1)* phi(2 ,0)* phi(3 ,0) +
                            phi(0, 0)* phi(1 ,0)*dphi(xj,2)* phi(3 ,0) +
                            phi(0, 0)* phi(1 ,0)* phi(2 ,0)*dphi(xj,3);
                        }
                    }
                        break;
                        
                    default:
                        DebugStop();
                }
            }

            REAL mult[] = {1.,1.,1.,1.,4.,4.,4.,4.,4.,4.,27.,27.,27.,27.,54.};
            for(is=4;is<NSides; is++)
            {
                phi(is,0) *= mult[is];
                dphi(0,is) *= mult[is];
                dphi(1,is) *= mult[is];
                dphi(2,is) *= mult[is];
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
		 * @brief Compute the internal functions of the tetrahedral shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dTetraInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index (also on the faces)
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            if(order < 4) return;
            int ord = order-3;
            
            TPZFNMatrix<100, REAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
            TPZShapeLinear::fOrthogonal(2.*x[0]-1.,ord,phi0,dphi0);
            TPZShapeLinear::fOrthogonal(2.*x[1]-1.,ord,phi1,dphi1);
            TPZShapeLinear::fOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);
            
            
            int index = 0;
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        if( i+j+k < ord ) {
                            //int index = ord*(ord*i+j)+k;
                            phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) =  2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =  2.* phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) =  2.* phi0(i,0)* phi1(j,0)*dphi2(0,k);
                            index++;
                        }
                    }
                }
            }
        }
        static void ShapeInternal(TPZVec<FADREAL> &x, int order,TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi)
        {
            if(order < 4) return;
            int ord = order-3;
            
            TPZFNMatrix<100, FADREAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
            TPZShapeLinear::FADfOrthogonal(2.*x[0]-1.,ord,phi0,dphi0);
            TPZShapeLinear::FADfOrthogonal(2.*x[1]-1.,ord,phi1,dphi1);
            TPZShapeLinear::FADfOrthogonal(2.*x[2]-1.,ord,phi2,dphi2);
            
            
            int index = 0;
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        if( i+j+k < ord ) {
                            //int index = ord*(ord*i+j)+k;
                            phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
                            dphi(0,index) =  2.*dphi0(0,i)* phi1(j,0)* phi2(k,0);
                            dphi(1,index) =  2.* phi0(i,0)*dphi1(0,j)* phi2(k,0);
                            dphi(2,index) =  2.* phi0(i,0)* phi1(j,0)*dphi2(0,k);
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
		static void ProjectPoint3dTetraToRib(int rib, TPZVec<REAL> &in, REAL &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param side side to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint3dTetrSide(int side, TPZVec<REAL> &in, REAL &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dTetraToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinates of the point on the face
		 */
		static void ProjectPoint3dTetrFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToTetra(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param face face index to which the point should be projected
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToTetra(int face,int num,TPZFMatrix<REAL> &dphi);
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gVet1dTetr[6];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceTrans3dTetr2d[4][2][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gVet2dTetr[4][2];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibTrans3dTetr1d[6][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceSum3dTetra2d[4][2];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gFaceTrans3dTetra2d[4][2][3];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibSum3dTetra1d[6];
		
		/** @brief Data structure which defines the tetrahedral transformations and topology */
		static REAL gRibTrans3dTetra1d[6][3];
		
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
        
        // temporarily maintaining the "old" signature
        static void CornerShape(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            if (side < 4 || side > 14) {
                DebugStop();
            }
            
            switch (side) {
                    
                case 4:
                case 5:
                case 6:
                case 7:
                case 8:
                case 9:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                    
                case 10:
                case 11:
                case 12:
                case 13:
                    
                {
                    pzshape::TPZShapeTriang::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                    
                case 14:
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
