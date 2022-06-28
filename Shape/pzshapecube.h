/**
 * @file
 * @brief Contains TPZShapeCube class which implements the shape functions of a hexaedral element.
 */

#ifndef SHAPECUBEHPP
#define SHAPECUBEHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "tpzcube.h"
#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshtmat.h"

#include "fadType.h"

/** @brief Groups all classes dedicated to the computation of shape functions */
namespace pzshape {
	
	/**
	 * @brief Implements the shape functions of a hexahedral (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/** 
	 * The range of the master element is [-1 ,1]
	 */
	class TPZShapeCube : public pztopology::TPZCube {
		
	public:
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a hexahedral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (19 connects in this case)
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
        static void ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders);//, TPZVec<int64_t> &sides;
        
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the internal shapefunctions of a side
         * @param sides is a vector with copy of side as much as needed, it depends on the order
         */
        static void SideShapeOrder(const int side,  const TPZVec<int64_t> &id, const int order, TPZGenMatrix<int> &shapeorders);
        
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a hexahedral element
		 * @param point (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions and derivatives
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
		 */
		static void ShapeCube(TPZVec<REAL> &point, TPZVec<int64_t> &id, TPZVec<int> &order, TPZVec<FADREAL> &phi);

		/**
		 * @brief Computes the corner shape functions for a hexahedral element
		 * @param pt (input) point where the shape function is computed already setup with derivatives
		 * @param phi (output) value of the (8) shape functions and derivatives
		 */
		static void ShapeCornerCube(TPZVec<FADREAL> &pt, TPZVec<FADREAL> &phi);
		
		/**
		 * @brief Compute the internal functions of the hexahedral shape function at a point
		 * @param x coordinate of the point (with derivatives already setup)
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values (and derivatives)
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner shape functions\n
		 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation \n
		 * determined by the transformation index
		 */
		static void Shape3dCubeInternal(TPZVec<FADREAL> &x, int order,TPZVec<FADREAL> &phi);//,int quad_transformation_index

		/**
		 * @brief Computes the corner shape functions for a hexahedral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (8) shape functions
		 * @param dphi (output) value of the derivatives of the (8) shape functions holding the derivatives in a column
		 */
//        static void ShapeCorner(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
//        {
//
//            REAL x[2],dx[2],y[2],dy[2],z[2],dz[2];
//            x[0]  = (1.-pt[0])/2.;
//            x[1]  = (1.+pt[0])/2.;
//            dx[0] = -0.5;
//            dx[1] =  0.5;
//            y[0]  = (1.-pt[1])/2.;
//            y[1]  = (1.+pt[1])/2.;
//            dy[0] = -0.5;
//            dy[1] =  0.5;
//            z[0]  = (1.-pt[2])/2.;
//            z[1]  = (1.+pt[2])/2.;
//            dz[0] = -0.5;
//            dz[1] =  0.5;
//
//            phi(0,0)  = x[0]*y[0]*z[0];
//            phi(1,0)  = x[1]*y[0]*z[0];
//            phi(2,0)  = x[1]*y[1]*z[0];
//            phi(3,0)  = x[0]*y[1]*z[0];
//            phi(4,0)  = x[0]*y[0]*z[1];
//            phi(5,0)  = x[1]*y[0]*z[1];
//            phi(6,0)  = x[1]*y[1]*z[1];
//            phi(7,0)  = x[0]*y[1]*z[1];
//            dphi(0,0) = dx[0]*y[0]*z[0];
//            dphi(1,0) = x[0]*dy[0]*z[0];
//            dphi(2,0) = x[0]*y[0]*dz[0];
//            dphi(0,1) = dx[1]*y[0]*z[0];
//            dphi(1,1) = x[1]*dy[0]*z[0];
//            dphi(2,1) = x[1]*y[0]*dz[0];
//            dphi(0,2) = dx[1]*y[1]*z[0];
//            dphi(1,2) = x[1]*dy[1]*z[0];
//            dphi(2,2) = x[1]*y[1]*dz[0];
//            dphi(0,3) = dx[0]*y[1]*z[0];
//            dphi(1,3) = x[0]*dy[1]*z[0];
//            dphi(2,3) = x[0]*y[1]*dz[0];
//            dphi(0,4) = dx[0]*y[0]*z[1];
//            dphi(1,4) = x[0]*dy[0]*z[1];
//            dphi(2,4) = x[0]*y[0]*dz[1];
//            dphi(0,5) = dx[1]*y[0]*z[1];
//            dphi(1,5) = x[1]*dy[0]*z[1];
//            dphi(2,5) = x[1]*y[0]*dz[1];
//            dphi(0,6) = dx[1]*y[1]*z[1];
//            dphi(1,6) = x[1]*dy[1]*z[1];
//            dphi(2,6) = x[1]*y[1]*dz[1];
//            dphi(0,7) = dx[0]*y[1]*z[1];
//            dphi(1,7) = x[0]*dy[1]*z[1];
//            dphi(2,7) = x[0]*y[1]*dz[1];
//        }
        template<class T>
         static void ShapeCorner(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
         {
             
             T x[2],dx[2],y[2],dy[2],z[2],dz[2];
             x[0]  = (1.-pt[0])/2.;
             x[1]  = (1.+pt[0])/2.;
             T half;
             if constexpr (std::is_same_v<FADREAL, T>)
             {
                 const int dim = pt[0].size();
                 half = FADREAL(dim,0.5);
             }
             else
             {
                 half = 0.5;
             }
             dx[0] = -half;
             dx[1] = half;
             y[0]  = (1.-pt[1])/2.;
             y[1]  = (1.+pt[1])/2.;
             dy[0] = -half;
             dy[1] = half;
             z[0]  = (1.-pt[2])/2.;
             z[1]  = (1.+pt[2])/2.;
             dz[0] = -half;
             dz[1] = half;
             
             phi(0,0)  = x[0]*y[0]*z[0];
             phi(1,0)  = x[1]*y[0]*z[0];
             phi(2,0)  = x[1]*y[1]*z[0];
             phi(3,0)  = x[0]*y[1]*z[0];
             phi(4,0)  = x[0]*y[0]*z[1];
             phi(5,0)  = x[1]*y[0]*z[1];
             phi(6,0)  = x[1]*y[1]*z[1];
             phi(7,0)  = x[0]*y[1]*z[1];
             dphi(0,0) = dx[0]*y[0]*z[0];
             dphi(1,0) = x[0]*dy[0]*z[0];
             dphi(2,0) = x[0]*y[0]*dz[0];
             dphi(0,1) = dx[1]*y[0]*z[0];
             dphi(1,1) = x[1]*dy[0]*z[0];
             dphi(2,1) = x[1]*y[0]*dz[0];
             dphi(0,2) = dx[1]*y[1]*z[0];
             dphi(1,2) = x[1]*dy[1]*z[0];
             dphi(2,2) = x[1]*y[1]*dz[0];
             dphi(0,3) = dx[0]*y[1]*z[0];
             dphi(1,3) = x[0]*dy[1]*z[0];
             dphi(2,3) = x[0]*y[1]*dz[0];
             dphi(0,4) = dx[0]*y[0]*z[1];
             dphi(1,4) = x[0]*dy[0]*z[1];
             dphi(2,4) = x[0]*y[0]*dz[1];
             dphi(0,5) = dx[1]*y[0]*z[1];
             dphi(1,5) = x[1]*dy[0]*z[1];
             dphi(2,5) = x[1]*y[0]*dz[1];
             dphi(0,6) = dx[1]*y[1]*z[1];
             dphi(1,6) = x[1]*dy[1]*z[1];
             dphi(2,6) = x[1]*y[1]*dz[1];
             dphi(0,7) = dx[0]*y[1]*z[1];
             dphi(1,7) = x[0]*dy[1]*z[1];
             dphi(2,7) = x[0]*y[1]*dz[1];
         }

        /**
         * @brief Computes the generating shape functions for a hexahedral element
         * @param pt (input) point where the shape function is computed
         * @param phi (input) value of the (8) corner shape functions
         * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        template<class T>
        static void ShapeGenerating(const TPZVec<T> &pt, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            int is;
            // contribute the ribs
            for(is=8; is<20; is++)
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
            for(is=20; is<26; is++)
            {
                int is1,is2;
                is1 = ContainedSideLocId(is,0);
                is2 = ContainedSideLocId(is,2);
                phi(is,0) = phi(is1,0)*phi(is2,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
            }
            // contribution of the volume
            for(is=26; is<27; is++)
            {
                int is1,is2;
                is1 = 0;
                is2 = 6;
                phi(is,0) = phi(is1,0)*phi(is2,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
                dphi(2,is) = dphi(2,is1)*phi(is2,0)+phi(is1,0)*dphi(2,is2);
            }

            // Make the generating shape functions linear and unitary
            for(is=8; is<27; is++)
            {
                TPZStack<int> highsides;
                HigherDimensionSides(is,highsides);
                int h, nh = highsides.NElements();
                for(h=0; h<nh; h++)
                {
                    int hs = highsides[h];
                    phi(is,0) += phi(hs,0);
                    dphi(0,is) += dphi(0,hs);
                    dphi(1,is) += dphi(1,hs);
                    dphi(2,is) += dphi(2,hs);
                }
                int dim = SideDimension(is);
                int mult = (dim == 1) ? 4 : (dim == 2) ? 16 : (dim ==3) ? 64 : 0;
                phi(is,0) *= mult;
                dphi(0,is) *= mult;
                dphi(1,is) *= mult;
                dphi(2,is) *= mult;
            }

        }
        
        /**
         * @brief Computes the generating shape functions for a hexahedral element
         * @param pt (input) point where the shape function is computed
         * @param phi (input) value of the (8) corner shape functions
         * @param dphi (input) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
        
	private:
		
	

		/**
		 * @brief Compute the internal functions of the hexahedral shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation \n
		 * determined by the transformation index
		 */
//		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
//								  TPZFMatrix<REAL> &dphi)
//        {
//            if((order-1) < 1) return;
//            int ord = order - 1;//fSideOrder[18]-1;
//            int nshape = ord*ord*ord;
//            phi.Resize(nshape,1);
//            dphi.Resize(3,nshape);
//            REAL store1[20],store2[20],store3[20],store4[20],store5[20],store6[20];
//            TPZFNMatrix<20, REAL> phi0(ord,1),phi1(ord,1),phi2(ord,1),
//            dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
//            TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);
//            TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);
//            TPZShapeLinear::fOrthogonal(x[2],ord,phi2,dphi2);
//            for (int i=0;i<ord;i++) {
//                for (int j=0;j<ord;j++) {
//                    for (int k=0;k<ord;k++) {
//                        int index = ord*(ord*i+j)+k;
//                        phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
//                        dphi(0,index) = dphi0(0,i)* phi1(j,0)* phi2(k,0);
//                        dphi(1,index) =  phi0(i,0)*dphi1(0,j)* phi2(k,0);
//                        dphi(2,index) =  phi0(i,0)* phi1(j,0)*dphi2(0,k);
//                    }
//                }
//            }
//        }
        template<class T>
        static void ShapeInternal(TPZVec<T> &x, int order,TPZFMatrix<T> &phi,
                                  TPZFMatrix<T> &dphi)
        {
            if((order-1) < 1) return;
            int ord = order - 1;//fSideOrder[18]-1;
            int nshape = ord*ord*ord;
            phi.Resize(nshape,1);
            dphi.Resize(3,nshape);
            TPZFNMatrix<20, T> phi0(ord,1),phi1(ord,1),phi2(ord,1),
            dphi0(1,ord),dphi1(1,ord),dphi2(1,ord);
            if constexpr (std::is_same_v<FADREAL, T>)
            {
                TPZShapeLinear::FADfOrthogonal(x[0],ord,phi0,dphi0);
                TPZShapeLinear::FADfOrthogonal(x[1],ord,phi1,dphi1);
                TPZShapeLinear::FADfOrthogonal(x[2],ord,phi2,dphi2);
            }
            else{
                TPZShapeLinear::fOrthogonal(x[0],ord,phi0,dphi0);
                TPZShapeLinear::fOrthogonal(x[1],ord,phi1,dphi1);
                TPZShapeLinear::fOrthogonal(x[2],ord,phi2,dphi2);
            }
            for (int i=0;i<ord;i++) {
                for (int j=0;j<ord;j++) {
                    for (int k=0;k<ord;k++) {
                        int index = ord*(ord*i+j)+k;
                        phi(index,0) =  phi0(i,0)* phi1(j,0)* phi2(k,0);
                        dphi(0,index) = dphi0(0,i)* phi1(j,0)* phi2(k,0);
                        dphi(1,index) =  phi0(i,0)*dphi1(0,j)* phi2(k,0);
                        dphi(2,index) =  phi0(i,0)* phi1(j,0)*dphi2(0,k);
                    }
                }
            }
        }
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeToRib(int face, TPZVec<REAL> &in, REAL &outval);
		
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element already setup with derivatives
		 * @param outval coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeToRib(int face, TPZVec<FADREAL> &in, FADREAL &outval);

		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element (with derivatives)
		 * @param outval coordinates of the point on the face (with derivatives)
		 */
		static void ProjectPoint3dCubeToFace(int face, TPZVec<FADREAL> &in, TPZVec<FADREAL> &outval);

		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param face rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint3dCubeSide(int face, TPZVec<REAL> &in, REAL &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinates of the point on the face
		 */
		static void ProjectPoint3dCubeFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Projects a point from the interior of the element to a face
		 * @param face face index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param outval coordinates of the point on the face
		 */
		static void ProjectPoint3dCubeToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
		
		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
		 * of the function with respect to the element. \n The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromRibToCube(int rib,int num,TPZFMatrix<REAL> &dphi);

		/**
		 * @brief Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
		 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphe shapefunctions (modified in place)
		 */
		static void TransformDerivativeFromFaceToCube(int rib,int num,TPZFMatrix<REAL> &dphi);

		/** @brief Data structure which defines the hexahedral transformations */
		static REAL gFaceTrans3dCube2d[6][2][3];
		/** @brief Data structure which defines the hexahedral transformations */
		static REAL gRibTrans3dCube1d[12][3];
		
	public:
		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param order  interpolation order associated with the side
		 * @return number of shape functions
		 */
		static int NConnectShapeF(int side, int order);
		
		/**
		 * @brief Total number of shapefunctions, considering the order of interpolation of the element
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NShapeF(const TPZVec<int> &order);
        
        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order, TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            if (side < 8 || side > 26) {
                DebugStop();
            }
            
            switch (side) {
                    
                case 8:
                case 9:
                case 10:
                case 11:
                case 12:
                case 13:
                case 14:
                case 15:
                case 16:
                case 17:
                case 18:
                case 19:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 20:
                case 21:
                case 22:
                case 23:
                case 24:
                case 25:
                {
                    pzshape::TPZShapeQuad::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 26:
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
