/**
 * @file
 * @brief Contains TPZShapeTriang class which implements the shape functions of a triangular element.
 */

#ifndef SHAPETRIANGHPP
#define SHAPETRIANGHPP

#include "tpztriangle.h"
#include "pzshtmat.h"
#include "pzshapelinear.h"
#include "fadType.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {
	
	/** 
	 * @brief Implements the shape functions of a triangular (2D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/**
	 * The triangular shape functions are also used in 3D elements \n
	 * The range of the master element is 0,1
	 */
	class TPZShapeTriang : public pztopology::TPZTriangle {
		
	public:

		/**
		 * @brief Computes the values of the shape functions and their derivatives for a triangular element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (4 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions		 
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear element for its implementation
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
		 * @brief Computes the corner shape functions for a triangular element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (3) shape functions
		 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
		 */
		static void ShapeCorner(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            phi(0,0) =  1.-pt[0]-pt[1];
            phi(1,0) =  pt[0];
            phi(2,0) =  pt[1];
            dphi(0,0) = -1.;
            dphi(1,0) = -1.;
            dphi(0,1) =  1.;
            dphi(1,1) =  0.;
            dphi(0,2) =  0.;
            dphi(1,2) =  1.;
        }
        
        static void ShapeCorner(const TPZVec<FADREAL> &pt, TPZFMatrix<FADREAL> &phi, TPZFMatrix<FADREAL> &dphi)
        {
            const int dim = pt[0].size();
            phi(0,0) =  1.-pt[0]-pt[1];
            phi(1,0) =  pt[0];
            phi(2,0) =  pt[1];
            dphi(0,0) = FADREAL(dim,-1.);
            dphi(1,0) = FADREAL(dim,-1.);
            dphi(0,1) = FADREAL(dim, 1.);
            dphi(1,1) = FADREAL(dim, 0.);
            dphi(0,2) = FADREAL(dim, 0.);
            dphi(1,2) = FADREAL(dim, 1.);
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
            for(is=3; is<6; is++)
            {
                int is1 = is%3;
                int is2 = (is+1)%3;
                phi(is,0) = phi(is1,0)*phi(is2,0);
                dphi(0,is) = dphi(0,is1)*phi(is2,0)+phi(is1,0)*dphi(0,is2);
                dphi(1,is) = dphi(1,is1)*phi(is2,0)+phi(is1,0)*dphi(1,is2);
            }
            int is1 = 0;
            int is2 = 1;
            int is3 = 2;
            phi(is,0) = phi(is1,0)*phi(is2,0)*phi(is3,0);
            dphi(0,is) = dphi(0,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(0,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(0,is3);
            dphi(1,is) = dphi(1,is1)*phi(is2,0)*phi(is3,0)+phi(is1,0)*dphi(1,is2)*phi(is3,0)+phi(is1,0)*phi(is2,0)*dphi(1,is3);

            // Make the generating shape functions linear and unitary
            REAL mult[] = {1.,1.,1.,4.,4.,4.,27.};
            for(is=3;is<NSides; is++)
            {
                phi(is,0) *= mult[is];
                dphi(0,is) *= mult[is];
                dphi(1,is) *= mult[is];
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
		 * @brief Compute the internal functions of the triangle shape function at a point
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 * @param triangle_transformation_index determines the transformation applied to the internal shape
		 * functions. \n This parameter is computed by the GetTransformId2dT method
		 * @see GetTransformId2dQ
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions. \n
		 * Shape2dTriangleInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,int triangle_transformation_index);
        
        
        static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {
            if((order - 2 ) <= 0) return;
            int numshape = ((order-2)*(order-1))/2;
            
            TPZManVector<REAL,2> out(2,0.0);
            out[0] = 2.*x[0]-1.;
            out[1] = 2.*x[1]-1.;
            
            if (phi.Rows() < numshape || dphi.Cols() < numshape) {
                PZError << "\nTPZCompEl::Shape2dTriangleInternal phi or dphi resized\n";
                phi.Resize(numshape,1);
                dphi.Resize(dphi.Rows(),numshape);
            }
            
            TPZFNMatrix<50,REAL> phi0(numshape,1),phi1(numshape,1);
            TPZFNMatrix<100,REAL> dphi0(1,numshape),dphi1(1,numshape);
            
            TPZShapeLinear::fOrthogonal(out[0],numshape,phi0,dphi0);
            TPZShapeLinear::fOrthogonal(out[1],numshape,phi1,dphi1);
            int index = 0;
            int i;
           
            for (int iplusj=0;iplusj<(order - 2);iplusj++) {
                for (int j=0;j<=iplusj;j++) {
                    i = iplusj-j;
                    phi(index,0) = phi0(i,0)*phi1(j,0);
                    dphi(0,index) = 2.0*dphi0(0,i)*phi1(j,0);
                    dphi(1,index) = 2.0*phi0(i,0)*dphi1(0,j);
                    index++;
                }
            }
        }
        static void ShapeInternal(TPZVec<FADREAL> &x, int order,TPZFMatrix<FADREAL> &phi,TPZFMatrix<FADREAL> &dphi)
        {
            if((order - 2 ) <= 0) return;
            int numshape = ((order-2)*(order-1))/2;
            
            TPZManVector<FADREAL,2> out(2,0.0);
            out[0] = 2.*x[0]-1.;
            out[1] = 2.*x[1]-1.;
            
            if (phi.Rows() < numshape || dphi.Cols() < numshape) {
                PZError << "\nTPZCompEl::Shape2dTriangleInternal phi or dphi resized\n";
                phi.Resize(numshape,1);
                dphi.Resize(dphi.Rows(),numshape);
            }
            
            TPZFNMatrix<50,FADREAL> phi0(numshape,1),phi1(numshape,1);
            TPZFNMatrix<100,FADREAL> dphi0(1,numshape),dphi1(1,numshape);
            
            TPZShapeLinear::FADfOrthogonal(out[0],numshape,phi0,dphi0);
            TPZShapeLinear::FADfOrthogonal(out[1],numshape,phi1,dphi1);
            int index = 0;
            int i;
           
            for (int iplusj=0;iplusj<(order - 2);iplusj++) {
                for (int j=0;j<=iplusj;j++) {
                    i = iplusj-j;
                    phi(index,0) = phi0(i,0)*phi1(j,0);
                    dphi(0,index) = 2.0*dphi0(0,i)*phi1(j,0);
                    dphi(1,index) = 2.0*phi0(i,0)*dphi1(0,j);
                    index++;
                }
            }
        }

        
        
		/**
		 * @brief Projects a point from the interior of the element to a rib
		 * @param rib rib index to which the point should be projected
		 * @param in coordinate of the point at the interior of the element
		 * @param out coordinate of the point on the rib
		 */
		static void ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out);
        
      

		/**
		 * @brief Transforms the derivative of a shapefunction computed on the rib into the two dimensional derivative \n
		 * of the function with respect to the element. The parameter dphi should be dimensioned (2,num), at least
		 * @param rib rib index along which the shapefunction is defined
		 * @param num number of shapefunction derivatives which need to be transformed
		 * @param dphi values of the derivatives of the shapefunctions
		 */
		static void TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Method which identifies the triangular transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point
		 */
		static int GetTransformId2dT(TPZVec<int64_t> &id);
		
		/**
		 * @brief Transform the coordinates of the point in the space of the triangle
		 * master element based on the transformation id
		 * @param transid identifier of the transformation of the element as obtained by the GetTransformId2dT method
		 * @param in coordinates of the variational parameter
		 * @param out coordinates of the transformed parameter
		 */
		static void TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out);
		
		/**
		 * @brief Transform the derivatives of num shapefunctions in place for a triangle
		 * @param transid identifier of the transformation of the triangle element as obtained by the GetTransformId2dT method
		 * @param num number of shapefunctions needed to transform
		 * @param in matrix containing the values of the derivatives of the shapefunctions as a row vector \n
		 * the values of the derivatives contained in this matrix are modified upon return
		 */
		static void TransformDerivative2dT(int transid, int num, TPZFMatrix<REAL> &in);
		
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gTrans2dT[6][2][2];
		
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gVet2dT[6][2];
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gRibTrans2dT1d[3][2];
		/** @brief Data structure which defines the triangle transformations*/
		static REAL gVet1dT[3];
		
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
        static TPZTransform<REAL>  ParametricTransform(int trans_id);
        
        // compute the internal shape function for the side. The point x is already projected in the side coordinates
        template<class T>
        static void ShapeInternal(int side, TPZVec<T> &x, int order,
                                           TPZFMatrix<T> &phi, TPZFMatrix<T> &dphi)
        {
            if (side < 3 || side > 6) {
                DebugStop();
            }
            
            switch (side) {
                    
                case 3:
                case 4:
                case 5:
                {
                    pzshape::TPZShapeLinear::ShapeInternal(x, order, phi, dphi);
                }
                    break;
                case 6:
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
