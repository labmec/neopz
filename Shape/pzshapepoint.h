/**
 * @file
 * @brief Contains TPZShapePoint class which implements the shape function associated with a point.
 */

#ifndef PZSHAPEPOINT
#define PZSHAPEPOINT

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "tpzpoint.h"
#include "pzshtmat.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape{
	
	/**
	 * @brief Compute the single shape function associated with a point. \ref shape "Shape"
	 * @ingroup shape
	 */
	class TPZShapePoint  : public pztopology::TPZPoint  {
	public:
		
		/** @brief Temporary storage to accelerate the computation of shape functions. To point it isn't necessary. */
		struct TMem  
		{
		};
		typedef pztopology::TPZPoint Top;
		
		/**
		 * @brief Computes the values of the shape functions and their derivatives for a quadrilateral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
			phi(0,0) = 1.;
		}
		
        static void ShapeCorner(const TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {
            phi(0,0) = 1.;
        }
        /**
         * @brief Computes the generating shape functions for a quadrilateral element
         * @param pt (input) point where the shape function is computed
         * @param phi (input/output) value of the (4) shape functions
         * @param dphi (input/output) value of the derivatives of the (4) shape functions holding the derivatives in a column
         */
        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZVec<int> &nshape, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){}

        static void ShapeGenerating(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi){}

        static void ShapeInternal(int side, TPZVec<REAL> &x, int order, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
        {
            
        }

        
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int64_t> &id, TPZVec<int> &order,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi) {
			if(side == 0) Shape(pt,id,order,phi,dphi);
		}
		
        /**
         * @brief returns the polynomial order in the natural ksi, eta of the side associated with each shapefunction
         */
        static void ShapeOrder(const TPZVec<int64_t> &id, const TPZVec<int> &order, TPZGenMatrix<int> &shapeorders)//, TPZVec<int64_t> &sides
        {
            shapeorders(0,0) = 0;
        }
        

		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NConnectShapeF(int side, int order) { return 1; }
		
		/**
		 * @brief Total number of shapefunctions, considering the order
		 * of interpolation of the element
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NShapeF(TPZVec<int> &order) { return 1; }
		
		/**
		 * @brief Compute the permutation of the connects of the sides such that the order of the \n 
		 * shape functions becomes independent of the element orientation
		 */
		static void PermuteSides(int side, TPZVec<int64_t> &id, TPZVec<int> &permutegather)
		{
			permutegather.Resize(1);
			permutegather[0] = 0;
		}
		
	};
	
};

#endif
