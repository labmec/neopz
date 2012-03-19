/**
 * @file
 * @brief Contains SPr class which implements the shape functions of a hexahedral three-dimensional element.
 */
// $Id: pzshapeextend.h,v 1.1 2008-10-08 02:44:18 phil Exp $
#ifndef SHAPEEXTENDHPP
#define SHAPEEXTENDHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
//#include "pzgeoelside.h"


/** @brief Groups all classes dedicated to the computation of shape functions */
namespace pzshape {
	
	/**
	 * @brief Implements the shape functions of a hexahedral (3D) element. \ref shape "Shape"
	 * @ingroup shape
	 */
	/** The range of the master element is \f$ [-1,1] \f$ */
	template<class TFather>
	class SPr : public TFather {
		
	public:
		
		
		typedef typename TFather::TMem FatMem;
		typedef typename TFather::Top FatTop;
		
		typedef  pztopology::Pr<FatTop> Top;

		/** @brief Temporary storage to accelerate the computation of shape functions. */
		class TMem
		{
			/** @brief Retained values of the upper shape functions */
			FatMem fUpper;
			/** @brief Retained values of the lower shape functions */
			FatMem fLower;
			/** @brief Retained values of the extension shape functions */
			FatMem fExtension;
			/** @brief Interpolation order of the extension sides */
			int fExtSideOrders[TFather::Top::NSides];
			/** @brief Sign of the shape functions (positive or negative) */
			TPZManVector<int,TFather::Top::NSides> fSign;
			/** @brief Number of shape functions before "my" shapefunctions (only for the extension sides) */
			TPZManVector<int,TFather::Top::NSides> fNShapeBefore;
			/** @brief Number of shape functions after "my" shapefunctions (only for the extension sides) */
			TPZManVector<int,TFather::Top::NSides> fNShapeAfter;
			/** @brief Values of the shape functions associated with the extension sides */
			TPZManVector<REAL> fExtShapes[TFather::Top::NSides];
			/** @brief Number of shape functions for each extension side */
			TPZManVector<int,TFather::Top::NSides> fNShape;
			/** @brief Value of maximum order of the extension sides */
			int fMaxOrder;
			/** @brief Values of the last computed shape functions */
			TPZManVector<REAL> fShape;
			
		};

		/**
		 * @brief Computes the values of the shape functions and their derivatives for a hexahedral element
		 * @param pt (input) point where the shape functions are computed
		 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
		 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
		 * @param phi (output) values of the shape functions
		 * @param dphi (output) values of the derivatives of the shapefunctions		 
		 * @param memory auxiliar structure to accellerate computing
		 */
		/**
		 * These values depend on the point, the order of interpolation and ids of the corner points
		 * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
		 */
		static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi, TMem &memory);
		static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		
		/**
		 * @brief Computes the corner shape functions for a hexahedral element
		 * @param pt (input) point where the shape function is computed
		 * @param phi (output) value of the (8) shape functions
		 * @param dphi (output) value of the derivatives of the (8) shape functions holding the derivatives in a column
		 */
		static void ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
		
		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param order order associated with the side
		 * @return number of shape functions
		 */
		static int NSideShapeF(int side, int order);
		
		/**
		 * Total number of shapefunctions, considering the order
		 * of interpolation of the element
		 * @ param order vector of integers indicating the interpolation order of the element
		 * @ return number of shape functions
		 */
		//static int NShapeF(TPZVec<int> &order);
		
		/**
		 * @brief Number of shapefunctions of the connect associated with the side, considering the order
		 * of interpolation of the element
		 * @param side associated side
		 * @param mem auxiliar structure to accellerate computing
		 * @return number of shape functions
		 */
		static int NSideShapeF(int side, const TMem &mem);
		
		/**
		 * @brief Total number of shapefunctions, considering the order
		 * of interpolation of the element
		 * @param order vector of integers indicating the interpolation order of the element
		 * @return number of shape functions
		 */
		static int NShapeF(TPZVec<int> &order);
		/**
		 * @brief Compute the internal functions of the hexahedral shape function at a point\n
		 * @param x coordinate of the point
		 * @param order maximum order of shape functions to be computed
		 * @param phi shapefunction values
		 * @param dphi values of the derivatives of the shape functions
		 */
		/**
		 * The internal shape functions are the shapefunctions before being multiplied by the corner
		 * shape functions\n
		 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation
		 * determined by the transformation index
		 */
		static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix<REAL> &phi,
								  TPZFMatrix<REAL> &dphi);//,int quad_transformation_index

	};
	
	template<class TFather>
	inline void SPr<TFather>::Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix<REAL> &shape,TPZFMatrix<REAL> &dshape, TMem &mem)
	{
		if(!mem.IsInitialized())
		{
			mem.Initialize(id,order);
		}
		else if(pt == mem.fpt)
		{
			shape = mem.fphi;
			dshape = mem.fdphi;
			return;
		}
		TPZFNMatrix<100> shapeupper,dshapeupper;
		TPZFNMatrix<100> shapelower,dshapelower;
		TPZFNMatrix<100> shapeextend,dshapeextend;
		// compute (or copy) the shape functions of the upper and lower sides
		TFather::Shape(pt,shapeupper,dshapeupper,TMem::fUpper);
		TFather::Shape(pt,shapelower,dshapelower,TMem::fLower);
		// compute (or copy) the shape function values of the shapefunctions which will be extended
		TFather::Shape(pt,shapeextend,dshapeextend,TMem::fExtension);
		
		REAL ExtendPoint = pt[Top::Dimension-1];
		TPZFNMatrix<100> shape1d(mem.fMaxOrder,1),dshape1d(1,mem.fMaxOrder);
		if(ExtendPoint == mem.fExtendPoint)
		{
			shape1d = mem.fshape1d;
			dshape1d = mem.fdshape1d;
		}
		else
		{
			// compute the onedimensional shapefunctions
		}
		// compute the shapefunctions of the lower and upper sides
		int i,j;
		int imax0 = shapelower.Rows();
		REAL fac = (1-ExtendPoint)/2.;
		for(i=0; i<imax0; i++)
		{
			shape(i,0) = shapelower(i,0)*fac;
			for(j=0; j< Top::Dimension-1; j++)
			{
				dshape(j,i) = dshapelower(j,i)*fac;
			}
			dshape(j,i) = -shapelower(i,0)/2;
		}
		int imax1 = imax0+shapeupper.Rows();
		fac = (1+ExtendPoint)/2.;
		for(i=imax0; i<imax1; i++)
		{
			shape(i,0) = shapeupper(i-imax0,0)*fac;
			for(j=0; j< Top::Dimension-1; j++)
			{
				dshape(j,i) = dshapelower(j,i-imax0)*fac;
			}
			dshape(j,i) = shapelower(i-imax0,0)/2;
		}
		// compute the contributions of the extended side shape functions
		
		REAL bubble = (1.+ExtendPoint)*(1.-ExtendPoint);
		// for each side
        int is;
		for(is = 0; is < TFather::Top::NSides; is++)
		{
			int ishape, ishapebef, ishapeaft, ishapeext;
			int next = mem.fExtShapes[is].NElements();
			for(ishape=0; ishape<mem.fNShape[is]; ishape++)
			{
				ishapebef = ishape%mem.fNShapeBefore[is];
				ishapeext = ishape/mem.fNShapeBefore[is];
				ishapeaft = ishape/(mem.fNShapeBefore[is]*next);
				int orgind = ishapebef+ishapeaft*mem.fNShapeBefore[is];
				// this could be lower too
				shape(index,0) = mem.fExtShapes[is][orgind]*shapelower[ishapeext];
			}
		}
		// project the point to the side
		// compute the internal shape functions
		// compute the extended derivatives of the internal shape functions
		// compute the extension function for the side
		// multiply both
	}
	
	/**
	 * Computes the corner shape functions for a hexahedral element
	 */
	template<class TFather>
    inline void SPr<TFather>::ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi)
	{
		TPZManVector<REAL,FatTop::NCornerNodes> phifather;
		TPZFNMatrix<TFather::TOP::Dimension*TFather::TOP::NCornerNodes> dphifather(FatTop::Dimension,FatTop::NCornerNodes);
		TFather::ShapeCorner(pt,phifather,dphifather);
		phi.Resize(Top::NCornerNodes,1);
		dphi.Resize(Top::Dimension,Top::NCornerNodes);
		REAL faclower = -(pt[Top::Dimension-1]-1.)/2.;
		REAL facupper = (pt[Top::Dimension-1]+1.)/2.;
		int i;
		for(i=0; i< FatTop::NCornerNodes; i++)
		{
			phi(i,0) = faclower*phifather(i,0);
			phi(i+FatTop::NCornerNodes,0) = facupper*phifather(i,0);
			dphi(Top::Dimension-1,i) = -0.5*phifather(i,0);
			dphi(Top::Dimension-1,i+FatTop::NCornerNodes) = 0.5*phifather(i,0);
			int j;
			for(j=0; j<FatTop::Dimension; j++)
			{
				dphi(j,i) = faclower*dphifather(j,i);
				dphi(j,i+FatTop::NCornerNodes) = facupper*dphifather(j,i);
			}
		}
	}
	
};

#endif
