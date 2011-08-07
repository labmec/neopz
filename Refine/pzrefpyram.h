/**
 * @file
 * @brief Contains the TPZRefPyramid class which implements the uniform refinement of a geometric hexahedral element.
 */
/* class that defines the default refinement of the tetrahedra element */

#ifndef TPZREFPYRAMIDH
#define TPZREFPYRAMIDH

#include "pzreal.h"
#include "pzstack.h"
class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric hexahedral element. \ref refine "Refine"
	 */
	class TPZRefPyramid{
		
	public:
		
		enum{NSubEl = 10};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(TPZGeoEl *gel,int side,int &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int &index);
		static void GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform GetTransform(int side,int son);
		static int FatherSide(int side,int son);
		static void MidCoordSide(int side,TPZVec<REAL> &coord);
		//static int NSubElements();
	};
	
};

#endif
