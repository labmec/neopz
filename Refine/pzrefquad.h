/**
 * @file
 * @brief Contains the TPZRefQuad class which implements the uniform refinement of a geometric quadrilateral element.
 */

#ifndef TPZREFQUADH
#define TPZREFQUADH

#include "pzreal.h"
#include "pzstack.h"
class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric quadrilateral element. \ref refine "Refine"
	 */
	class TPZRefQuad {
		
	public:
		
		enum{NSubEl = 4};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(const TPZGeoEl *gel,int side,long &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,long &index);
		static void GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform GetTransform(int side,int son);
		static int FatherSide(int side,int son);
	};
	
};

#endif
