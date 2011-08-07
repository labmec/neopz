/**
 * @file
 * @brief Contains the TPZRefLinear class which implements the uniform refinement of a geometric linear element.
 */
/* class that defines the default refinement of the hexaedral element */

#ifndef TPZREFLINEARH
#define TPZREFLINEARH

#include "pzstack.h"
class TPZGeoEl;
class TPZTransform;
class TPZGeoElSide;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric linear element. \ref refine "Refine"
	 */
	class TPZRefLinear {
		
	public:
		
		enum{NSubEl = 2};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(TPZGeoEl *gel,int side,int &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int &index);
		static void GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform GetTransform(int side,int son);
		static int FatherSide(int side,int son);
	};
	
};

#endif
