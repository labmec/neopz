/**
 * @file
 * @brief Contains the TPZRefPoint class which implements the uniform refinement of a geometric point element.
 */
/* class that defines the default refinement of the hexaedral element */

#ifndef TPZREFPOINTH
#define TPZREFPOINTH

#include "pzstack.h"

class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

/** 
 * @brief Groups all classes which model the h refinement \n
 * These classes are used as template arguments of @see TPZGeoElement
 */
namespace pzrefine {
	
	/**
	 * @brief Implements the uniform refinement of a geometric point element. \ref refine "Refine"
	 * @ingroup refine
	 * @note Objects of this class implement the uniform refinement of an element 
	*/
	class TPZRefPoint {
		
	public:
		
		enum{NSubEl = 1};
		
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
