/**
 * @file
 * @brief Contains the TPZRefPoint class which implements the uniform refinement of a geometric point element.
 */

#ifndef TPZREFPOINTH
#define TPZREFPOINTH

#include "pzstack.h"
#include "TPZSavable.h"

class TPZGeoEl;
class TPZGeoElSide;
template<class T>
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
	class TPZRefPoint : public TPZSavable {
		
	public:
		
		enum{NSubEl = 1};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index);
		static void GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform<REAL> GetTransform(int side,int son);
		static int FatherSide(int side,int son);
                public:
int ClassId() const override;

	};
	
};

#endif
