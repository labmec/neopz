/**
 * @file
 * @brief Contains the TPZRefLinear class which implements the uniform refinement of a geometric linear element.
 */
/* class that defines the default refinement of the hexaedral element */

#ifndef TPZREFLINEARH
#define TPZREFLINEARH

#include "pzstack.h"
#include "TPZSavable.h"
class TPZGeoEl;
template<class T>
class TPZTransform;
class TPZGeoElSide;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric linear element. \ref refine "Refine"
	 */
	class TPZRefLinear : public TPZSavable {
		
	public:
		
		enum{NSubEl = 2};
		
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
