/**
 * @file
 * @brief Contains the TPZRefPrism class which implements the uniform refinement of a geometric prism element.
 */

#ifndef TPZREFPRISMH
#define TPZREFPRISMH

#include "pzstack.h"
#include "TPZSavable.h"
class TPZGeoEl;
class TPZGeoElSide;
template<class T>
class TPZTransform;

namespace pzrefine {
	
	/** 
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric prism element. \ref refine "Refine"
	 */
	class TPZRefPrism : public TPZSavable{
		
	public:
		
		enum{NSubEl = 8};
		
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
