/**
 * @file
 * @brief Contains the TPZRefPyramid class which implements the uniform refinement of a geometric hexahedral element.
 */

#ifndef TPZREFPYRAMIDH
#define TPZREFPYRAMIDH

#include "pzreal.h"
#include "pzstack.h"
#include "TPZSavable.h"
class TPZGeoEl;
class TPZGeoElSide;
template<class T>
class TPZTransform;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric hexahedral element. \ref refine "Refine"
	 */
	class TPZRefPyramid : public TPZSavable{
		
	public:
		
		enum{NSubEl = 10, NMidSideNodes = 9};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index);
		static void GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform<REAL> GetTransform(int side,int son);
		static int FatherSide(int side,int son);
		static void MidCoordSide(int side,TPZVec<REAL> &coord);
                public:
int ClassId() const override;

	};
	
};

#endif
