/**
 * @file
 * @brief Contains the TPZRefCube class which implements the uniform refinement of a geometric hexahedral element.
 */
/* class that defines the default refinement of the hexaedral element */

#ifndef TPZREFCUBEH
#define TPZREFCUBEH

#include "pzstack.h"
class TPZGeoEl;
class TPZGeoElSide;
template<class T>
class TPZTransform;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric hexahedral element. \ref refine "Refine"
	 */
	class TPZRefCube : public TPZSavable {
		
	public:
		
		enum{NSubEl = 8};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index);
		static void GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		//static int NSideSubElements(int side);
		static TPZTransform<> GetTransform(int side,int son);
		static int FatherSide(int side,int son);
		//static int NSubElements();
                public:
int ClassId() const override;

	};
	
};

#endif
