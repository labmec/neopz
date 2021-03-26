/**
 * @file
 * @brief Contains the TPZRefTetrahedra class which implements the uniform refinement of a geometric tetrahedral element.
 */

#ifndef TPZREFTETRAHEDRAH
#define TPZREFTETRAHEDRAH

class TPZGeoEl;
class TPZGeoElSide;
template<class T>
class TPZTransform;

template<class T>
class TPZVec;
template<class T, int N>
class TPZStack;

namespace pzrefine {
	
	/**
	 * @ingroup refine
	 * @brief Implements the uniform refinement of a geometric tetrahedral element. \ref refine "Refine"
	 */
	class TPZRefTetrahedra : public TPZSavable {
		
	public:
		
		enum{NSubEl = 6};
		
		static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
		static void MidSideNodeIndex(const TPZGeoEl *gel,int side,int64_t &index);
		static void NewMidSideNode(TPZGeoEl *gel,int side,int64_t &index);
		static void GetSubElements(const TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
		static int NSideSubElements(int side);
		static TPZTransform<> GetTransform(int side,int son);
		static int FatherSide(int side,int son);
                public:
int ClassId() const override;

	};
	
};

#endif
