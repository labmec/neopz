// _*_ c++ _*_
/* class that defines the default refinement of the tetrahedra element */


#ifndef TPZREFTETRAHEDRAH
#define TPZREFTETRAHEDRAH

class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

template<class T>
class TPZVec;
template<class T, int N>
class TPZStack;

namespace pzrefine {

/// implements the uniform refinement of a geometric tetrahedral element
class TPZRefTetrahedra {

public:

	enum{NSubEl = 6};

	static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
	static void MidSideNodeIndex(TPZGeoEl *gel,int side,int &index);
	static void NewMidSideNode(TPZGeoEl *gel,int side,int &index);
	static void GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
	static int NSideSubElements(int side);
	static TPZTransform GetTransform(int side,int son);
	static int FatherSide(int side,int son);
	//static int NSubElements();
};

};
#endif
