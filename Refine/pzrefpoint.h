/* class that defines the default refinement of the hexaedral element */
// -*- c++ -*-

#ifndef TPZREFPOINTH
#define TPZREFPOINTH

#include "pzstack.h"

class TPZGeoEl;
class TPZGeoElSide;
class TPZTransform;

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
#endif
