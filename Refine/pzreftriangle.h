/* class that defines the default refinement of the hexaedral element */
// -*- c++ -*-

#ifndef TPZREFTRIANGH
#define TPZREFTRIANGH

#include "pzstack.h"

class TPZGeoEl;
class TPZTransform;

class TPZRefTriangle{
public:

	enum{NSubEl = 4};

	static void Divide(TPZGeoEl *geo,TPZVec<TPZGeoEl *> &SubElVec);
	static void MidSideNodeIndex(TPZGeoEl *gel,int side,int &index);
	static void NewMidSideNode(TPZGeoEl *gel,int side,int &index);
	static void GetSubElements(TPZGeoEl *father,int side, TPZStack<TPZGeoElSide> &subel);
	static int NSideSubElements(int side);
	//static int NSideSubElements(int side);
	static TPZTransform GetTransform(int side,int son);
	static int FatherSide(int side,int son);
	//static int NSubElements();
};
#endif
