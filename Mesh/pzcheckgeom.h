//$Id: pzcheckgeom.h,v 1.3 2005-02-28 22:08:50 phil Exp $

#ifndef TPZCHECKGEOMH
#define TPZCHECKGEOMH

#include "pzgeoel.h"
#include "pzgmesh.h"

/// This class performs a series of consistency tests on geometric transformations between elements
class TPZCheckGeom {

	TPZGeoMesh *fMesh;

public:
	TPZCheckGeom();

	int PerformCheck();
	int CheckElement(TPZGeoEl *gel);
	int CheckRefinement(TPZGeoEl *gel);
	int CheckSideTransform(TPZGeoEl *gel, int sidefrom, int sideto);
	int CheckSubFatherTransform(TPZGeoEl *subel, int sidesub);
	void CreateMesh();
	static int main();
//template <class TShape>
//	TPZTransform PrintHighDimTransforms(int side, TPZGeoEl *gel, ostream &out);

};

template<class TShape>
void BuildHigherDimensionSides(TPZStack<int> &highdim, int side);
template<class TShape>
void PrintHigherDimensionSides(ostream &out);
template<class TShape>
void PrintHighDimTransforms(int side, ostream &out);
template<class TShape>
void PrintHighDimTransforms(ostream &out);


#endif

