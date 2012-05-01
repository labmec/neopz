/**
 * @file
 * @brief Contains declaration of TPZCheckGeom class which performs a series of consistency tests on geometric transformation.
 */

#ifndef TPZCHECKGEOMH
#define TPZCHECKGEOMH

#include "pzgeoel.h"
#include "pzgmesh.h"

/**
 * @ingroup geometry
 * @brief This class performs a series of consistency tests on geometric transformations between elements. \ref geometry "Geometry"
 */
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
	
};

#endif

