/**
 * @file
 * @brief Contains declaration of TPZCheckGeom class which performs a series of consistency tests on geometric transformation.
 */
//$Id: pzcheckgeom.h,v 1.4 2005-04-25 02:31:46 phil Exp $

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
	//template <class TShape>
	//	TPZTransform PrintHighDimTransforms(int side, TPZGeoEl *gel, ostream &out);
	
};

/**
template<class TShape>
void BuildHigherDimensionSides(TPZStack<int> &highdim, int side);
template<class TShape>
void PrintHigherDimensionSides(std::ostream &out);
template<class TShape>
void PrintHighDimTransforms(int side, std::ostream &out);
template<class TShape>
void PrintHighDimTransforms(std::ostream &out);
*/

#endif

