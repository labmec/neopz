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
	TPZCheckGeom(TPZGeoMesh *gmesh = NULL);
	
    /// verify compatibility between elements and their father and between elements and their neighbours
	int PerformCheck();
    
    /// divide all elements and call PerformCheck
    int DivideandCheck();
    
	int CheckElement(TPZGeoEl *gel);
    
    /// check the internal side transformations
    int CheckInternalTransforms(TPZGeoEl *);
    
    /// check the maps between the element and its father
	int CheckRefinement(TPZGeoEl *gel);
    
    /// verify if the mapping between neighbouring elements is conforming
    int CheckNeighbourMap(TPZGeoEl *gel);
    
	int CheckSideTransform(TPZGeoEl *gel, int sidefrom, int sideto);
	int CheckSubFatherTransform(TPZGeoEl *subel, int sidesub);
	void CreateMesh();
	static int main();
	
};

#endif

