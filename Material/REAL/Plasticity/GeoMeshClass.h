#ifndef GEOMESHCLASSH
#define GEOMESHCLASSH

/*
 *  GeoMeshClass.h
 *  PZ
 *
 *  Created by Diogo Cecilio on 10/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

//#include "TPZTensor.h"
class TPZGeoMesh;

class GeoMeshClass {

public:
	static TPZGeoMesh * Talude();
    static void WellBore2d(TPZGeoMesh *gMesh);
	
};

#endif