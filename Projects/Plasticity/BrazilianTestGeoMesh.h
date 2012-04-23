/*
 *  BrazilianTestGeoMesh.h
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 10/5/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzgmesh.h"
#include <iostream>
#include <string>
#include "pzlog.h"
#include "pzmatrix.h"
#include "pzlog.h"
//#include "TPZTensor.h"
#include "TPZGeoCube.h"
#include "tpzgeoblend.h"


#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"


#include "TPZCompElDisc.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"

#include <iostream>
#include <cstdlib>
#include <math.h>
#include "pzelasmat.h" 
#include "pzelast3d.h"
#include "pzgeopoint.h"
#include "TPZVTKGeoMesh.h"
//#include "TPZTensor.h"


class BrazilianTestGeoMesh {
public:
	
	BrazilianTestGeoMesh();
	~BrazilianTestGeoMesh();
	static TPZGeoMesh *GeoMesh(int h);
	static TPZGeoMesh *GeoBlendMesh(int h);
	static void TransformBlendToLinearMesh(TPZGeoMesh *newlinearmesh,int h);
	static void TransformBlendToLinearMesh2(TPZGeoMesh *newlinearmesh,int h);
	static void ReadMesh(TPZGeoMesh &mesh);
	static TPZGeoMesh * TwoDMesh(int h, int dir);
	static TPZGeoMesh * TwoDMeshSlopeStability(int h, int dir);
	static TPZGeoMesh * TwoDMeshSlopeStability45(int h, int dir);
	static TPZGeoMesh * TwoDMeshSlopeStability452(int h, int dir);
	static TPZGeoMesh * MisesPressure(int h, int dir);
	static TPZGeoMesh * MisesPressure2(int h, int dir);
    static TPZGeoMesh * MalhaPredio();
    static TPZGeoMesh * TwoDMeshII(int h, int dir);
	
	
};
