/*
 *  TPZVTKGeoMesh.h
 *  Crack
 *
 *  Created by Cesar Lucci on 16/08/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZVTKGEOMESHH
#define TPZVTKGEOMESHH 

/*
 *  TPZVTKGeoMesh.h
 *
 *  Created by Cesar Lucci on 17/03/10.
 *  Copyright 2010 LabMeC. All rights reserved.
 *
 */

using namespace std;

#include <set>
#include "pzgeoel.h"

class TPZVTKGeoMesh
{
	
public:
	
	TPZVTKGeoMesh();
	~TPZVTKGeoMesh();
	
	/**
	 * Generate an output of all geomesh to VTK
	 */
	static void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, bool matColor = false);
	
	/**
	 * Generate an output of all geomesh to VTK, associating to each one the given data
	 */
	static void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file, TPZVec<int> &elData);
	
	/**
	 * Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKneighbour_material(TPZGeoMesh * gmesh, std::ofstream &file, int neighMaterial, bool matColor = false);
	
	/**
	 * Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKmy_material(TPZGeoMesh * gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor = false);
	
	static int GetVTK_ElType(TPZGeoEl * gel);
};

#endif
