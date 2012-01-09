/**
 * @file
 * @brief Contains the TPZVTKGeoMesh class which implements the graphical mesh to VTK environment to geometric mesh.
 */
/*
 *  TPZVTKGeoMesh.h
 *  Crack
 *
 *  Created by Cesar Lucci on 16/08/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 */

#ifndef TPZVTKGEOMESHH
#define TPZVTKGEOMESHH 

#include <set>
#include "pzgeoel.h"
#include "pzcompel.h"

/**
 * @ingroup post
 * @author Cesar Lucci
 * @since 16/08/10
 * @brief To export a graphical mesh to VTK environment to geometric mesh. \ref post "Post processing"
 */
class TPZVTKGeoMesh
{
	
public:
	
	TPZVTKGeoMesh();
	~TPZVTKGeoMesh();
	
	/** @brief Generate an output of all geomesh to VTK */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = false);
	
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, TPZVec<int> &elData);
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data (int), creates a file with filename given */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename, TPZChunkVector<int> &elData);
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data. Print the values of the variable var */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename,int var);
	
	/**
	 * @brief Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKneighbour_material(TPZGeoMesh *gmesh, std::ofstream &file, int neighMaterial, bool matColor = false);
    
    /**
     * @brief Print the elements that surround a givel geoel
     */
    static void PrintGMeshVTKneighbourhood(TPZGeoMesh * gmesh, int elId, std::ofstream &file);
    static void SetMaterial(TPZGeoEl * gel, int mat);
	
	/**
	 * @brief Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file
	 */
	static void PrintGMeshVTKmy_material(TPZGeoMesh *gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor = false);
	
	static int GetVTK_ElType(TPZGeoEl *gel);
};

#endif
