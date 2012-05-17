/*
 *  MeshGeneration.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/15/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "MeshGeneration.h"
#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzgeoelside.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "TPZVTKGeoMesh.h"
#include <pzgengrid.h>

using namespace std;


/**
 * Void Constructor
 */	
MeshGeneration::MeshGeneration()
{	
}

/**
 * Destrutor
 */
MeshGeneration::~MeshGeneration()
{
}


void MeshGeneration::Setdimensions(REAL xdimension, REAL ydimension, REAL zthickness)
{
	fxLength = xdimension;
	fyLength = ydimension;
	fThickness = zthickness;	
}

void MeshGeneration::Setdimensions(REAL rdimension, REAL zthickness)
{
	frLength = rdimension;
	fThickness = zthickness;	
}


void MeshGeneration::GeometricMesh2DValidation()
{
	// 2D Cylindrical Domain boundaries
	const int arc1 = -1;
	const int arc2 = -2;
	const int arc3 = -3;
	const int arc4 = -4;
	const int matId = 1;
	const int WellPoint = -5;
	
	gRefDBase.InitializeRefPatterns();
	gRefDBase.InitializeAllUniformRefPatterns();
	ofstream RefinElPatterns("RefElPatterns.txt");
	gRefDBase.WriteRefPatternDBase(RefinElPatterns);
	gRefDBase.ReadRefPatternDBase("RefElPatterns.txt");
	
	int nodenumber = 9;
	REAL ModelRadius = this->frLength;
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	gmesh->NodeVec().Resize(nodenumber);
	
	// Setting node coordantes for Arc3D 1
	int id = 0;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
	gmesh->NodeVec()[id].SetCoord(1,ModelRadius);//coord Y
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
	gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y	
	id++;	
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
	gmesh->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y	
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
	gmesh->NodeVec()[id].SetCoord(1,sqrt(2)*ModelRadius/2.);//coord Y	
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
	id++;
	gmesh->NodeVec()[id].SetNodeId(id);
	gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
	id++;
	
	int elementid = 0;
	// Create Geometrical Arc #1
	// Definition of Arc coordenates
	TPZVec < int > nodeindex(3,0.0);
	nodeindex[0] = 1;	
	nodeindex[1] = 2;
	nodeindex[2] = 5;
	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc1 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1,*gmesh);
	elementid++;
	
	// Create Geometrical Arc #2
	nodeindex[0] = 2;	
	nodeindex[1] = 3;
	nodeindex[2] = 6;			
	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc2 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2,*gmesh);
	elementid++;
	
	// Create Geometrical Arc #3
	nodeindex[0] = 3;	
	nodeindex[1] = 4;
	nodeindex[2] = 7;			
	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc3 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc3,*gmesh);
	elementid++;
	
	// Create Geometrical Arc #4
	nodeindex[0] = 4;	
	nodeindex[1] = 1;
	nodeindex[2] = 8;			
	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc4 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc4,*gmesh); 
	elementid++;
	
	// Create Geometrical Point for fluid injection or Production #1	
	nodeindex.resize(1);
	nodeindex[0] = 0;	
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, WellPoint,*gmesh); 
	elementid++;
	
//	// Create Geometrical triangle #1	
//	nodeindex.resize(3);
//	nodeindex[0] = 0;
//	nodeindex[1] = 1;
//	nodeindex[2] = 2;	
//	TPZGeoElRefPattern< TPZGeoBlend < pzgeom::TPZGeoTriangle > > *triangle1 = new TPZGeoElRefPattern< TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical triangle #2		
//	nodeindex[0] = 0;
//	nodeindex[1] = 2;
//	nodeindex[2] = 3;		
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle2 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical triangle #3		
//	nodeindex[0] = 0;
//	nodeindex[1] = 3;
//	nodeindex[2] = 4;		
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle3 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical triangle #4		
//	nodeindex[0] = 0;
//	nodeindex[1] = 4;
//	nodeindex[2] = 1;		
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle4 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
	
	gmesh->BuildConnectivity();
	
}


