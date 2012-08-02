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
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "tpzcube.h"
#include "pzgeopoint.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"

#include "TPZVTKGeoMesh.h"
#include <pzgengrid.h>

using namespace std;
using namespace pzgeom;


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

void MeshGeneration::Setdimensions()
{
	fxLength = 0.0;
	fyLength = 0.0;
	fThickness = 0.0;	
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


TPZGeoMesh	* MeshGeneration::MalhaGeom(TPZVec < int > matIdlist)
{

//	Rectangular Mesh without GenGrid 
	int Qnodes = 9;
	int matId = matIdlist[0];
	int bcBottom = matIdlist[1];
	int bcRight = matIdlist[2];
	int bcTop = matIdlist[3];
	int bcLeft = matIdlist[4];
	int WellPoint = matIdlist[5];
	int xfixedPoints = matIdlist[6];
	int yfixedPoints = matIdlist[7];	
	
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	TPZVec<TPZGeoNode> Node(Qnodes);	
	gmesh->SetMaxNodeId(Qnodes-1);
	gmesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	//indice dos nos
	int idN = 0;
	for(int xi = 0; xi < (Qnodes-5)/2; xi++)
	{
		Node[idN].SetNodeId(idN);
		Node[idN].SetCoord(0 ,xi*fxLength);//coord X
		Node[idN].SetCoord(1 ,0. );//coord Y
		gmesh->NodeVec()[idN] = Node[idN];
		idN++;		
	}
	
	for(int yi = 0; yi < (Qnodes-5)/2; yi++)
	{
		Node[idN].SetNodeId(idN);
		Node[idN].SetCoord(0 ,(1-yi)*fxLength );//coord X
		Node[idN].SetCoord(1 ,fyLength );//coord Y
		gmesh->NodeVec()[idN] = Node[idN];
		idN++;
	}
	REAL Lmid = fxLength/2.0;
	REAL hmid = fyLength/2.0;	
	
	Node[idN].SetNodeId(idN);
	Node[idN].SetCoord(0,Lmid);//coord X
	Node[idN].SetCoord(1,0.0);//coord Y
	gmesh->NodeVec()[idN] = Node[idN];
	idN++;	
	
	Node[idN].SetNodeId(idN);
	Node[idN].SetCoord(0,fxLength);//coord X
	Node[idN].SetCoord(1,hmid);//coord Y
	gmesh->NodeVec()[idN] = Node[idN];
	idN++;
	
	Node[idN].SetNodeId(idN);
	Node[idN].SetCoord(0,Lmid);//coord X
	Node[idN].SetCoord(1,fyLength);//coord Y
	gmesh->NodeVec()[idN] = Node[idN];		
	idN++;	
	
	Node[idN].SetNodeId(idN);
	Node[idN].SetCoord(0,0.0);//coord X
	Node[idN].SetCoord(1,hmid);//coord Y
	gmesh->NodeVec()[idN] = Node[idN];		
	idN++;
	
	//	Production or injection point
	Node[idN].SetNodeId(idN);
	Node[idN].SetCoord(0,Lmid);//coord X
	Node[idN].SetCoord(1,hmid);//coord Y
	gmesh->NodeVec()[idN] = Node[idN];		
	
	
	
	//	Element Index
	
	int idEl = 0;
	TopolLine[0] = 0;
	TopolLine[1] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcBottom,*gmesh);
	idEl++;
	
	TopolLine[0] = 4;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcBottom,*gmesh);
	idEl++;	
	
	TopolLine[0] = 1;
	TopolLine[1] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcRight,*gmesh);
	idEl++;
	
	TopolLine[0] = 5;
	TopolLine[1] = 2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcRight,*gmesh);
	idEl++;	
	
	TopolLine[0] = 2;
	TopolLine[1] = 6;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcTop,*gmesh);
	idEl++;
	
	TopolLine[0] = 6;
	TopolLine[1] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcTop,*gmesh);
	idEl++;	
	
	TopolLine[0] = 3;
	TopolLine[1] = 7;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcLeft,*gmesh);
	idEl++;
	
	TopolLine[0] = 7;
	TopolLine[1] = 0;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcLeft,*gmesh);
	idEl++;	
	
	TopolQuad[0] = 0;
	TopolQuad[1] = 4;
	TopolQuad[2] = 8;
	TopolQuad[3] = 7;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (idEl,TopolQuad,matId,*gmesh);
	idEl++;
	
	TopolQuad[0] = 4;
	TopolQuad[1] = 1;
	TopolQuad[2] = 5;
	TopolQuad[3] = 8;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (idEl,TopolQuad,matId,*gmesh);
	idEl++;
	
	TopolQuad[0] = 8;
	TopolQuad[1] = 5;
	TopolQuad[2] = 2;
	TopolQuad[3] = 6;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (idEl,TopolQuad,matId,*gmesh);
	idEl++;
	
	TopolQuad[0] = 7;
	TopolQuad[1] = 8;
	TopolQuad[2] = 6;
	TopolQuad[3] = 3;
	new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (idEl,TopolQuad,matId,*gmesh);
	idEl++;
	
	TopolPoint[0] = 6;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idEl,TopolPoint,xfixedPoints,*gmesh);
	idEl++;
	
	TopolPoint[0] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idEl,TopolPoint,xfixedPoints,*gmesh);
	idEl++;
	
	TopolPoint[0] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idEl,TopolPoint,yfixedPoints,*gmesh);
	idEl++;
	
	TopolPoint[0] = 7;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idEl,TopolPoint,yfixedPoints,*gmesh);
	idEl++;
	
	TopolPoint[0] = 8;
	new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idEl,TopolPoint,WellPoint,*gmesh);
	idEl++;	
	
	gmesh->BuildConnectivity();
	
	ofstream arg("gmesh.txt");
	gmesh->Print(arg);
	
	return gmesh;	
	
}


// description of Geometry and application
TPZGeoMesh * MeshGeneration::GeometricMesh2DValidation(TPZVec < int > matIdlist)
{
	// 2D Cylindrical Domain boundaries
	int matId = matIdlist[0];
	int arc1 = matIdlist[1];
	int arc2 = matIdlist[2];
	int arc3 = matIdlist[3];
	int arc4 = matIdlist[4];
	int Point1 = matIdlist[5];
	int Point2 = matIdlist[6];
	int Point3 = matIdlist[7];
	int Point4 = matIdlist[8];	
	int WellPoint = matIdlist[9];
	
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
	
	// Create Geometrical Point #1	
	nodeindex.resize(1);
	nodeindex[0] = 1;
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint1 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point1,*gmesh); 
	elementid++;
	
	// Create Geometrical Point #2
	nodeindex.resize(1);
	nodeindex[0] = 3;	
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint2 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point2,*gmesh); 
	elementid++;
	
	// Create Geometrical Point #3
	nodeindex.resize(1);
	nodeindex[0] = 2;	
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint3 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point3,*gmesh); 
	elementid++;
	
	// Create Geometrical Point #4
	nodeindex.resize(1);
	nodeindex[0] = 4;	
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint4 = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, Point4,*gmesh); 
	elementid++;	
	
	// Create Geometrical Point for fluid injection or Production #1	
	nodeindex.resize(1);
	nodeindex[0] = 0;	
	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, WellPoint,*gmesh); 
	elementid++;	
	
	// Create Geometrical triangle #1	
	nodeindex.resize(3);
	nodeindex[0] = 0;
	nodeindex[1] = 1;
	nodeindex[2] = 2;	
	TPZGeoElRefPattern< TPZGeoBlend < pzgeom::TPZGeoTriangle > > *triangle1 = new TPZGeoElRefPattern< TPZGeoBlend < pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
	elementid++;
	
	// Create Geometrical triangle #2		
	nodeindex[0] = 0;
	nodeindex[1] = 2;
	nodeindex[2] = 3;		
	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle2 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
	elementid++;
	
	// Create Geometrical triangle #3		
	nodeindex[0] = 0;
	nodeindex[1] = 3;
	nodeindex[2] = 4;		
	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle3 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
	elementid++;
	
	// Create Geometrical triangle #4		
	nodeindex[0] = 0;
	nodeindex[1] = 4;
	nodeindex[2] = 1;		
	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > *triangle4 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle > > (elementid,nodeindex, matId,*gmesh);
	elementid++;
	
	gmesh->BuildConnectivity();
	
	ofstream arg("CricularGeoMesh.txt");
	gmesh->Print(arg);	
	
	return gmesh;
	
}

TPZGeoMesh * MeshGeneration::MalhaGeoGravenobj(int nLayers, REAL LlengthFootFault, REAL DipFaultAngleleft, REAL DipFaultAngleright, REAL WellFaultlength, TPZVec <bool> Productionlayer, bool InterfaceElement)
{
	// This method created a graven geometry
	
	// Rock Blocks bondaries
	// Right Block 
	int RBright = -23;
	int RBleft = 43;
	int RBBot = -13;
	int RBTop = -33;
	
	// Left Block 
	int LBright = 22;
	int LBleft = -42;
	int LBBot = -12;
	int LBTop = -32;
	
	// Graven Block 
	int GBright = -21;
	int GBleft = -41;
	int GBBot = -11;
	int GBTop = -31;
	
	int WellLine = -5;	
	
	bool mergeleft = true;
	bool mergeright = true;
	
	
	// Creating Geometrical Mesh Objects
	
	int nProdLayesr = 0.0;
	REAL ndivideleft = 1.0;
	REAL ndivideright = 1.0;
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{	
		if (Productionlayer[ilayer] && Productionlayer[ilayer+1]) {
			nProdLayesr++;
		}
		
	}
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;
	int TotalNodes = 8+4*(nLayers-1)+nProdLayesr*(2*((ndivideleft+ndivideright-1))+3);
	gmesh->SetMaxNodeId(TotalNodes-1);
	gmesh->NodeVec().Resize(TotalNodes);
	
	
	
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolTria(3);	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolPoint(1);
	TPZVec<TPZGeoNode> Node(TotalNodes);
	
	//	Constan h for all layers
	REAL hlayerj = 100.0;
	TPZVec<REAL> hlayers(nLayers+1,hlayerj);
	//	Horizon at level zero
	hlayers[0] = 0.0;
	//	Constant Material ID for all layers
	TPZVec<int> layersIdRB(nLayers,3);
	TPZVec<int> layersIdLB(nLayers,2);
	TPZVec<int> layersIdGB(nLayers,1);
	
	//	for(int i = 0; i < nLayers; i++)
	//	{
	//		layersIdLB[i]=i+1;
	//		layersIdRB[i]=i+2;
	//		layersIdGB[i]=i+3;
	//	}
	
	REAL X = 0.0;
	REAL Y = 0.0;
	REAL PI = atan(1.)*4.;	
	
	//	Distance between foot faults
	REAL GravenBaseLength = 1.0*LlengthFootFault;
	
	
	// Node index
	int idN = 0;
	int idE = 0;
	int lasSequencenodId = 0;
	int contlayers = 0;
	
	for(int cont = 0; cont < (4 + 2*(nLayers-1)); cont++)
	{
		
		idN = cont;
		if (cont%2==0) 
		{
			Y = 0.0;
			// Calculate the sum of all thickness
			for(int l = 0; l <= contlayers; l++)
			{
				Y += hlayers[l];
			}
			X = 0.0;
			contlayers++;
		}
		else 
		{ 
			X = LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleleft))/(sin((PI/180.0)*DipFaultAngleleft)));
		}
		
		Node[idN].SetNodeId(idN);
		Node[idN].SetCoord(0,X);//coord X
		Node[idN].SetCoord(1,Y);//coord Y
		gmesh->NodeVec()[idN] = Node[idN];
		
	}
	idN++;
	
	
	//	Left Rock block Element Index
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{
		TopolQuad[0] = 0+2*ilayer;
		TopolQuad[1] = 1+2*ilayer;
		TopolQuad[2] = 3+2*ilayer;
		TopolQuad[3] = 2+2*ilayer;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdLB[ilayer],*gmesh);
		idE++;
	}
	
	
	//	Bottom bondary
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,LBBot,*gmesh);
	idE++;
	
	
	for(int inode = 0; inode < (idN - 2); inode++)
	{
		//	Left bondary
		if (inode%2==0) 
		{
			TopolLine[0] = inode;
			TopolLine[1] = inode+2;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,LBleft,*gmesh);
			idE++;		
		}
		//	Right bondary
		else 	
		{
			TopolLine[0] = inode;
			TopolLine[1] = inode+2;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,LBright,*gmesh);
			idE++;
		}
		
		
	}
	
	//	Top bondary	
	TopolLine[0] = 0+2*nLayers;
	TopolLine[1] = 1+2*nLayers;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,LBTop,*gmesh);
	idE++;
	
	
	//	Second block with the same Dip Fault Angle
	
	// Node index
	contlayers = 0;
	
	for(int cont = idN; cont < (8 + 4*(nLayers-1)); cont++)
	{
		
		idN = cont;
		if (cont%2==0) 
		{
			Y = 0.0;
			// Calculate the sum of all thickness
			for(int l = 0; l <= contlayers; l++)
			{
				Y += hlayers[l];
			}
			X = GravenBaseLength+2*LlengthFootFault;
			contlayers++;
		}
		else 
		{ 
			X = GravenBaseLength+2*LlengthFootFault - (LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleright))/(sin((PI/180.0)*DipFaultAngleright))));
		}
		
		Node[idN].SetNodeId(idN);
		Node[idN].SetCoord(0,X);//coord X
		Node[idN].SetCoord(1,Y);//coord Y
		gmesh->NodeVec()[idN] = Node[idN];
		
	}
	
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{
		TopolQuad[0] = 0+2*nLayers+2+2*ilayer;
		TopolQuad[1] = 1+2*nLayers+2+2*ilayer;
		TopolQuad[2] = 3+2*nLayers+2+2*ilayer;
		TopolQuad[3] = 2+2*nLayers+2+2*ilayer;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdRB[ilayer],*gmesh);
		idE++;
	}
	
	
	//	Bottom bondary
	TopolLine[0] = 0+2*nLayers+2;
	TopolLine[1] = 1+2*nLayers+2;	
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBBot,*gmesh);
	idE++;
	
	
	for(int inode = 0+2*nLayers+2; inode <= (idN - 2); inode++)
	{
		//	Left bondary
		if (inode%2!=0) 
		{
			TopolLine[0] = inode;
			TopolLine[1] = inode+2;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBleft,*gmesh);
			idE++;		
		}
		//	Right bondary
		else 	
		{
			TopolLine[0] = inode;
			TopolLine[1] = inode+2;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBright,*gmesh);
			idE++;
		}
		
		
	}
	
	//	Top bondary	
	TopolLine[0] = 0+4*nLayers+2;
	TopolLine[1] = 1+4*nLayers+2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBTop,*gmesh);
	idE++;
	
	// inserting injections/productions points
	REAL Xleft = 0.0;
	REAL Xright = 0.0;
	REAL deltaxleft = 0.0;
	REAL deltaxright = 0.0;	
	idN++;// current node
	
	//	Insering right and left boundaries
	
	lasSequencenodId = idN;	
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{		
		
		if (Productionlayer[ilayer]==true) 
		{
			
			Y = 0.0;
			// Calculate the sum of all thickness
			for(int l = 0; l <= ilayer; l++)
			{
				Y += hlayers[l];
			}
			Xleft = LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleleft))/(sin((PI/180.0)*DipFaultAngleleft)));	
			Xright = GravenBaseLength+2*LlengthFootFault - (LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleright))/(sin((PI/180.0)*DipFaultAngleright))));
			deltaxleft = abs(LlengthFootFault+WellFaultlength-Xleft)/ndivideleft;				
			deltaxright = abs(Xright-LlengthFootFault-WellFaultlength)/ndivideright;
			
			for(int cont = 0; cont < ((ndivideleft+ndivideright)-1); cont++)
			{
				
				if (cont < ndivideleft) 
				{	
					Xleft += deltaxleft;
					Node[idN].SetNodeId(idN);
					Node[idN].SetCoord(0,Xleft);//coord X
					Node[idN].SetCoord(1,Y);//coord Y
					gmesh->NodeVec()[idN] = Node[idN];
					idN++;					
				}
				else 
				{
					Xleft += deltaxright;
					Node[idN].SetNodeId(idN);
					Node[idN].SetCoord(0,Xleft);//coord X
					Node[idN].SetCoord(1,Y);//coord Y
					gmesh->NodeVec()[idN] = Node[idN];
					idN++;					
				}
				
				
				
			}
			
			if (Productionlayer[ilayer] && Productionlayer[ilayer+1]) {
				Y = 0.0;
				// Calculate the sum of all thickness in the middle of the layer
				for(int l = 0; l <= ilayer; l++)
				{
					Y += hlayers[l];
				}
				Y += hlayers[ilayer+1]/2;
				Xleft = LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleleft))/(sin((PI/180.0)*DipFaultAngleleft)));	
				Xright = GravenBaseLength+2*LlengthFootFault - (LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleright))/(sin((PI/180.0)*DipFaultAngleright))));
				deltaxleft = abs(LlengthFootFault+WellFaultlength-Xleft)/ndivideleft;				
				deltaxright = abs(Xright-LlengthFootFault-WellFaultlength)/ndivideright;
				
				Xleft = LlengthFootFault+WellFaultlength-deltaxleft;
				Node[idN].SetNodeId(idN);
				Node[idN].SetCoord(0,Xleft);//coord X
				Node[idN].SetCoord(1,Y);//coord Y
				gmesh->NodeVec()[idN] = Node[idN];					
				idN++;					
				
				Xleft = LlengthFootFault+WellFaultlength+deltaxright;
				Node[idN].SetNodeId(idN);
				Node[idN].SetCoord(0,Xleft);//coord X
				Node[idN].SetCoord(1,Y);//coord Y
				gmesh->NodeVec()[idN] = Node[idN];
				idN++;					
				
				Xleft = LlengthFootFault+WellFaultlength;
				Node[idN].SetNodeId(idN);
				Node[idN].SetCoord(0,Xleft);//coord X
				Node[idN].SetCoord(1,Y);//coord Y
				gmesh->NodeVec()[idN] = Node[idN];
				TopolPoint[0] = idN;
				new TPZGeoElRefPattern< pzgeom::TPZGeoPoint > (idE,TopolPoint,WellLine,*gmesh);
				idE++;					
				idN++;						
			}				
			
			
		}
		
		if (ilayer+1 == nLayers && Productionlayer[nLayers]==true ) 
		{
			
			Y = 0.0;
			// Calculate the sum of all thickness
			for(int l = 0; l <= nLayers; l++)
			{
				Y += hlayers[l];
			}
			Xleft = LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleleft))/(sin((PI/180.0)*DipFaultAngleleft)));	
			Xright = GravenBaseLength+2*LlengthFootFault - (LlengthFootFault - Y*((cos((PI/180.0)*DipFaultAngleright))/(sin((PI/180.0)*DipFaultAngleright))));
			deltaxleft = abs(LlengthFootFault+WellFaultlength-Xleft)/ndivideleft;				
			deltaxright = abs(Xright-LlengthFootFault-WellFaultlength)/ndivideright;
			
			for(int cont = 0; cont < ((ndivideleft+ndivideright)-1); cont++)
			{
				
				if (cont < ndivideleft) 
				{	
					Xleft += deltaxleft;
					Node[idN].SetNodeId(idN);
					Node[idN].SetCoord(0,Xleft);//coord X
					Node[idN].SetCoord(1,Y);//coord Y
					gmesh->NodeVec()[idN] = Node[idN];
					idN++;					
				}
				else 
				{
					Xleft += deltaxright;
					Node[idN].SetNodeId(idN);
					Node[idN].SetCoord(0,Xleft);//coord X
					Node[idN].SetCoord(1,Y);//coord Y
					gmesh->NodeVec()[idN] = Node[idN];
					idN++;					
				}
			}
		}
		
	}

	
	//	Bottom bondary Graven Block
	TopolLine[0] = 1;
	TopolLine[1] = 1+2*nLayers+2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,GBBot,*gmesh);
	idE++;	
		
	//	Graven Block with inyection point
	
	//	Insertion of additional nodes
	
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{	
		//	Next sequence of nodes ids
		idN = lasSequencenodId;	
		
		if (Productionlayer[ilayer] && Productionlayer[ilayer+1]) {
			
			// First lef Quad in the i-layer
			TopolQuad[0] = idN;
			TopolQuad[1] = 1+2*ilayer;
			TopolQuad[2] = 3+2*ilayer;
			TopolQuad[3] = idN+(ndivideleft+ndivideright-1)+3;
			idN++;
			new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
			idE++;		
			
			// intern left Quads in the i-layer
			while (idN < (ndivideleft+lasSequencenodId-1)) 
			{
				TopolQuad[0] = idN;
				TopolQuad[1] = idN-1;
				TopolQuad[2] = idN-1+(ndivideleft+ndivideright-1)+3;
				TopolQuad[3] = idN+(ndivideleft+ndivideright-1)+3;
				idN++;			
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;				
			}
			
			// intern four well neighborhood Quads in the i-layer
			if (idN == (ndivideleft+lasSequencenodId-1)) {
				TopolQuad[0] = idN;
				TopolQuad[1] = idN-1;
				TopolQuad[2] = ndivideleft+ndivideright+lasSequencenodId-1;
				TopolQuad[3] = ndivideleft+ndivideright+lasSequencenodId+1;		
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;
				
				TopolQuad[0] = ndivideleft+ndivideright+lasSequencenodId+1;
				TopolQuad[1] = ndivideleft+ndivideright+lasSequencenodId-1;
				TopolQuad[2] = idN-1+(ndivideleft+ndivideright-1)+3;
				TopolQuad[3] = idN+(ndivideleft+ndivideright-1)+3;		
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;
				
				TopolQuad[0] = idN+1;
				TopolQuad[1] = idN;
				TopolQuad[2] = ndivideleft+ndivideright+lasSequencenodId+1;
				TopolQuad[3] = ndivideleft+ndivideright+lasSequencenodId;		
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;
				
				TopolQuad[0] = ndivideleft+ndivideright+lasSequencenodId;
				TopolQuad[1] = ndivideleft+ndivideright+lasSequencenodId+1;
				TopolQuad[2] = idN+(ndivideleft+ndivideright-1)+3;
				TopolQuad[3] = idN+1+(ndivideleft+ndivideright-1)+3;	
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;		
			}
			idN++;
			
			// intern right Quads in the i-layer
			while (idN > (ndivideleft+lasSequencenodId-1) && idN < (ndivideleft+ndivideright+lasSequencenodId-2)) 
			{
				TopolQuad[0] = idN+1;
				TopolQuad[1] = idN;
				TopolQuad[2] = idN+(ndivideleft+ndivideright-1)+3;
				TopolQuad[3] = idN+1+(ndivideleft+ndivideright-1)+3;
				idN++;			
				new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
				idE++;				
			}			
			
			
			// Last right Quad in the i-layer
			TopolQuad[0] = 1+2*nLayers+2+2*ilayer;
			TopolQuad[1] = idN;
			TopolQuad[2] = idN+(ndivideleft+ndivideright-1)+3;
			TopolQuad[3] = 3+2*nLayers+2+2*ilayer;	
			new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);		
			idE++;
			
			//	updating the next i-layer node sequence
			lasSequencenodId += (ndivideleft+ndivideright-1)+3;
		}
		else {
			
			TopolQuad[0] = 1+2*nLayers+2+2*ilayer;
			TopolQuad[1] = 1+2*ilayer;
			TopolQuad[2] = 3+2*ilayer;
			TopolQuad[3] = 3+2*nLayers+2+2*ilayer;
			new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdGB[ilayer],*gmesh);
			idE++;
			
		}
		
		
	}
	
	//	Top bondary Graven Block	
	TopolLine[0] = 1+2*nLayers;
	TopolLine[1] = 1+4*nLayers+2;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,GBTop,*gmesh);
	idE++;		
	
	//---------  Interface Elements  --------------------------------------------
	if (InterfaceElement==true) 
	{
		gmesh->AddInterfaceMaterial(1,2,22);
		gmesh->AddInterfaceMaterial(2,1,22);
		gmesh->AddInterfaceMaterial(3,1,43);
		gmesh->AddInterfaceMaterial(1,3,43);
	}
	
	gmesh->BuildConnectivity();
	
	ofstream arg("GeologicalGravenmesh.txt");
	gmesh->Print(arg);	
	
	return gmesh;
}



// Create a Mesh based on GID Mesh Dump file 
// You must define the boundaries IDs
// Your must to defined Material IDs
// You must to defined interface elements and its IDs
// Into the GID mesh creation process

TPZGeoMesh * MeshGeneration::GeometricGIDMesh(std::string FiletoRead)

{
	// File to read
	
	string FileName;
	std::string stringTemp;
	FileName = FiletoRead;

	
	int nMats = 0;
	int numnodes=0;
	int numelements=0;
	int elements3DT=0;
	int elements3DH=0;	
	int elements2DT=0;
	int elements2DQ=0;	
	int elements1D=0;
	int elements0D=0;
	
	
//	Scanning for total Number of Nodes and differents Dimension Elements
	int NumEntitiestoRead;
	std::vector <std::string> SentinelString;
	{
		
		// reading a general mesh information by filter
		ifstream read (FileName.c_str());
		std::string FlagString;
		int cont = 0;
		int flag = -1;		

		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			flag = str.find("---");
			
			if (flag >= 0)
			{
				if(str != "--- ELEMENTS ---") 
				{
					SentinelString.push_back(str);
				}
			}
			if(str == "LINEAR") 
			{
				SentinelString.push_back(str);
			}
			if(str == "TRIANGLE") 
			{
				SentinelString.push_back(str);
			}
			if(str == "QUADRILATERAL") 
			{
				SentinelString.push_back(str);
			}
			if(str == "TETRAHEDRA") 
			{
				SentinelString.push_back(str);
			}
			if(str == "HEXAHEDRA") 
			{
				SentinelString.push_back(str);
			}			
			
		}
		
		FlagString = "EndReading";
		SentinelString.push_back(FlagString);
	}
	

	
	
	NumEntitiestoRead = SentinelString.size();
	std::vector <int> GeneralData(NumEntitiestoRead,0);
	std::vector <int> DataToProcess(NumEntitiestoRead,-1);		
	
	{
		
		// reading a general mesh information by filter
		ifstream read (FileName.c_str());
		std::string FlagString;
		int cont = 0;
		
		
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			if(SentinelString[cont] == "") 
			{
				cont++;	
			}			
						
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if(str != "" && FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) {
					char buftemp[1024];
					read.getline(buftemp, 1024);
					std::string strtemp(buftemp);
					GeneralData[cont]++;
					if(strtemp == "") 
					{
						FlagString = "";
						GeneralData[cont]--;
						cont++;
						std::cout << "Scanning General Data -> done!" << std::endl;
						break;
					}
				}
				
			}	
		
			
		}
	}
	
	for (int i = 0 ; i < NumEntitiestoRead; i++ )
	{
		if(SentinelString[i] == "--- USED MATERIALS ---") 
		{
			nMats=GeneralData[i];
			DataToProcess[i]=0;
		}
		if(SentinelString[i] == "--- CONDITIONS OVER NODES ---") 
		{
			if(GeneralData[i] !=0)
			{
			GeneralData[i]--;	
			elements0D=GeneralData[i];
			}
			else 
			{

			}

			DataToProcess[i]=1;			
		}	
		if(SentinelString[i] == "--- NODES ---") 
		{
			numnodes=GeneralData[i];
			DataToProcess[i]=2;			
		}		
		if(SentinelString[i] == "LINEAR") 
		{
			elements1D=GeneralData[i];
			DataToProcess[i]=3;		
		}
		if(SentinelString[i] == "TRIANGLE") 
		{
			elements2DT=GeneralData[i];
			DataToProcess[i]=4;		
		}
		if(SentinelString[i] == "QUADRILATERAL") 
		{
			elements2DQ=GeneralData[i];
			DataToProcess[i]=5;			
		}
		if(SentinelString[i] == "TETRAHEDRA") 
		{
			elements3DT=GeneralData[i];
			DataToProcess[i]=6;			
		}
		if(SentinelString[i] == "HEXAHEDRA") 
		{
			elements3DH=GeneralData[i];
			DataToProcess[i]=7;			
		}		
		
	}

		
	numelements=elements3DT+elements3DH+elements2DT+elements2DQ+elements1D+elements0D;

	
//  Mesh Creation
	
	TPZGeoMesh * gmesh = new TPZGeoMesh;	
	gmesh -> NodeVec().Resize(numnodes);
	// needed for node insertion
	const int Tnodes = numnodes;
	TPZVec <TPZGeoNode> Node(Tnodes);
	int nodeId;
	double nodecoordX , nodecoordY , nodecoordZ ;
	// needed for triangular elements insertion
	int elementId = 0;
	int matElId = 0;
	int Layerid = 0; // data not used here
	int ContMats = 0;
	int ContNode = 0;
	int ContPoint = 0;
	int ContLine = 0;
	int ContTrian = 0;
	int ContQuad = 0;
	int ContTet = 0;
	int ContHex = 0;	

	TPZManVector <int> TopolPoint(1);	
	TPZManVector <int> TopolLine(2);
	TPZManVector <int> TopolTriangle(3);
	TPZManVector <int> TopolQuad(4);
	TPZManVector <int> TopolTet(4);
	TPZManVector <int> TopolHex(8);	


	{
		
		// reading a general mesh information by filter
		ifstream read (FileName.c_str());
		std::string FlagString;
		int cont = 0;
		int dim = 0;
		int flag = 0;
		while(read)
		{
			char buf[1024];	
			read.getline(buf, 1024);
			std::string str(buf);
			std::string strtemp="InitialState";			
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			
			if(SentinelString[cont] == "") 
			{
				cont++;	
			}			
			
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if(str != "" && FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) {
					
			switch (DataToProcess[cont]) {
				case 0:
				{
					 //"--- USED MATERIALS ---"
					if (GeneralData[cont] != 0)
					{
						read.getline(buf, 1024);
						ContMats++;
					}
					if(ContMats == nMats)
					{
						strtemp = "";	
					}					
				}
					break;
				case 1:
				{
					//"--- CONDITIONS OVER NODES ---"
					// 0D Elements
					if (GeneralData[cont] != 0)
					{
						if(flag == 0)
						{
						read.getline(buf, 1024);
						flag++;
						}
					read >> TopolPoint[0]; //node 1
					read >> nodecoordZ;												
					TopolPoint[0]--;
					ContPoint++;						
					TPZGeoEl *Point = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (numelements - elements0D + ContPoint, TopolPoint, nMats+ContPoint,*gmesh); 						

					}
					if(ContPoint == elements0D)
					{
						strtemp = "";	
					}					
				}
					break;
				case 2:
				{
					if (GeneralData[cont] != 0)
					{
					//"--- NODES ---"

					if (!dim)
					{							
						read.getline(buf, 1024);
						int cont2 = 0;
						char *p=buf, *q;
						while (p) {
							q = strchr(p,' ');
							if(!q)
							{
								break;
							}
							*q = 0;
							if(cont2==0) nodeId = atoi(p);
							else if(cont2==1) Layerid = atoi(p);
							else if(cont2==2) nodecoordX = atof(p);
							else if(cont2==3) nodecoordY = atof(p);
//							else if(cont2==4) nodecoordZ = atof(p);
							p = q+1;							
							while(q && *p==' ')
								p++;
							cont2++;
						}
						if (cont2 == 3)
						{
							dim = 2;
							nodecoordY = atof(p);							
							nodecoordZ = 0.0;							
						}
						else 
						{
							dim = 3;
							nodecoordZ = atof(p);				
						}

					}
					else 
					{
						read >> nodeId;
						read >> Layerid;						
						read >> nodecoordX;
						read >> nodecoordY;
						if (dim == 2) {
							nodecoordZ = 0.0;
						}		
						else 
						{
							read >> nodecoordZ;// 2D needed to correct
						}
					}

					Node[nodeId-1].SetNodeId(nodeId);
					Node[nodeId-1].SetCoord(0,nodecoordX);
					Node[nodeId-1].SetCoord(1,nodecoordY);
					Node[nodeId-1].SetCoord(2,nodecoordZ);
					gmesh->NodeVec()[nodeId-1] = Node[nodeId-1];
					ContNode++;
					}
					if(ContNode == numnodes)
					{
						strtemp = "";	
					}					
				}
					break;
				case 3:
				{
					//"LINEAR"
					if (GeneralData[cont] != 0)
					{
					read >> elementId;
					read >> matElId;  // Material ID 
					read >> Layerid;						
					read >> TopolLine[0]; //node 2
					read >> TopolLine[1]; //node 3
					elementId--;
					TopolLine[0]--;
					TopolLine[1]--;
					TPZGeoEl * Line = new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (elementId,TopolLine,matElId,*gmesh);
					ContLine++;					
					}
					if(ContLine == elements1D)
					{
						strtemp = "";	
					}						
				}
					break;
				case 4:
				{
					if (GeneralData[cont] != 0)
					{
					//"TRIANGLE"
					read >> elementId;
					read >> matElId;  // Material ID
					read >> Layerid;						
					read >> TopolTriangle[0]; //node 1
					read >> TopolTriangle[1]; //node 2
					read >> TopolTriangle[2]; //node 3
					elementId--;						
					TopolTriangle[0]--;
					TopolTriangle[1]--;
					TopolTriangle[2]--;					
					TPZGeoEl * triangle = new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle> (elementId, TopolTriangle, matElId, *gmesh);
					ContTrian++;							
					}
					if(ContTrian == elements2DT)
					{
						strtemp = "";	
					}					
				}
					break;
				case 5:
				{
					//"QUADRILATERAL"
					if (GeneralData[cont] != 0)
					{
					read >> elementId;
					read >> matElId;  // Material ID
					read >> Layerid;						
					read >> TopolQuad[0]; //node 1
					read >> TopolQuad[1]; //node 2
					read >> TopolQuad[2]; //node 3
					read >> TopolQuad[3]; //node 4						
					elementId--;						
					TopolQuad[0]--;
					TopolQuad[1]--;
					TopolQuad[2]--;	
					TopolQuad[3]--;							
					TPZGeoEl * Quadrilateral = new TPZGeoElRefPattern< pzgeom::TPZGeoQuad> (elementId, TopolQuad, matElId, *gmesh);
					ContQuad++;						
					}
					if(ContQuad == elements2DQ)
					{
						strtemp = "";	
					}						
				}
					break;
				case 6:
				{
					//"TETRAHEDRA"
					if (GeneralData[cont] != 0)
					{
					read >> elementId;
					read >> matElId;  // Material ID
					read >> Layerid;						
					read >> TopolTet[0]; //node 1
					read >> TopolTet[1]; //node 2
					read >> TopolTet[2]; //node 3
					read >> TopolTet[3]; //node 4						
					elementId--;						
					TopolTet[0]--;
					TopolTet[1]--;
					TopolTet[2]--;	
					TopolTet[3]--;	
					TPZGeoEl * tetra = new TPZGeoElRefPattern< pzgeom::TPZGeoTetrahedra> (elementId, TopolTet, matElId, *gmesh);
					ContTet++;
					}
					if(ContTet == elements3DT)
					{
						strtemp = "";	
					}
				}
					break;
				case 7:
				{
					//"HEXAHEDRA"
					if (GeneralData[cont] != 0)
					{
					read >> elementId;
					read >> matElId;  // Material ID
					read >> Layerid;						
					read >> TopolHex[0]; //node 1
					read >> TopolHex[1]; //node 2
					read >> TopolHex[2]; //node 3
					read >> TopolHex[3]; //node 4
					read >> TopolHex[4]; //node 5
					read >> TopolHex[5]; //node 6
					read >> TopolHex[6]; //node 7
					read >> TopolHex[7]; //node 8						
					elementId--;						
					TopolHex[0]--;
					TopolHex[1]--;
					TopolHex[2]--;	
					TopolHex[3]--;
					TopolHex[4]--;
					TopolHex[5]--;
					TopolHex[6]--;	
					TopolHex[7]--;						
					TPZGeoEl * Hex = new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (elementId, TopolHex, matElId, *gmesh);					
					ContHex++;
					}
					if(ContHex == elements3DH)
					{
						strtemp = "";	
					}
				}
					break;					
				default:
				{
						strtemp = "";
				}
					break;
			}		
					
					if(strtemp == "") 
					{
						FlagString = "";
						cont++;
						std::cout << "Reading Data -> done!" << std::endl;
						break;
					}
				}
				
			}	
			
			
		}
	}	
	
	gmesh->BuildConnectivity();
		
	
//TPZGeoElBC(TPZGeoEl *el,int side,int matid, TPZGeoMesh &mesh);
//TPZGeoElBC(TPZGeoElSide &elside,int matid, TPZGeoMesh &mesh);

	
//	ofstream arg("malhaPZ1BC.txt");
//	gMesh->Print(arg);
//	
//	std::ofstream out("MyGrid.vtk");
//	TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
	
	return gmesh;
//	
//}


}

// Rectangular geometric mesh using TPZGenGrid 

//	TPZVec < int > nx(2);
//	TPZVec < REAL > corx0(2);
//	TPZVec < REAL > corx1(2);
//	int  	numlayer = 1; // Layers Numbers
//	REAL  	rotation = 0.5; // For testing purpose 
//	
//	// refinement level
//	nx[0] = 2;
//	nx[1] = 2;
//	//	x0	lower left coordinate
//	corx0[0] = 0.0;
//	corx0[1] = 0.0;	
//	//	x1	upper right coordinate 
//	corx1[0] = 50000.0;
//	corx1[1] = 50000.0;	
//	
//	TPZGenGrid geomesh(nx,corx0,corx1,numlayer);
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	geomesh.Read(*gmesh);
//	
//	// Setting BC conditions
//	geomesh.SetBC(gmesh,0,bcBottom);
//	geomesh.SetBC(gmesh,1,bcRight);
//	geomesh.SetBC(gmesh,2,bcTop);
//	geomesh.SetBC(gmesh,3,bcLeft);	
//	
//	TPZVec < REAL > PointSourceCor(3);
//	// Injection point in the center of the model
//	PointSourceCor[0]=25000.0;
//	PointSourceCor[1]=25000.0;
//	PointSourceCor[2]=0.0;	
//	geomesh.SetPointBC(gmesh,PointSourceCor, pointsource);	
//	
//	gmesh->BuildConnectivity();	

// End Rectangular geometric mesh using TPZGenGrid



// Begin Half Sapce for Flamant Problem Semi-Circular geometric mesh using tpzarc3d

//	int nodenumber = 6;
//	REAL ModelRadius = 30000.0;
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	gmesh->NodeVec().Resize(nodenumber);
//	
//	// Setting node coordantes for Arc3D 1
//	int id = 0;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,-ModelRadius );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,0.0 );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-ModelRadius);//coord Y	
//	id++;	
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,ModelRadius );//coord X
//	gmesh->NodeVec()[id].SetCoord(1,0.0);//coord Y
//	id++;	
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,-sqrt(2)*ModelRadius/2.);//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
//	id++;
//	gmesh->NodeVec()[id].SetNodeId(id);
//	gmesh->NodeVec()[id].SetCoord(0,sqrt(2)*ModelRadius/2.);//coord X
//	gmesh->NodeVec()[id].SetCoord(1,-sqrt(2)*ModelRadius/2.);//coord Y	
//	id++;
//	
//	int elementid = 0;
//	// Create Geometrical Arc #1
//	// Definition of Arc coordenates
//	TPZVec < int > nodeindex(3,0.0);
//	nodeindex[0] = 1;	
//	nodeindex[1] = 4;
//	nodeindex[2] = 2;
//	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc1 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc1,*gmesh);
//	elementid++;
//	
//	// Create Geometrical Arc #2
//	nodeindex[0] = 2;	
//	nodeindex[1] = 5;
//	nodeindex[2] = 3;			
//	TPZGeoElRefPattern < pzgeom::TPZArc3D > * elarc2 = new TPZGeoElRefPattern < pzgeom::TPZArc3D > (elementid,nodeindex, arc2,*gmesh);
//	elementid++;
//
//	// Create Geometrical triangle #1	
//	nodeindex.resize(3);
//	nodeindex[0] = 0;
//	nodeindex[1] = 1;
//	nodeindex[2] = 2;	
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > *triangle1 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical triangle #2		
//	nodeindex[0] = 0;
//	nodeindex[1] = 2;
//	nodeindex[2] = 3;		
//	TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > *triangle2 = new TPZGeoElRefPattern<TPZGeoBlend< pzgeom::TPZGeoTriangle> > (elementid,nodeindex, matId,*gmesh);
//	elementid++;
//	
//	// Create Geometrical Line #1	
//	nodeindex.resize(2);
//	nodeindex[0] = 3;
//	nodeindex[0] = 0;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoLinear > * elline1 = new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc3,*gmesh); 
//	elementid++;
//	
//	// Create Geometrical Line #2	
//	nodeindex.resize(2);
//	nodeindex[0] = 0;
//	nodeindex[0] = 1;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoLinear > * elline2 = new TPZGeoElRefPattern < pzgeom::TPZGeoLinear > (elementid,nodeindex, arc4,*gmesh); 
//	elementid++;	
//	
//	// Create Geometrical Point for fluid injection or Production #1	
//	nodeindex.resize(1);
//	nodeindex[0] = 0;	
//	TPZGeoElRefPattern < pzgeom::TPZGeoPoint > * elpoint = new TPZGeoElRefPattern < pzgeom::TPZGeoPoint > (elementid,nodeindex, WellPoint,*gmesh); 
//	elementid++;

//	gmesh->BuildConnectivity();

//	TPZVec < int > nx(2);
//	TPZVec < REAL > corx0(2);
//	TPZVec < REAL > corx1(2);
//	int  	numlayer = 1; // Layers Numbers
//	REAL  	rotation = 0.5; // For testing purpose 
//	
//	// refinement level
//	nx[0] = 2;
//	nx[1] = 2;
//	//	x0	lower left coordinate
//	corx0[0] = 0.0;
//	corx0[1] = 0.0;	
//	//	x1	upper right coordinate 
//	corx1[0] = 100000.0;
//	corx1[1] = 50000.0;	
//	
//	TPZGenGrid geomesh(nx,corx0,corx1,numlayer);
//	TPZGeoMesh * gmesh = new TPZGeoMesh;
//	geomesh.Read(*gmesh);
//	
//	// Setting BC conditions
//	geomesh.SetBC(gmesh,0,bcBottom);
//	geomesh.SetBC(gmesh,1,bcRight);
//	geomesh.SetBC(gmesh,2,bcTop);
//	geomesh.SetBC(gmesh,3,bcLeft);	
//	
//	TPZVec < REAL > PointSourceCor(3);
//	// Injection point in the center of the model
//	PointSourceCor[0]=50000.0;
//	PointSourceCor[1]=50000.0;
//	PointSourceCor[2]=0.0;	
//	geomesh.SetPointBC(gmesh,PointSourceCor, pointsource);	
//
//	gmesh->BuildConnectivity();		

//	
// End Half Sapce for Flamant Problem Semi-Circular geometric mesh using tpzarc3d	

// Use this for irregular mesh created with GID format
// Not implemented	
// End Use this for irregular mesh created with GID format	


