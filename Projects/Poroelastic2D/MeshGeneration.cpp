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
	
	TopolLine[0] = 1;
	TopolLine[1] = 4;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcBottom,*gmesh);
	idEl++;	
	
	TopolLine[0] = 1;
	TopolLine[1] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcRight,*gmesh);
	idEl++;
	
	TopolLine[0] = 2;
	TopolLine[1] = 5;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcRight,*gmesh);
	idEl++;	
	
	TopolLine[0] = 2;
	TopolLine[1] = 6;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcTop,*gmesh);
	idEl++;
	
	TopolLine[0] = 3;
	TopolLine[1] = 6;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcTop,*gmesh);
	idEl++;	
	
	TopolLine[0] = 3;
	TopolLine[1] = 7;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idEl,TopolLine,bcLeft,*gmesh);
	idEl++;
	
	TopolLine[0] = 0;
	TopolLine[1] = 7;
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
	int WellPoint = matIdlist[5];
	
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
	int RBright = -2;
	int RBleft = -4;
	int RBBot = -1;
	int RBTop = -3;
	
	// Left Block 
	int LBright = -2;
	int LBleft = -4;
	int LBBot = -1;
	int LBTop = -3;
	
	// Graven Block 
	int GBright = -2;
	int GBleft = -4;
	int GBBot = -1;
	int GBTop = -3;
	
	int WellLine = -5;	
	
	
	// Creating Geometrical Mesh Objects
	
	int nProdLayesr = 0.0;
	REAL ndivideleft = 3.0;
	REAL ndivideright = 3.0;
	
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
	TPZVec<int> layersIdRB(nLayers,1);
	TPZVec<int> layersIdLB(nLayers,1);
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
	int idE = 0;
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{
		TopolQuad[0] = 0+2*ilayer;
		TopolQuad[1] = 1+2*ilayer;
		TopolQuad[2] = 3+2*ilayer;
		TopolQuad[3] = 2+2*ilayer;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdRB[ilayer],*gmesh);
		idE++;
	}
	
	
	//	Bottom bondary
	TopolLine[0] = 0;
	TopolLine[1] = 1;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBBot,*gmesh);
	idE++;
	
	
	for(int inode = 0; inode < (idN - 2); inode++)
	{
		//	Left bondary
		if (inode%2==0) 
		{
			TopolLine[0] = inode;
			TopolLine[1] = inode+2;
			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBleft,*gmesh);
			idE++;		
		}
		//	Right bondary
		else 	
		{
			//			TopolLine[0] = inode;
			//			TopolLine[1] = inode+2;
			//			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBright,*gmesh);
			//			idE++;
		}
		
		
	}
	
	//	Top bondary	
	TopolLine[0] = 0+2*nLayers;
	TopolLine[1] = 1+2*nLayers;
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBTop,*gmesh);
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
	idN++;
	int lasSequencenodId = idN;	
	
	for(int ilayer = 0; ilayer < nLayers; ilayer++)
	{
		TopolQuad[0] = 0+2*nLayers+2+2*ilayer;
		TopolQuad[1] = 1+2*nLayers+2+2*ilayer;
		TopolQuad[2] = 3+2*nLayers+2+2*ilayer;
		TopolQuad[3] = 2+2*nLayers+2+2*ilayer;
		new TPZGeoElRefPattern< pzgeom::TPZGeoBlend< pzgeom::TPZGeoQuad> > (idE,TopolQuad,layersIdLB[ilayer],*gmesh);
		idE++;
	}
	
	
	//	Bottom bondary
	TopolLine[0] = 0+2*nLayers+2;
	TopolLine[1] = 1+2*nLayers+2;	
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBBot,*gmesh);
	idE++;
	
	
	for(int inode = 0+2*nLayers+2; inode < (idN - 2); inode++)
	{
		//	Left bondary
		if (inode%2!=0) 
		{
			//			TopolLine[0] = inode;
			//			TopolLine[1] = inode+2;
			//			new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBleft,*gmesh);
			//			idE++;		
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
	
	// inserting injections points
	REAL Xleft = 0.0;
	REAL Xright = 0.0;
	REAL deltaxleft = 0.0;
	REAL deltaxright = 0.0;	
	
	
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
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBBot,*gmesh);
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
	new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (idE,TopolLine,RBTop,*gmesh);
	idE++;		
	
	//---------  Interface Elements  --------------------------------------------
	if (InterfaceElement==true) 
	{
		gmesh->AddInterfaceMaterial(1,1,1);
		gmesh->AddInterfaceMaterial(1,1,1);
	}
	
	gmesh->BuildConnectivity();
	
	ofstream arg("GeologicalGravenmesh.txt");
	gmesh->Print(arg);	
	
	return gmesh;
}
