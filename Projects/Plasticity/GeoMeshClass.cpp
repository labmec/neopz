/*
 *  GeoMeshClass.cpp
 *  PZ
 *
 *  Created by Diogo Cecilio on 10/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "GeoMeshClass.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "TPZRefPatternTools.h"

using namespace pzgeom;


TPZGeoMesh * GeoMeshClass::Beangm()
{
	REAL height =0.30;
	REAL l = 2.;
	int h=2;
	int Qnodes = 4;	
	
	TPZGeoMesh *gMesh = new TPZGeoMesh;
	
	gRefDBase.InitializeRefPatterns();
	gMesh->NodeVec().Resize(Qnodes);
	
	TPZVec <int> TopolLine(2);
	TPZVec <int> TopolQuad(4);
	TPZVec <int> TopolPoint(1);
	
	TPZVec<TPZGeoNode> Node(Qnodes);
	/////ZERO////////
	Node[0].SetNodeId(0);
	Node[0].SetCoord(0 , 0.);//coord X
	Node[0].SetCoord(1 , 0.);//coord Y
	Node[0].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[0] = Node[0];
	
	Node[1].SetNodeId(1);
	Node[1].SetCoord(0 , l);//coord X
	Node[1].SetCoord(1 , 0.);//coord Y
	Node[1].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[1] = Node[1];
	
	Node[2].SetNodeId(2);
	Node[2].SetCoord(0 , l);//coord X
	Node[2].SetCoord(1 , height);//coord Y
	Node[2].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[2] = Node[2];
	
	Node[3].SetNodeId(3);
	Node[3].SetCoord(0 , 0.);//coord X
	Node[3].SetCoord(1 , height);//coord Y
	Node[3].SetCoord(2,  0.);//coord Z
	gMesh->NodeVec()[3] = Node[3];
	
	
	
	//APOIO1 ID -1
	TopolPoint[0] = 0;
	new TPZGeoElRefPattern< TPZGeoPoint  > (0,TopolPoint,-1,*gMesh);
	
	//APOIO2 ID -2
	TopolPoint[0] = 1;
	new TPZGeoElRefPattern< TPZGeoPoint  > (1,TopolPoint,-2,*gMesh);
	
	//CARGA DISTRIBUIDA ID -3
	TopolLine[0] = 0;	TopolLine[1] = 1;
	new TPZGeoElRefPattern< TPZGeoLinear  > (3,TopolLine,-3,*gMesh);
	
	TopolQuad[0] = 0;	TopolQuad[1] = 1;
	TopolQuad[2] = 2; 	TopolQuad[3] = 3;  
	new TPZGeoElRefPattern< TPZGeoQuad > (2,TopolQuad,1,*gMesh);
	
	
	
	gMesh->BuildConnectivity();
	
	for(int ref = 0; ref < h; ref++)
	{
		TPZVec<TPZGeoEl *> tatara;
		int n = gMesh->NElements();
		for(int i = 0; i < n; i++)
		{
			TPZGeoEl * gel = gMesh->ElementVec()[i];
			gel->Divide(tatara);
		}
	}
	
	
	return gMesh;
}
