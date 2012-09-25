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


TPZGeoMesh * GeoMeshClass::Talude()
{
    
	int numnodes;//=-1;
	int numelements;//=-1;
	
	string FileName;
	FileName = "taludemeshgid.txt";
    ifstream read (FileName.c_str());
    
    // gRefDBase.InitializeRefPatterns();
    
    
    int nodeId = 0, elementId = 0, matElId = 1;
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    read >> numnodes;
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    const int Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    for(int in=0; in<numnodes; in++)
    {
        read >> nodeId;
        read >> nodecoordX;
        read >> nodecoordY;
        read >> nodecoordZ;
        Node[nodeId-1].SetNodeId(nodeId);
        Node[nodeId-1].SetCoord(0,nodecoordX);
        Node[nodeId-1].SetCoord(1,nodecoordY);
        Node[nodeId-1].SetCoord(2,nodecoordZ);
        gMesh->NodeVec()[nodeId-1] = Node[nodeId-1];
        
    }
    
    //string FileName2;
	//FileName2 = "gidmeshelementsrefinado.txt";
    //ifstream read2 (FileName2.c_str());
    //read2 >> numelements;
    read>> numelements;
    TPZVec <int> TopolTriangle(3);
    TPZVec <int> TopolLine(2);
	
    for(int el=0; el<numelements; el++)
    {

        int topol1,topol2,topol3;
        read >> elementId;
        read >> topol1; //node 1
        read >> topol2; //node 2
        read >> topol3; //node 3
        
        TopolTriangle[0]=topol1;
        TopolTriangle[1]=topol2;
        TopolTriangle[2]=topol3;
        
        TopolTriangle[0]--;
        TopolTriangle[1]--;
        TopolTriangle[2]--;
        
        
        new TPZGeoElRefPattern< TPZGeoTriangle  > (elementId,TopolTriangle,1,*gMesh);
        
        
    }
    
    
    //LINHA Inferior
    TopolLine[0] = 903;	TopolLine[1] =0;
	new TPZGeoElRefPattern<TPZGeoLinear> (960,TopolLine,-1,*gMesh);
    TopolLine[0] = 958;	TopolLine[1] =903;
    new TPZGeoElRefPattern<TPZGeoLinear> (961,TopolLine,-2,*gMesh);
    TopolLine[0] = 0;	TopolLine[1] =216;
    new TPZGeoElRefPattern<TPZGeoLinear> (962,TopolLine,-3,*gMesh);
    
    
    
    gMesh->BuildConnectivity();
    
    ofstream arg("taludemeshoutpz.txt");
    gMesh->Print(arg);
    ofstream predio("taludemeshgid.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true);
    
    //	
	return gMesh;
}
