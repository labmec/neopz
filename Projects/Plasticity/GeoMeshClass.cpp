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
#include "pzgeoelbc.h"

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

#include "tpzarc3d.h"


TPZGeoMesh * GeoMeshClass::WellBore2d()
{
    
	int numnodes;
	int numelements;
	
	string FileName;
	FileName = "wellcil.txt";
    ifstream read (FileName.c_str());
    
    // gRefDBase.InitializeRefPatterns();
    
    
    int nodeId = 0, elementId = 0;
    
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
    
    read>> numelements;
    TPZVec <int> TopoQuad(4);
    TPZVec <int> TopoTri(3);
    TPZVec <int> TopolLine(2);
    TPZVec <int> arc(3);
	
    for(int el=0; el<numelements; el++)
    {
        
        int topol1,topol2,topol3,topol4;
        read >> elementId;
        read >> topol1; //node 1
        read >> topol2; //node 2
        read >> topol3; //node 3
        read >> topol4; //node 3
        
        TopoQuad[0]=topol1;
        TopoQuad[1]=topol2;
        TopoQuad[2]=topol3;
        TopoQuad[3]=topol4;
        
        TopoQuad[0]--;
        TopoQuad[1]--;
        TopoQuad[2]--;
        TopoQuad[3]--;
        
        
        new TPZGeoElRefPattern< TPZGeoQuad  > (elementId,TopoQuad,1,*gMesh);
        
        
    }
    
    gMesh->BuildConnectivity();
    
    for (int iel=0; iel<numelements; iel++) {
        TPZGeoEl *gel = gMesh->ElementVec()[iel];
        if (!gel) {
            DebugStop();
        }
        int ncorner = gel->NCornerNodes();
        int nsides = gel->NSides();
        for (int side=ncorner; side<nsides-1; side++) {
            TPZGeoElSide gelside(gel,side);
            if (gelside.Dimension() != 1) {
                continue;
            }
            TPZGeoElSide neighbour = gelside.Neighbour();
            if (neighbour != gelside) {
                continue;
            }
            int node1 = gelside.SideNodeIndex(0);
            int node2 = gelside.SideNodeIndex(1);
            TPZManVector<REAL,3> co1(3),co2(3);
            gMesh->NodeVec()[node1].GetCoordinates(co1);
            gMesh->NodeVec()[node2].GetCoordinates(co2);
            REAL radius1 = sqrt(co1[0]*co1[0]+co1[1]*co1[1]);
            REAL radius2 = sqrt(co2[0]*co2[0]+co2[1]*co2[1]);
            bool inner1 = fabs(radius1-0.1) < 0.02;
            bool inner2 = fabs(radius2-0.1) < 0.02;
            bool left1 = fabs(co1[0]) < 0.01;
            bool left2 = fabs(co2[0]) < 0.01;
            bool bottom1 = fabs(co1[1]) < 0.01;
            bool bottom2 = fabs(co2[1]) < 0.01;
            bool outer1 = fabs(radius1-1.) < 0.03;
            bool outer2 = fabs(radius2-1.) < 0.03;
            if (bottom1 && bottom2) {
                // material -4
                TPZGeoElBC(gel,side,-4);
                continue;
                cout << "\n "<<TopolLine;
            }
            if (left1 && left2) {
                // material -5
                TPZGeoElBC(gel,side,-5);
                continue;
            }
            if (inner1 && inner2) {
                // material -2
                TPZGeoElBC(gel,side,-2);
                continue;
            }
            if (outer1 && outer2) {
                // material -3
                TPZGeoElBC(gel,side,-3);
                continue;
            }
            DebugStop();
        }
    }
    
    
    gMesh->BuildConnectivity();
    
    ofstream arg("wellboremeshOut.txt");
    gMesh->Print(arg);
    ofstream predio("wellboremesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true);
    
    //
	return gMesh;
}

