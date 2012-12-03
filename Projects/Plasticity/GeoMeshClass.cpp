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
    int val=1000;
    REAL tol=1./val;
    int ndiv= val;
    REAL delta=(M_PI/2)/val;
    vector<int> ids,ids2,ids3,ids4;
    for(int j=0;j<=ndiv;j++)
    {
        TPZVec<REAL> co(3),co2(3);
        co[0]= 1*cos(delta*j);
        co[1]= 1*sin(delta*j);
        co2[0]=0.1 *cos(delta*j);
        co2[1]=0.1*sin(delta*j);
        for(int i=0;i<numnodes;i++)
        {
            TPZVec<REAL> cord(3);
            gMesh->NodeVec()[i].GetCoordinates(cord);
            REAL val1,val2,val3,val4;
            val1=fabs(cord[0]-co[0]);
            val2=fabs(cord[1]-co[1]);
            val3=fabs(cord[0]-co2[0]);
            val4=fabs(cord[1]-co2[1]);
            //cout << "\n coord = "<< cord<<endl;
            //cout << "co = "<< co<<endl;
            if(val1<tol && val2 < tol)
            {
                ids.push_back(i);
            }
            if(val3<tol && val4 < tol)
            {
                ids2.push_back(i);
            }
    
        }
    }
    
    tol=1./ndiv;
    delta=1./ndiv;
    for(int j=0;j<=ndiv;j++)
    {
        
        
        TPZVec<REAL> co3(3),co4(3);
        co3[0]= delta*j;
        co3[1]= 0.;
        co4[0]=0.;
        co4[1]=delta*j;
        
        for(int i=0;i<numnodes;i++)
        {
            TPZVec<REAL> cord(3);
            gMesh->NodeVec()[i].GetCoordinates(cord);
            REAL val1,val2,val3,val4;
            
            
            val1=fabs(cord[0]-co3[0]);
            val2=fabs(cord[1]-co3[1]);
            
            val3=fabs(cord[0]-co4[0]);
            val4=fabs(cord[1]-co4[1]);
            
            //cout << "\n coord = "<< cord<<endl;
            //cout << "co = "<< co3 <<endl;
            //cout << "co = "<< co4 <<endl;
            if(val1<tol && val2<tol)
            {
                ids3.push_back(i);
            }
            if(val3 < tol && val4 <tol)
            {
                ids4.push_back(i);
            }
            
        }
    }
    
    //TPZVec<int> arc(3);
    int id=0;
    //poco
    
    for(int i=0;i<(ids2.size())-1;i++)
    {
        
        TopolLine[0] =ids2[i];	TopolLine[1] = ids2[i+1];
        if(TopolLine[0] != TopolLine[1])
        {
            new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-2,*gMesh);
            cout << "\n "<<TopolLine;
        }
        
    }
    //for field

    for(int i=0;i<(ids.size())-1;i++)
    {
            TopolLine[0] =ids[i] ;	TopolLine[1] =ids[i+1];
            if(TopolLine[0] != TopolLine[1])
            {
                new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-3,*gMesh);
            }

    }
    
    
    //linha inferior
    int sz = ids3.size()-1;
    
    TopolLine[0] =ids3[0] ;	TopolLine[1] =ids3[sz];
    new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-4,*gMesh);
    
//    for(int i=0;i<(ids3.size())-1;i++)
//    {
//        cout << "\n "<<ids3[i];
//        TopolLine[0] =ids3[i+1] ;	TopolLine[1] =ids3[i];
//        if(TopolLine[0] != TopolLine[1])
//        {
//            new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-4,*gMesh);
//        }
//        
//    }

    
    
    
    sz = ids4.size()-1;
    
    TopolLine[0] =ids4[0] ;	TopolLine[1] =ids4[sz];
    new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-5,*gMesh);
//    
//    
//    for(int i=0;i<(ids4.size())-1;i++)
//    {
//        //cout << "\n "<<ids4[i];
//        TopolLine[0] =ids4[i+1] ;	TopolLine[1] =ids4[i];
//        if(TopolLine[0] != TopolLine[1])
//        {
//            new TPZGeoElRefPattern<TPZGeoLinear> (id,TopolLine,-5,*gMesh);
//        }
//        
//    }
    
    
    gMesh->BuildConnectivity();
    
    ofstream arg("wellboremeshOut.txt");
    gMesh->Print(arg);
    ofstream predio("wellboremesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gMesh,predio,true);
    
    //
	return gMesh;
}

