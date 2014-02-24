/*
 *  GeoMeshClass.cpp
 *  PZ
 *
 *  Created by Diogo Cecilio on 10/6/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

/*
#include "pzmatrix.h"
#include "pzlog.h"
#include "TPZGeoCube.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"
#include "TPZRefPattern.h"
#include "pzgeopoint.h"
#include "pzmatrix.h"
#include "pzlog.h"
#include "TPZGeoCube.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"

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
//#include "tpzgeolinear.h"
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
#include "pzelast3d.h"
#include "pzgeopoint.h"
*/
#include "GeoMeshClass.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "TPZRefPatternTools.h"
#include "pzgeoelbc.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "TPZVTKGeoMesh.h"

#include "pzelasmat.h" 

using namespace pzgeom;


TPZGeoMesh * GeoMeshClass::Talude()
{
    
	long numnodes;//=-1;
	long numelements;//=-1;
	
	string FileName;
	FileName = "taludemeshgid.txt";
    ifstream read (FileName.c_str());
    
    // gRefDBase.InitializeRefPatterns();
    
    
    long nodeId = 0, elementId = 0;
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    read >> numnodes;
    
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    
    gMesh -> NodeVec().Resize(numnodes);
    
    const int Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    for(long in=0; in<numnodes; in++)
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
    TPZVec <long> TopolTriangle(3);
    TPZVec <long> TopolLine(2);
	
    for(long el=0; el<numelements; el++)
    {

        long topol1,topol2,topol3;
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


void GeoMeshClass::WellBore2d(TPZGeoMesh *gMesh)
{
    
	long numnodes;
	long numelements;
	
	string FileName;
	FileName = "../wellcil.txt";
    ifstream read (FileName.c_str());
    
    // gRefDBase.InitializeRefPatterns();
    
    
    long nodeId = 0, elementId = 0;
    
    double nodecoordX , nodecoordY , nodecoordZ ;
    read >> numnodes;
    
    
    gMesh -> NodeVec().Resize(numnodes);
    
    const int Qnodes = numnodes;
    TPZVec <TPZGeoNode> Node(Qnodes);
    
    for(long in=0; in<numnodes; in++)
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
    TPZVec <long> TopoQuad(4);
    TPZVec <long> TopoTri(3);
    TPZVec <long> TopolLine(2);
    TPZVec <long> arc(3);
	
    for(long el=0; el<numelements; el++)
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
    
    /*
    TPZGeoEl *gel = gMesh->ElementVec()[0];
    TPZGeoElBC(gel,4,-4);
    TPZGeoElBC(gel, 5,-2);
    TPZGeoElBC(gel,6,-2);
    TPZGeoElBC(gel,7,-5);
    */
    
    for (long iel=0; iel<numelements; iel++) {
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
            long node1 = gelside.SideNodeIndex(0);
            long node2 = gelside.SideNodeIndex(1);
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
	//return gMesh;
}

