/*
 *  TPZCurve.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 03/27/15.
 *  Copyright 2015 __Labmec__. All rights reserved.
 *
 */


#include "TPZSurface.h"

TPZSurface::TPZSurface()
{
    fradius = 1.0;
    fdimension = 3;
    Pi = M_PI;
    fIsclosed = false;
    fgeometricmesh = NULL;
}

TPZSurface::~TPZSurface()
{
    
}

void TPZSurface::MakeCube()
{
    fgeometricmesh = new TPZGeoMesh;
    fgeometricmesh->SetDimension(fdimension);
    
    int nodes =  8;
    fgeometricmesh->SetMaxNodeId(nodes-1);
    fgeometricmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,2> TopolQuad(4);
    TPZManVector<REAL,3> coord(3,0.);
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    
    coord = ParametricSphere(Pi-cphi,Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,-Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi-cphi,-Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    coord = ParametricSphere(Pi-cphi,3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,-3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi-cphi,-3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    int id = 0;
    int matid = 1;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;

    TopolQuad[0] = 4;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 5;
    TopolQuad[3] = 1;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;
    
    TopolQuad[0] = 3;
    TopolQuad[1] = 7;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 7;
    TopolQuad[3] = 3;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (id,TopolQuad,matid,*fgeometricmesh);
    id++;
    
    
    fgeometricmesh->BuildConnectivity();
    
}

void TPZSurface::MakeRhombohedron()
{
    fgeometricmesh = new TPZGeoMesh;
    fgeometricmesh->SetDimension(fdimension);
    
    int nodes =  6;
    fgeometricmesh->SetMaxNodeId(nodes-1);
    fgeometricmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,2> TopolTriangle(3);
    TPZManVector<REAL,3> coord(3,0.);
    
    int nodeindex = 0;
    
    coord = ParametricSphere(Pi/2.0,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,Pi/2.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,Pi);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,-Pi/2.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    coord = ParametricSphere(0.0,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;

    
    
    int id = 0;
    int matid = 1;
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 1;
    TopolTriangle[1] = 2;
    TopolTriangle[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 2;
    TopolTriangle[1] = 3;
    TopolTriangle[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 3;
    TopolTriangle[1] = 0;
    TopolTriangle[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 1;
    TopolTriangle[1] = 2;
    TopolTriangle[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 2;
    TopolTriangle[1] = 3;
    TopolTriangle[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;
    
    TopolTriangle[0] = 3;
    TopolTriangle[1] = 0;
    TopolTriangle[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZGeoTriangle > (id,TopolTriangle,matid,*fgeometricmesh);
    id++;

    
    
    fgeometricmesh->BuildConnectivity();
    
}

void TPZSurface::MakeSphereFromTriangle()
{
    fgeometricmesh = new TPZGeoMesh;
    fgeometricmesh->SetDimension(fdimension);
    
    int nodes =  6;
    fgeometricmesh->SetMaxNodeId(nodes-1);
    fgeometricmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,2> TopolTriangle(3);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    int nodeindex = 0;
    
    coord = ParametricSphere(Pi/2.0,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,Pi/2.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,Pi);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi/2.0,-Pi/2.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    coord = ParametricSphere(0.0,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi,0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    
    int id = 0;
    int matid = 1;
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 4;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad1->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 1;
    TopolTriangle[1] = 2;
    TopolTriangle[2] = 4;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad2->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 2;
    TopolTriangle[1] = 3;
    TopolTriangle[2] = 4;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad3->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 3;
    TopolTriangle[1] = 0;
    TopolTriangle[2] = 4;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad4->Geom().SetData(fradius, xc);
    id++;
    
    
    TopolTriangle[0] = 0;
    TopolTriangle[1] = 1;
    TopolTriangle[2] = 5;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad5->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 1;
    TopolTriangle[1] = 2;
    TopolTriangle[2] = 5;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad6->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 2;
    TopolTriangle[1] = 3;
    TopolTriangle[2] = 5;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad7 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad7->Geom().SetData(fradius, xc);
    id++;
    
    TopolTriangle[0] = 3;
    TopolTriangle[1] = 0;
    TopolTriangle[2] = 5;
    TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > * quad8 = new TPZGeoElRefPattern< pzgeom::TPZTriangleSphere< pzgeom::TPZGeoTriangle > > (id,TopolTriangle,matid,*fgeometricmesh);
    quad8->Geom().SetData(fradius, xc);
    id++;
    
    
    
    
    fgeometricmesh->BuildConnectivity();
    
}

void TPZSurface::MakeSphereFromQuadrilateral()
{
    fgeometricmesh = new TPZGeoMesh;
    fgeometricmesh->SetDimension(fdimension);
    
    int nodes =  8;
    fgeometricmesh->SetMaxNodeId(nodes-1);
    fgeometricmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,4> Node(nodes);
    
    TPZManVector<int64_t,2> TopolQuad(4);
    TPZManVector<REAL,3> coord(3,0.);
    TPZVec<REAL> xc(3,0.);
    
    REAL cphi = atan(sqrt(2.0));
    
    int nodeindex = 0;
    
    coord = ParametricSphere(Pi-cphi,Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,-Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi-cphi,-Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    coord = ParametricSphere(Pi-cphi,3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(cphi,-3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricSphere(Pi-cphi,-3.0*Pi/4.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    
    int id = 0;
    int matid = 1;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 1;
    TopolQuad[2] = 2;
    TopolQuad[3] = 3;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad1 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad1->Geom().SetData(fradius, xc);
    id++;

    TopolQuad[0] = 4;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 7;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad2 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad2->Geom().SetData(fradius, xc);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 5;
    TopolQuad[3] = 1;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad3 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad3->Geom().SetData(fradius, xc);
    id++;
    
    TopolQuad[0] = 3;
    TopolQuad[1] = 7;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad4 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad4->Geom().SetData(fradius, xc);
    id++;
    
    TopolQuad[0] = 0;
    TopolQuad[1] = 4;
    TopolQuad[2] = 7;
    TopolQuad[3] = 3;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad5 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad5->Geom().SetData(fradius, xc);
    id++;
    
    TopolQuad[0] = 1;
    TopolQuad[1] = 5;
    TopolQuad[2] = 6;
    TopolQuad[3] = 2;
    TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > * quad6 = new TPZGeoElRefPattern< pzgeom::TPZQuadSphere< pzgeom::TPZGeoQuad > > (id,TopolQuad,matid,*fgeometricmesh);
    quad6->Geom().SetData(fradius, xc);
    id++;
    
    
    fgeometricmesh->BuildConnectivity();
    
}

//void TPZSurface::MakeSphereFromTriangle()
//{
//    fgeometricmesh = new TPZGeoMesh;
//    fgeometricmesh->SetDimension(fdimension);
//    
//    int nodes =  8;
//    fgeometricmesh->SetMaxNodeId(nodes-1);
//    fgeometricmesh->NodeVec().Resize(nodes);
//    TPZManVector<TPZGeoNode,8> Node(nodes);
//
//    fgeometricmesh->BuildConnectivity();
//    
//}



void TPZSurface::PrintMe()
{
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fgeometricmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fgeometricmesh,Dummyfile, true);
}

void TPZSurface::RefineMe(int i)
{
    int nref = i;
    TPZVec<TPZGeoEl *> sons;
    for (int iref = 0; iref < nref; iref++)
    {
        int nel = fgeometricmesh->NElements();
        for (int iel = 0; iel < nel; iel++)
        {
            TPZGeoEl *gel = fgeometricmesh->ElementVec()[iel];
            if (gel->HasSubElement()) {
                continue;
            }
            gel->Divide(sons);
        }
    }
}

TPZManVector<REAL,3> TPZSurface::ParametricCircle(REAL t)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = fradius * cos(t);
    xcoor[1] = fradius * sin(t);
    xcoor[2] = 0.0;
    return xcoor;
}

TPZManVector<REAL,3> TPZSurface::ParametricSphere(REAL phi,REAL theta)
{
    TPZManVector<REAL,3> xcoor(3,0.0);
    xcoor[0] = fradius * cos(theta) * sin(phi);
    xcoor[1] = fradius * sin(theta) * sin(phi);
    xcoor[2] = fradius * cos(phi) ;
    return xcoor;
}
