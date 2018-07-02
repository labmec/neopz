/*
 *  TPZCurve.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 03/27/15.
 *  Copyright 2015 __Labmec__. All rights reserved.
 *
 */


#include "TPZCurve.h"

TPZCurve::TPZCurve()
{
  fradius = 1.0;
  fdimension = 3;
  Pi = M_PI;
  fbegin_point.Resize(3);
  fbegin_point.Fill(0.0);
  fend_point.Resize(3);
  fend_point.Fill(1.0);
  fIsclosed = false;
  fgeometricmesh = NULL; 
//   if(fend_point == fbegin_point) {fIsclosed = true;}
}

TPZCurve::~TPZCurve()
{

}

void TPZCurve::MakeRhombus()
{
  fgeometricmesh = new TPZGeoMesh;
  fgeometricmesh->SetDimension(fdimension);

  int nodes =  4;
  fgeometricmesh->SetMaxNodeId(nodes-1);
  fgeometricmesh->NodeVec().Resize(nodes);
  TPZManVector<TPZGeoNode,4> Node(nodes);

  TPZManVector<int64_t,2> TopolLine(2);
  TPZManVector<REAL,3> coord(3,0.);

  int nodeindex = 0;

  coord = ParametricCircle(0.0);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.5*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.5*Pi + Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  int id = 0;
  int matid = 1;

  TopolLine[0] = 0;
  TopolLine[1] = 1;
  new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matid,*fgeometricmesh);
  id++;

  TopolLine[0] = 1;
  TopolLine[1] = 2;
  new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matid,*fgeometricmesh);
  id++;

  TopolLine[0] = 2;
  TopolLine[1] = 3;
  new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matid,*fgeometricmesh);
  id++;

  TopolLine[0] = 3;
  TopolLine[1] = 0;
  new TPZGeoElRefPattern< pzgeom::TPZGeoLinear > (id,TopolLine,matid,*fgeometricmesh);
  id++;

  fgeometricmesh->BuildConnectivity();

}

void TPZCurve::MakeCircleFromArc()
{
  fgeometricmesh = new TPZGeoMesh;
  fgeometricmesh->SetDimension(fdimension);

  int nodes =  8;
  fgeometricmesh->SetMaxNodeId(nodes-1);
  fgeometricmesh->NodeVec().Resize(nodes);
  TPZManVector<TPZGeoNode,8> Node(nodes);

  TPZManVector<int64_t,2> TopolLine(2);
  TPZManVector<int64_t,3> TopolArc(3);
  TPZManVector<REAL,3> coord(3,0.);

  int nodeindex = 0;

  coord = ParametricCircle(0.0);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.5*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.5*Pi + Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.25*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.75*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(1.25*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++; 
  
  coord = ParametricCircle(1.75*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++; 
  
  int id = 0;
  int matid = 1;

  TopolArc[0] = 0;
  TopolArc[1] = 1;
  TopolArc[2] = 4;
  new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
  id++;

  TopolArc[0] = 1;
  TopolArc[1] = 2;
  TopolArc[2] = 5;
  new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
  id++;

  TopolArc[0] = 2;
  TopolArc[1] = 3;
  TopolArc[2] = 6;
  new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
  id++;

  TopolArc[0] = 3;
  TopolArc[1] = 0;
  TopolArc[2] = 7;
  new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
  id++;

  fgeometricmesh->BuildConnectivity();

}

void TPZCurve::MakeCircleQuadratic()
{
    fgeometricmesh = new TPZGeoMesh;
    fgeometricmesh->SetDimension(fdimension);
    
    int nodes =  8;
    fgeometricmesh->SetMaxNodeId(nodes-1);
    fgeometricmesh->NodeVec().Resize(nodes);
    TPZManVector<TPZGeoNode,8> Node(nodes);
    REAL epsilon = 0.0;
    
    TPZManVector<int64_t,2> TopolLine(2);
    TPZManVector<int64_t,3> TopolArc(3);
    TPZManVector<REAL,3> coord(3,0.);
    
    int nodeindex = 0;
    
    coord = ParametricCircle(0.0);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(0.5*Pi);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(Pi);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(0.5*Pi + Pi);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(0.25*Pi + epsilon);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(0.75*Pi +epsilon);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(1.25*Pi + epsilon);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    coord = ParametricCircle(1.75*Pi + epsilon);
    fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
    fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
    nodeindex++;
    
    int id = 0;
    int matid = 1;
    
    TopolArc[0] = 0;
    TopolArc[1] = 1;
    TopolArc[2] = 4;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine > (id,TopolArc,matid,*fgeometricmesh);
    id++;
    
    TopolArc[0] = 1;
    TopolArc[1] = 2;
    TopolArc[2] = 5;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine > (id,TopolArc,matid,*fgeometricmesh);
    id++;
    
    TopolArc[0] = 2;
    TopolArc[1] = 3;
    TopolArc[2] = 6;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine > (id,TopolArc,matid,*fgeometricmesh);
    id++;
    
    TopolArc[0] = 3;
    TopolArc[1] = 0;
    TopolArc[2] = 7;
    new TPZGeoElRefPattern< pzgeom::TPZQuadraticLine > (id,TopolArc,matid,*fgeometricmesh);
    id++;
    
    fgeometricmesh->BuildConnectivity();
    
}


void TPZCurve::MakeCircleWave()
{

  fgeometricmesh = new TPZGeoMesh;
  fgeometricmesh->SetDimension(fdimension);

  int nodes =  4;
  fgeometricmesh->SetMaxNodeId(nodes-1);
  fgeometricmesh->NodeVec().Resize(nodes);
  TPZManVector<TPZGeoNode,4> Node(nodes);

  TPZManVector<int64_t,2> TopolLine(2);
  TPZManVector<REAL,3> coord(3,0.);

  int nodeindex = 0;

// coord = ParametricCircle(0.0);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

  coord = ParametricCircle(0.0*Pi);
  fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
  fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
  nodeindex++;

//   coord = ParametricCircle(Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(0.5*Pi + Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;

  int id = 0;
  int matid = 1;
  int numwaves;
  TPZVec<REAL> wavedir(3,0.0);  

  TopolLine[0] = 0;
  TopolLine[1] = 1;
  TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wave1 = new TPZGeoElRefPattern < pzgeom::TPZWavyLine > (id,TopolLine,matid,*fgeometricmesh);
  wavedir[0] = 0.0;
  wavedir[1] = 1.0;
  numwaves = 4;
  wave1->Geom().SetData(wavedir,numwaves);
  id++;

//   TopolLine[0] = 1;
//   TopolLine[1] = 2;
//   TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wave2 = new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (id,TopolLine,matid,*fgeometricmesh);
//   wavedir[0] = -1.0;
//   wavedir[1] = -1.0;
//   numwaves = 10;
//   wave2->Geom().SetData(wavedir,numwaves);
//   id++;
// 
//   TopolLine[0] = 2;
//   TopolLine[1] = 3;
//   TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wave3 = new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (id,TopolLine,matid,*fgeometricmesh);
//   wavedir[0] = 1.0;
//   wavedir[1] = -1.0;
//   numwaves = 10;
//   wave3->Geom().SetData(wavedir,numwaves);
//   id++;
// 
//   TopolLine[0] = 3;
//   TopolLine[1] = 0;
//   TPZGeoElRefPattern< pzgeom::TPZWavyLine > * wave4 = new TPZGeoElRefPattern< pzgeom::TPZWavyLine > (id,TopolLine,matid,*fgeometricmesh);
//   wavedir[0] = 1.0;
//   wavedir[1] = 1.0;
//   numwaves = 10;
//   wave4->Geom().SetData(wavedir,numwaves);
//   id++;

  fgeometricmesh->BuildConnectivity();

  
//   fgeometricmesh = new TPZGeoMesh;
//   fgeometricmesh->SetDimension(fdimension);
// 
//   int nodes =  8;
//   fgeometricmesh->SetMaxNodeId(nodes-1);
//   fgeometricmesh->NodeVec().Resize(nodes);
//   TPZManVector<TPZGeoNode,8> Node(nodes);
// 
//   TPZManVector<int64_t,2> TopolLine(2);
//   TPZManVector<int64_t,3> TopolArc(3);
//   TPZManVector<REAL,3> coord(3,0.);
// 
//   int nodeindex = 0;
// 
//   coord = ParametricCircle(0.0);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(0.5*Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(0.5*Pi + Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(0.25*Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(0.75*Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(1.25*Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   coord = ParametricCircle(1.75*Pi);
//   fgeometricmesh->NodeVec()[nodeindex].SetCoord(coord);
//   fgeometricmesh->NodeVec()[nodeindex].SetNodeId(nodeindex);
//   nodeindex++;
// 
//   int id = 0;
//   int matid = 1;
// 
//   TopolArc[0] = 0;
//   TopolArc[1] = 1;
//   TopolArc[2] = 4;
//   new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
//   id++;
// 
//   TopolArc[0] = 1;
//   TopolArc[1] = 2;
//   TopolArc[2] = 5;
//   new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
//   id++;
// 
//   TopolArc[0] = 2;
//   TopolArc[1] = 3;
//   TopolArc[2] = 6;
//   new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
//   id++;
// 
//   TopolArc[0] = 3;
//   TopolArc[1] = 0;
//   TopolArc[2] = 7;
//   new TPZGeoElRefPattern< pzgeom::TPZArc3D > (id,TopolArc,matid,*fgeometricmesh);
//   id++;
// 
//   TopolLine[0] = 0;
//   TopolLine[1] = 1;
//   new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoLinear > > (id,TopolLine,matid,*fgeometricmesh);
//   id++;
// 
//   TopolLine[0] = 1;
//   TopolLine[1] = 2;
//   new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoLinear > > (id,TopolLine,matid,*fgeometricmesh);
//   id++;
// 
//   TopolLine[0] = 2;
//   TopolLine[1] = 3;
//   new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoLinear > > (id,TopolLine,matid,*fgeometricmesh);
//   id++;
// 
//   TopolLine[0] = 3;
//   TopolLine[1] = 0;
//   new TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoLinear > > (id,TopolLine,matid,*fgeometricmesh);
//   id++;  
// 
// //   TPZGeoElRefPattern< pzgeom::TPZGeoBlend < pzgeom::TPZGeoQuad > >
// 
//   fgeometricmesh->BuildConnectivity();

}

void TPZCurve::PrintMe()
{
    //  Print Geometrical Base Mesh
    std::ofstream argument("GeometicMesh.txt");
    fgeometricmesh->Print(argument);
    std::ofstream Dummyfile("GeometricMesh.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(fgeometricmesh,Dummyfile, true);
}

void TPZCurve::RefineMe(int i)
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

TPZManVector<REAL,3> TPZCurve::ParametricCircle(REAL t)
{
  TPZManVector<REAL,3> xcoor(3,0.0);
  xcoor[0] = fradius * cos(t);
  xcoor[1] = fradius * sin(t);
  xcoor[2] = 0.0;
  return xcoor;
}
