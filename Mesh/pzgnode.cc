//$Id: pzgnode.cc,v 1.2 2003-11-05 16:02:21 tiago Exp $

//METHODS DEFINITION FOR CLASS NODE


#include <stdio.h>
#include "pzgnode.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgmesh.h"


TPZGeoNode::TPZGeoNode(int id,TPZVec<REAL> &coord,TPZGeoMesh &mesh) {
  mesh.SetNodeIdUsed(id);
  fId = id;
  int i,dim=coord.NElements();
  if(dim > 3) dim = 3;
  for(i=0;i<dim;i++) fCoord[i] = coord[i];
  for(;i<3;i++) fCoord[i] = 0.;
}
TPZGeoNode::TPZGeoNode() {
  fId = -1;
  for(int i=0;i<3;i++) fCoord[i] = 0.;
}
TPZGeoNode::TPZGeoNode(TPZGeoNode &node) {
  fId = node.Id();
  for(int i=0;i<3;i++) fCoord[i] = node.Coord(i);
}
void TPZGeoNode::Initialize(TPZVec<REAL> &coord,TPZGeoMesh &mesh) {
  fId = mesh.CreateUniqueNodeId();
  int i,dim = coord.NElements();
  if(dim > 3) dim = 3;
  for(i=0;i<dim;i++) fCoord[i]=coord[i];
  for(;i<3;i++) fCoord[i]=0.;
}
void TPZGeoNode::Initialize(int id,TPZVec<REAL> &coord,TPZGeoMesh &mesh) {
  fId = id;
  mesh.SetNodeIdUsed(id);
  int i,dim = coord.NElements();
  if(dim > 3) dim = 3;
  for(i=0;i<dim;i++) fCoord[i]=coord[i];
  for(;i<3;i++) fCoord[i]=0.;
}

REAL TPZGeoNode::Coord(int i) {
  if(i > 2 || i < 0) {
    PZError << "Not exist (TPZGeoNode) coordinate " << i << endl;
    return 1.e12;
  }
  return fCoord[i];
}

void TPZGeoNode::SetCoord(double *x,int dim) {
  if(dim > 3 || dim < 1) {
    PZError << "TPZGeoNode::SetCoord with bad parameter dim." << endl;
    return;
  }
  int i;
  for(i=0;i<dim;i++) fCoord[i] = x[i];
  for(;i<3;i++) fCoord[i] = 0.;
}

void TPZGeoNode::SetCoord(int i,double coord) {
  if(i > 2 || i < 0) {
    PZError << "TPZGeoNode::SetCoord with bad parameter i-th coordinate." << endl;
    return;
  }
  fCoord[i] = coord;
}

void TPZGeoNode::Print(ostream & out) {
  out << "Node : fId = " << fId;
  out << endl << "Coordinates";
  for(int i=0;i<3;i++) out << "\t" << fCoord[i];
  out << "\n";
}


