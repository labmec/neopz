/**
 * @file
 * @brief Contains the implementation of the TPZGeoNode methods.
 */
//$Id: pzgnode.cpp,v 1.10 2011-03-11 13:27:49 fortiago Exp $

//METHODS DEFINITION FOR CLASS NODE
#include <stdio.h>
#include "pzgnode.h"
#include "pzerror.h"
#include "pzvec.h"
#include "pzgmesh.h"

using namespace std;
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
TPZGeoNode::TPZGeoNode(const TPZGeoNode &node) {
	fId = node.Id();
	for(int i=0;i<3;i++) fCoord[i] = node.Coord(i);
}

TPZGeoNode & TPZGeoNode::operator=(const TPZGeoNode &node){
	fId = node.Id();
	for(int i=0;i<3;i++) fCoord[i] = node.Coord(i);
	return *this;
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
void TPZGeoNode::Initialize(const TPZGeoNode &node,TPZGeoMesh &mesh) {
	fId = node.fId;
	mesh.SetNodeIdUsed(fId);
	int i;
	for(i=0;i<3;i++) fCoord[i]=node.fCoord[i];
}

void TPZGeoNode::SetCoord(const TPZVec<REAL> &x){
	const int dim = x.NElements();
#ifdef DEBUG
	if(dim > 3 || dim < 1) {
		PZError << "TPZGeoNode::SetCoord with bad parameter dim." << endl;
		DebugStop();
		return;
	}
#endif
	int i;
	for(i=0;i<dim;i++) fCoord[i] = x[i];
	for(;i<3;i++) fCoord[i] = 0.;
}

void TPZGeoNode::SetCoord(int i,REAL coord) {
	if(i > 2 || i < 0) {
		PZError << "TPZGeoNode::SetCoord with bad parameter i-th coordinate." << endl;
		return;
	}
	fCoord[i] = coord;
}

// fill the coordinates of the node
void TPZGeoNode::GetCoordinates(TPZVec<REAL> &co)
{
	for(int i=0; i<3; i++) co[i] = fCoord[i];
}

void TPZGeoNode::Print(ostream & out) {
	out << "Node : fId = " << fId;
	out << endl << "Coordinates";
	for(int i=0;i<3;i++) out << "\t" << fCoord[i];
	out << "\n";
}

// return the id of the class (used for writing reading the object)
int TPZGeoNode::ClassId() const
{
	return TPZGEONODEID;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZGeoNode,TPZGEONODEID>;
#endif

