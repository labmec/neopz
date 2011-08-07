/**
 * @file
 * @brief Contains the implementation of the TPZGraphNode methods. 
 */
#include "pzgraphnode.h"
#include "pzgraphel.h"

using namespace std;

TPZGraphNode::TPZGraphNode(TPZConnect *cn, TPZGraphMesh *gm)
{
	fConnect = cn;
	fGraphMesh = gm;
	fGraphEl = NULL;
	if(cn) fSequenceNumber = cn->SequenceNumber();
	else fSequenceNumber = -1;
}

TPZGraphNode::TPZGraphNode()
{
	fConnect = 0;
	fGraphMesh = 0;
	fGraphEl = NULL;
	fSequenceNumber = -1;
}


TPZGraphNode::~TPZGraphNode(void)
{
}

long TPZGraphNode::FirstPoint(void)
{
	return(fPointNum);
}

void TPZGraphNode::SetElement(TPZGraphEl *gel)
{
	fGraphEl = gel;
}

void TPZGraphNode::SetPointNumber(long num)
{
	fPointNum = num;
}

void TPZGraphNode::SetConnect(TPZConnect *connect) {
	fConnect = connect;
}

void TPZGraphNode::SetGraphMesh(TPZGraphMesh *mesh) {
	fGraphMesh = mesh;
}

int TPZGraphNode::NPoints()
{
	if(fGraphEl) return(fGraphEl->NPoints(this));
	return 0;
}

void TPZGraphNode::DrawCo(TPZDrawStyle st)
{
	if(fGraphEl) fGraphEl->DrawCo(this,st);
}

void TPZGraphNode::DrawSolution(int solind, TPZDrawStyle st) {
	if(fGraphEl) fGraphEl->DrawSolution(this,solind,st);
}

void TPZGraphNode::DrawSolution(TPZVec<int> &solind, TPZDrawStyle st) {
	if(fGraphEl) fGraphEl->DrawSolution(this,solind,st);
}

void TPZGraphNode::DrawSolution(TPZBlock &bl, TPZDrawStyle st) {
	if(fGraphEl) fGraphEl->DrawSolution(this,bl,st);
}

void TPZGraphNode::Print(ostream &out) {
	out << "Connect Sequence number = " << fSequenceNumber <<
	" First Point = " << fPointNum << endl;
}
