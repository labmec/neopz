/**
 * @file
 * @brief Contains the implementation of the TPZGraphMesh methods. 
 */

#include "pzgraphmesh.h"
#include "pzcmesh.h"
#include "pzsubcmesh.h"
#include "pzintel.h"
#include "pzgmesh.h"
#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzgraphel.h"
#include "pztrigraph.h"
#include "pzgraphnode.h"
#include "pzmaterial.h"
#include "TPZCompElDisc.h"

#ifndef STATE_COMPLEX
	#include "TPZAgglomerateEl.h"
#endif

using namespace std;

TPZGraphMesh::TPZGraphMesh(TPZCompMesh *cm, int dimension, TPZMaterial * mat)
{
	long nel,i;
	fElementList.Resize(0);
	fElementList.CompactDataStructure(1);
	fNodeMap.Resize(0);
	fNodeMap.CompactDataStructure(1);
	fMaterial = mat;
	fCompMesh = cm;
	fDimension = dimension;
	
	TPZAdmChunkVector<TPZCompEl *> &celvec = fCompMesh->ElementVec();
	TPZCompEl *ce;
	nel = celvec.NElements();
	for(i=0;i<nel;i++) {
		ce = (TPZCompEl *) celvec[i];
		if(!ce) continue;
		ce->CreateGraphicalElement(*this, dimension);
	}
	
	fScalarNames = "";
	fVecNames = "";
}

TPZGraphMesh::~TPZGraphMesh(void)
{
	long nel = fElementList.NElements();
	TPZGraphEl *el;
	for(long i=0;i<nel;i++) {
		el = fElementList[i]; 
		if(!el) continue;
		if(el) delete el;
	}
}

static TPZGraphNode gn;

TPZGraphNode &TPZGraphMesh::FindNode(long sid)
{
	
	long nnod = fNodeMap.NElements();
	for(long index=0;index<nnod;index++) {
        TPZGraphNode *node = &fNodeMap[index];
        if(node && node->SequenceNumber() != -1 && node->SequenceNumber()==sid) return fNodeMap[index];
	}
	return gn;
}

TPZGraphEl *TPZGraphMesh::FindElement(long sid)
{
	long nelem = fElementList.NElements();
	if(sid > nelem-1) {
		cout << "TPZGraphMesh::FindNode, sid out of range sid = " << sid << endl;
		return 0;
	}
	for(long index=0;index<nelem;index++) {
        TPZGraphEl *el = fElementList[index];
        if(el && el->Id()==sid) return fElementList[index];
	}
	return NULL;
}


TPZAdmChunkVector<TPZGraphEl *> &TPZGraphMesh::ElementList() {//MElementType type
	
	TPZAdmChunkVector<TPZGraphEl *> *list = &fElementList;
	if(!list) {
        list->Resize(0);
	}
	return fElementList;
}

TPZAdmChunkVector<TPZGraphNode> &TPZGraphMesh::NodeMap() {
	return fNodeMap;
}

void TPZGraphMesh::SequenceNodes(){
	
	long nnod = fNodeMap.NElements();
	long seq = 0;
	for(long i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) {
			int numnod = n->NPoints();
			n->SetPointNumber(seq);
			seq += numnod;
		}
	}
	long el, nelem = fElementList.NElements();
	long firstel=0;
	for(el=0; el<nelem; el++) {
		TPZGraphEl *grel = fElementList[el];
		if(!grel) continue;
		grel->SetId(firstel);
		firstel += grel->NElements();
	}
}

void TPZGraphMesh::SetFileName(const std::string &filename)
{
	fFileName = filename;
	if(fOutFile.is_open())
	{
		fOutFile.close();
	}
}

void TPZGraphMesh::DrawNodes()
{
	long nnod = fNodeMap.NElements();
	for(long i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawCo(fStyle);
	}
}

void TPZGraphMesh::DrawMesh(int /*numcases*/){
}

void TPZGraphMesh::DrawSolution(int /*step*/, REAL /*time*/){
	cout << "TPZGraphMesh::DrawSolution called\n";
}

void TPZGraphMesh::DrawConnectivity(MElementType type)
{
	TPZAdmChunkVector<TPZGraphEl *> *list;
	list = (TPZAdmChunkVector<TPZGraphEl *> *) &fElementList;
	if(!list) return;
	long nel = fElementList.NElements();
	TPZGraphEl *e;
	for(long i=0;i<nel;i++) {
		e = (TPZGraphEl *) fElementList[i];
		if(e && e->Type() == type) e->Connectivity(fStyle);
	}
}

long TPZGraphMesh::NPoints() {
	long nn = 0;
	long i;
	long nnod = fNodeMap.NElements();
	for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) nn += n->NPoints();
	}
	return nn;
}

long TPZGraphMesh::NElements(MElementType type) {
	
	long numel = 0;
	long nel = fElementList.NElements();
	for(long j=0;j<nel;j++) {
		TPZGraphEl *el = (TPZGraphEl *) fElementList[j];
		if(el->Type() == type) numel += el->NElements();
	}
	return numel;
}

void TPZGraphMesh::Print(ostream &out) {
	
	long i;
	long nnod = fNodeMap.NElements();
	for(i=0;i<nnod;i++) {
		TPZGraphNode *nod = &fNodeMap[i];
		if(nod) nod->Print(out);
	}
	long nel = fElementList.NElements();
	for(i=0;i<nel;i++) {
		TPZGraphEl *el = (TPZGraphEl *) fElementList[i];
		if(el) el->Print(out);
	}
}

void TPZGraphMesh::SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames) {
	fScalarNames = scalarnames;
	fVecNames = vecnames;
}

TPZMaterial * TPZGraphMesh::Material() {
	return fMaterial;
}

void TPZGraphMesh::SetCompMesh(TPZCompMesh *mesh, TPZMaterial * &mat){
	if(fCompMesh == mesh && mat == fMaterial) return;
	long i;
	fCompMesh = mesh;
	fMaterial = mat;
	long nel = fElementList.NElements();
	TPZGraphEl *el;
	for(i=0;i<nel;i++) {
		el = fElementList[i]; 
		if(!el) continue;
		if(el) delete el;
	}
	fElementList.Resize(0);
	fNodeMap.Resize(0);
	TPZAdmChunkVector<TPZCompEl *> &celvec = fCompMesh->ElementVec();
	TPZCompEl *ce;
	nel = celvec.NElements();
	for(i=0;i<nel;i++) {
		ce = (TPZCompEl *) celvec[i];
		if(ce)  ce->CreateGraphicalElement(*this, fDimension);
	}	
}

TPZCompEl *TPZGraphMesh::FindFirstInterpolatedElement(TPZCompMesh *mesh, int dim) {
	long nel = mesh->NElements();
	TPZCompEl *cel;
	long iel;
	for(iel=0; iel<nel; iel++) {
		cel = mesh->ElementVec()[iel];
		if(!cel) continue;
		int type = cel->Type();
		if(type == EAgglomerate){
#ifndef STATE_COMPLEX
			if(!cel->Reference()) continue;
			TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
			if(agg && agg->Dimension() == dim) return agg;
#else
			DebugStop();
#endif
		}
		if(type == EDiscontinuous){
			TPZCompElDisc *disc = dynamic_cast<TPZCompElDisc *>(cel);
			if(disc && disc->Reference()->Dimension() == dim) return disc;
		}
		TPZCompEl *intel = dynamic_cast<TPZInterpolatedElement *>(cel);
		if(intel && intel->Reference()->Dimension() == dim) return intel;
		TPZSubCompMesh *subcmesh = dynamic_cast<TPZSubCompMesh *> (cel);
		if(subcmesh) {
			intel = FindFirstInterpolatedElement(subcmesh,dim);
			if(intel) return intel;
		}
	}
	return 0;
}
