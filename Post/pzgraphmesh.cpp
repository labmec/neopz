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
#include "TPZAgglomerateEl.h"

TPZGraphMesh::TPZGraphMesh(TPZCompMesh *cm, int dimension, TPZMaterial *mat)
{
  int nel,i;
  fElementList.Resize(0);
  fElementList.CompactDataStructure(1);
  fNodeMap.Resize(0);
  fNodeMap.CompactDataStructure(1);
  fMaterial = mat;
  fCompMesh = cm;
  fDimension = dimension;
  TPZGeoMesh *geomesh = fCompMesh->Reference();
  TPZAdmChunkVector<TPZGeoEl *> &gelvec = geomesh->ElementVec();
  TPZGeoEl *ge;
  nel = gelvec.NElements();
  for(i=0;i<nel;i++) {
    ge = gelvec[i];
    if(!ge || !ge->Reference() || ge->Reference()->Type() == EInterface || !ge->Reference()->IsInterpolated()) continue;
//     TPZCompEl *celtype = ge->Reference();
//     if(celtype->Type() == EAgglomerate){
//       TPZStack<TPZCompEl *> discvec;
//       dynamic_cast<TPZAgglomerateElement *>(celtype)->ListOfDiscEl(discvec);
//       int nd = discvec.NElements(),el;
//       for(el=0;el<nd;el++)
// 	discvec[el]->CreateGraphicalElement(*this, dimension);
//     } else {
      ge->Reference()->CreateGraphicalElement(*this, dimension);
//     }
  }
/*
  TPZAdmChunkVector<TPZCompEl *> &celvec = fCompMesh->ElementVec();
  TPZCompEl *ce;
  nel = celvec.NElements();
  for(i=0;i<nel;i++) {
    ce = (TPZCompEl *) celvec[i];
    if(!ce || !ce->IsInterpolated()) continue;
    ce->CreateGraphicalElement(*this, dimension);
  }
*/
  fOutFile = 0;
  fScalarNames = 0;
  fVecNames = 0;
}

TPZGraphMesh::~TPZGraphMesh(void)
{
   int nel = fElementList.NElements();
   TPZGraphEl *el;
   for(int i=0;i<nel;i++) {
   	el = fElementList[i]; 
      if(!el) continue;
      if(el) delete el;
	}
   if(fOutFile) {
   	delete fOutFile;
      fOutFile = 0;
   }
}

static TPZGraphNode gn;

TPZGraphNode &TPZGraphMesh::FindNode(int sid)
{

   int nnod = fNodeMap.NElements();
   for(int index=0;index<nnod;index++) {
        TPZGraphNode *node = &fNodeMap[index];
        if(node && node->SequenceNumber() != -1 && node->SequenceNumber()==sid) return fNodeMap[index];
   }
   return gn;
}

TPZGraphEl *TPZGraphMesh::FindElement(int sid)
{
   int nelem = fElementList.NElements();
   if(sid > nelem-1) {
   	cout << "TPZGraphMesh::FindNode, sid out of range sid = " << sid << endl;
      return 0;
   }
   for(int index=0;index<nelem;index++) {
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

   int nnod = fNodeMap.NElements();
   long seq = 0;
   for(int i=0;i<nnod;i++) {
      TPZGraphNode *n = &fNodeMap[i];
      if(n) {
           int numnod = n->NPoints();
           n->SetPointNumber(seq);
           seq += numnod;
      }
   }
   int el, nelem = fElementList.NElements();
   int firstel=0;
   for(el=0; el<nelem; el++) {
	   TPZGraphEl *grel = fElementList[el];
	   if(!grel) continue;
	   grel->SetId(firstel);
	   firstel += grel->NElements();
   }
}

void TPZGraphMesh::SetOutFile(ostream &out)
{
	fOutFile = &out;
}

ostream *TPZGraphMesh::Out()
{
	return fOutFile;
}


void TPZGraphMesh::DrawNodes()
{
	int nnod = fNodeMap.NElements();
	for(int i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) n->DrawCo(fStyle);
	}
}

void TPZGraphMesh::DrawMesh(int /*numcases*/){
}


void TPZGraphMesh::DrawSolution(char * /*var*/)
{
}

void TPZGraphMesh::DrawSolution(TPZBlock &/*bl*/)
{
}

void TPZGraphMesh::DrawSolution(int /*step*/, REAL /*time*/){
  cout << "TPZGraphMesh::DrawSolution called\n";
}



void TPZGraphMesh::DrawConnectivity(MElementType type)
{
  TPZAdmChunkVector<TPZGraphEl *> *list;
  list = (TPZAdmChunkVector<TPZGraphEl *> *) &fElementList;
  if(!list) return;
  int nel = fElementList.NElements();
  TPZGraphEl *e;
  for(int i=0;i<nel;i++) {
    e = (TPZGraphEl *) fElementList[i];
    if(e && e->Type() == type) e->Connectivity(fStyle);
  }
}

long TPZGraphMesh::NPoints() {
	long nn = 0;
   int i;
   int nnod = fNodeMap.NElements();
	for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) nn += n->NPoints();
	}
	return nn;
}

long TPZGraphMesh::NElements(MElementType type) {

  long numel = 0;
  int nel = fElementList.NElements();
  for(int j=0;j<nel;j++) {
    TPZGraphEl *el = (TPZGraphEl *) fElementList[j];
    if(el->Type() == type) numel += el->NElements();
  }
  return numel;
}

void TPZGraphMesh::Print(ostream &out) {

	int i;
   int nnod = fNodeMap.NElements();
   for(i=0;i<nnod;i++) {
      TPZGraphNode *nod = &fNodeMap[i];
      if(nod) nod->Print(out);
   }
   int nel = fElementList.NElements();
   for(i=0;i<nel;i++) {
	   TPZGraphEl *el = (TPZGraphEl *) fElementList[i];
			if(el) el->Print(out);
	}
}

void TPZGraphMesh::SetNames(TPZVec<char *>&scalarnames, TPZVec<char *>&vecnames) {
	fScalarNames = scalarnames;
	fVecNames = vecnames;
}


TPZMaterial *TPZGraphMesh::Material() {
   return fMaterial;
}

void TPZGraphMesh::SetCompMesh(TPZCompMesh *mesh, TPZMaterial *mat){
  if(fCompMesh == mesh && mat == fMaterial) return;
	int i;
  fCompMesh = mesh;
  fMaterial = mat;
  int nel = fElementList.NElements();
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

TPZCompEl *TPZGraphMesh::FindFirstInterpolatedElement(TPZCompMesh *mesh, int dim){

  int nel = mesh->NElements();
  TPZCompEl *cel;
  int iel;
  for(iel=0; iel<nel; iel++) {
    cel = mesh->ElementVec()[iel];
    if(!cel) continue;
    int type = cel->Type();
    if(type == EAgglomerate){
      if(!cel->Reference()) continue;
      TPZAgglomerateElement *agg = dynamic_cast<TPZAgglomerateElement *>(cel);
      if(agg && agg->Dimension() == dim) return agg;
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
