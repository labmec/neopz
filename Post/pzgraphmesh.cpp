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
#include "TPZMaterial.h"
#include "TPZCompElDisc.h"

#ifndef STATE_COMPLEX
	#include "TPZAgglomerateEl.h"
#endif

using namespace std;

TPZGraphMesh::TPZGraphMesh(TPZCompMesh *cm, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames) :
    fCompMesh(cm), fDimension(dimension), fMaterialIds(matids), fScalarNames(scalarnames), fVecNames(vecnames), fTensorNames()
{
	int64_t nel,i;
	fElementList.Resize(0);
	fElementList.CompactDataStructure(fElementList.NOW);
	fNodeMap.Resize(0);
	fNodeMap.CompactDataStructure(fNodeMap.NOW);
	
	TPZAdmChunkVector<TPZCompEl *> &celvec = fCompMesh->ElementVec();
	TPZCompEl *ce;
	nel = celvec.NElements();
	for(i=0;i<nel;i++) {
		ce = (TPZCompEl *) celvec[i];
		if(!ce) continue;
		ce->CreateGraphicalElement(*this, dimension);
	}
	
}

TPZGraphMesh::TPZGraphMesh(TPZCompMesh *cm, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensornames) :
fCompMesh(cm), fDimension(dimension), fMaterialIds(matids), fScalarNames(scalarnames), fVecNames(vecnames), fTensorNames(tensornames)
{
    int64_t nel,i;
    fElementList.Resize(0);
    fElementList.CompactDataStructure(fElementList.NOW);
    fNodeMap.Resize(0);
    fNodeMap.CompactDataStructure(fNodeMap.NOW);
    
    TPZAdmChunkVector<TPZCompEl *> &celvec = fCompMesh->ElementVec();
    TPZCompEl *ce;
    nel = celvec.NElements();
    for(i=0;i<nel;i++) {
        ce = (TPZCompEl *) celvec[i];
        if(!ce) continue;
        ce->CreateGraphicalElement(*this, dimension);
    }
    
}

TPZGraphMesh::~TPZGraphMesh(void)
{
	int64_t nel = fElementList.NElements();
	TPZGraphEl *el;
	for(int64_t i=0;i<nel;i++) {
		el = fElementList[i]; 
		if(!el) continue;
		if(el) delete el;
	}
}

static TPZGraphNode gn;

TPZGraphNode &TPZGraphMesh::FindNode(int64_t sid)
{
	
	int64_t nnod = fNodeMap.NElements();
	for(int64_t index=0;index<nnod;index++) {
        TPZGraphNode *node = &fNodeMap[index];
        if(node && node->SequenceNumber() != -1 && node->SequenceNumber()==sid) return fNodeMap[index];
	}
	return gn;
}

TPZGraphEl *TPZGraphMesh::FindElement(int64_t sid)
{
	int64_t nelem = fElementList.NElements();
	if(sid > nelem-1) {
		cout << "TPZGraphMesh::FindNode, sid out of range sid = " << sid << endl;
		return 0;
	}
	for(int64_t index=0;index<nelem;index++) {
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
	
	int64_t nnod = fNodeMap.NElements();
	int64_t seq = 0;
	for(int64_t i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) {
			int numnod = n->NPoints();
			n->SetPointNumber(seq);
			seq += numnod;
		}
	}
	int64_t el, nelem = fElementList.NElements();
	int64_t firstel=0;
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
	int64_t nnod = fNodeMap.NElements();
	for(int64_t i=0;i<nnod;i++) {
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
	int64_t nel = fElementList.NElements();
	TPZGraphEl *e;
	for(int64_t i=0;i<nel;i++) {
		e = (TPZGraphEl *) fElementList[i];
		if(e && e->Type() == type) e->Connectivity(fStyle);
	}
}

int64_t TPZGraphMesh::NPoints() {
	int64_t nn = 0;
	int64_t i;
	int64_t nnod = fNodeMap.NElements();
	for(i=0;i<nnod;i++) {
		TPZGraphNode *n = &fNodeMap[i];
		if(n) nn += n->NPoints();
	}
	return nn;
}

int64_t TPZGraphMesh::NElements(MElementType type) {
	
	int64_t numel = 0;
	int64_t nel = fElementList.NElements();
	for(int64_t j=0;j<nel;j++) {
		TPZGraphEl *el = (TPZGraphEl *) fElementList[j];
		if(el->Type() == type) numel += el->NElements();
	}
	return numel;
}

void TPZGraphMesh::Print(ostream &out) {
	
	int64_t i;
	int64_t nnod = fNodeMap.NElements();
	for(i=0;i<nnod;i++) {
		TPZGraphNode *nod = &fNodeMap[i];
		if(nod) nod->Print(out);
	}
	int64_t nel = fElementList.NElements();
	for(i=0;i<nel;i++) {
		TPZGraphEl *el = (TPZGraphEl *) fElementList[i];
		if(el) el->Print(out);
	}
}

void TPZGraphMesh::SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames) {
	fScalarNames = scalarnames;
	fVecNames = vecnames;
    fTensorNames.resize(0);
}

void TPZGraphMesh::SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames, const TPZVec<std::string>& tensornames) {
    fScalarNames = scalarnames;
    fVecNames = vecnames;
    fTensorNames = tensornames;
}


void TPZGraphMesh::SetMaterialIds(const std::set<int> & matids){
    SetCompMesh(fCompMesh, matids);
}

std::set<int> TPZGraphMesh::MaterialIds(){
    return fMaterialIds;
}

bool TPZGraphMesh::Material_Is_PostProcessed(int matid){
    return fMaterialIds.find(matid) != fMaterialIds.end();
}

void TPZGraphMesh::SetCompMesh(TPZCompMesh *mesh, const std::set<int> & matids){
	if(fCompMesh == mesh && matids == fMaterialIds) return;
	int64_t i;
	fCompMesh = mesh;
	fMaterialIds = matids;
	int64_t nel = fElementList.NElements();
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
	int64_t nel = mesh->NElements();
	TPZCompEl *cel;
	int64_t iel;
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

int TPZGraphMesh::ClassId() const {
    return Hash("TPZGraphMesh");
}

void TPZGraphMesh::Read(TPZStream& buf, void* context) {
    fCompMesh = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
    fGeoMesh = dynamic_cast<TPZGeoMesh *>(TPZPersistenceManager::GetInstance(&buf));
    TPZManVector<int> mat_ids;
    buf.Read(mat_ids);
    for (auto matid: mat_ids) {
        fMaterialIds.insert(matid);
    }
    buf.Write(&fDimension);
    buf.Read(&fDimension);
    buf.ReadPointers(fElementList);
    fNodeMap.Read(buf, context);
    buf.Read(&fResolution);
    int fStyleInt;
    buf.Read(&fStyleInt);
    fStyle = TPZDrawStyle(fStyleInt);
    buf.Read(&fFileName);
    this->SetFileName(fFileName); ///Forcing to close the previously open file, if any.
    buf.Read(fScalarNames);
    buf.Read(fVecNames);
    buf.Read(fTensorNames);
}

void TPZGraphMesh::Write(TPZStream& buf, int withclassid) const {
    TPZPersistenceManager::WritePointer(fCompMesh, &buf);
    TPZPersistenceManager::WritePointer(fGeoMesh, &buf);
    TPZManVector<int> mat_ids(fMaterialIds.size());
    auto it = fMaterialIds.begin();
    for (int i = 0; i < fMaterialIds.size(); i++, it++ ) {
        mat_ids[i] = *it;
    }
    buf.Write(mat_ids);
    buf.Write(&fDimension);
    buf.WritePointers(fElementList);
    fNodeMap.Write(buf, withclassid);
    buf.Write(&fResolution);
    int fStyleInt = as_integer(fStyle);
    buf.Write(&fStyleInt);
    buf.Write(&fFileName);
    buf.Write(fScalarNames);
    buf.Write(fVecNames);
    buf.Write(fTensorNames);
}
