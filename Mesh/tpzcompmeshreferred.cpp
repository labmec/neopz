/**
 * @file
 * @brief Contains the implementation of the TPZCompMeshReferred methods.
 */

#include "tpzcompmeshreferred.h"
#include "pzgmesh.h"
#include "pzcompel.h"

void TPZCompMeshReferred::Print(std::ostream & out) const {
	out << __PRETTY_FUNCTION__ << "\n";
	TPZCompMesh::Print(out);
	out << "ReferredMesh = " << this->ReferredMesh() << "\n";
}//void

TPZCompMeshReferred::TPZCompMeshReferred()
: TPZRegisterClassId(&TPZCompMeshReferred::ClassId),
TPZCompMesh(0), fReferredIndices(0), fReferred(0)
{
}


TPZCompMeshReferred::TPZCompMeshReferred(TPZGeoMesh *gmesh)
: TPZRegisterClassId(&TPZCompMeshReferred::ClassId),
TPZCompMesh(gmesh), fReferredIndices(0), fReferred(0)
{
}

TPZCompMeshReferred::TPZCompMeshReferred(const TPZCompMeshReferred &copy)
: TPZRegisterClassId(&TPZCompMeshReferred::ClassId),
TPZCompMesh(copy), fReferredIndices(copy.fReferredIndices), fReferred(copy.fReferred)
{
}

TPZCompMeshReferred::~TPZCompMeshReferred()
{
}


void TPZCompMeshReferred::LoadReferred(TPZCompMesh *mesh)
{
	fReferredIndices.Resize(this->NElements());
	fReferred = mesh;
	if(!mesh) return;
	TPZGeoMesh *gmesh = mesh->Reference();
	gmesh->ResetReference();
	mesh->LoadReferences();
	int64_t iel,nel = NElements();
	for(iel=0; iel<nel; iel++)
	{
		TPZCompEl *cel = fElementVec[iel];
		fReferredIndices[iel] = -1;
		if (!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel)continue;
		TPZCompEl *cel2 = gel->Reference();
		if(!cel2) continue;
		fReferredIndices[iel] = cel2->Index();
	}
}

void TPZCompMeshReferred::ResetReferred()
{
	fReferredIndices.Resize(0);
	fReferred = 0;
}

TPZCompEl *TPZCompMeshReferred::ReferredEl(int64_t index)
{
	if(!fReferred) return 0;
	int64_t celindex = fReferredIndices[index];
	if(celindex <0) return 0;
	TPZCompEl *cel = fReferred->ElementVec()[celindex];
	return cel;
}

void TPZCompMeshReferred::DivideReferredEl(TPZVec<TPZCompEl *> WhichRefine, TPZCompMesh * cmesh){
	TPZCompMeshReferred * me = dynamic_cast<TPZCompMeshReferred*>(cmesh);
	TPZCompMesh * other = NULL;
	if (me) other = me->ReferredMesh();
	
	const int64_t nel2ref = WhichRefine.NElements();
	TPZVec<TPZCompEl *> Other2Refine(nel2ref, NULL);
	if (other){
		for(int64_t i = 0; i < nel2ref; i++){
			TPZCompEl * cel = WhichRefine[i];
			if (!cel) continue;
			Other2Refine[i] = me->ReferredEl(cel->Index());
		}
	}
	
	TPZGeoMesh * gmesh = me->Reference();
	gmesh->ResetReference();
	me->ResetReferred();
	me->LoadReferences();
	TPZVec<int64_t> filhos;
	for ( int64_t iref = 0; iref < nel2ref; iref++ )
	{
		TPZCompEl * cel = WhichRefine[iref];
		if (!cel) continue;
		cel->Divide ( cel->Index(), filhos );
	}
	me->ExpandSolution();
	
	if (other){
		TPZCompMeshReferred::DivideReferredEl(Other2Refine, other);
		me->LoadReferred(other);
	}
	
}

/** @brief Returns the unique identifier for reading/writing objects to streams */
int TPZCompMeshReferred::ClassId() const{
    return Hash("TPZCompMeshReferred") ^ TPZCompMesh::ClassId() << 1;
}
/** @brief Save the element data to a stream */
void TPZCompMeshReferred::Write(TPZStream &buf, int withclassid) const {
    TPZCompMesh::Write(buf, withclassid);
    buf.Write(fReferredIndices);
    TPZPersistenceManager::WritePointer(fReferred, &buf);
}

/** @brief Read the element data from a stream */
void TPZCompMeshReferred::Read(TPZStream &buf, void *context) {
    TPZCompMesh::Read(buf, context);
    buf.Read(fReferredIndices);
    fReferred = dynamic_cast<TPZCompMesh *>(TPZPersistenceManager::GetInstance(&buf));
}

template class TPZRestoreClass<TPZCompMeshReferred>;

