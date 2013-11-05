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
: TPZCompMesh(0), fReferredIndices(0), fReferred(0)
{
}


TPZCompMeshReferred::TPZCompMeshReferred(TPZGeoMesh *gmesh)
: TPZCompMesh(gmesh), fReferredIndices(0), fReferred(0)
{
}

TPZCompMeshReferred::TPZCompMeshReferred(const TPZCompMeshReferred &copy)
: TPZCompMesh(copy), fReferredIndices(copy.fReferredIndices), fReferred(copy.fReferred)
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
	long iel,nel = NElements();
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

TPZCompEl *TPZCompMeshReferred::ReferredEl(long index)
{
	if(!fReferred) return 0;
	long celindex = fReferredIndices[index];
	if(celindex <0) return 0;
	TPZCompEl *cel = fReferred->ElementVec()[celindex];
	return cel;
}

void TPZCompMeshReferred::DivideReferredEl(TPZVec<TPZCompEl *> WhichRefine, TPZCompMesh * cmesh){
	TPZCompMeshReferred * me = dynamic_cast<TPZCompMeshReferred*>(cmesh);
	TPZCompMesh * other = NULL;
	if (me) other = me->ReferredMesh();
	
	const long nel2ref = WhichRefine.NElements();
	TPZVec<TPZCompEl *> Other2Refine(nel2ref, NULL);
	if (other){
		for(long i = 0; i < nel2ref; i++){
			TPZCompEl * cel = WhichRefine[i];
			if (!cel) continue;
			Other2Refine[i] = me->ReferredEl(cel->Index());
		}
	}
	
	TPZGeoMesh * gmesh = me->Reference();
	gmesh->ResetReference();
	me->ResetReferred();
	me->LoadReferences();
	TPZVec<long> filhos;
	for ( long iref = 0; iref < nel2ref; iref++ )
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
int TPZCompMeshReferred::ClassId() const
{
    return TPZCOMPMESHREFERREDID;
}
/** @brief Save the element data to a stream */
void TPZCompMeshReferred::Write(TPZStream &buf, int withclassid)
{
    TPZCompMesh::Write(buf, withclassid);
    TPZSaveable::WriteObjects(buf, this->fReferredIndices);
}

/** @brief Read the element data from a stream */
void TPZCompMeshReferred::Read(TPZStream &buf, void *context)
{
    fReferred = (TPZCompMesh *) context;
    context = fReferred->Reference();
    TPZCompMesh::Read(buf, context);
    TPZSaveable::ReadObjects(buf, this->fReferredIndices);
}

template class TPZRestoreClass<TPZCompMeshReferred,TPZCOMPMESHREFERREDID>;

