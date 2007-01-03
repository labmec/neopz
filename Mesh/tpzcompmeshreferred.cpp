//
// C++ Implementation: tpzcompmeshreferred
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzcompmeshreferred.h"
#include "pzgmesh.h"
#include "pzcompel.h"

TPZCompMeshReferred::TPZCompMeshReferred(TPZGeoMesh *gmesh)
 : TPZCompMesh(gmesh), fReferredIndices(0), fReferred(0)
{
}

TPZCompMeshReferred::TPZCompMeshReferred(const TPZCompMeshReferred &copy)
 : TPZCompMesh(*this), fReferredIndices(copy.fReferredIndices), fReferred(copy.fReferred)
{
}

TPZCompMeshReferred::~TPZCompMeshReferred()
{
}


void TPZCompMeshReferred::LoadReferred(TPZCompMesh *mesh)
{
  fReferredIndices.resize(this->NElements());
  fReferred = mesh;
  if(!mesh) return;
  TPZGeoMesh *gmesh = mesh->Reference();
  gmesh->ResetReference();
  mesh->LoadReferences();
  int iel,nel = NElements();
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
  fReferredIndices.resize(0);
  fReferred = 0;
}

TPZCompEl *TPZCompMeshReferred::ReferredEl(int index)
{
  if(!fReferred) return 0;
  int celindex = fReferredIndices[index];
  if(celindex <0) return 0;
  TPZCompEl *cel = fReferred->ElementVec()[celindex];
  return cel;
}

