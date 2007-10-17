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

void TPZCompMeshReferred::Print(std::ostream & out){
  out << __PRETTY_FUNCTION__ << "\n";
  TPZCompMesh::Print(out);
  out << "ReferredMesh = " << this->ReferredMesh() << "\n";
}//void

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
  fReferredIndices.Resize(0);
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

void TPZCompMeshReferred::DivideReferredEl(TPZVec<TPZCompEl *> WhichRefine, TPZCompMesh * cmesh){
  TPZCompMeshReferred * me = dynamic_cast<TPZCompMeshReferred*>(cmesh);
  TPZCompMesh * other = NULL;
  if (me) other = me->ReferredMesh();

  const int nel2ref = WhichRefine.NElements();
  TPZVec<TPZCompEl *> Other2Refine(nel2ref, NULL);
  if (other){
    for(int i = 0; i < nel2ref; i++){
      TPZCompEl * cel = WhichRefine[i];
      if (!cel) continue;
      Other2Refine[i] = me->ReferredEl(cel->Index());
    }
  }

  TPZGeoMesh * gmesh = me->Reference();
  gmesh->ResetReference();
  me->ResetReferred();
  me->LoadReferences();
  TPZVec<int> filhos;
  for ( int iref = 0; iref < nel2ref; iref++ )
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

