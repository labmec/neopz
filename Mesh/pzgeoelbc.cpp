//$Id: pzgeoelbc.cpp,v 1.3 2004-04-26 14:27:03 phil Exp $

#include "pzgeoelbc.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgmesh.h"
#include "pzmeshid.h"
#include "pzstream.h"


TPZGeoElBC::TPZGeoElBC() {
    fElement = 0;
    fBCElement = 0;
    fSide = -1;
    fId = -1;
}

TPZGeoElBC::TPZGeoElBC(TPZGeoEl *el,int side,int id, TPZGeoMesh &mesh) {
    fElement = el;
    fBCElement = 0;
    fSide = side;
    fId = id;
    int index = mesh.BCElementVec().AllocateNewElement();
    mesh.BCElementVec()[index] = *this;
}

TPZGeoElBC::TPZGeoElBC(TPZGeoElSide &elside,int id, TPZGeoMesh &mesh) {
    fElement = elside.Element();
    fBCElement = 0;
    fSide = elside.Side();
    fId = id;
    int index = mesh.BCElementVec().AllocateNewElement();
    mesh.BCElementVec()[index] = *this;
}

void TPZGeoElBC::Print(ostream &out)
{
	out << "TPZGeoElBC\n";
	if(fElement) out << "Id of the associated element " << fElement->Id() << endl;
	else out << "No Element associated\n";
	out << "Side of the element " << fSide << endl;
	out << "Id of the boundary condition " << fId << endl;
	if(fBCElement) out << "Id of element associated with the boundary condition " << fBCElement->Id();
	else out << "Associated boundary element not created yet\n";

}

int TPZGeoElBC::ClassId() const {
  return TPZGEOELBCID;
}

void TPZGeoElBC::Write(TPZStream &buf, int withclassid) {
  TPZSaveable::Write(buf,withclassid);
  int bcelindex = -1;
  int elemindex = -1;
  if(fBCElement) bcelindex = fBCElement->Index();
  if(fElement) elemindex = fElement->Index();
  buf.Write(&bcelindex,1);
  buf.Write(&elemindex,1);
  buf.Write(&fId,1);
  buf.Write(&fSide,1);
}

void TPZGeoElBC::Read(TPZStream &buf, void *context) {
  TPZSaveable::Read(buf,context);  
  TPZGeoMesh *gmesh = (TPZGeoMesh *) context;
  int bcelindex = -1;
  int elemindex = -1;
  buf.Read(&bcelindex,1);
  buf.Read(&elemindex,1);
  buf.Read(&fId,1);
  buf.Read(&fSide,1);
  if(bcelindex != -1) fBCElement = gmesh->ElementVec()[bcelindex];
  if(elemindex != -1) fElement = gmesh->ElementVec()[elemindex];
}
