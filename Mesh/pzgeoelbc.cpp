#include "pzgeoelbc.h"
#include "pzgeoel.h"
#include "pzgeoelside.h"
#include "pzgmesh.h"


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

