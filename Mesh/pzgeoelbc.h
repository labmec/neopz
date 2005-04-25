//$Id: pzgeoelbc.h,v 1.5 2005-04-25 02:31:48 phil Exp $

#ifndef PZGEOELBCH
#define PZGEOELBCH

#include <iostream>
#include "pzsave.h"


class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;

/*******       TPZGeoElBC       *******/

/// Associates an geometric element side with a boundary condition
/**
within the pz environment specific geometric elements represent the boundary conditions
This class simplifies the creation of these boundary elements
The constructor of the class automatically creates a copy of the object in the mesh object which is passed as parameter
@ingroup geometry
*/
struct TPZGeoElBC : public TPZSaveable {
  TPZGeoEl		*fElement;
  TPZGeoEl		*fBCElement;
  int			fSide;
  int			fId;

  TPZGeoElBC();

  TPZGeoElBC(TPZGeoEl *el,int side,int id, TPZGeoMesh &mesh);

  TPZGeoElBC(TPZGeoElSide &elside,int id, TPZGeoMesh &mesh);

  void Print(std::ostream &out = std::cout);
 
virtual int ClassId() const;

virtual void Write(TPZStream &buf, int withclassid);

virtual void Read(TPZStream &buf, void *context); 
};

#endif
