//$Id: pzgeoelbc.h,v 1.3 2004-04-26 14:27:03 phil Exp $

#ifndef PZGEOELBCH
#define PZGEOELBCH

#include <iostream>
#include "pzsave.h"

using namespace std;

class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;

/*******       TPZGeoElBC       *******/

struct TPZGeoElBC : public TPZSaveable {
  TPZGeoEl		*fElement;
  TPZGeoEl		*fBCElement;
  int			fSide;
  int			fId;

  TPZGeoElBC();

  TPZGeoElBC(TPZGeoEl *el,int side,int id, TPZGeoMesh &mesh);

  TPZGeoElBC(TPZGeoElSide &elside,int id, TPZGeoMesh &mesh);

  void Print(ostream &out = cout);
 
virtual int ClassId() const;

virtual void Write(TPZStream &buf, int withclassid);

virtual void Read(TPZStream &buf, void *context); 
};

#endif
