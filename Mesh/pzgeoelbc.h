#ifndef PZGEOELBCH
#define PZGEOELBCH

#include <iostream>

using namespace std;

class TPZGeoMesh;
class TPZGeoEl;
class TPZGeoElSide;

/*******       TPZGeoElBC       *******/

struct TPZGeoElBC {
  TPZGeoEl		*fElement;
  TPZGeoEl		*fBCElement;
  int			fSide;
  int			fId;

  TPZGeoElBC();

  TPZGeoElBC(TPZGeoEl *el,int side,int id, TPZGeoMesh &mesh);

  TPZGeoElBC(TPZGeoElSide &elside,int id, TPZGeoMesh &mesh);
	
public:
	void Print(ostream &out = cout);
};

#endif
