#ifndef ELGQ2DCYL
#define ELGQ2DCYL

#include "pzreal.h"
#include "pzelgq2d.h"
#include "pzcosys.h"

/// Quadrilateral element defined in a cylindrical coordinate system
/**
THIS CLASS NEEDS TO BE REWRITTEN
*/
class TPZGeoElQ2dCyl : public TPZGeoElQ2d {

public:
  /**Constructors. Parameters: id - element id, nodeindexes - vector containing node indexes,
     matind - material index, refind - index of the node which indicates the reference direction*/
  TPZGeoElQ2dCyl(int id,TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex);
  TPZGeoElQ2dCyl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh, int cosysindex);

  //Destructor
  ~TPZGeoElQ2dCyl();
virtual  TPZGeoElQ2d* CreateGeoEl(TPZVec<int> &nodeindexes,int matind,TPZGeoMesh &mesh);


private :
  REAL fX[4][3];
  int fCosysIndex;
public:
  TPZCosys *fCosys;
  void Jacobian(TPZVec<REAL> &param,TPZFMatrix &jacobian,TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);
  void X(TPZVec<REAL> & loc,TPZVec<REAL> &result);
  void NormalVector(int side,TPZVec<REAL> &param,TPZVec<REAL> &normal,TPZFMatrix &axes,TPZFMatrix &jac1d);
  void VerifyTheta();
  void Print(ostream & out);
};

#endif
