//---------------------------------------------------------------------------

#ifndef SmallerGMeshGeneratorH
#define SmallerGMeshGeneratorH

#include "pzvec.h"
//---------------------------------------------------------------------------

class TPZGeoMesh;

class SmallerGMeshGenerator
{
public:
  SmallerGMeshGenerator();
  ~SmallerGMeshGenerator();

  static TPZGeoMesh * GetGMesh(TPZVec<REAL> & XYZanchor,
                               double Lx, double Ly, double Lz, double Lcharact);
};


#endif



