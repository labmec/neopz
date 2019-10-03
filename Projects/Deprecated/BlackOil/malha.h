//$Id: malha.h,v 1.3 2009-02-18 11:54:03 fortiago Exp $

#include "pzcmesh.h"

/** Gera malha tridimensional para problema com solucao unidimensional.
  * O problema eh injecao de agua em reservatorio saturado de oleo
  */
TPZCompMesh *Unidimensional(int h, double deltaT);

/** Gera malha tridimensional para problema com solucao unidimensional.
  * O problema eh agua em cima e oleo em baixo que devem migrar.
  */
TPZCompMesh *UnidimensionalGravidade(int h, double deltaT);

TPZGeoMesh * QuarterFiveSpot(int ndiv_xy, int ndiv_z, double lxy, double lz, double wellDiam, double factor);

TPZGeoMesh * QuarterFiveSpotReg(double lxy, double lz);

TPZCompMesh * QuarterFiveSpot(TPZGeoMesh * gmesh, double deltaT);

double AreaPoco(TPZGeoMesh * gmesh, int matid);

void DivideMalha(TPZGeoMesh * gmesh);

void DivideTornoPocos(TPZGeoMesh * gmesh);

void ScaleVec(TPZVec<REAL> &NodeIni, TPZVec<REAL> &NodeFin, double Norm, TPZVec<REAL> &OriginVec, TPZVec<REAL> &OutputVec);

