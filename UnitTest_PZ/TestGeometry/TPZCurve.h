#ifndef TPZCurveH
#define TPZCurveH

/*
 *  TPZCurve.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 03/27/15.
 *  Copyright 2015 __Labmec__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

#include "TPZVTKGeoMesh.h"


#include "TPZGeoLinear.h"
#include "tpzarc3d.h"
#include "tpzgeoblend.h"
#include "tpzquadraticline.h"
#include "TPZWavyLine.h"

#include "pzgmesh.h"
#include "tpzgeoelrefpattern.h"


class TPZCurve {
  
private:  

  int fdimension;
  REAL fradius;
  TPZManVector<REAL,3> fbegin_point;
  TPZManVector<REAL,3> fend_point;
  bool fIsclosed;
  
  REAL Pi;
  
  TPZGeoMesh * fgeometricmesh;
  
public:
  
  TPZCurve();

  ~TPZCurve();

  void SetBeginPoint(TPZManVector<REAL,3> &begin_point) { fbegin_point = begin_point; }
  TPZManVector<REAL,3> GetBeginPoint() {  return fbegin_point; }

  void SetEndPoint(TPZManVector<REAL,3> &end_point) { fend_point = end_point; }
  TPZManVector<REAL,3> GetEndPoint() {  return fend_point; }

  void SetRadius(REAL r) { fradius = r; }
  REAL GetRadius() {  return fradius; }

  TPZGeoMesh * GetGeometry() {  return fgeometricmesh; }
  
  void MakeRhombus();
  void MakeCircleWave();
  void MakeCircleFromArc();
  void MakeCircleQuadratic();
  void PrintMe();
  void RefineMe(int i);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TPZManVector<REAL,3> ParametricCircle(REAL t);
  
};


#endif
