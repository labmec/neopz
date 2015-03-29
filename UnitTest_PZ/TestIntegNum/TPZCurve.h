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


class TPZCurve {
  
private:  

  int fdimension;
  STATE fradius;
  TPZManVector<STATE,3> fbegin_point;
  TPZManVector<STATE,3> fend_point;
  bool fIsclosed;
  
  STATE Pi;
  
  TPZGeoMesh * fgeometricmesh;
  
public:
  
  TPZCurve();

  ~TPZCurve();

  void SetBeginPoint(TPZManVector<STATE,3> &begin_point) { fbegin_point = begin_point; }
  TPZManVector<STATE,3> GetBeginPoint() {  return fbegin_point; }

  void SetEndPoint(TPZManVector<STATE,3> &end_point) { fend_point = end_point; }
  TPZManVector<STATE,3> GetEndPoint() {  return fend_point; }

  void SetRadius(STATE r) { fradius = r; }
  STATE GetRadius() {  return fradius; }

  TPZGeoMesh * GetGeometry() {  return fgeometricmesh; }
  
  void MakeRhombus();
  void MakeCircleWave();
  void MakeCircleFromArc();
  void PrintMe();
  void RefineMe(int i);

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  TPZManVector<STATE,3> ParametricCircle(STATE t);
  
};


#endif