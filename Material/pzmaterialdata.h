//$Id: pzmaterialdata.h,v 1.3 2007-04-16 13:50:14 tiago Exp $

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

#include "pzmanvector.h"
#include "pzfmatrix.h"

/**
This class implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods. It request to the material which attributes must be computed by the computational element and trigger their computation. Attributes are solution and its derivatives, X coordinate, etc.

@since April 10, 2007
*/

class TPZMaterialData{

public:

/** Flags indicating whether some attributes shall be computed or not */
  bool fNeedsSol, fNeedsNeighborSol, fNeedsHSize;

/** Attributes to be computed in CalcStiff */
  TPZFNMatrix<220> phi, phil, phir;
  TPZFNMatrix<660> dphix, dphixl, dphixr;
  TPZFNMatrix<9> axes, axesleft, axesright;
  TPZFNMatrix<9> jacobian, leftjac, rightjac;
  TPZFNMatrix<9> jacinv, leftjacinv, rightjacinv;
  TPZManVector<REAL,3> normal;
  TPZManVector<REAL,3> x;
  int p, leftp, rightp;
  TPZManVector<REAL,10> sol, soll, solr;
  TPZFNMatrix<30> dsol, dsoll, dsolr;
  REAL HSize;
  REAL detjac, leftdetjac, rightdetjac;

/** Class constructor */
  TPZMaterialData();

/** Class destructor */
  ~TPZMaterialData();

/** Set all flags at once */
  void SetAllRequirements(bool set);

};

#endif
