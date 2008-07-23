//$Id: pzmaterialdata.h,v 1.7 2008-07-23 21:38:26 erick Exp $

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

/** Index of the current integration point being evaluated **/
/** Needed for materials with memory **/

  int intPtIndex;

/** Class constructor */
  TPZMaterialData();

/** Copy constructor */
  TPZMaterialData( const TPZMaterialData &cp );

/** Class destructor */
  ~TPZMaterialData();

/** Set all flags at once */
  void SetAllRequirements(bool set);

  void InvertLeftRightData();

  TPZMaterialData &operator= (const TPZMaterialData &cp );

};

#endif
