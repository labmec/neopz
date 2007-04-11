//$Id: pzmaterialdata.h,v 1.1 2007-04-11 14:26:23 tiago Exp $

#ifndef PZMATERIALDATA_H
#define PZMATERIALDATA_H

class TPZMaterial;
class TPZCompEl;
class TPZElementMatrix;
#include "pzmanvector.h"
#include "pzfmatrix.h"

/**
This class implements an interface between TPZCompEl::CalcStiff and TPZMaterial::Contribute methods. It request to the material which attributes must be computed by the computational element and trigger their computation. Attributes are solution and its derivatives, X coordinate, etc.

@since April 10, 2007
*/
class TPZMaterialData{

public:

  bool fNeedsSol, fNeedsX, fNeedsNeighborSol, fNeedsPOrder, fNeedsHSize;

private:

  TPZManVector<REAL> sol, leftsol, rightsol;
  TPZFNMatrix<100> dsol, leftdsol, rightdsol;
  TPZManVector<REAL,3> X;
  TPZFNMatrix<9> jacobian, jacinv, axes, leftaxes, rightaxes;
  TPZFNMatrix<220> phi;
  TPZFNMatrix<660> dphix;
  TPZManVector<REAL,3> normal;

  REAL Hsize, leftHsize, rightHsize;
  int POrder, leftPOrder, rightPOrder;

private:

  TPZFNMatrix<660> dphi;

private:

  TPZCompEl * fEl;

  TPZMaterial * fMat;

  void FillData(TPZVec<REAL> &qsi, REAL &weight);
  int AveragePOrder(TPZCompEl * cel);

public:

  TPZMaterialData(TPZCompEl &cel, TPZMaterial &material);

  ~TPZMaterialData();

  void Contribute(TPZVec<REAL> &qsi, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
  void ContributeInterface(TPZVec<REAL> &qsi, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);

  void SetAllRequirements(bool set){
    this->fNeedsSol = set;
    this->fNeedsX = set;
    this->fNeedsNeighborSol = set;
    this->fNeedsPOrder = set;
    this->fNeedsHSize = set;
  }

  void InitializeAttributes(TPZElementMatrix & ek, TPZElementMatrix & ef);

};

#endif
