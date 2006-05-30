//$Id: pzreferredcompel.h,v 1.2 2006-05-30 17:51:24 tiago Exp $

// -*- c++ -*-
#ifndef PZSPECIAL
#define PZSPECIAL

class TPZElementMatrix;
class TPZCompEl;
class TPZInterpolatedElement;
class TPZCompElDisc;
class TPZGeoEl;
class TPZCompMesh;
class TPZFMatrix;
#include "pzvec.h"
#include "pzmanvector.h"

template<class TCOMPEL>
class TPZReferredCompEl : public TCOMPEL {
  public:
  
   TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  ~TPZReferredCompEl();
  
 /**
  * Computes solution and its derivatives in local coordinate qsi
  * @param qsi master element coordinate
  * @param phi matrix containing shape functions compute in qsi point
  * @param dphix matrix containing the derivatives of shape functions with respect of global coordinates: D[phi,x], D[phi,y], D[phi,z]
  * @param sol finite element solution
  * @param dsol solution derivatives
  */
  virtual void ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix, TPZVec<REAL> &sol, TPZFMatrix &dsol);
  
};

  TPZCompEl *CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif
