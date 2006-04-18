//$Id: pzreferredcompel.h,v 1.1 2006-04-18 20:39:54 tiago Exp $

// -*- c++ -*-
#ifndef PZSPECIAL
#define PZSPECIAL

class TPZElementMatrix;
class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;
class TPZFMatrix;
#include "pzquad.h"

template<class TCOMPEL>
class TPZReferredCompEl : public TCOMPEL {
  public:
  
   TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  ~TPZReferredCompEl();
  
  /**
   * Compute the element stifness matrix
   * @param ek element stiffness matrix
   * @param ef element loads matrix
   */
  virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);  
  
//  static TPZCompEl *CreateFunction(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

  protected:
  
    void ComputeSolutionInOtherMesh(TPZFMatrix & phix, TPZFMatrix & dphix, TPZVec<REAL> &sol, TPZFMatrix &dsol);

};

  TPZCompEl *CreateSpecialDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
  TPZCompEl *CreateSpecialTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif