// -*- c++ -*-
#ifndef PZELCTEMPH
#define PZELCTEMPH

#include "pzintel.h"

template<class TGEO, class TSHAPE>
class TPZIntelGen : public TPZInterpolatedElement {

  int fConnectIndexes[TSHAPE::NSides];
  //  int fSideOrder[TSHAPE::NSides-TSHAPE::NNodes];
  int fPreferredSideOrder;
  typename TGEO::IntruleType fIntRule;

public:

  TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

virtual ~TPZIntelGen();

  virtual MElementType Type();

  virtual int NConnects() {
    return TSHAPE::NSides;
  }

  virtual void SetConnectIndex(int i, int connectindex);

  virtual int NConnectShapeF(int connect);

  virtual int Dimension() {
    return TSHAPE::Dimension;
  }

  virtual int NCornerConnects() {
    return TSHAPE::NNodes;
  }

  virtual int NSideConnects(int side);

  virtual int SideConnectLocId(int node, int side);

  virtual void SetIntegrationRule(int ord);

  /**Sets the interpolation order for the interior of the element*/
  virtual void SetInterpolationOrder(int order);

  /**Identifies the interpolation order on the interior of the element*/
  virtual void GetInterpolationOrder(TPZVec<int> &ord);

  /**return the preferred order of the polynomial along side iside*/
  virtual int PreferredSideOrder(int iside);

  /**Sets the preferred interpolation order along a side
  This method only updates the datastructure of the element
  In order to change the interpolation order of an element, use the method PRefine*/
  virtual void SetPreferredSideOrder(int side, int order);

  /**sets the interpolation order of side to order*/
  virtual void SetSideOrder(int side, int order);

  /**returns the actual interpolation order of the polynomial along the side*/
  virtual int SideOrder(int side);

  /**transform a point in the parameter space of the side into a point in the space
     of the master element*/
  //  virtual void SideParameterToElement(int side,TPZVec<REAL> &par,TPZVec<REAL> &point);

  /**transform a point in the parameter space of the master element into a point in the
     space of the side*/
  //  virtual void ElementToSideParameter(int side, TPZVec<REAL> &point, TPZVec<REAL> &par);


  /**compute the values of the shape function of the side*/
  virtual void SideShapeFunction(int side,TPZVec<REAL> &point,TPZFMatrix &phi,TPZFMatrix &dphi);

  void Shape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);

  void CreateGraphicalElement(TPZGraphMesh &grafgrid, int dimension);

virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

  /** Jorge 09/06/2001
   * Returns the transformation which transform a point from the side to the interior of the element
   */
  TPZTransform TransformSideToElement(int side);



};

#endif
