//$Id: pzelctemp.h,v 1.15 2009-06-17 22:02:15 fortiago Exp $

// -*- c++ -*-
#ifndef PZELCTEMPH
#define PZELCTEMPH

#include "pzintel.h"
#include "pzquad.h"

/// This class implements a "generic" computational element
/**
By varying the classes passed as template arguments, the complete family of computational elements are implemented
@ingroup CompElement
*/
template<class TSHAPE>
class TPZIntelGen : public TPZInterpolatedElement {

protected:
	
  int fConnectIndexes[TSHAPE::NSides];

  //int fPreferredSideOrder;

  typename TSHAPE::IntruleType fIntRule;

public:

  TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

	TPZIntelGen(TPZCompMesh &mesh, TPZGeoEl *gel, int &index, int nocreate);

  TPZIntelGen(TPZCompMesh &mesh, const TPZIntelGen<TSHAPE> &copy);

  /**
   * used to generate patch mesh... generates a map of connect index from
   * global mesh to clone mesh
   */
  TPZIntelGen(TPZCompMesh &mesh,
              const TPZIntelGen<TSHAPE> &copy,
              std::map<int,int> & gl2lcConMap,
              std::map<int,int> & gl2lcElMap);

  TPZIntelGen();

  virtual ~TPZIntelGen();

  virtual TPZCompEl *Clone(TPZCompMesh &mesh) const {
    return new TPZIntelGen<TSHAPE> (mesh, *this);
  }

  /**
   * Create a copy of the given element. The clone copy have the connect indexes
   * mapped to the local clone connects by the given map
   * @param mesh Patch clone mesh
   * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
   * @param gl2lcElMap map the indexes of the elements between the original element and the patch element
   */
  virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,std::map<int,int> & gl2lcConMap,std::map<int,int>&gl2lcElMap) const
  {
    return new TPZIntelGen<TSHAPE> (mesh, *this, gl2lcConMap, gl2lcElMap);
  }


  virtual MElementType Type();

  virtual int NConnects() const{
    return TSHAPE::NSides;
  }

  virtual void SetConnectIndex(int i, int connectindex);

  virtual int NConnectShapeF(int connect);

  virtual int Dimension() const {
    return TSHAPE::Dimension;
  }

  virtual int NCornerConnects() {
    return TSHAPE::NCornerNodes;
  }

  virtual int NSideConnects(int side);

  virtual int SideConnectLocId(int node, int side);

  virtual int ConnectIndex(int node) const;

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
  virtual void SetPreferredOrder(int order);

  /**sets the interpolation order of side to order*/
  virtual void SetSideOrder(int side, int order);

  /**returns the actual interpolation order of the polynomial along the side*/
  virtual int SideOrder(int side);
 /**returns the actual interpolation order of the polynomial for a connect*/
	virtual int ConnectOrder(int connect);

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

  //virtual void Solution(TPZVec<REAL> &qsi,int var,TPZManVector<REAL> &sol);

  /** Jorge 09/06/2001
   * Returns the transformation which transform a point from the side to the interior of the element
   */
  TPZTransform TransformSideToElement(int side);

  virtual TPZIntPoints &GetIntegrationRule() {
    return fIntRule;
  }

  /**
  * returns the unique identifier for reading/writing objects to streams
  */
  virtual int ClassId() const;
  /**
  Save the element data to a stream
  */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  virtual void Read(TPZStream &buf, void *context);

};

#endif
