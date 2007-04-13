//$Id: pzinterpolationspace.h,v 1.1 2007-04-13 13:14:44 tiago Exp $

#ifndef PZINTERPOLATIONSPACE_H
#define PZINTERPOLATIONSPACE_H

#include <pzcompel.h>
class TPZMaterialData;

/**
This class implements the interfaces for TPZCompElDisc, TPZInterfaceElement and TPZInterpolatedElement.

@since April 11, 2007
*/
class TPZInterpolationSpace : public TPZCompEl
{
public:

  /**
   * Simple Constructor
   */
  TPZInterpolationSpace();

  /**
   * Simple destructor
   */
  virtual ~TPZInterpolationSpace();

  /**
   * put a copy of the element in the referred mesh
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy);

  /**
   * put a copy of the element in the patch mesh
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int,int> &gl2lcElMap);

  /**
   * copy of the element in the new mesh whit alocated index
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int &index);

  /**
   * Create a computational element within mesh
   * Inserts the element within the data structure of the mesh
   * @param mesh mesh wher will be created the element
   * @param index new elemen index
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  /**
   * it returns the shapes number of the element
   */
  virtual int NShapeF() = 0;

  /**returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int inod) = 0;

  /** Returns the max order of interpolation. */
  int MaxOrder();

  /**computes the shape function set at the point x. This method uses the order of interpolation
   * of the element along the sides to compute the number of shapefunctions
   * @param qsi point in master element coordinates
   * @param phi vector of values of shapefunctions, dimension (numshape,1)
   * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
   */
  virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix &phi,TPZFMatrix &dphi) = 0;

  /** Compute shape functions based on master element in the classical FEM manner.
   */
  virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes, 
                            REAL &detjac, TPZFMatrix &jacinv, TPZFMatrix &phi, TPZFMatrix &dphix);

  /** Initialize a material data and its attributes based on element dimension, number
   * of state variables and material definitions
   */
  void InitMaterialData(TPZMaterialData &data);

  /** Compute and fill data with requested attributes */
  void ComputeRequiredData(TPZMaterialData &data,
                           TPZVec<REAL> &qsi);

  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**returns a reference to an integration rule suitable for integrating
     the interior of the element */
  virtual TPZIntPoints &GetIntegrationRule() = 0;

  /**
   * Returns the inner radius value.
   */
  virtual REAL InnerRadius();

};

#endif
