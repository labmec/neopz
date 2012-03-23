// -*- c++ -*-

#ifndef PZSHAPEPOINT
#define PZSHAPEPOINT

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "tpzpoint.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape{

/// Computes the shpae functions associated with a point
/**
 Compute the single shape function associated with a point
 @ingroup shape
*/
class TPZShapePoint  : public pztopology::TPZPoint  {
public:
  
  struct TMem  
  {
  };
  typedef pztopology::TPZPoint Top;

/**
 * Computes the values of the shape functions and their derivatives for a quadrilateral element
 * These values depend on the point, the order of interpolation and ids of the corner points
 * The shapefunction computation uses the shape functions of the linear element for its implementation
 * @param pt (input) point where the shape functions are computed
 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
 * @param phi (output) values of the shape functions
 * @param dphi (output) values of the derivatives of the shapefunctions

*/
  static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
		    TPZFMatrix &phi,TPZFMatrix &dphi) {
    phi(0,0) = 1.;
  }

  static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,TPZFMatrix &phi,TPZFMatrix &dphi) {
    if(side == 0) Shape(pt,id,order,phi,dphi);
  }


  /**
 * Number of shapefunctions of the connect associated with the side, considering the order
 * of interpolation of the element
 * @param side associated side
 * @order vector of integers indicating the interpolation order of the element
 * @return number of shape functions
 */
  static int NConnectShapeF(int side, int order) { return 1;}


};

};
#endif
