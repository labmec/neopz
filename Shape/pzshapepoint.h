// -*- c++ -*-

#ifndef PZSHAPEPOINT
#define PZSHAPEPOINT

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"

//template<class T>
//class TPZVec<REAL>;
//class TPZVec<int>;

class TPZShapePoint {
public:
  enum {NNodes = 1, NSides = 1, Dimension = 0};

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

};


#endif
