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

  /**
   * returns the local node number of the node "node" along side "side"
   */
  static int SideNodeLocId(int side, int node) {
    return 0;
  }

 /**volume of the master element*/
  static REAL RefElVolume(){
    return 0.;
  }

 /**
  * return the number of nodes (not connectivities) associated with a side
  */
  static int NSideNodes(int side) {
    return 1;
  }
 /**
  * returns all sides whose closure contains side
  * @param side smaller dimension side
  * @param high vector which will contain all sides whose closure contain sidefrom
  */
  static void HigherDimensionSides(int side, TPZStack<int> &high) {
  }

 /**
  * return the number of nodes (not connectivities) associated with a side
  */
  static int NSideConnects(int side) {
    return 1;
  }
 /**
  * returns the local connect number of the connect "c" along side "side"
  */
  static int SideConnectLocId(int side, int c) {
    return 0;
  }
 /**
  * returns the barycentric coordinates in the master element space of the original element
  */

  static void CenterPoint(int side, TPZVec<REAL> &center) {
  }

 /**
  * returns the transformation which takes a point from the side sidefrom ot
  * the side sideto
  * @param sidefrom side where the point resides
  * @param sideto side whose closure contains sidefrom
  */
  static TPZTransform SideToSideTransform(int sidefrom, int sideto) {
    TPZTransform result(0,0);
    return result;
  }
/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=2)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
  static TPZTransform TransformSideToElement(int side) {
    TPZTransform result(0,0);
    return result;
  }

 /**
  * returns the dimension of the side
  */
  static int SideDimension(int side) {
    return 0;
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


#endif
