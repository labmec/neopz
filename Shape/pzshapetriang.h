// -*- c++ -*-
// $ Id: $
#ifndef SHAPETRIANGHPP
#define SHAPETRIANGHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"

/** 
 *
 * @brief Implements the shape functions of a triangular (2D) element

 * The triangular shape functions are also used in 3D elements
 * The range of the master element is 0,1
 * @ingroup shape
 */
class TPZShapeTriang {

 public:
	 enum {NSides = 7, NNodes= 3, Dimension = 2};


/**
 * Computes the values of the shape functions and their derivatives for a triangular element
 * These values depend on the point, the order of interpolation and ids of the corner points
 * The shapefunction computation uses the shape functions of the linear element for its implementation
 * @param pt (input) point where the shape functions are computed
 * @param in (input) indexes of the corner points which determine the orientation of the shape functions
 * @param order (input) order of the side connects different from the corner connects (4 connects in this case)
 * @param phi (output) values of the shape functions
 * @param dphi (output) values of the derivatives of the shapefunctions

 */
static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
             			   TPZFMatrix &phi,TPZFMatrix &dphi);
static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,
             			   TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * Computes the corner shape functions for a triangular element
 * @param pt (input) point where the shape function is computed
 * @param phi (output) value of the (3) shape functions
 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
 */
static void ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);

/**compute the internal functions of the triangle shape function at a point
 * the internal shape functions are the shapefunctions before being multiplied by the corner
 * shape functions\n
 * Shape2dTriangleInternal is basically a call to the orthogonal shapefunction with the transformation
 * determined by the transformation index
 * @param x coordinate of the point
 * @param order maximum order of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 * @param triangle_transformation_index determines the transformation applied to the internal shape
 * functions. This parameter is computed by the GetTransformId2dT method
 * @see GetTransformId2dQ
*/
static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,TPZFMatrix &dphi,int triangle_transformation_index);


/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinate of the point on the rib
 */
static void ProjectPoint2dTriangToRib(int rib, TPZVec<REAL> &in, REAL &out);

/**
 * Transforms the derivative of a shapefunction computed on the rib into the two dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (2,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @num number of shapefunction derivatives which need to be transformed
 * @dphi values of the derivatives of the shapefunctions
 */
static void TransformDerivativeFromRibToTriang(int rib,int num,TPZFMatrix &dphi);

/**
 * Method which identifies the triangular transformation based on the IDs
 * of the corner nodes
 * @param id indexes of the corner nodes
 * @return index of the transformation of the point
 */
static int GetTransformId2dT(TPZVec<int> &id);

/**
 * Transform the coordinates of the point in the space of the triangle
 * master element based on the transformation id
 * @param transid identifier of the transformation of the element as obtained by the GetTransformId2dT method
 * @param in coordinates of the variational parameter
 * @param out coordinates of the transformed parameter
 */
static void TransformPoint2dT(int transid, TPZVec<REAL> &in, TPZVec<REAL> &out);

/**
 * Transform the derivatives of num shapefunctions in place for a triangle
 * @param transid identifier of the transformation of the triangle element as obtained by the GetTransformId2dT method
 * @param num number of shapefunctions needed to transform
 * @in matrix containing the values of the derivatives of the shapefunctions as a row vector
 * the values of the derivatives contained in this matrix are modified upon return
 */
static void TransformDerivative2dT(int transid, int num, TPZFMatrix &in);

/**
 * Returns the transformation which transform a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=8)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);

/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=2)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);

/** Data structure which defines the triangle transformations*/
static REAL gTrans2dT[6][2][2];

/** Data structure which defines the triangle transformations*/
static REAL gVet2dT[6][2];
/** Data structure which defines the triangle transformations*/
static REAL gRibTrans2dT1d[3][2];
/** Data structure which defines the triangle transformations*/
static REAL gVet1dT[3];

/**
 * Number of connects of the element (7)
 * @return number of connects of the element
 */
static int NConnects();

/**
 * Number of shapefunctions of the connect associated with the side, considering the order
 * of interpolation of the element
 * @param side associated side
 * @order vector of integers indicating the interpolation order of the element
 * @return number of shape functions
 */
static int NConnectShapeF(int side, int order);

/**
 * Total number of shapefunctions, considering the order
 * of interpolation of the element
 * @order vector of integers indicating the interpolation order of the element
 * @return number of shape functions
 */
static int NShapeF(TPZVec<int> &order);

 /**
  * returns the dimension of the side
  */
static int SideDimension(int side);
 /**
  * returns the transformation which takes a point from the side sidefrom ot
  * the side sideto
  * @param sidefrom side where the point resides
  * @param sideto side whose closure contains sidefrom
  */
static TPZTransform SideToSideTransform(int sidefrom, int sideto);
 /**
  * returns all sides whose closure contains side
  * @param side smaller dimension side
  * @param high vector which will contain all sides whose closure contain sidefrom
  */
static void HigherDimensionSides(int side, TPZStack<int> &high);
 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideNodes(int side);
 /**
  * returns the local node number of the node "node" along side "side"
  */
static int SideNodeLocId(int side, int node);
 /**
  * return the number of nodes (not connectivities) associated with a side
  */
static int NSideConnects(int side);
 /**
  * returns the local connect number of the connect "c" along side "side"
  */
static int SideConnectLocId(int side, int c);

/**
 * returns the barycentric coordinates in the master element space of the original element
 */

 static void CenterPoint(int side, TPZVec<REAL> &center);

 /**volume of the master element*/
static REAL RefElVolume(){return 0.5;}
};
#endif
