// -*- c++ -*-
#ifndef SHAPEPIRAMHPP
#define SHAPEPIRAMHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"


/** 
 *
 * @brief Implements the shape functions of a pyramid (3D) element

 * The range of the master element is -1,1 for the base and [0,1] for the height
 * @ingroup shape
 */
class TPZShapePiram {

public:

	enum {NSides = 19, NNodes = 5, Dimension = 3};
  /**
   * Computes the values of the shape functions and their derivatives for a pyramid element
   * These values depend on the point, the order of interpolation and ids of the corner points
   * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
   * @param pt (input) point where the shape functions are computed
   * @param id (input) indexes of the corner points which determine the orientation of the shape functions
   * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
   * @param phi (output) values of the shape functions
   * @param dphi (output) values of the derivatives of the shapefunctions   
   */
  static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi);
  static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi);
  
  /**
   * Computes the corner shape functions of the element
   * @param pt (input) point where the shape function is computed
   * @param phi (output) value of the (5) shape functions
   * @param dphi (output) value of the derivatives of the (5) shape functions holding the derivatives in a column
   */
  static void CornerShape(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);
  
  /** 
   * Compute the internal functions of the pyramid shape function at a point\n
   * the internal shape functions are the shapefunctions before being multiplied by the corner
   * shape functions\n
   * @param x coordinate of the point
   * @param order maximum order of shape functions to be computed
   * @param phi shapefunction values
   * @param dphi values of the derivatives of the shape functions
   */
  static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                   TPZFMatrix &dphi);

  /**
   * Projects a point from the interior of the element to a rib
   * @param rib rib index to which the point should be projected
   * @param in coordinate of the point at the interior of the element
   * @param outval coordinate of the point on the rib
   */
  static void ProjectPoint3dPiramToRib(int rib, TPZVec<REAL> &in, REAL &outval);

/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinates of the point on the face
 */
  static void ProjectPoint3dPiramToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);

  /**
   * Transforms a point on the face by the corresponding transformation
   * This method applies two dimensional transformations which are implemented by the classes
   * associated with triangles and quadrilaterals
   * @param transid transformation index of the face
   * @param face face to which the coordinates refer
   * @param in coordinate of the point on the face before transformation
   * @param out coordinate of the point after transformation
   */
  static void TransformPoint3dPiramFace(int transid, int face, TPZVec<REAL> &in, TPZVec<REAL> &out);

/**
 * Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @param num number of shapefunction derivatives which need to be transformed
 * @param dphi values of the derivatives of the shapefunctions (modified in place)
 */
  static void TransformDerivativeFromRibToPiram(int rib,int num,TPZFMatrix &dphi);

/**
 * Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @param num number of shapefunction derivatives which need to be transformed
 * @param dphi values of the derivatives of the shapefunctions (modified in place)
 */
  static void TransformDerivativeFromFaceToPiram(int face,int num,TPZFMatrix &dphi);

  /**
   * Transform the derivatives of the shapefunction on the shape (i.e. two dimensional derivative) to acount
   * for the transformation on the corresponding face
   * @param transid id of the transformation which needs to be applied to the face
   * @param face face to which the coordinates of the point refers
   * @param num number of derivatives of shapefunctions which need to be transformed
   * @param in matrix of derivatives of dimension (2,num)
   */
  static void TransformDerivativeFace3dPiram(int transid, int face, int num, TPZFMatrix &in);

/**
 * Returns the transformation which projects a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=20)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
  static TPZTransform TransformElementToSide(int side);

/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=20)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
  static TPZTransform TransformSideToElement(int side);

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static REAL gFaceTrans3dPiram2d[5][2][3];

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static REAL gRibTrans3dPiram1d[8][3];

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static REAL gFaceSum3dPiram2d[5][2];

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int SideNodes[8][2];

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int FaceNodes[5][4];

/** 
 * Data structure which defines the pyramid transformations and topology
 */
static int ShapeFaceId[5][4];

/**
 * Number of connects of the element (21)
 * @return number of connects of the element
 */
static int NConnects();

/**
 * Number of shapefunctions of the connect associated with the side, considering the order
 * of interpolation of the element
 * @param side associated side
 * @param order vector of integers indicating the interpolation order of the element
 * @return number of shape functions
 */
static int NConnectShapeF(int side, TPZVec<int> &order);

/**
 * Total number of shapefunctions, considering the order
 * of interpolation of the element
 * @param order vector of integers indicating the interpolation order of the element
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
static REAL RefElVolume(){return (4./3.);}

};
#endif
