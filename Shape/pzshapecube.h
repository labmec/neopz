// -*- c++ -*-
// $Id: pzshapecube.h,v 1.5 2004-10-06 19:12:23 phil Exp $
#ifndef SHAPECUBEHPP
#define SHAPECUBEHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
//#include "pzgeoelside.h"


#ifdef _AUTODIFF
#include "fadType.h"
#endif

/** 
 *
 * @brief Implements the shape functions of a hexahedral (3D) element

 * The range of the master element is -1,1
 * @ingroup shape
 */
class TPZShapeCube {

 public:
	 enum {NSides = 27, NNodes = 8, Dimension = 3};

/**
 * Computes the values of the shape functions and their derivatives for a hexahedral element
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

#ifdef _AUTODIFF
/**
 * Computes the values of the shape functions and their derivatives for a hexahedral element
 * These values depend on the point, the order of interpolation and ids of the corner points
 * The shapefunction computation uses the shape functions of the linear and quadrilateral element for its implementation
 * @param point (input) point where the shape functions are computed
 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
 * @param order (input) order of the side connects different from the corner connects (5 connects in this case)
 * @param phi (output) values of the shape functions and derivatives
 */
static void ShapeCube(TPZVec<REAL> &point, TPZVec<int> &id, TPZVec<int> &order, TPZVec<FADREAL> &phi);
#endif

/**
 * Computes the corner shape functions for a hexahedral element
 * @param pt (input) point where the shape function is computed
 * @param phi (output) value of the (8) shape functions
 * @param dphi (output) value of the derivatives of the (8) shape functions holding the derivatives in a column
 */
static void ShapeCorner(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);

#ifdef _AUTODIFF
/**
 * Computes the corner shape functions for a hexahedral element
 * @param pt (input) point where the shape function is computed already setup with derivatives
 * @param phi (output) value of the (8) shape functions and derivatives
 */
static void ShapeCornerCube(TPZVec<FADREAL> &pt, TPZVec<FADREAL> &phi);
#endif

/**
 * Compute the internal functions of the hexahedral shape function at a point\n
 * the internal shape functions are the shapefunctions before being multiplied by the corner
 * shape functions\n
 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation
 * determined by the transformation index
 * @param x coordinate of the point
 * @param order maximum order of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 */
static void ShapeInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                   TPZFMatrix &dphi);//,int quad_transformation_index
#ifdef _AUTODIFF
/**
 * Compute the internal functions of the hexahedral shape function at a point\n
 * the internal shape functions are the shapefunctions before being multiplied by the corner
 * shape functions\n
 * Shape3dCubeInternal is basically a call to the orthogonal shapefunction with the transformation
 * determined by the transformation index
 * @param x coordinate of the point (with derivatives already setup)
 * @param order maximum order of shape functions to be computed
 * @param phi shapefunction values (and derivatives)
 */
static void Shape3dCubeInternal(TPZVec<FADREAL> &x, int order,TPZVec<FADREAL> &phi);//,int quad_transformation_index
#endif
/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param outval coordinate of the point on the rib
 */
static void ProjectPoint3dCubeToRib(int side, TPZVec<REAL> &in, REAL &outval);

#ifdef _AUTODIFF
/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element already setup with derivatives
 * @param outval coordinate of the point on the rib
 */
static void ProjectPoint3dCubeToRib(int side, TPZVec<FADREAL> &in, FADREAL &outval);
#endif

/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param outval coordinate of the point on the rib
 */
static void ProjectPoint3dCubeSide(int side, TPZVec<REAL> &in, REAL &out);

/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinates of the point on the face
 */
static void ProjectPoint3dCubeFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);

/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinates of the point on the face
 */
static void ProjectPoint3dCubeToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);
#ifdef _AUTODIFF
/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element (with derivatives)
 * @param out coordinates of the point on the face (with derivatives)
 */
static void ProjectPoint3dCubeToFace(int face, TPZVec<FADREAL> &in, TPZVec<FADREAL> &outval);
#endif
/**
 * Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @num number of shapefunction derivatives which need to be transformed
 * @dphi values of the derivatives of the shapefunctions (modified in place)
 */
static void TransformDerivativeFromRibToCube(int rib,int num,TPZFMatrix &dphi);
#ifdef _AUTODIFF
/**
 * Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @num number of shapefunction derivatives which need to be transformed
 * @phi values of the derivatives of the shapefunctions (modified in place)
 */
//static void TransformDerivativeFromRibToCube(int rib,int num,TPZVec<FADREAL> &phi);
#endif
/**
 * Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @param num number of shapefunction derivatives which need to be transformed
 * @param dphi values of the derivatives of the shapefunctions (modified in place)
 */
static void TransformDerivativeFromFaceToCube(int rib,int num,TPZFMatrix &dphi);

#ifdef _AUTODIFF
/**
 * Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @param num number of shapefunction derivatives which need to be transformed
 * @param phi values of the shapefunctions (modified in place) with derivatives.
 */
//static void TransformDerivativeFromFaceToCube(int rib,int num,TPZVec<FADREAL> &phi);
#endif
/**
 * Returns the transformation which projects a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=26)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformElementToSide(int side);

/**
 * Returns the transformation which transform a point from the side to the interior of the element
 * @param side side from which the point will be tranformed (0<=side<=26)
 * @return TPZTransform object
 * @see the class TPZTransform
 */
static TPZTransform TransformSideToElement(int side);

/** 
 * Data structure which defines the hexahedral transformations
 */
static REAL gFaceTrans3dCube2d[6][2][3];
/** 
 * Data structure which defines the hexahedral transformations
 */
static REAL gRibTrans3dCube1d[12][3];

/** 
 * Data structure which defines the hexahedral transformations
 */
static int FaceNodes[6][4];

/** 
 * Data structure which defines the hexahedral transformations
 */
static int SideNodes[12][2];

/** 
 * Data structure which defines the hexahedral transformations
 */
static int ShapeFaceId[6][2];
/**
 * Number of connects of the element (27)
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
static int NConnectShapeF(int side, int order);

/**
 * Total number of shapefunctions, considering the order
 * of interpolation of the element
 * @param order vector of integers indicating the interpolation order of the element
 * @return number of shape functions
 */
static int NShapeF(TPZVec<int> &order);
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
  * return the number from corner nodes of the element
  */
static int NCornerNodes() { return 8;}
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
 * it returns the sides from lesser dimension associates to the side of the element
 */
static void LowerDimensionSides(int side,TPZStack<int> &smallsides);

static int SideDimension(int side);

 /**
  * returns the barycentric coordinates in the master element space of the original element
  */
 static void CenterPoint(int side, TPZVec<REAL> &center);

 /**volume of the master element*/
static REAL RefElVolume(){return 8.0;}
};

#endif
