#ifndef SHAPETETRAHPP
#define SHAPETETRAHPP

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"

/** 
 *
 * @brief Implements the shape functions of a tetrahedral (3D) element
 
 * The range of the master element is 0,1
 * @ingroup shape
 */
class TPZShapeTetra {

public:
	enum {NSides = 15, NNodes = 4, Dimension = 3};
  
  
/**
 * Computes the values of the shape functions and their derivatives for a tetrahedral element
 * These values depend on the point, the order of interpolation and ids of the corner points
 * The shapefunction computation uses the shape functions of the linear and triangular element for its implementation
 * @param pt (input) point where the shape functions are computed
 * @param id (input) indexes of the corner points which determine the orientation of the shape functions
 * @param order (input) order of the side connects different from the corner connects (11 connects in this case)
 * @param phi (output) values of the shape functions
 * @param dphi (output) values of the derivatives of the shapefunctions
 */
  static void ShapeTetra(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order, TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * Computes the corner shape functions for a tetrahedral element
 * @param pt (input) point where the shape function is computed
 * @param phi (output) value of the (4) shape functions
 * @param dphi (output) value of the derivatives of the (4) shape functions holding the derivatives in a column
 */
static void CornerShapeTetraedro(TPZVec<REAL> &pt, TPZFMatrix &phi, TPZFMatrix &dphi);

/** 
 * Compute the internal functions of the tetrahedral shape function at a point\n
 * the internal shape functions are the shapefunctions before being multiplied by the corner
 * shape functions\n
 * Shape3dTetraInternal is basically a call to the orthogonal shapefunction with the transformation
 * determined by the transformation index (also on the faces)
 * @param x coordinate of the point
 * @param order maximum order of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 */
static void Shape3dTetraInternal(TPZVec<REAL> &x, int order,TPZFMatrix &phi,
				                    TPZFMatrix &dphi);

/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param outval coordinate of the point on the rib
 */
static void ProjectPoint3dTetraToRib(int rib, TPZVec<REAL> &in, REAL &outval);

/**
 * Projects a point from the interior of the element to a rib
 * @param rib rib index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param outval coordinate of the point on the rib
 */
static void ProjectPoint3dTetrSide(int side, TPZVec<REAL> &in, REAL &out);

/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinates of the point on the face
 */
static void ProjectPoint3dTetraToFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &outval);

/**
 * Projects a point from the interior of the element to a face
 * @param face face index to which the point should be projected
 * @param in coordinate of the point at the interior of the element
 * @param out coordinates of the point on the face
 */
static void ProjectPoint3dTetrFace(int face, TPZVec<REAL> &in, TPZVec<REAL> &out);

/**
 * Transforms the derivative of a shapefunction computed on the rib into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @num number of shapefunction derivatives which need to be transformed
 * @dphi values of the derivatives of the shapefunctions (modified in place)
 */
static void TransformDerivativeFromRibToTetra(int rib,int num,TPZFMatrix &dphi);

/**
 * Transforms the derivative of a shapefunction computed on the face into the three dimensional derivative
 * of the function with respect to the element. The parameter dphi should be dimensioned (3,num), at least
 * @param rib rib index along which the shapefunction is defined
 * @param num number of shapefunction derivatives which need to be transformed
 * @param dphi values of the derivatives of the shapefunctions (modified in place)
 */
static void TransformDerivativeFromFaceToTetra(int face,int num,TPZFMatrix &dphi);

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
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gVet1dTetr[6];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gFaceTrans3dTetr2d[4][2][3];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gVet2dTetr[4][2];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gRibTrans3dTetr1d[6][3];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gFaceSum3dTetra2d[4][2];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gFaceTrans3dTetra2d[4][2][3];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gRibSum3dTetra1d[6];

/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static REAL gRibTrans3dTetra1d[6][3];


/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static int FaceNodes[4][3];


/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static int SideNodes[6][2];


/** 
 * Data structure which defines the tetrahedral transformations and topology
 */
static int ShapeFaceId[4][3];

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
static REAL RefElVolume(){return (1./6.);}
};
#endif
