// -*- c++ -*-
// $ Id: $
#ifndef SHAPELINEARHPP
#define SHAPELINEARHPP

//#include "..\..\MULTIGRID\pztrnsform.h"	// Added by ClassView
#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/**
 *
 * @brief Implements the shape functions of a linear (1D) element

 * The linear shape functions form the basis of all other shape function computations
 * The range of the master element is -1,1
 * The orthogonal function which generates the linear shape functions can be modified
 * by changing the function pointer fOrthogonal
 * all static tables and functions concerning one-d elements will be grouped in this class
 * @ingroup shape
 */
class TPZShapeLinear {

public:
	enum {NNodes = 2, NSides = 3, Dimension = 1};
/**
 *	pointer to function which returns num orthogonal functions at the point x
 * @param x coordinate of the point
 * @param num number of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 */
static void (*fOrthogonal)(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi);


#ifdef _AUTODIFF
/**
 *	pointer to function which returns num orthogonal functions at the point x
 * @param x coordinate of the point with derivatives already setup
 * @param num number of shape functions to be computed
 * @param phi shapefunction values with derivatives
 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
 */
static void (*FADfOrthogonal)(FADREAL&,int ,TPZVec<FADREAL> &);
#endif
/**
 * Chebyshev orthogonal function, computes num orthogonal functions at the point x
 * @param x coordinate of the point
 * @param num number of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 */
static void Chebyshev(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi);

/**
 * Legendre orthogonal function, computes num orthogonal functions at the point x
 * @param x coordinate of the point
 * @param num number of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 */
static void Legendre(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi);

static REAL fJacobiAlfa;
static REAL fJacobiBeta;

static void Jacobi(REAL x,int num,TPZFMatrix & phi,TPZFMatrix & dphi);

#ifdef _AUTODIFF
/**
 * Chebyshev orthogonal function, computes num orthogonal functions at the point x
 * @param x coordinate of the point (with derivative already setup)
 * @param num number of shape functions to be computed
 * @param phi shapefunction values and derivatives
 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
 */
static void Chebyshev(FADREAL & x,int num,TPZVec<FADREAL> &phi);
#endif
 public:

    //TPZShapeLinear();
     //~TPZShapeLinear() {};

/**
 * Functions which computes the shapefunctions of a one-d element
 * The orientation of the shapefunctions depend on the order of the id parameters
 * if id[0] < id[1] the shapefunctions are unchanged
 * if id[0] > id[1] the odd ordered shapefunctions are inverted
 * @param x coordinate of the point
 * @param order order of the shape functions to be computed 0<= order
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 * @param id determines the orientation of the shape functions
 */
static void Shape(TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,TPZFMatrix &phi,TPZFMatrix &dphi);

 static void SideShape(int side, TPZVec<REAL> &pt, TPZVec<int> &id, TPZVec<int> &order,TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * Computes the values of the orthogonal shapefunctions before multiplying them by the
 * corner shapefunctions
 * The shape1dInternal function is extensively used by the shapefunction computation of
 * the other elements
 * @param x coordinate of the point
 * @param num number of shape functions to be computed
 * @param phi shapefunction values
 * @param dphi values of the derivatives of the shape functions
 * @param transformation_index determines the transformation applied to the internal shape
 * functions. This parameter is computed by the GetTransformId1d method
 * @see GetTransformId1d
 */
static void ShapeInternal(TPZVec<REAL> &x,int ord,TPZFMatrix &phi,TPZFMatrix &dphi,int transformation_index);

#ifdef _AUTODIFF
/**
 * Computes the values of the orthogonal shapefunctions before multiplying them by the
 * corner shapefunctions
 * The shape1dInternal function is extensively used by the shapefunction computation of
 * the other elements
 * @param x coordinate of the point (with derivative already setup)
 * @param num number of shape functions to be computed
 * @param phi shapefunction values and derivatives
 * @param transformation_index determines the transformation applied to the internal shape
 * functions. This parameter is computed by the GetTransformId1d method
 * @see GetTransformId1d
 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
 */
static void ShapeInternal(FADREAL & x,int num,TPZVec<FADREAL> & phi,int transformation_index);

#endif
/**
 * Computes the transformation applied to the variational parameter of the one-d element
 * The transformation applied to compensate for odd/even brokerage between elements can also
 * be viewed by the transformation of a variational parameter.
 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
 * @param in coordinate of the variational parameter
 * @param out transformed parameter
 */
static void TransformPoint1d(int transid,REAL in,REAL &out);
#ifdef _AUTODIFF
/**
 * Computes the transformation applied to the variational parameter of the one-d element
 * The transformation applied to compensate for odd/even brokerage between elements can also
 * be viewed by the transformation of a variational parameter.
 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
 * @param in coordinate of the variational parameter (with derivatives)
 * @param out transformed parameter (with derivatives)
 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
 */
static void TransformPoint1d(int transid,FADREAL & in,FADREAL &out);
#endif
/**
 * Applies the transformation on the values of the derivatives of the shape functions of the
 * internal shape functions
 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
 * @param num number of shapefunctions needed to transform
 * @in matrix containing the values of the derivatives of the shapefunctions as a row vector
 * the values of the derivatives contained in this matrix are modified upon return
 */
static void TransformDerivative1d(int transid,int num,TPZFMatrix &in);
#ifdef _AUTODIFF
/**
 * Applies the transformation on the values of the derivatives of the shape functions of the
 * internal shape functions
 * @param transid identifier of the transformation of the one-d element as obtained by the GetTransformId1d method
 * @param num number of shapefunctions needed to transform
 * @param in TPZVec of a differentiable type containing the values of the derivatives of the shapefunctions.
 * REMARK: The Derivative classes MUST store at least only one derivative - 1d problem
 */
//static void TransformDerivative1d(int transid,int num,TPZVec<FADREAL> &in);
#endif
/**
 * Computes the id of the transformation which will be applied on the parameter of the element\n
 * the return value is used in several methods of this class
 * @param id contains two distinct integer numbers which determine the orientation of the element
 * @return index of the tranformation
 */
static int GetTransformId1d(TPZVec<int> &id);
/**
 * Returns the transformation which transform a point from the interior of the element to the side
 * @param side side to which the point will be tranformed (0<=side<=2)
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


/**
 * Number of connects of the element (3)
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
static REAL RefElVolume(){return 2.0;}


};
#endif

