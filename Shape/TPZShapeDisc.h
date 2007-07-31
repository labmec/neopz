// -*- c++ -*-

// $Id: TPZShapeDisc.h,v 1.9 2007-07-31 23:16:45 tiago Exp $
#ifndef SHAPEDISCHPP
#define SHAPEDISCHPP

#include "pzfmatrix.h"
#include "pzvec.h"

/// groups all classes dedicated to the computation of shape functions
namespace pzshape {


/** 
 *
 * @brief Implements the shape functions discontinuous elements

 * This class computes the shape functions for n-dimensional elements
 * the shape functions can be tensor based or interpolation order based
 * The Shape2dFull method also computes higher order derivatives
 * @ingroup shape
 */
class TPZShapeDisc {

public:

/**
 * @param C:      normalizing constant, it provided by the discontinous element
 * @param x0:     interior point of the deformed element, it can be, for example, the center of mass of the element
 * @param x:      coordinate of the point
 * @param degree: degree of interpolation of the element
 * @param phi:    shapefunction values
 * @param dphi:   values of the derivatives of the shape functions
 * @param n:      number of derivatives to be computed
 */
static void Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi, int n = 1);

static void PolynomialWithoutScale(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi, int n = 1);

static void Legendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi, int n = 1);

/**
 * UseOrthoShape = 1 means it will be used Legendre polynomial as shape function. 
 * @since Mar 31, 2004
 */
static void (*fOrthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix & phi, TPZFMatrix & dphi, int n);

public:

enum MShapeType {ETensorial, EOrdemTotal};

TPZShapeDisc();

TPZShapeDisc(const TPZShapeDisc &copy);

~TPZShapeDisc();
  
protected:

/**
 * Number of shapefunctions dependent on the dimension and order of interpolation
 */
static int NDiscShapeF(int degree, int dimension, MShapeType type);
/**
 * discontinous polynomials of the line element
 */
static void Shape0D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * discontinous polynomials of the line element
 */
static void Shape1D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * discontinous bases of the two-dimensional elements 
 */
static void Shape2D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type);

/**
 * discontinous bases of the two-dimensional elements with many derivatives!
 */
static void Shape2DFull(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type);

/**
 * discontinous bases of the three-dimensional elements 
 */
static void Shape3D(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi, MShapeType type);

/**
 * @param x:      coordinate of the point
 * @param SingularPoint: singular point coordinate
 * @param phi:    shapefunction values
 * @param dphi:   values of the derivatives of the shape functions
 * @param n:      number of derivatives to be computed
 * Remember when set this attribute to also set fNumberOfSingularFunctions.
 */  
  void (*fSingularFunctions)(const TPZVec<REAL>& x, const TPZVec<REAL> &SingularPoint, TPZFMatrix & phi, TPZFMatrix & dphi, int n);

/** Number of singular functions computed in void *fSingularFunctions.
 * It must be according to void *fSingularFunctions.
 */
  int fNumberOfSingularFunctions;
  
/** Stores the singular side of the TPZCompElDisc element */  
  int fSingularSide;

  /** Computes the distance between two points */
  static REAL Distance(const TPZVec<REAL>& A, const TPZVec<REAL> &B);  
  
public:

/** Singular function of type phi = Sqrt[r] */
static void SqrtFunction(const TPZVec<REAL>& pt, const TPZVec<REAL> &SingularPoint, TPZFMatrix & phi, TPZFMatrix & dphi, int n = 1);  
  
/**
 * Number of shapefunctions dependent on the dimension and order of interpolation
 */
  int NShapeF(int degree, int dimension, TPZShapeDisc::MShapeType type);

/**
 * Discontinous base functions
 */
  void Shape(const TPZVec<REAL> &SingularPoint, REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree, int dim, TPZFMatrix &phi,TPZFMatrix &dphi, TPZShapeDisc::MShapeType type);
  void Shape(REAL C,TPZVec<REAL> &X0,TPZVec<REAL> &X,int degree, int dim, TPZFMatrix &phi,TPZFMatrix &dphi, TPZShapeDisc::MShapeType type);

/** Defines singular functions */
  void SetSingularShapeFunction(void (*f)(const TPZVec<REAL>& x, const TPZVec<REAL> &SingularPoint, TPZFMatrix & phi, TPZFMatrix & dphi, int n),
                                int NumberOfSingularFunctions,
                                int SingularSide);

/** Returns existence of singular shape function */
  bool HasSingularFunction(){
    return (fSingularFunctions != NULL);
  }
  
  /** Returns attribue fSingularSide */
  int SingularSide(){
    return this->fSingularSide;
  }

};

};

#endif
