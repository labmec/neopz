// -*- c++ -*-

// $Id: TPZShapeDisc.h,v 1.7 2005-02-28 22:11:26 phil Exp $
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

static void Legendre(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi, int n = 1);

/**
 * UseOrthoShape = 1 means it will be used Legendre polynomial as shape function. 
 * @since Mar 31, 2004
 */
static void (*fOrthogonal)(REAL C, REAL x0, REAL x,int degree, TPZFMatrix & phi, TPZFMatrix & dphi, int n);

public:

enum MShapeType {ETensorial, EOrdemTotal};

TPZShapeDisc();

~TPZShapeDisc() {};
  
/**
 * Number of shapefunctions dependent on the dimension and order of interpolation
 */
static int NShapeF(int degree, int dimension, MShapeType type);
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

};

};

#endif
