
#ifndef SHAPEDISCHPP
#define SHAPEDISCHPP

#include "pzfmatrix.h"
#include "pzvec.h"


class TPZShapeDisc {

public:

/**
 * @param C:      normalizing constant
 * @param x0:     interior point of the deformed element
 * @param x:      coordinate of the point
 * @param degree: degree of interpolation of the element
 * @param phi:    shapefunction values
 * @param dphi:   values of the derivatives of the shape functions
 */
static void Polynomial(REAL C,REAL x0,REAL x,int degree,TPZFMatrix & phi,TPZFMatrix & dphi);

 public:

    TPZShapeDisc();
    ~TPZShapeDisc() {};

/**
 * discontinous polynomials of the line element
 */
static void Shape1D(REAL C,TPZVec<REAL> X0,TPZVec<REAL> X,int degree,TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * discontinous bases of the two-dimensional elements 
 */
static void Shape2D(REAL C,TPZVec<REAL> x0,TPZVec<REAL> x,int degree,TPZFMatrix &phi,TPZFMatrix &dphi);

/**
 * discontinous bases of the three-dimensional elements 
 */
static void Shape3D(REAL C,TPZVec<REAL> x0,TPZVec<REAL> x,int degree,TPZFMatrix &phi,TPZFMatrix &dphi);

};
#endif
