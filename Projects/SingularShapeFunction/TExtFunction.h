//$Id: TExtFunction.h,v 1.1 2008-02-06 18:18:44 tiago Exp $

#ifndef TEXTFUNCTION_H
#define TEXTFUNCTION_H

#include "pzsave.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfunction.h"

/**
  * Implements a function. Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
  * @since August 01, 2007
  */
class TExtFunction : public TPZFunction{
public:

  /**
   * Class constructor
   */
  TExtFunction();

  /**
   * Class destructor
   */
  ~TExtFunction();
    
  /**
   * Performs function computation
   * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
   * @param f function values
   * @param df function derivatives
   */
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df);
  
 /** Returns number of functions.
  */ 
  virtual int NFunctions();
  
  /** Polynomial order of this function. In case of non-polynomial
   * function it can be a reasonable approximation order.
   */
  virtual int PolynomialOrder();

};

#endif

