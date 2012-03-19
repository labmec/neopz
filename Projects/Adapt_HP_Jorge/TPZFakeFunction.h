//$Id: TPZFakeFunction.h,v 1.1 2009-08-28 19:32:09 fortiago Exp $

#ifndef TPZFAKEFUNCTION_H
#define TPZFAKEFUNCTION_H

#include "pzsave.h"
#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzfunction.h"

/**
  * Implements a function. Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
  * @since August 01, 2007
  */
class TPZFakeFunction : public TPZFunction{
public:

  /**
   * Class constructor
   */
   TPZFakeFunction(){
   }

  /**
   * Class destructor
   */
  ~TPZFakeFunction(){
  }
    
  /**
   * Performs function computation
   * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
   * @param f function values
   * @param df function derivatives
   */
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df){
    f.Resize(3);
    f.Fill(0.);
    df.Resize(3,3);
    df.Zero();
  }
  
 /** Returns number of functions.
  */ 
  virtual int NFunctions(){ return 3; }
  
  /** Polynomial order of this function. In case of non-polynomial
   * function it can be a reasonable approximation order.
   */
  virtual int PolynomialOrder(){ return 0; }

};

#endif

