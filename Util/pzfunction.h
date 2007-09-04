//$Id: pzfunction.h,v 1.1 2007-09-04 12:35:22 tiago Exp $

#ifndef PZFUNCTION_H
#define PZFUNCTION_H

#include "pzsave.h"
#include "pzvec.h"
#include "pzfmatrix.h"

const int TPZFUNCTIONID = 9000;

/**
  * Implements a function. Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
  * @since August 01, 2007
  */
class TPZFunction : public TPZSaveable{
public:

  /**
   * Class constructor
   */
  TPZFunction();

  /**
   * Class destructor
   */
  ~TPZFunction();
    
  /**
   * Performs function computation
   * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
   * @param f function values
   * @param df function derivatives
   */
  virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix &df) = 0;
  
  /**
   * Unique identifier for serialization purposes
   */
  virtual int ClassId() const;

  /**
   * Save the element data to a stream
   */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  virtual void Read(TPZStream &buf, void *context);

};

#endif
