//HEADER FILE FOR CLASS CYLINSYS

#ifndef CYLINSYSTHPP
#define CYLINSYSTHPP

#include "pzcosys.h"
#include "pzreal.h"

/**
 * @ingroup cosys
 * Defines the cylindrical coordinate system
 */

class  TPZCylinsys : public TPZCosys {

public:

  /**
   * Default empty constructor for a TPZCylinsys object
   */
  TPZCylinsys();

  /**
   * Create object from one reference object and new origin
   * @param num index of the coordinate system
   * @param ref reference coordinate system
   */
  TPZCylinsys(int num, TPZCartsys* ref = NULL);

  /**
   * Default Destructor 
   */
  ~TPZCylinsys() {;}

  /**
   * Converts point given in cylindric coordinate system to cartesian reference system.
   * The reference cartesian coordinate system has the same origin and "z" axe of the
   * cylindric system 
   * @param point point to transfer to reference coordinate system   
   */
  void ToReference(TPZVec<REAL> &point);

  /**
   * Converts point given in cartesian reference system to cylindric coordinate system.
   * The reference cartesian coordinate system has the same origin and "z" axe of the
   * cylindric system 
   * @param point point to transfer from reference coordinate system
   */
  void FromReference(TPZVec<REAL> &point);

  /**
   * Method to return the current coordinate system type. 
   */
  int Type() { return cylindric;}

  /**
   * Verify if the difference between two nodes are greater than PI
   * @param points points to verify the angle
   */
  void VerifyRange(TPZFMatrix &points);

  /**
   * Calculates the transformation gradient Gradx given Gradient of  the point X
   * @param X point that want to transform gradient
   * @param GradX gradient of point X
   * @param x will receive point X in terms of the destintation coordinate system
   * @param Gradx will receive the tranformed gradient
   * @param dest coordinate system to transform the gradient
   */
  void TransformGradient(TPZVec<REAL> &X, TPZFMatrix &GradX, TPZVec<REAL> &x, TPZFMatrix &Gradx, TPZCosys *dest = 0);


};

#endif

