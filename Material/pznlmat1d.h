//$Id: pznlmat1d.h,v 1.2 2007-01-27 14:49:27 phil Exp $
// -*- c++ -*-

#ifndef TPZNLMAT1D_H
#define TPZNLMAT1D_H

#include <pzmaterial.h>
//#include <pzelasmat.h>

/**
Virtual class that implements the whole structure for evaluta non linear truss elements.
The theory can be found at section 3.1 of the M. A. Crisfield book: Non-Linear Finite Element Analysis of Solids and Structures: Volume 1 Essentials
 
@author Edimar Cesar Rylo
*/
class TPZNLMat1d : public TPZMaterial//TPZElasticityMaterial
{
public:
  /**
   * Constructor
   */
  TPZNLMat1d(int id);

  /**
   * Destructor
   */
  ~TPZNLMat1d();

  /**
    * Returns the name of the material
    */
  virtual char *Name() { return "nonlinear_1dMaterial"; }

  /**
    * Returns the integrable dimension of the material:
    * Material is 1d
    */
  virtual int Dimension() {return  1;}

  /**
    * Returns the number of state variables associated with the material:
    * Only w?
    */
  virtual int NStateVariables() {return  1;}

  /**
    * Print out the data associated with the material
    */
  virtual void Print(std::ostream &out = std::cout);

  /**
    * Returns the variable index associated with the name
    */
  virtual int VariableIndex(char *name);

  /**
    * Returns the number of variables associated with the variable\
    * indexed by var.  var is obtained by calling VariableIndex
    */
  virtual int NSolutionVariables(int var);

  /**
    * Returns the solution associated with the var index based on\
    * the finite element approximation
    */
  virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
                        TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);


  /**
    * Compute contribution to the stiffness matrix and right hand\
    * side at an integration point
    */
  virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                          TPZVec<REAL> &sol, TPZFMatrix &dsol,
                          REAL weight,TPZFMatrix &axes,
                          TPZFMatrix &phi, TPZFMatrix &dphi,
                          TPZFMatrix &ek, TPZFMatrix &ef);

  /**
    * Compute contribution to the stiffness matrix and right hand\
    * side at the integration point of a boundary
    */
  virtual void ContributeBC(TPZVec<REAL> &x, TPZVec<REAL> &sol,
                            REAL weight, TPZFMatrix &axes,
                            TPZFMatrix &phi, TPZFMatrix &ek,
                            TPZFMatrix &ef, TPZBndCond &bc);

  /**
    * To create another material of the same type
    */
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();

  /**
    * Read data of the material from a istream (file data)
    */
  virtual void SetData(std::istream &data);

  /**
    * Compute contribution to the stiffness matrix and right hand\
    * side at an integration point
    */
  virtual void Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
                          TPZVec<REAL> &sol, TPZFMatrix &dsol, REAL weight,
                          TPZFMatrix &axes, TPZFMatrix &phi,
                          TPZFMatrix &dphi, TPZFMatrix &ef);


  /**
    * Save the element data to a stream
    */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
    * Read the element data from a stream
    */
  virtual void Read(TPZStream &buf, void *context);

  virtual int ClassId() const;

  virtual REAL Eps(TPZVec<REAL> &sol,TPZFMatrix &axes,TPZFMatrix &dphi) = 0;

protected:
  /**
   * Cross Section Area
   */
   double fArea;

  /**
   * Young's modulus
   */
   double fE;
};

#endif
