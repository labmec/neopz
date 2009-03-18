//$Id: pzinterpolationspace.h,v 1.12 2009-03-18 13:56:43 fortiago Exp $

#ifndef PZINTERPOLATIONSPACE_H
#define PZINTERPOLATIONSPACE_H

#include "pzcompel.h"
class TPZMaterialData;

/**
This class implements the interfaces for TPZCompElDisc, TPZInterfaceElement and TPZInterpolatedElement.

@since April 11, 2007
*/
class TPZInterpolationSpace : public TPZCompEl
{
public:

  /**
   * Simple Constructor
   */
  TPZInterpolationSpace();

  /**
   * Simple destructor
   */
  virtual ~TPZInterpolationSpace();

  /**
   * put a copy of the element in the referred mesh
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy);

  /**
   * put a copy of the element in the patch mesh
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int,int> &gl2lcElMap);

  /**
   * copy of the element in the new mesh whit alocated index
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int &index);

  /**
   * Create a computational element within mesh
   * Inserts the element within the data structure of the mesh
   * @param mesh mesh wher will be created the element
   * @param index new elemen index
   */
  TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);

  /**
   * it returns the shapes number of the element
   */
  virtual int NShapeF() = 0;

  /**returns the number of shapefunctions associated with a connect*/
  virtual int NConnectShapeF(int inod) = 0;

  /** Returns the max order of interpolation. */
  virtual int MaxOrder();

  /**computes the shape function set at the point x. This method uses the order of interpolation
   * of the element along the sides to compute the number of shapefunctions
   * @param qsi point in master element coordinates
   * @param phi vector of values of shapefunctions, dimension (numshape,1)
   * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
   */
  virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix &phi,TPZFMatrix &dphi) = 0;

  /** Compute shape functions based on master element in the classical FEM manner.
   */
  virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes,
                            REAL &detjac, TPZFMatrix &jacinv, TPZFMatrix &phi, TPZFMatrix &dphix);

  /** Initialize a material data and its attributes based on element dimension, number
   * of state variables and material definitions
   */
  void InitMaterialData(TPZMaterialData &data);

  /** Compute and fill data with requested attributes */
  virtual void ComputeRequiredData(TPZMaterialData &data,
                           TPZVec<REAL> &qsi);

  /**
   * CalcStiff computes the element stiffness matrix and right hand side
   * @param ek element matrix
   * @param ef element right hand side
   */
  virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /**
   * CalcResidual only computes the element residual
   * @param ef element residual
   */
  virtual void CalcResidual(TPZElementMatrix &ef);

  /** Initialize element matrix in which is computed CalcStiff
  */
  void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);

  /** Initialize element matrix in which is computed in CalcResidual
  */
  void InitializeElementMatrix(TPZElementMatrix &ef);

  /** Returns minimum and maximum values for each state variable.
   * It is not a cheap method because it computes solution for
   * all integration points ( with intrule.MaxOrder() )
   */
  void MinMaxSolutionValues(TPZVec<REAL> &min, TPZVec<REAL> &max);

  /**returns a reference to an integration rule suitable for integrating
     the interior of the element */
  virtual TPZIntPoints &GetIntegrationRule() = 0;

  /**
   * Returns the inner radius value.
   */
  virtual REAL InnerRadius();

  /**
   * Post processing method which computes the solution for the var post processed variable. The var index is obtained
   * by calling the TPZMaterial::VariableIndex method with a post processing name
   * @param qsi coordinate of the point in master element space where the solution will be evaluated
   * @param var variable which will be computed
   * @param sol (output) solution computed at the given point
   * @see TPZMaterial::VariableIndex
   * @see TPZMaterial::NSolutionVariables
   * @see TPZMaterial::Solution
   */
  virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);

  /**
   * Interpolates the solution into the degrees of freedom nodes from the degrees
   * of freedom nodes from the coarse element
   */
  void InterpolateSolution(TPZInterpolationSpace &coarsel);

  /** Create interfaces between this and its neighbours.
   * Param BetweenContinuous allows to create interface between two elements that are not TPZCompElDisc.
   * If param is false, it is necessary to have at least one TPZCompElDisc.
   */
  void CreateInterfaces(bool BetweenContinuous = false);

  /** Create an interface between this and the neighbour by side side.
   * @param side : side where interface must be created
   * @param BetweenContinuous allows to create interface between two elements that are not TPZCompElDisc. If param is false, it is necessary to have at least one TPZCompElDisc.
   * Returns the interface created.
   */
  TPZInterfaceElement * CreateInterface(int side, bool BetweenContinuous = false);

  /** Verify existence of interface */
  int ExistsInterface(TPZGeoElSide geosd);

  /** Remove interfaces connected to this element */
  void RemoveInterfaces();

  /** Remove interface which is neighbour from side side */
  void RemoveInterface(int side);

  /**
   * Perform an error estimate on the elemen
   * @param fp function pointer which computes the exact solution
   * @param true_error (output)  the true error of the solution
   * @param L2_error (output) the L2 norm of the error of the solution
   * @param flux (input) value of the interpolated flux values
   * @param estimate (output) estimated error based on the implemented criterium
   */
  virtual void EvaluateError(  void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
                               TPZVec<REAL> &errors,TPZBlock * flux );
                               
  /**
   * ComputeError computes the element error estimator
   */
  virtual void ComputeError(int errorid, TPZVec<REAL> &error);

  /**
   * Integrate a variable over the element.
   */
   virtual void Integrate(int variable, TPZVec<REAL> & value);

   /** 
    * Integrate the solution over the element
    */
   virtual void IntegrateSolution(TPZVec<REAL> & value);

  /**
   * Will project the flux associated with the variational statement onto the finite element interpolation space
   * The ek matrix corresponds to an L2 (scalar) projection, the ef matrix contains multiple right hand sides, one
   * for each component of the flux
   * @param ek projection matrix
   * @param ef inner product of the flux with the finite element interpolation space
   */
  void ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef);

protected:

  int fPreferredOrder;

public:

  /**
   *  Define the desired order for entire element.
   */
  virtual void SetPreferredOrder ( int order ) = 0;

  /**
   * Return the prefered order for the element
   */
  virtual int GetPreferredOrder () { return fPreferredOrder; }

  /**
   * Change the preferred order for the element and proceed the
   * adjust of the aproximation space taking in acount the type
   * of formulation and the neighbours of the element
   */
  virtual void PRefine ( int order ) = 0;

public:

  /**
   * Save the element data to a stream
   */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  virtual void Read(TPZStream &buf, void *context);

  /**
   * Accumulates the transfer coefficients between the current element and the
   * coarse element into the transfer matrix, using the transformation t
   * This method forms the basis for the multigrid method
   * @param coarsel larger element with respect to which the transfer matrix is computed
   * @param t transformation which maps the master element space of the current element into the master element space of the coarse element
   * @param transfer transfer matrix mapping the solution of the coarse mesh into the fine mesh
   */
  void BuildTransferMatrix(TPZInterpolationSpace &coarsel, TPZTransform &t, TPZTransfer &transfer);

  protected:

  /**
   * Auxiliary method to expand a vector of shapefunctions and their derivatives to acount for constraints
   * As input the regular values of the shapefunctions are given and their derivatives\n
   * if these shapefunctions are dependent upon other shapefunctions (because of constraints) then the vectors
   * are expanded to include the value of the independent shapefunctions and their derivatives as well
   * @param (input) connectlist vector of all connects to which the element will contribute
   * @param (input) dependencyorder vector of indices which indicate the order in which the connects will be processed
   * @param blocksizes (output) number of shapefunctions associated with each connect
   * @param phi (input/output) values of the shapefunctions
   * @param dphi (input/output) values of the derivatives of the shapefunctions
   */
    void ExpandShapeFunctions(TPZVec<int> &connectlist, TPZVec<int> &dependencyorder, TPZVec<int> &blocksizes, TPZFMatrix &phi, TPZFMatrix &dphi);

};

#endif
