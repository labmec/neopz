/**
 * @file
 * @brief Contains declaration of TPZInterpolationSpace class which implements the interface for interpolated computational and interface elements.
 */

#ifndef PZINTERPOLATIONSPACE_H
#define PZINTERPOLATIONSPACE_H

#include "pzcompel.h"
class TPZMaterialData;

/**
 * @brief Implements the interfaces for TPZCompElDisc, TPZInterfaceElement and TPZInterpolatedElement. \ref CompElement "Computational element"
 * @since April 11, 2007
 * @ingroup CompElement
 */
class TPZInterpolationSpace : public TPZCompEl
{
public:
	
	/** @brief Default constructor */
	TPZInterpolationSpace();
	
	/** @brief Default destructor */
	virtual ~TPZInterpolationSpace();
	
	/** @brief Puts a copy of the element in the referred mesh */
	TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy);
	
	/** @brief Puts a copy of the element in the patch mesh */
	TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, std::map<int,int> &gl2lcElMap);
	
	/** @brief Copy of the element in the new mesh whit alocated index */
	TPZInterpolationSpace(TPZCompMesh &mesh, const TPZInterpolationSpace &copy, int &index);
	
	/**
	 * @brief Create a computational element within mesh
	 * @param mesh mesh wher will be created the element
	 * @param gel geometrical element to insert
	 * @param index new elemen index
	 */
	/** Inserts the element within the data structure of the mesh */
	TPZInterpolationSpace(TPZCompMesh &mesh, TPZGeoEl *gel, int &index);
	
	/** @brief It returns the shapes number of the element */
	virtual int NShapeF() const = 0;
	
	/** @brief Returns the number of shapefunctions associated with a connect*/
	virtual int NConnectShapeF(int inod) const = 0;
	
	/** @brief Returns the max order of interpolation. */
	virtual int MaxOrder();
	
	/** 
	 * @brief Computes the shape function set at the point x. 
	 * @param qsi point in master element coordinates
	 * @param phi vector of values of shapefunctions, dimension (numshape,1)
	 * @param dphi matrix of derivatives of shapefunctions, dimension (dim,numshape)
	 */
	/**
	 * This method uses the order of interpolation
	 * of the element along the sides to compute the number of shapefunctions
	 */
	virtual void Shape(TPZVec<REAL> &qsi,TPZFMatrix &phi,TPZFMatrix &dphi) = 0;
	
	/** @brief Compute shape functions based on master element in the classical FEM manner. */
	virtual void ComputeShape(TPZVec<REAL> &intpoint, TPZVec<REAL> &X, TPZFMatrix &jacobian, TPZFMatrix &axes,
							  REAL &detjac, TPZFMatrix &jacinv, TPZFMatrix &phi, TPZFMatrix &dphix);
	
	/** 
	 * @brief Initialize a material data and its attributes based on element dimension, number
	 * of state variables and material definitions
	 */
	virtual void InitMaterialData(TPZMaterialData &data);
	
	/** @brief Compute and fill data with requested attributes */
	virtual void ComputeRequiredData(TPZMaterialData &data,
									 TPZVec<REAL> &qsi);
	
	/** @brief Computes the proper normal vector towards the neighbour element */
	virtual void ComputeNormal(TPZMaterialData & data);
	
	/** @brief Computes the vectorial product of two vectors and normalize the result if unitary is set to true */
	void VectorialProd(TPZVec<REAL> & ivec, TPZVec<REAL> & jvec, TPZVec<REAL> & kvec, bool unitary = false);
	
	/**
	 * @brief Computes the element stiffness matrix and right hand side
	 * @param ek element matrix
	 * @param ef element right hand side
	 */
	virtual void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/**
	 * @brief Only computes the element residual
	 * @param ef element residual
	 */
	virtual void CalcResidual(TPZElementMatrix &ef);
	
	/** @brief Initialize element matrix in which is computed CalcStiff */
	void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
	/** @brief Initialize element matrix in which is computed in CalcResidual */
	void InitializeElementMatrix(TPZElementMatrix &ef);
	
	/** @brief Returns minimum and maximum values for each state variable */
	/** 
	 * It is not a cheap method because it computes solution for
	 * all integration points ( with intrule.MaxOrder() )
	 */
	void MinMaxSolutionValues(TPZVec<REAL> &min, TPZVec<REAL> &max);
	
	/** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
	virtual const TPZIntPoints &GetIntegrationRule() const = 0;
	
	/** @brief Returns the inner radius value. */
	virtual REAL InnerRadius();
	
	/**
	 * @brief Post processing method which computes the solution for the var post processed variable.
	 * @param qsi coordinate of the point in master element space where the solution will be evaluated
	 * @param var variable which will be computed
	 * @param sol (output) solution computed at the given point
	 * @see TPZMaterial::VariableIndex
	 * @see TPZMaterial::NSolutionVariables
	 * @see TPZMaterial::Solution
	 */
	/** The var index is obtained by calling the TPZMaterial::VariableIndex method with a post processing name */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<REAL> &sol);
	
	/**
	 * @brief Interpolates the solution into the degrees of freedom nodes from the degrees
	 * of freedom nodes from the coarse element
	 */
	void InterpolateSolution(TPZInterpolationSpace &coarsel);
	
	/**
	 * @brief Create interfaces between this and its neighbours.
	 * @param BetweenContinuous allows to create interface between two elements that are not TPZCompElDisc.
	 */
	/** If param is false, it is necessary to have at least one TPZCompElDisc. */
	void CreateInterfaces(bool BetweenContinuous = false);
	
	/** @brief Create an interface between this and the neighbour by side side.
	 * @param side : side where interface must be created
	 * @param BetweenContinuous allows to create interface between two elements that are not TPZCompElDisc. If param is false, it is necessary to have at least one TPZCompElDisc.
	 * Returns the interface created.
	 */
	TPZInterfaceElement * CreateInterface(int side, bool BetweenContinuous = false);
	
	/** @brief Verify existence of interface */
	int ExistsInterface(TPZGeoElSide geosd);
	
	/** @brief Remove interfaces connected to this element */
	void RemoveInterfaces();
	
	/** @brief Remove interface which is neighbour from side side */
	void RemoveInterface(int side);
	
	/**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors (output) the L2 norm or true error of the error of the solution
	 * @param flux (input) value of the interpolated flux values
	 */
	virtual void EvaluateError(  void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &val,TPZFMatrix &deriv),
                               TPZVec<REAL> &errors,TPZBlock * flux );
	
	/** @brief Computes the element error estimator */
	virtual void ComputeError(int errorid, TPZVec<REAL> &error);
	
	/** @brief Integrate a variable over the element. */
	virtual void Integrate(int variable, TPZVec<REAL> & value);
	
	/** @brief Integrate the solution over the element */
	virtual void IntegrateSolution(TPZVec<REAL> & value);
	
	/**
	 * @brief Will project the flux associated with the variational statement onto the finite element interpolation space
	 * @param ek projection matrix
	 * @param ef inner product of the flux with the finite element interpolation space
	 */
	/**
	 * The ek matrix corresponds to an L2 (scalar) projection, the ef matrix contains multiple right hand sides, one
	 * for each component of the flux
	 */
	void ProjectFlux(TPZElementMatrix &ek, TPZElementMatrix &ef);
	
protected:
	
	int fPreferredOrder;
	
public:
	
	/**  @brief Defines the desired order for entire element. */
	virtual void SetPreferredOrder ( int order ) = 0;
	
	/** @brief Returns the prefered order for the element */
	virtual int GetPreferredOrder () { return fPreferredOrder; }
	
	/**
	 * @brief Change the preferred order for the element and proceed the
	 * adjust of the aproximation space \n taking in acount the type
	 * of formulation and the neighbours of the element
	 */
	virtual void PRefine ( int order ) = 0;
	
public:
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	/**
	 * @brief Accumulates the transfer coefficients between the current element and the
	 * coarse element \n into the transfer matrix, using the transformation t
	 * @param coarsel larger element with respect to which the transfer matrix is computed
	 * @param t transformation which maps the master element space of the current element into the master element space of the coarse element
	 * @param transfer transfer matrix mapping the solution of the coarse mesh into the fine mesh
	 */
	/**
	 * This method forms the basis for the multigrid method
	 */
	void BuildTransferMatrix(TPZInterpolationSpace &coarsel, TPZTransform &t, TPZTransfer &transfer);
	
protected:
	
	/**
	 * @brief Auxiliary method to expand a vector of shapefunctions and their derivatives to acount for constraints
	 * @param connectlist (input) vector of all connects to which the element will contribute
	 * @param dependencyorder (input) vector of indices which indicate the order in which the connects will be processed
	 * @param blocksizes (output) number of shapefunctions associated with each connect
	 * @param phi (input/output) values of the shapefunctions
	 * @param dphi (input/output) values of the derivatives of the shapefunctions
	 */
	/**
	 * As input the regular values of the shapefunctions are given and their derivatives\n
	 * if these shapefunctions are dependent upon other shapefunctions (because of constraints) then the vectors
	 * are expanded to include the value of the independent shapefunctions and their derivatives as well
	 */
    void ExpandShapeFunctions(TPZVec<int> &connectlist, TPZVec<int> &dependencyorder, TPZVec<int> &blocksizes, TPZFMatrix &phi, TPZFMatrix &dphi);
	
};

#endif
