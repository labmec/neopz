/*
 *  pznondarcyflow.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzdiscgal.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "TPZAutoPointer.h"
#include "problemdata.h"


#ifndef TPZNONDARCYFLOW
#define TPZNONDARCYFLOW


#ifdef _AUTODIFF
#include "fad.h"



typedef TFad<3, REAL> VarFad;


class TPZNonDarcyFlow : public TPZMaterial {

	
	
private:
	
	
	
public:
	
    	
	/**
	 * Empty Constructor
	 */
	TPZNonDarcyFlow();
	
	/** Creates a material object and inserts it in the vector of
	 *  material pointers of the mesh.
	 */
    TPZNonDarcyFlow(int matid);
	
	
	/** Creates a material object based on the referred object and
	 *  inserts it in the vector of material pointers of the mesh.
	 */
	TPZNonDarcyFlow(const TPZNonDarcyFlow &mat);
	
	/**
	 * Destructor
	 */
	virtual ~TPZNonDarcyFlow();
	
	/** Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered
	 * as necessary. Each derived class may optimize performance by selecting
	 * only the necessary data.
	 * @since April 10, 2007
	 */
	virtual void FillDataRequirements(TPZMaterialData &data);
	
	virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data);	
	
	/** returns the name of the material */
	virtual std::string Name() {
		return "TPZNonDarcyFlow";
	}
	
	/** returns the integrable dimension of the material */
	virtual int Dimension() const {return 2;}
	
	/** returns the number of state variables associated with the material */
	virtual int NStateVariables() {return 1;} // for hdiv are 3
	
	/** print out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout);
	
	/** returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	/** returns the number of variables associated with the variable
	 indexed by var.  var is obtained by calling VariableIndex */
	virtual int NSolutionVariables(int var);
	
	/** returns the solution associated with the var index based on
	 * the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
	///Metodos Contribute
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data[in] stores all input data
	 * @param weight[in] is the weight of the integration rule
	 * @param ek[out] is the stiffness matrix
	 * @param ef[out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data[in] stores all input data
	 * @param weight[in] is the weight of the integration rule
	 * @param ek[out] is the stiffness matrix
	 * @param ef[out] is the load vector
	 * @param bc[in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	/**
	 * Unique identifier for serialization purposes
	 */
	public:
virtual int ClassId() const;

	
	/**
	 * Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/**
	 * Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
	/**
	 * Return the viscosity based on p and if FAD also the correct derivative
	 */
	template<class T>
	void ReturnVisc(T &p, T &visc);
	
	template<class T>
	void ReturnRho(T &p, T &rho);
	
	template<class T>
	void ReturnPor(T &p, T &por);
	
};

#endif

#endif
