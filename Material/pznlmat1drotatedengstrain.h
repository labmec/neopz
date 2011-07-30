/**
 * \file
 * @brief Contains the TPZNLMat1dRotatedEngStrain class which implements a non linear 1d material.
 */
//$Id: pznlmat1drotatedengstrain.h,v 1.7 2009-11-16 18:41:59 diogo Exp $
// -*- c++ -*-

#ifndef TPZNLMAT1DROTATEDENGSTRAIN_H
#define TPZNLMAT1DROTATEDENGSTRAIN_H

#include "pznlmat1d.h"
#include <fstream>
#include <iostream>

/**
 * @ingroup material
 * @brief Implements a non linear 1d material based on a rotated engineering strain \
 measurement.
 @author Edimar Cesar Rylo
 @since May, 2006
 */
/**
 * The implementation is based on the section 3.1.1 of the book Non-linear \n
 Finite Element Analysis of Solids and Structures - Volume 1: Essentials \n
 of M. A. Crisfield
 */
class TPZNLMat1dRotatedEngStrain : public TPZNLMat1d
{
public:
    /**
     * Simple Constructor
     */
    TPZNLMat1dRotatedEngStrain(int id);
	
    /**
	 * Destructor
     */
    virtual ~TPZNLMat1dRotatedEngStrain();
	
    /**
     * Returns the name of the material
     */
    virtual std::string Name() { return "nonlinear_RotatedEnginneringStrain"; }
	
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
	virtual int VariableIndex(const std::string &name);
	
    /**
     * Returns the number of variables associated with the variable\
     * indexed by var.  var is obtained by calling VariableIndex
     */
    virtual int NSolutionVariables(int var);
	
protected:
	/**
	 * Returns the solution associated with the var index based on\
	 * the finite element approximation
	 */
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
						  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		Solution(data.sol,data.dsol,data.axes,var,Solout);
    }
	
	
    /**
     * Compute contribution to the stiffness matrix and right hand\
     * side at an integration point
     */
    virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
							TPZFMatrix &ef);
	
    /**
     * Compute contribution to the stiffness matrix and right hand\
     * side at the integration point of a boundary
     */
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
	{
		TPZNLMat1d::ContributeBC(data,weight,ef,bc);
	}
	
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
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix &ef);
	
	
    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid);
	
    /**
     * Read the element data from a stream
     */
    virtual void Read(TPZStream &buf, void *context);
	
    virtual REAL Eps(TPZVec<REAL> &sol,TPZFMatrix &axes,TPZFMatrix &dphi) = 0;
};

#endif
