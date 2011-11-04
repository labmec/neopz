/**
 * \file
 * @brief Contains the TPZNLMat1d class which implements the structure to evaluate non linear elements.
 */
//$Id: pznlmat1d.h,v 1.7 2009-09-01 19:44:47 phil Exp $
// -*- c++ -*-

#ifndef TPZNLMAT1D_H
#define TPZNLMAT1D_H

#include "pzmaterial.h"
//#include <pzelasmat.h>

/**
 * @ingroup material
 * @brief Virtual class that implements the whole structure for evaluta non linear truss elements.
 * @author Edimar Cesar Rylo
 */
/**
 * The theory can be found at section 3.1 of the M. A. Crisfield book: \n
 * Non-Linear Finite Element Analysis of Solids and Structures: Volume 1 Essentials
 */
class TPZNLMat1d : public TPZMaterial//TPZElasticityMaterial
{
public:
	/** @brief Simple constructor */
	TPZNLMat1d(int id);
	
	/** @brief Default destructor */
	~TPZNLMat1d();
	
	/** @brief Returns the name of the material */
	virtual std::string Name() { return "nonlinear_1dMaterial"; }
	
	/** @brief Returns the integrable dimension of the material: Material is 1d */
	virtual int Dimension() {return  1;}
	
	/** @brief Returns the number of state variables associated with the material: Only w? */
	virtual int NStateVariables() {return  1;}
	
	/** @brief Prints out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout);
	
	/** Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix &DSol,
						  TPZFMatrix &axes, int var, TPZVec<REAL> &Solout);
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at an integration point */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ek,
                            TPZFMatrix &ef);
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
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
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	/** @brief To create another material of the same type */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief Reads data of the material from a istream (file data) */
	virtual void SetData(std::istream &data);
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at an integration point */
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef);

	virtual void Write(TPZStream &buf, int withclassid);
	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual int ClassId() const;
	
	virtual REAL Eps(TPZVec<REAL> &sol,TPZFMatrix &axes,TPZFMatrix &dphi) = 0;
	
protected:
	/** @brief Cross Section Area */
	REAL fArea;
	
	/** @brief Young's modulus */
	REAL fE;
};

#endif
