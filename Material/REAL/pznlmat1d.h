/**
 * @file
 * @brief Contains the TPZNLMat1d class which implements the structure to evaluate non linear elements.
 */

#ifndef TPZNLMAT1D_H
#define TPZNLMAT1D_H

#include "TPZMaterial.h"

/**
 * @ingroup material
 * @brief Virtual class that implements the whole structure for evaluta non linear truss elements.
 * @author Edimar Cesar Rylo
 */
/**
 * The theory can be found at section 3.1 of the M. A. Crisfield book: \n
 * Non-Linear Finite Element Analysis of Solids and Structures: Volume 1 Essentials
 */
class TPZNLMat1d : public TPZMaterial
{
public:
	/** @brief Simple constructor */
	TPZNLMat1d(int id);
	
	/** @brief Default destructor */
	virtual ~TPZNLMat1d();
	
	/** @brief Returns the name of the material */
	virtual std::string Name()  override { return "nonlinear_1dMaterial"; }
	
	/** @brief Returns the integrable dimension of the material: Material is 1d */
	virtual int Dimension() const  override {return  1;}
	
	/** @brief Returns the number of state variables associated with the material: Only w? */
	virtual int NStateVariables() const override {return  1;}
	
	/** @brief Prints out the data associated with the material */
	virtual void Print(std::ostream &out = std::cout) override;
	
//	/** Returns the variable index associated with the name */
//	virtual int VariableIndex(const std::string &name);
	
//	virtual int NSolutionVariables(int var);
	
protected:
//	/** @brief Returns the solution associated with the var index based on the finite element approximation */
//	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol,
//						  TPZFMatrix<REAL> &axes, int var, TPZVec<REAL> &Solout);
public:
	
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	/** @brief Computes contribution to the stiffness matrix and right hand side at an integration point */
	virtual void Contribute(TPZMaterialData &data, REAL weight,
                            TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
//	/** @brief Computes contribution to the right hand side at an integration point */
//	virtual void Contribute(TPZMaterialData &data, REAL weight,
//                            TPZFMatrix<STATE> &ef);

	/** @brief Computes contribution to the stiffness matrix and right hand side at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
	/** @brief Computes contribution to load vector (right hand side) at the integration point of a boundary */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
//	/** @brief To create another material of the same type */
//	virtual TPZMaterial * NewMaterial();
	
//	/** @brief Reads data of the material from a istream (file data) */
//	virtual void SetData(std::istream &data);
		
//	void Write(TPZStream &buf, int withclassid) const override;
	
//	void Read(TPZStream &buf, void *context) override;
	
//	public:
virtual int ClassId() const override;

	
	virtual STATE Eps(TPZVec<STATE> &sol,TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &dphi) = 0;
	
protected:
	/** @brief Cross Section Area */
	REAL fArea;
	
	/** @brief Young's modulus */
	REAL fE;
};

#endif
