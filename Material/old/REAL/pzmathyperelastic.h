/**
 * @file
 * @brief Contains the TPZMatHyperElastic class which implements a hyper elasticity material.
 */

#ifndef MATHYPERELASTICHPP
#define MATHYPERELASTICHPP

#include "TPZMaterial.h"
#include "pzfmatrix.h"

#include "fadType.h"

/**
 * @ingroup material
 * @brief Implements a hyper elasticity material.
 */
class TPZMatHyperElastic : public TPZMaterial {
	
    STATE fXf[3];

	STATE fE1[3],fE5[3],fE9[3];

    
	STATE fLambda,fNu,fE,fMu;
	STATE fCoef1,fCoef2,fCoef3;
	
	public :
	
	TPZMatHyperElastic(int nummat,STATE e,STATE mu,STATE nu=-1.,STATE lambda=-1.,STATE coef1=-1.,STATE coef2=-1.,STATE coef3=-1.);
	
	virtual ~TPZMatHyperElastic();
	
	void SetMaterial(TPZFMatrix<STATE> &xfin){
		fXf[0] = xfin(0,0);
		fXf[1] = xfin(1,0);
		fXf[2] = xfin(2,0);
	}
	
	int Dimension() const  override { return 3;}
	
	virtual int NStateVariables() const override;
	
	virtual void Print(std::ostream & out) override;
	
	std::string Name()  override { return "TPZMatHyperElastic"; }
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	

	/** @brief Computes contribution to the energy at an integration point*/
	virtual void ContributeEnergy(TPZVec<REAL> &x,
								  TPZVec<FADFADREAL> &sol,
								  TPZVec<FADFADREAL> &dsol,
								  FADFADREAL &U,
								  REAL weight);
	
	static void ComputeEnergy(STATE lambda, STATE mu,  TPZFMatrix<STATE> &dsol, TFad<9,TFad<9,STATE> > &energy);
	
	/** @brief Computes contribution of BC to the Energy*/
	virtual void ContributeBCEnergy(TPZVec<REAL> & x,
									TPZVec<FADFADREAL> & sol, FADFADREAL &U,
									REAL weight, TPZBndCond &bc);
	

	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
		        TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
    public:
virtual int ClassId() const override;
 
};

#endif
