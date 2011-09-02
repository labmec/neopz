/**
 * \file
 * @brief Contains the TPZMatHyperElastic class which implements a hyper elasticity material.
 */
#ifndef MATHYPERELASTICHPP
#define MATHYPERELASTICHPP

#include "pzmaterial.h"
#include "pzfmatrix.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/**
 * @ingroup material
 * @brief Implements a hyper elasticity material.
 */
class TPZMatHyperElastic : public TPZMaterial {
	
	REAL fK2[3][3],fK3[3][3],fK4[3][3],fK6[3][3],fK7[3][3],fK8[3][3],fXf[3];
	REAL fL1[3][3],fL2[3][3],fL3[3][3],fL4[3][3],fL5[3][3],fL6[3][3],fL7[3][3],fL8[3][3],fL9[3][3],fGradtrC[3][3];
	REAL fE1[3],fE5[3],fE9[3],fGradDetF[3][3];
	REAL fLambda,fNu,fE,fMu;
	REAL fCoef1,fCoef2,fCoef3;
	
	public :
	
	TPZMatHyperElastic(int nummat,REAL e,REAL mu,REAL nu=-1.,REAL lambda=-1.,REAL coef1=-1.,REAL coef2=-1.,REAL coef3=-1.);
	
	virtual ~TPZMatHyperElastic();
	
	void SetMaterial(TPZFMatrix &xfin){
		fXf[0] = xfin(0,0);
		fXf[1] = xfin(1,0);
		fXf[2] = xfin(2,0);
	}
	
	int Dimension() { return 3;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	std::string Name() { return "TPZMatHyperElastic"; }
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef);
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef, TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix &ef, TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
#ifdef _AUTODIFF
	
	/** @brief Computes contribution to the energy at an integration point*/
	virtual void ContributeEnergy(TPZVec<REAL> &x,
								  TPZVec<FADFADREAL> &sol,
								  TPZVec<FADFADREAL> &dsol,
								  FADFADREAL &U,
								  REAL weight);
	
	static void ComputeEnergy(REAL lambda, REAL mu,  TPZFMatrix &dsol, TFad<9,TFad<9,REAL> > &energy);
	
	/** @brief Computes contribution of BC to the Energy*/
	virtual void ContributeBCEnergy(TPZVec<REAL> & x,
									TPZVec<FADFADREAL> & sol, FADFADREAL &U,
									REAL weight, TPZBndCond &bc);
	
#endif
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 9;}
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif
