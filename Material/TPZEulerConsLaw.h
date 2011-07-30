/**
 * \file
 * @brief Contains the TPZEulerConsLaw class which implements a Euler equation as conservation law formulation.
 */
/*
 Dissipa�
 substantivo feminino
 1. acto ou efeito de dissipar ou dissipar-se;
 2. desaparecimento; desvanecimento;
 3. desperd�io de meios; gasto exagerado de dinheiro;
 4. devassid�;
 (Do lat. dissipati�e-, id.)
 */

#ifndef EULERCONSLAWHPP
#define EULERCONSLAWHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "TPZConservationLaw.h"

/**
 * @ingroup material
 * @brief 
 */
class TPZEulerConsLaw  : public TPZConservationLaw {
	
	/**
	 * ratio between specific heat is constant and the specific heat the constant
	 * volume of a polytropic gas
	 */
	REAL fGamma;
	
	/**
	 * Termo que adiciona estabilidade ao m�odo num�ico de aproxima�o
	 * SUPG
	 * LS
	 * Bornhaus
	 */
	std::string fArtificialDiffusion;
	
	//int fIntegrationDegree;//grau de integra� da solu� inicial:opcional
	
	public :
	
	TPZEulerConsLaw(int nummat,REAL delta_t,REAL gamma,int dim,const std::string &artdiff);
	
	/**copy constructor*/
	TPZEulerConsLaw(TPZEulerConsLaw & copy);
	
	/**To create another material of the same type*/
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	~TPZEulerConsLaw();
	
	/**
	 * compute the boundary condition left solution
	 */
	virtual void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);
	
	/**
	 * compute the boundary condition right solution
	 */
	virtual void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
	
	/**
	 * termodinamic pressure determined by the law of an ideal gas
	 */
	virtual REAL Pressure(TPZVec<REAL> &U);
	
	virtual REAL Gamma(){return fGamma;}
	
	virtual void SetDeltaTime(REAL maxveloc,REAL deltax,int degree);
	
	/**
	 * tensor of the three-dimensional flux of Euler
	 */
	void Flux(TPZVec<REAL> &U,TPZVec<REAL> &Fx,TPZVec<REAL> &Fy,TPZVec<REAL> &Fz);
	
	//virtual void SetIntegDegree(int degree){fIntegrationDegree = degree;}
	
	//virtual int IntegrationDegree(){return fIntegrationDegree;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZEulerConsLaw"; }
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef);
	
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
		TPZConservationLaw::ContributeBC(data,weight,ef,bc);
	}
	
	virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix &ef);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return Dimension();}
protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<REAL> &Solout)
	{
		Solution(data.sol,data.dsol,data.axes,var,Solout);
	}
	
	/**compute the value of the flux function to be used by ZZ error estimator*/
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes,
					  TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
	// PARA TESTE PARA TESTE PARA PARA TESTE
	void ContributeTESTE(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
						 REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
						 TPZFMatrix &ek,TPZFMatrix &ef);// PARA TESTE PARA TESTE PARA PARA TESTE
	
	void TestOfRoeFlux(REAL &tetainit,REAL &tetamax,REAL &tol,REAL &increment);
};

#endif
