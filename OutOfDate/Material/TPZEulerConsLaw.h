/**
 * \file
 * @brief DEPRECATED FILE. This file contains the TPZEulerConsLaw class which implements the weak statement \n
 * of the compressible euler equations as conservation law formulation.
 */
/*
 Dissipacao
 substantivo feminino
 1. acto ou efeito de dissipar ou dissipar-se;
 2. desaparecimento; desvanecimento;
 3. desperdizio de meios; gasto exagerado de dinheiro;
 4. devassidao;
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
 * @deprecated DEPRECATED Euler equation material CLASS.
 * @brief Implements the weak statement of the compressible euler equations.
 */
class TPZEulerConsLawDEP  : public TPZConservationLawDEP {
	
	/** @brief Ratio between specific heat is constant and the specific heat the constant volume of a polytropic gas */
	REAL fGamma;
	
	/** @brief Term that adds stability to the numerical method of approach: SUPG, LS, Bornhaus */
	std::string fArtificialDiffusion;
		
	public :
	
	TPZEulerConsLawDEP(int nummat,REAL delta_t,REAL gamma,int dim,const std::string &artdiff);
	
	/** @brief Copy constructor */
	TPZEulerConsLawDEP(TPZEulerConsLawDEP & copy);
	
	/** @brief To create another material of the same type */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	~TPZEulerConsLawDEP();
	
	/** @brief Computes the boundary condition left solution */
	virtual void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);
	
	/** @brief Computes the boundary condition right solution */
	virtual void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
	
	/** @brief Termodinamic pressure determined by the law of an ideal gas */
	virtual REAL Pressure(TPZVec<REAL> &U);
	
	virtual REAL Gamma(){return fGamma;}
	
	virtual void SetDeltaTime(REAL maxveloc,REAL deltax,int degree);
	
	/** @brief tensor of the three-dimensional flux of Euler */
	void Flux(TPZVec<REAL> &U,TPZVec<REAL> &Fx,TPZVec<REAL> &Fy,TPZVec<REAL> &Fz);
	
	//virtual void SetIntegDegree(int degree){fIntegrationDegree = degree;}
	
	//virtual int IntegrationDegree(){return fIntegrationDegree;}
	
	int NStateVariables();
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZEulerConsLawDEP"; }
	
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
		TPZConservationLawDEP::ContributeBC(data,weight,ef,bc);
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
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes,
					  TPZVec<REAL> &flux);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
	void ContributeTESTE(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
						 REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
						 TPZFMatrix &ek,TPZFMatrix &ef);
	
	void TestOfRoeFlux(REAL &tetainit,REAL &tetamax,REAL &tol,REAL &increment);
};

#endif
