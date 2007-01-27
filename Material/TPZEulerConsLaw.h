/*
  Dissipa¢ão
  substantivo feminino
  1. acto ou efeito de dissipar ou dissipar-se;
  2. desaparecimento; desvanecimento;
  3. desperdício de meios; gasto exagerado de dinheiro;
  4. devassidão;
  (Do lat. dissipatióne-, «id.»)
*/

#ifndef EULERCONSLAWHPP
#define EULERCONSLAWHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "TPZConservationLaw.h"

class TPZEulerConsLaw  : public TPZConservationLaw {

  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas 
   */
  REAL fGamma;
  
  /**
   * Termo que adiciona estabilidade ao método numérico de aproximação
   * SUPG
   * LS
   * Bornhaus
   */
  char *fArtificialDiffusion;

  //int fIntegrationDegree;//grau de integra¢ão da solu¢ão inicial:opcional
  
  public :
  
  TPZEulerConsLaw(int nummat,REAL delta_t,REAL gamma,int dim,char *artdiff);

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
  
  char *Name() { return "TPZEulerConsLaw"; }
  
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,TPZFMatrix &axes,
			    TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ef);
  
  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  virtual int NFluxes(){ return Dimension();}
  
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
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
