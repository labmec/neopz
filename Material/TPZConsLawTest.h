#ifndef CONSLAWTESTHPP
#define CONSLAWTESTHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "TPZConservationLaw.h"

class TPZConsLawTest  : public TPZConservationLaw {
  
  TPZFMatrix fXf;//fonte
  TPZVec<REAL> fB;
  int fArtificialDiffusion;
  //int fIntegrationDegree;//integra¢ão da solu¢ão inicial:opcional
  int fTest;
  
  public :
  
  TPZConsLawTest(int nummat, TPZVec<REAL> B,int artdiff,REAL delta_t,int dim,REAL delta,int test=0);
  
  virtual ~TPZConsLawTest();
  
  void SetMaterial(TPZFMatrix &xfin){
    fXf = xfin;
  }

  REAL DeltaOtimo();

  REAL CFL(int degree);

  REAL B(int i,TPZVec<REAL> &x);

  //  virtual void SetIntegDegree(int degree){fIntegrationDegree = degree;}

  //virtual int IntegrationDegree(){return fIntegrationDegree;}
  
  REAL Delta();

  REAL T(int jn,TPZVec<REAL> &x);
  
  int NStateVariables();
  
  virtual void Print(ostream & out);
  
  char *Name() { return "TPZConsLawTest"; }
  
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
  
  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  virtual int NFluxes(){ return 1;}
  
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
  
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);

  void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
  void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);
};

#endif
