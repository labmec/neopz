#ifndef CONSERVATIONLAWHPP
#define CONSERVATIONLAWHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"


class TPZConservationLaw  : public TPZMaterial {

  int fDim;
  REAL fTimeStep;
  
  public :
  
  TPZConservationLaw(int nummat,REAL delta_t,int dim);
  
  virtual ~TPZConservationLaw();
  
  REAL Delta();

  REAL TimeStep(){return fTimeStep;}

  int Dimension() { return fDim;}
  
  int NStateVariables();
  
  virtual void Print(ostream & out);
  
  char *Name() { return "TPZConservationLaw"; }
  
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZFMatrix &axes,TPZFMatrix &phiL,TPZFMatrix &phiR,
				   TPZFMatrix &dphiL,TPZFMatrix &dphiR,TPZFMatrix &ek,TPZFMatrix &ef) = 0;
  
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
};

#endif
