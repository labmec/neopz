#ifndef MATPOISSON3DHPP
#define MATPOISSON3DHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"

//#include "pzmanvector.h"


class TPZMatPoisson3d : public TPZMaterial {
  
  TPZFMatrix fXf;//fonte
  int fDim;
  
  public :
    
    static int problema;
  
  TPZMatPoisson3d(int nummat, int dim);
  
  virtual ~TPZMatPoisson3d();
  
  void SetMaterial(TPZFMatrix &xfin){
    fXf = xfin;
  }
  
  int Dimension() { return fDim;}
  
  int NStateVariables();
  
  virtual void Print(ostream & out);
  
  char *Name() { return "TPZMatPoisson3d"; }
  
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  
  
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
  
  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  virtual int NFluxes(){ return 3;}
  
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
  
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif
