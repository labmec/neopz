#ifndef HIBRIDAH
#define HIBRIDAH


#include "pzmaterial.h"
#include "pzfmatrix.h"

class TPZMatHybrid : public TPZMaterial {

  REAL fNumMat;//material id
  TPZFMatrix fXf;//fonte

public :

TPZMatHybrid(int nummat);

virtual ~TPZMatHybrid();

void SetMaterial(TPZFMatrix &xkin){
   fXf = xkin;
}

int Dimension() { return fXf.Rows();}

int NStateVariables();

virtual void Print(ostream & out);

char *Name() { return "TPZMatHybrid"; }

virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);


virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

virtual int VariableIndex(char *name);

virtual int NSolutionVariables(int var);

virtual int NFluxes(){ return 3;}

virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

/**compute the value of the flux function to be used by ZZ error estimator*/
 virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

 void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	    TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};

#endif
