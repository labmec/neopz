#ifndef MATTESTHPP
#define MATTESTHPP


#include "pzmaterial.h"
#include "pzfmatrix.h"

class TPZMaterialTest : public TPZMaterial {

	REAL fNumMat;//material id
	REAL fAlfa;//intensity of singularity
	REAL fX0;//singularity
   TPZFMatrix fXf;//fonte

public :

TPZMaterialTest(int nummat, REAL alfa, REAL x0);

virtual ~TPZMaterialTest();

REAL Alfa() {return fAlfa;}

REAL   X0() {return fX0;}

void SetMaterial(TPZFMatrix &xkin){
   fXf = xkin;
}

int Dimension() { return 2;}

int NStateVariables();

virtual void Print(ostream & out);

char *Name() { return "TPZMaterialTest"; }

virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);


virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

virtual int VariableIndex(char *name);

virtual int NSolutionVariables(int var);

virtual int NFluxes(){ return 2;}

virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);

/**compute the value of the flux function to be used by ZZ error estimator*/
virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
				  TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		        TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif
