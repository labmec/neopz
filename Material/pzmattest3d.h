#ifndef MATTEST3DHPP
#define MATTEST3DHPP


#include "pzmaterial.h"
#include "pzfmatrix.h"

//#include "pzmanvector.h"


class TPZMaterialTest3D : public TPZMaterial {

   TPZFMatrix fXf;//fonte

public :

static int eq3;//Cedric : para testes no programa main 3dmaterial.c

TPZMaterialTest3D(int nummat);

virtual ~TPZMaterialTest3D();

void SetMaterial(TPZFMatrix &xfin){
   fXf = xfin;
}

int Dimension() { return 3;}

int NStateVariables();

virtual void Print(std::ostream & out);

char *Name() { return "TPZMaterialTest3D"; }

virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);


virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
  
  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  virtual int NFluxes(){ return 3;}
  
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();
  
  //virtual void Clone(TPZAdmChunkVector<TPZMaterial *> &matvec);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
  
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric
};

#endif
