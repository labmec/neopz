#ifndef MAT1DLINHPP
#define MAT1DLINHPP

#include <iostream>

using namespace std;

#include "pzfmatrix.h"
#include "pzmaterial.h"
#include "pzvec.h"

struct TPZElementMatrix;
class TPZBndCond;
template<class T>
class TPZVec;


class TPZMat1dLin : public TPZMaterial{

  TPZFMatrix		fXk,fXc,fXb,fXf;

  public :


    TPZMat1dLin(int num) : TPZMaterial(num) , fXk(), fXc(), fXb(), fXf() {
  }

  virtual int NStateVariables() { return fXk.Rows(); }

  int Dimension() { return 1;}

  void Print(ostream & out);

  void SetMaterial(TPZFMatrix &xkin,TPZFMatrix &xcin,TPZFMatrix &xbin,TPZFMatrix &xfin){
    fXk = xkin;
    fXc = xcin;
    fXb = xbin;
    fXf = xfin;
  }

  virtual char *Name() { return "TPZMat1dLin"; }

  int NFluxes() { return NStateVariables(); }

  void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &,REAL weight,REAL ,
		  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);

  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
		      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
};

#endif
