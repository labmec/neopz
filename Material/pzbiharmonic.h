//$Id: pzbiharmonic.h,v 1.1 2003-11-27 16:02:51 igor Exp $

#ifndef  TPZBIHARMONICHPP
#define TPZBIHARMONICHPP

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * This class implements discontinuous Galerkin formulation for the bi-harmonic equation.
 * @since Nov 27, 2003
 * @author Igor Mozolevski 
 */
class TPZBiharmonic : public TPZDiscontinuousGalerkin {

 private:
  TPZFMatrix fXf, fXk;

  //Problem dimention
  int fDim;
  
  public :

 static REAL gLambda1, gLambda2, gSigma, gL_alpha, gM_alpha, gL_betta, gM_betta;

 TPZBiharmonic(int nummat, int dim);
  
  virtual ~TPZBiharmonic();
  //  dimension of f and k is one!!!!  
  void SetMaterial(TPZFMatrix &xfin,TPZFMatrix  &xkin){
    fXf = xfin;
    fXk =xkin;
  }
  
  int Dimension() { return fDim;}

  // Returns one because of scalar problem 
  int NStateVariables(){
    return 1;
  };
  
  virtual void Print(ostream & out);
  
  char *Name() { return "TPZBiharmonic"; }
  
  //Implements integral over  element's volume
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
			  TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef);
  // Implements boundary conditions for continuous Galerkin
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
			    TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  //<!> ????  
virtual int NFluxes(){ return 3;}
  
 virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
  
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef);
  
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
			    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

};

#endif
