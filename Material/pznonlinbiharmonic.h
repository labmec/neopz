// -*- c++ -*-
//$Id: pznonlinbiharmonic.h,v 1.1 2005-02-04 12:21:30 paulo Exp $

#ifndef TPZNONLINBIHARMONICHPP
#define TPZNONLINBIHARMONICHPP

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * This class implements discontinuous Galerkin formulation for the non-linear bi-harmonic equation.
 * @since Jan 31, 2005
 * @author Igor Mozolevski e Paulo Bosing
 */
class TPZNonLinBiharmonic : public TPZDiscontinuousGalerkin {

private:
  REAL  fXf;
  
  //Problem dimention
  
public :

  static REAL gLambda1, gLambda2, gSigmaA,gSigmaB, gL_alpha, gM_alpha, gL_betta, 
              gM_betta, g_teta, Re;

  /**
   * Inicialisation of biharmonic material
   */
  TPZNonLinBiharmonic(int nummat, REAL f);
  
  virtual ~TPZNonLinBiharmonic();
 
  /**
   * Returns the number of norm errors. Default is 3: energy, L2,  H1, semi-norm H2 and H2.
   */
  virtual int NEvalErrors() {return 5;}

  void SetMaterial(REAL &xfin){
    fXf = xfin;
  }
  
  int Dimension() { return 2;}

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
  
  virtual int NFluxes(){ return 0;}
  
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

  virtual void ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
				   TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize);
  
  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
				     TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
				     TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize);

};

#endif
