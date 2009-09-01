// -*- c++ -*-
//$Id: pzbiharmonic.h,v 1.11 2009-09-01 19:44:46 phil Exp $

#ifndef  TPZBIHARMONICHPP
#define TPZBIHARMONICHPP

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

/**
 * This class implements discontinuous Galerkin formulation for the bi-harmonic equation.
 * @since Nov 27, 2003
 * @author Igor Mozolevski e Paulo Bosing
 */
class TPZBiharmonic : public TPZDiscontinuousGalerkin {

private:
  REAL  fXf;

  //Problem dimention

public :

  static REAL gLambda1, gLambda2, gSigmaA,gSigmaB, gL_alpha, gM_alpha, gL_betta, gM_betta;

  /**
   * Inicialisation of biharmonic material
   */
  TPZBiharmonic(int nummat, REAL f);

  virtual ~TPZBiharmonic();

  /**
   * Returns the number of norm errors. Default is 3: energy, L2,  H1, semi-norm H2 and H2.
   */
  virtual int NEvalErrors() {return 8;}

  void SetMaterial(REAL &xfin){
    fXf = xfin;
  }

  int Dimension() { return 2;}

  // Returns one because of scalar problem
  int NStateVariables(){
    return 1;
  };

  virtual void Print(std::ostream & out);

  virtual std::string Name() { return "TPZBiharmonic"; }

  //Implements integral over  element's volume
  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
			    TPZFMatrix &ek,
                            TPZFMatrix &ef);
  //Implements integral over  element's volume
  virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
							TPZFMatrix &ef)
  {
	   TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
  }
  // Implements boundary conditions for continuous Galerkin
  virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
				  TPZFMatrix &ek,
							  TPZFMatrix &ef,
							  TPZBndCond &bc);

 // Implements boundary conditions for continuous Galerkin
  virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix &ef,
							  TPZBndCond &bc)
  {
	   TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
  }

  virtual int VariableIndex(const std::string &name);

  virtual int NSolutionVariables(int var);

  virtual int NFluxes(){ return 0;}

protected:
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:
      /**returns the solution associated with the var index based on
       * the finite element approximation*/
	  virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	  {
			TPZDiscontinuousGalerkin::Solution(data,var,Solout);
      }


  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);


  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);//Cedric


  virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ek,
									 TPZFMatrix &ef);


  virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc);

  virtual void ContributeInterface(TPZMaterialData &data,
									 REAL weight,
									 TPZFMatrix &ef)
  {
		TPZDiscontinuousGalerkin::ContributeInterface(data,weight,ef);
  }


  virtual void ContributeBCInterface(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
  {
      TPZDiscontinuousGalerkin::ContributeBCInterface(data,weight,ef,bc);
  }


};

#endif
