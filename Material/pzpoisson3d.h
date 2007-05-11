// -*- c++ -*-

//$Id: pzpoisson3d.h,v 1.25 2007-05-11 19:15:18 joao Exp $

#ifndef MATPOISSON3DH
#define MATPOISSON3DH

#include <iostream>
#include "pzdiscgal.h"
#include "pzfmatrix.h"

//#include "pzmanvector.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

// -fK Laplac(u) + fC * div(fConvDir*u) = - fXf
class TPZMatPoisson3d : public TPZDiscontinuousGalerkin {
 
  protected :

  /** Forcing function value */
  REAL fXf;
  
  /** Problem dimension */
  int fDim;
  
  /** Coeficient which multiplies the Laplacian operator. */
  REAL fK;
  
  /** Coeficient which multiplies the Laplacian operator associated to right neighbour element of the interface.
   * The coefficient of left neighbour is fK.
   * It is for use with discontinuous Galerkin method. Default value for fRightK is fK.
   */
  REAL fRightK;
  
  /** Variable which multiplies the convection term of the equation */
  REAL fC;
  
  /** Direction of the convection operator */
  REAL fConvDir[3];
  
  /** Symmetry coefficient of elliptic term.
   * Symmetrical formulation - Global element method - has coefficient = -1.
   * Non-symmetrical formulation - Baumann's formulation - has coefficient = +1.
   */
  REAL fSymmetry;
  
  /**
   * multiplication value for the streamline diffusion term
   */
   REAL fSD;

public:

  /** Constant multiplyer of penalty term, when required.
   * Penalty terms are not used always. It's used when TPZInterfaceElement::CalcStiffPenalty 
   * is set.
   */
  REAL fPenaltyConstant;

  /** Usado em InterfaceErrors */
  static REAL gAlfa;
  
  TPZMatPoisson3d(int nummat, int dim);

  virtual ~TPZMatPoisson3d();

  TPZMatPoisson3d(const TPZMatPoisson3d &copy) : TPZDiscontinuousGalerkin(copy){
    fXf  = copy.fXf;
    fDim = copy.fDim;
    fK   = copy.fK;
    this->fRightK = copy.fRightK;
    fC   = copy.fC;
    for (int i = 0; i < 3; i++) fConvDir[i] = copy.fConvDir[i];
    fSymmetry = copy.fSymmetry;
    fSD = copy.fSD;
    fPenaltyConstant = copy.fPenaltyConstant;

  }
  
  /** Set material elliptic term as the global element method, i.e. the symmetrical formulation.
  */
  void SetSymmetric(){
    this->fSymmetry = -1.0;
  }
  
  /** Set material elliptic term as the Baumann's formulation, i.e. the non-symmetrical formulation.
  */  
  void SetNonSymmetric() {
    this->fSymmetry = +1.0; 
  }
  
  bool IsSymetric(){
    if (fSymmetry == -1.0) return true;
    if (fSymmetry == +1.0) return false;
    PZError << __PRETTY_FUNCTION__ << "\n Comparacao de numeros reais da errado\n";
    return false;
  }

  virtual TPZAutoPointer<TPZMaterial> NewMaterial(){
    return new TPZMatPoisson3d(*this);
  }
    
  int Dimension() { return fDim;}

  int NStateVariables();

  void SetParameters(REAL diff, REAL conv, TPZVec<REAL> &convdir);
  
  void GetParameters(REAL &diff, REAL &conv, TPZVec<REAL> &convdir);

  void SetInternalFlux(REAL flux)
  {
    fXf = flux;
  }
  
  void SetSD(REAL sd)
  {
    fSD = sd;
  }
  
  /** Define fK of right neighbour element. It is used with discontinuous Galerkin
   * on the computation of ContributeInterface methods.
   * Attention that method SetParameters override the modifications of this method. Then call it after SetParameters
   * and never before or it will have no effect.
   */
  void SetRightK(REAL rightK){
    this->fRightK = rightK;
  }
  REAL GetRightK(){
    return this->fRightK;    
  }
  
  virtual void Print(std::ostream & out);
  
  char *Name() { return "TPZMatPoisson3d"; }
  
  virtual void Contribute(TPZMaterialData &data,REAL weight,
			  TPZFMatrix &ek,TPZFMatrix &ef);
#ifdef _AUTODIFF
  /**Compute contribution to the energy at an integration point*/
  void ContributeEnergy(TPZVec<REAL> &x,
			      TPZVec<FADFADREAL> &sol,
			      TPZVec<FADFADREAL> &dsol,
			      FADFADREAL &U,
			      REAL weight);
#endif

  virtual void ContributeBC(TPZMaterialData &data,REAL weight,
			    TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

#ifdef _AUTODIFF

  virtual void ContributeBCEnergy(TPZVec<REAL> & x,
	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
	REAL weight, TPZBndCond &bc);

#endif

  virtual int VariableIndex(char *name);
  
  virtual int NSolutionVariables(int var);
  
  virtual int NFluxes(){ return 3;}
  
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
  
  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);
  
  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);


  virtual void ContributeBCInterface(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
			    
  virtual void ContributeInterface(TPZMaterialData &data,REAL weight,
				   TPZFMatrix &ek,TPZFMatrix &ef);


  void InterfaceErrors(TPZVec<REAL> &/*x*/,
		       TPZVec<REAL> &leftu, TPZFMatrix &leftdudx, /* TPZFMatrix &leftaxes,*/ 
		       TPZVec<REAL> &rightu, TPZFMatrix &rightdudx, /* TPZFMatrix &rightaxes,*/ 
		       TPZVec<REAL> &/*flux*/,
		       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values, 
		       TPZVec<REAL> normal, REAL elsize);
                       
  /** Compute interface jumps 
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[0] + values[1]
   * @since Feb 14, 2006
   */  
  virtual void InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu, TPZVec<REAL> &leftNormalDeriv,
                              TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                              TPZVec<REAL> &values);
                              
  /** Compute interface jumps from element to Dirichlet boundary condition
   * values[1] = (solleft - solright)^2
   * values[2] = (dsolleft - dsolright)^2
   * values[0] = values[0] + values[1]
   * @since Mar 08, 2006
   */  
  virtual void BCInterfaceJumps(TPZVec<REAL> &leftu,
                                TPZBndCond &bc,
                                TPZVec<REAL> &values);                              

  virtual int IsInterfaceConservative(){ return 1;}

};

#endif
