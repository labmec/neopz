// -*- c++ -*-

//$Id: pzpoisson3d.h,v 1.35 2011-05-11 02:22:56 phil Exp $

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

   /** Enumerate for penalty term definitions
    */
   enum EPenaltyType {ENoPenalty = 0, EFluxPenalty = 1, ESolutionPenalty, EBoth};

   /** Penalty term definition
    */
   EPenaltyType fPenaltyType;

public:

  /** Constant multiplyer of penalty term, when required.
   * is set.
   */
  REAL fPenaltyConstant;

  /** Define no penalty terms in ContributeInterface
   */
  void SetNoPenalty(){ this->fPenaltyType = ENoPenalty;}

  /** Define flux penalty terms in ContributeInterface
   */
  void SetFluxPenalty(){ this->fPenaltyType = EFluxPenalty; }

  /** Define solution penalty terms in ContributeInterface
   */
  void SetSolutionPenalty(){ this->fPenaltyType = ESolutionPenalty; }

  /** Define solution and flux penalty terms in ContributeInterface
   */
  void SetBothPenalty(){ this->fPenaltyType = EBoth; }

  /** Usado em InterfaceErrors */
  static REAL gAlfa;

  TPZMatPoisson3d(int nummat, int dim);

  TPZMatPoisson3d();

  TPZMatPoisson3d(const TPZMatPoisson3d &copy);

  virtual ~TPZMatPoisson3d();

 TPZMatPoisson3d &operator=(const TPZMatPoisson3d &copy);

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

  virtual std::string Name() { return "TPZMatPoisson3d"; }

  virtual void Contribute(TPZMaterialData &data,REAL weight,
			  TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void Contribute(TPZMaterialData &data,REAL weight,
			  TPZFMatrix &ef)
			  {
				  TPZDiscontinuousGalerkin::Contribute(data,weight,ef);
              }
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

  virtual void ContributeBC(TPZMaterialData &data,REAL weight,
				TPZFMatrix &ef,TPZBndCond &bc)
				{
					TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
                }

#ifdef _AUTODIFF

  virtual void ContributeBCEnergy(TPZVec<REAL> & x,
	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
	REAL weight, TPZBndCond &bc);

#endif

  virtual int VariableIndex(const std::string &name);

  virtual int NSolutionVariables(int var);

  virtual int NFluxes(){ return 3;}

protected:
  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout);
public:

   /**returns the solution associated with the var index based on
    * the finite element approximation*/
  virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);

  /**compute the value of the flux function to be used by ZZ error estimator*/
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux);

  void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);

  virtual int NEvalErrors() {return 5;}


  virtual void ContributeBCInterface(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);

  virtual void ContributeInterface(TPZMaterialData &data,REAL weight,
				   TPZFMatrix &ek,TPZFMatrix &ef);

  virtual void ContributeBCInterface(TPZMaterialData &data,REAL weight,TPZFMatrix &ef,TPZBndCond &bc)
  {
	  TPZDiscontinuousGalerkin::ContributeBCInterface(data,weight,ef,bc);
  }

  virtual void ContributeInterface(TPZMaterialData &data,REAL weight,
				   TPZFMatrix &ef)
				   {
                       TPZDiscontinuousGalerkin::ContributeInterface(data,weight,ef);
                   }

  /**
   * Compute square of residual of the differential equation at one integration point.
   * @param X is the point coordinate (x,y,z)
   * @param sol is the solution vector
   * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
   */   
  virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<REAL> &sol, TPZFMatrix &dsol);


  void InterfaceErrors(TPZVec<REAL> &/*x*/,
		       TPZVec<REAL> &leftu, TPZFMatrix &leftdudx, /* TPZFMatrix &leftaxes,*/
		       TPZVec<REAL> &rightu, TPZFMatrix &rightdudx, /* TPZFMatrix &rightaxes,*/
		       TPZVec<REAL> &/*flux*/,
		       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values,
		       TPZVec<REAL> normal, REAL elsize);

  /** Compute interface jump from element to Dirichlet boundary condition
   * Returns sol-u_dirichlet
   * @since Mar 08, 2006
   */
  virtual void BCInterfaceJump(TPZVec<REAL> &x, TPZVec<REAL> &leftu,TPZBndCond &bc,TPZVec<REAL> & jump);

  virtual int IsInterfaceConservative(){ return 1;}

  /**
   * Unique identifier for serialization purposes
   */
  virtual int ClassId() const;

  /**
   * Save the element data to a stream
   */
  virtual void Write(TPZStream &buf, int withclassid);

  /**
   * Read the element data from a stream
   */
  virtual void Read(TPZStream &buf, void *context);

};

#endif
