//$Id: pzconslaw.h,v 1.10 2004-01-21 00:23:32 erick Exp $

#ifndef PZCONSLAW_H
#define PZCONSLAW_H

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzdiscgal.h"


enum TPZTimeDiscr
{
   None_TD = -1,
   Explicit_TD = 0,
   ApproxImplicit_TD = 1,
   Implicit_TD = 2,
   Unknown_TD = 3
};

enum TPZContributeTime
{
   None_CT = -1,
   Last_CT = 0,
   Advanced_CT = 1
};

enum TPZResidualType
{
   None_RT = -1,
   Residual_RT = 0,
   Flux_RT = 1
};

class TPZConservationLaw2  : public TPZDiscontinuousGalerkin
{
public:

  TPZConservationLaw2(int nummat,REAL timeStep,int dim);

  virtual ~TPZConservationLaw2();

//------------------attributes and parameters

  /**
   * Returns the dimension of the problem
   */
  int Dimension();

  /**
   * Returns the value of the time step
   */
  REAL TimeStep();

  /**
   * Sets the time step used for time integration
   *
   * @param timeStep [in]
   */
  void SetTimeStep(REAL timeStep);

  /**
   * Sets the time step used for time integration
   *
   * @param maxveloc [in] maximal speed of flow inside the cell
   * @param deltax [in] greatest dimension
   * @param degree [in] interpolation degree
   */
  virtual void SetTimeStep(REAL maxveloc,REAL deltax,int degree)=0;

  /**
   * Returns the CFL number
   */
  REAL CFL();

  /**
   * Sets the CFL number
   *
   * @param CFL [in]
   */
  void SetCFL(REAL CFL);

  /**
   * Returns the value of delta for the artificial diffusion term
   */
//  REAL Delta();

  /**
   * Sets the value of the artificial diffusion delta parameter
   *
   * @param delta [in]
   */
//  void SetDelta(REAL delta);


  /**
   * Returns the value of Gamma (constant of gas)
   */
  REAL Gamma();

  /**
   * Sets the value of Gamma (constant of gas)
   *
   * @param gamma [in]
   */
  void SetGamma(int gamma);

  /**
   * Sets whether the contribution is advanced or referring to the last
   * state.
   *
   */
  void SetContributionTime(TPZContributeTime time);

  /**
   * Residual_RT for calculations and Flux_RT for convergence check.
   */
  void TPZConservationLaw2::SetResidualType(TPZResidualType type);

  /**
   * Number of state variables according to the dimension
   */
  virtual int NStateVariables() = 0;

  /**
   * termodynamic pressure determined by the law of an ideal gas
   *
   * @param U [in] vector of state variables (sol)
   */
  virtual REAL Pressure(TPZVec<REAL> &U)=0;

  /**
   * Prints the state of internal variables
   *
   * @param out [in]
   */
  virtual void Print(ostream & out);

  /**
   * Returns the material name
   */
  virtual char *Name()=0;

  /**
   * Returns the relative index of a variable according to its name
   *
   * @param name [in]
   */
  virtual int VariableIndex(char *name)=0;

  virtual int NSolutionVariables(int var)=0;

  /**
   * Returns the number of fluxes associated to this material
   */
  virtual int NFluxes();


//------------------solutions
  /**
   * compute the boundary condition left solution
   */
/*  virtual void ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,
			TPZVec<REAL> &normal,TPZBndCond *bcleft) = 0;
*/
  /**
   * compute the boundary condition right solution
   */
/*  virtual void ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,
			TPZVec<REAL> &normal,TPZBndCond *bcright)=0;
*/

  virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,
			TPZFMatrix &axes,int var,
			TPZVec<REAL> &Solout)=0;

  /**
   * computes the value of the flux function to be used
   * by ZZ error estimator
   */
  virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol,
			TPZFMatrix &DSol, TPZFMatrix &axes,
			TPZVec<REAL> &flux);

  virtual void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
	      TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &flux,
	      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,
	      TPZVec<REAL> &values);


//------------------contributions

  /**
   * Contributes to the residual vector and tangent matrix the
   * volume-based quantities.
   */
  virtual void Contribute(TPZVec<REAL> &x,TPZFMatrix &,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight, TPZFMatrix &axes,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)=0;

  /**
   * Contributes to the residual vector and tangent matrix the
   * face-based quantities.
   */
  virtual void ContributeInterface(TPZVec<REAL> &x,
			TPZVec<REAL> &solL,TPZVec<REAL> &solR,
			TPZFMatrix &dsolL,TPZFMatrix &dsolR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &dphiL,TPZFMatrix &dphiR,
			TPZFMatrix &ek,TPZFMatrix &ef)=0;
  /**
   * Contributes to the residual vector the boundary conditions
   *
   */
  virtual void ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,
			REAL weight, TPZFMatrix &axes,
			TPZFMatrix &phi,
			TPZFMatrix &ek,TPZFMatrix &ef,
			TPZBndCond &bc)=0;

//--------------------
  //virtual int IntegrationDegree() = 0;

  //virtual void SetIntegDegree(int degree) = 0;

//---------------------Attributes
protected:

  /**
   * Dimension of the problem
   */
  int fDim;

  /**
   * Time step used for time integration
   */
  REAL fTimeStep;

  /**
   * CFL number
   */
  REAL fCFL;

  /**
   * Delta coefficient related to the Artificial diffusive term
   */
//  REAL fDelta;

  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  REAL fGamma;

  /**
   * Variable indicating the context of the solution.
   * If advanced, then the implicit terms are to be
   * contributed. If last, then the explicit.
   */
  TPZContributeTime fContributionTime;

  /**
   * Variable to indicate the type of residual to be computed by Assemble.
   * A Flux Evaluation type is interesting for residual evaluation, and
   * a complete residual for the global invertion.
   */
  TPZResidualType fResidualType;

};

/*
inline void TPZConservationLaw2::SetDelta(REAL delta)
{
   fDelta = delta;
}

inline REAL TPZConservationLaw2::Delta()
{
   return fDelta;
}
*/
inline int TPZConservationLaw2::Dimension()
{
   return fDim;
}

inline REAL TPZConservationLaw2::CFL()
{
   return fCFL;
}

inline void TPZConservationLaw2::SetCFL(REAL CFL)
{
   fCFL = CFL;
}

inline void TPZConservationLaw2::SetGamma(int gamma)
{
   fGamma = gamma;
}

inline REAL TPZConservationLaw2::Gamma()
{
   return fGamma;
}

inline void TPZConservationLaw2::SetTimeStep(REAL timeStep)
{
   fTimeStep = timeStep;
}

inline REAL TPZConservationLaw2::TimeStep()
{
   if(fResidualType == Residual_RT)return fTimeStep;
   return 1.;
}
/*
inline char * TPZConservationLaw2::Name()
{
   return "TPZConservationLaw2";
}
*/
inline int TPZConservationLaw2::NFluxes()
{
   return 1;
}

inline void TPZConservationLaw2::SetContributionTime(TPZContributeTime time)
{
   fContributionTime = time;
}


inline void TPZConservationLaw2::SetResidualType(TPZResidualType type)
{
   fResidualType = type;
   /*
   if(fResidualType == Flux_RT)cout << "Flux Residual Type\n";
   if(fResidualType == Residual_RT)cout << "Residual Residual Type\n";*/
}

#endif
