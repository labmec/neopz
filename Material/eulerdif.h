/**
 * \file
 * @brief Contains the TEulerDiffusivity class which implements a numerical diffusivity coefficient for SUPG.
 */
#ifndef EULERDIFFUSIONHH
#define EULERDIFFUSIONHH

#include "pzvec.h" 
#include "pzfmatrix.h"

/** 
 * @brief Class which implements a numerical diffusivity coeficient for the SUPG method.
 * @ingroup analysis
 */
/**
 * This class is to be used as a template argument for a different class
 */
class TEulerDiffusivity {

static  REAL  fGamma;

 public:

static  REAL Pressure(TPZVec<REAL> &U);

static  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
static  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix &Ajacob,TPZFMatrix &Bjacob);
static	void JacobFlux(TPZVec<REAL> &U,TPZFMatrix &jacob,TPZVec<REAL> &normal);

static  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix &valjacob,TPZVec<REAL> &normal);

static  void MatrixDiff(TPZVec<REAL> &sol, TPZFMatrix &axes, TPZFMatrix &jacinv,TPZFMatrix
&ATauA,TPZFMatrix &ATauB,TPZFMatrix &BTauA,TPZFMatrix &BTauB);

static  void InvJacob2d(TPZFMatrix &axes,TPZFMatrix &jacinv);

static  void InverseJacob(TPZFMatrix &jac);

static int main();

};

#endif
