#ifndef PZARTDIFF_H
#define PZARTDIFF_H


#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzdiffmatrix.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/**
 * This classroom adds to the term of diffusion to the variacional formulation
 * of the differential equation partial compressible of Euler (hyperbolic). This
 * term is introduced to stabilize the numerical method of approach.
 *
 * @author Erick Slis
 * @author Cedric Ayala
 */
class TPZArtDiff {

 public:
  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  static REAL fGamma;

  /**
   * scalar coefficient of the element in the diffusion term
   */
  static REAL fDelta;

  /**
   * parámetro that it limits the condition of stability of the numerical approach
   */
  static REAL fCFL;

  static REAL OptimalDelta();

  static REAL Delta();

  static REAL CFL(int degree);

  /**
   * Jacobian of the tensor flux of Euler
   * @param dim [in]
   * @param U [in] dim+2 solutions at given point
   * @param Ai [out] vector of dim tensors of (dim+2)*(dim*2),
   * representing the derivatives of F with respect to the
   * dim spatial dimensions
   */

  template <class T>
  static void JacobFlux(int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai);

  /**
   * operation product point in the diffusion term
   * @param dphi [in] vector of scalars that make part of the ODot operation
   * @param M [in] vector of matrices to operate
   * @param Result [out]
   *
   */

  template <class T>
  static void ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZDiffMatrix<T> > &M, TPZDiffMatrix<T> &Result);

  /**
   * Computes the common values A B C and Tau vector of matrixes
   * for contributions
   *
   * @param dim [in]
   * @param U [in] vector of solutions at the point
   * @param Ai [out] matrixes A B C
   * @param Tau [out] diffusive vector of matrixes
   */
  template <class T>
  static void PrepareDiff(int dim, char * ArtDiff, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  /**
   * Computes the value of the diffusive parcell at a given
   * solution U with respect to the nthphi test funtion
   * @param dim
   * @param Ai vector of Matrixes A B C
   * @param Tau vector of tensors
   * @param U values of dim+2 solutions at the given point
   * @param dphi matrix of derivatives of shape functions
   * @param nthphi index of test function to use
   * @param diff diffusive term value
   * @param approxDiffDer an approximation to the diff term derivative
   */
/*  static void Diff(int dim,
		 TPZVec<TPZDiffMatrix<REAL> > & Ai,
		 TPZVec<TPZDiffMatrix<REAL> > & Tau, TPZVec<REAL> &U,
		 TPZFMatrix & dphi, const int nthphi,
		 TPZVec<REAL> & diff,
		 TPZDiffMatrix<REAL> &approxDiffDer);
*/
#ifdef _AUTODIFF
  /**
   * Computes the value of the diffusive parcell at a given
   * solution U with respect to the nthphi test funtion
   * @param dim
   * @param Ai vector of Matrixes A B C
   * @param Tau vector of tensors
   * @param U values of dim+2 solutions at the given point (with derivatives)
   * @param dU values of (dim+2)*dim derivatives of solutions at the given point (with derivatives)
   * @param nthphi index of test function to use
   * @param dphi matrix of derivatives of shape functions
   * @param diff diffusive term value (with exact derivatives
   */
/*  static void Diff(int dim,
		TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
		TPZVec<TPZDiffMatrix<FADREAL> > & Tau,
		TPZVec<FADREAL> &U,
		TPZVec<FADREAL> &dU, int nthphi,
		TPZVec<FADREAL> &diff);
*/
#endif

  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param ArtDiff [in] Name of the artificial diffusion to apply
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div
   * @param pTaudDiv [out] pointer of a vector to matrices to store the
   * approximate derivatives of Tau_i*Div. If Null, then no approximation
   * derivative is evaluated.
   *
   */

static void PrepareFastDiff(int dim, char * ArtDiff, TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv,
		 TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv = NULL);

#ifdef _AUTODIFF
  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param ArtDiff [in] Name of the artificial diffusion to apply
   * @param sol [in] solution of the dim+2 state functions setup with shape functions
   * @param dsol [in] derivatives of sol with respect to the dim dimensions with derivatives
   * *aligned as a vector)
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div with derivatives
   *
   */
static void PrepareFastDiff(int dim, char * ArtDiff, TPZVec<FADREAL> &sol,
                 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv);
#endif

  /**
   * Evaluates the divergent of F
   * @param dsol [in] vector of solution values derived with respect to the spatial variables
   * @param dphi [in] matrix containig the derivatives of shapefuntions
   * @param Ai [in] vector of matrices tangent to the flux components
   * @param Div [out] value of the divergent
   * @param dDiv [out] an apposimation to the matrix tangent to the divergent
   *
   */
static void Divergent(TPZFMatrix &dsol,
			   TPZFMatrix & dphi,
			   TPZVec<TPZDiffMatrix<REAL> > & Ai,
			   TPZVec<REAL> & Div,
			   TPZDiffMatrix<REAL> * dDiv);

#ifdef _AUTODIFF

  /**
   * Evaluates the divergent of F with FAD type
   * @param dsol [in] derivatives of the solution aligned in one vector
   * @param Ai [in] vector of matrices tangent to the flux components
   * @param Div [out] divergent with derivatives.
   *
   *
   */
static void Divergent(TPZVec<FADREAL> &dsol,
			   TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
			   TPZVec<FADREAL> & Div);
#endif

  /**
   * Computes the diffusive term according to the name
   * @param ArtDiff [in]  Name of the diffusive term
   * @param Ai [in] vector of tensors tangent to the directional fluxes
   * @param Tau [out] Diffusive tensors
   */
template <class T>
static void ComputeTau(int dim, char * ArtDiff,
		 TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau);
private:

  template <class T>
  static void SUPG(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  static void LS(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  static void Bornhaus(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

public:
  /**
   * Flux of Roe (MOUSE program)
   */
  template <class T>
  static void Roe_Flux(const T & rho_f, const T & rhou_f, const T & rhov_f, const T & rhow_f, const T & rhoE_f,
	       const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhow_t, const T & rhoE_t,
	       const REAL nx, const REAL ny, const REAL nz, const REAL gam, T & flux_rho, T & flux_rhou,
	       T & flux_rhov,T & flux_rhow, T & flux_rhoE);

  template <class T>
  static void Roe_Flux(const T & rho_f,const T & rhou_f,const T & rhov_f,const T & rhoE_f,
		   const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhoE_t,
		   const REAL nx, const REAL ny, const REAL gam,
		   T &flux_rho, T &flux_rhou,
		   T &flux_rhov, T &flux_rhoE);
};

#endif

