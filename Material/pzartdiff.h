#ifndef PZARTDIFF_H
#define PZARTDIFF_H


#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzdiffmatrix.h"
#include "pzstring.h"
#include "TPZCompElDisc.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif


/* Enumerate to define the possible types of
artificial diffusion term*/
enum TPZArtDiffType
{
  None_AD = -1,
  LeastSquares_AD = 1,
  Bornhaus_AD,
  SUPG_AD
};

/**
 * This classroom adds to the term of diffusion to the variacional formulation
 * of the differential equation partial compressible of Euler (hyperbolic). This
 * term is introduced to stabilize the numerical method of approach.
 *
 * @author Erick Slis
 * @author Cedric Ayala
 */
class TPZArtDiff
{
public:

   TPZArtDiff(TPZArtDiffType type, REAL gamma, REAL CFL = 0., REAL delta = 0.);

   ~TPZArtDiff();

//-----------------Attributes and parameters

   /**
    * Returns the type of artifical diffusion
    */
   TPZArtDiffType ArtDiffType();

   /**
    * Configures the type of artificial diffusion
    */
   void SetArtDiffType(TPZArtDiffType type);

   /**
    * Returns the name of diffusive term
    *
    * @param Type [in] type of diffusive term
    */
   TPZString DiffusionName();

   /**
    * Returns the best value for delta based on
    * the interpolation degree
    */
   REAL OptimalDelta();

   /**
    * Returns the stored value for delta
    */
   REAL Delta();

   /**
    * Sets the value for delta
    */
   void SetDelta(REAL delta);

   /**
    * Returns the best value for CFL based on the
    * interpolation degree
    */
   REAL OptimalCFL(int degree = TPZCompElDisc::gDegree);

//-------------------A B C matrices and operations

  /**
   * Jacobian of the tensor flux of Euler
   * @param dim [in]
   * @param U [in] dim+2 solutions at given point
   * @param Ai [out] vector of dim tensors of (dim+2)*(dim*2),
   * representing the derivatives of F with respect to the
   * dim spatial dimensions
   */

  template <class T>
  void JacobFlux(int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai);

  /**
   * operation product point in the diffusion term
   * @param dphi [in] vector of scalars that make part of the ODot operation
   * @param M [in] vector of matrices to operate
   * @param Result [out]
   *
   */

  template <class T>
  void ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZDiffMatrix<T> > &M, TPZDiffMatrix<T> &Result);

  /**
   * operation product point in the diffusion term
   * @param dphi [in] vector of scalars that make part of the ODot operation
   * @param TauDiv [in] vector of vectors to operate
   * @param Result [out]
   *
   */
  template <class T>
  void ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZVec<T> > &TauDiv, TPZVec<T> &Result);


  /**
   * Evaluates the divergent of F
   * @param dsol [in] vector of solution values derived with respect to the spatial variables
   * @param dphi [in] matrix containig the derivatives of shapefuntions
   * @param Ai [in] vector of matrices tangent to the flux components
   * @param Div [out] value of the divergent
   * @param dDiv [out] an apposimation to the matrix tangent to the divergent
   *
   */
  void Divergent(TPZFMatrix &dsol,
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
  void Divergent(TPZVec<FADREAL> &dsol,
			   TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
			   TPZVec<FADREAL> & Div);
#endif

//----------------------Tau tensor

  /**
   * Computes the diffusive term according to the name
   * @param dim [in] Spatial dimension of the problem
   * @param Ai [in] vector of tensors tangent to the directional fluxes
   * @param Tau [out] Diffusive tensors
   */
  template <class T>
  void ComputeTau(int dim,
		 TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau);

private:

  template <class T>
  void SUPG(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  void LS(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  void Bornhaus(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

public:

//------------------ Diff setup

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
  void PrepareDiff(int dim, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div
   * @param pTaudDiv [out] pointer of a vector to matrices to store the
   * approximate derivatives of Tau_i*Div. If Null, then no approximation
   * derivative is evaluated.
   *
   */

void PrepareFastDiff(int dim,TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv,
		 TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv = NULL);

#ifdef _AUTODIFF
  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param sol [in] solution of the dim+2 state functions setup with shape functions
   * @param dsol [in] derivatives of sol with respect to the dim dimensions with derivatives
   * *aligned as a vector)
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div with derivatives
   *
   */
void PrepareFastDiff(int dim, TPZVec<FADREAL> &sol,
                 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv);
#endif

//-----------------Contribute

  /**
   * Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
   * @param dim [in] dimension
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * is to be computed and contributed to ek.
   */
void ContributeApproxImplDiff(int dim,
                           TPZVec<REAL> &sol, TPZFMatrix &dsol,
			   TPZFMatrix &dphix,
			   TPZFMatrix &ek, TPZFMatrix &ef,
			   REAL weight, REAL timeStep);

  /**
   * Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
   * @param dim [in] dimension
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * is to be computed and contributed to ek.
   */
void ContributeExplDiff(int dim,
                           TPZVec<REAL> &sol, TPZFMatrix &dsol,
			   TPZFMatrix &dphix,
			   TPZFMatrix &ef,
			   REAL weight, REAL timeStep);

#ifdef _AUTODIFF

  /**
   * Contributes the diffusion term to the tangent matrix (ek) and residual vector (ef)
   * @param dim [in] dimension
   * @param sol [in] solution of the dim+2 state functions setup with derivatives
   * @param dsol [in] derivatives of U with respect to the dim dimensions setup with derivatives (xi components aligned as a vector)
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   */
void ContributeImplDiff(int dim,
                           TPZVec<FADREAL> &sol, TPZVec<FADREAL> &dsol,
			   TPZFMatrix &ek, TPZFMatrix &ef,
			   REAL weight, REAL timeStep);

#endif


private:

  /**
   * kind of artificial diffusion term to apply
   */
  TPZArtDiffType fArtDiffType;

  /**
   * ratio between specific heat is constant and the specific heat the constant
   * volume of a polytropic gas
   */
  REAL fGamma;

  /**
   * scalar coefficient of the element in the diffusion term
   */
  REAL fDelta;

  /**
   * CFL number
   */
  REAL fCFL;


};

#endif

