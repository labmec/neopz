#ifndef PZARTDIFF_H
#define PZARTDIFF_H


#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzdiffmatrix.h"
#include "pzstring.h"
#include "TPZCompElDisc.h"
#include "pzsave.h"
#include "pzmaterialid.h"

#ifdef _AUTODIFF
   #include "fadType.h"
#endif

/* Enumerate to define the possible types of
artificial diffusion term*/
enum TPZArtDiffType
{
  None_AD = -1,
  LeastSquares_AD = 1,
  Bornhaus_AD = 2,
  SUPG_AD = 3,
  TrnLeastSquares_AD = 4
};

/**
 * This classroom adds to the term of diffusion to the variacional formulation
 * of the differential equation partial compressible of Euler (hyperbolic). This
 * term is introduced to stabilize the numerical method of approach.
 *
 * @author Erick Slis
 * @author Cedric Marcelo Augusto Ayala Bravo
 */
class TPZArtDiff : public TPZSaveable
{
public:

   TPZArtDiff();

   TPZArtDiff(TPZArtDiffType type, REAL gamma, REAL CFL = 0., REAL delta = 0.);

   virtual ~TPZArtDiff();

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
    *
    * @para [in] deltaX
    * if fDelta is positive, then it is used as an overall
    * scale on the deltaX.
    * if fDelta is negative, then the absolute value of it
    * is used as the delta for all elements.
    * if fDelta is zero, then the routine OptimalDelta
    * is called, but this method is useful only for
    * the reflected shock problem.
    * Prefer to use fDelta = 1.
    */
   REAL Delta(double deltaX);

   /**
    * Sets the value for delta
    */
   void SetDelta(REAL delta);

   /**
    * Returns the best value for CFL based on the
    * interpolation degree
    */
   REAL OptimalCFL(int degree = TPZCompElDisc::gDegree);

     /**
  Save the element data to a stream
  */
  void Write(TPZStream &buf, int withclassid);

  /**
  Read the element data from a stream
  */
  void Read(TPZStream &buf, void *context);

  /**
  Class identificator
  */
  virtual int ClassId() const {
     return TPZARTDIFFID;
  }


   
   /**
    * pressure
    */
   template< class T >
   static void Pressure(REAL gamma, int dim, T& press, TPZVec<T> &U);

//-------------------A B C matrices and operations

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
   * @param Ai [in] vector of matrices tangent to the flux components.
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
   * Evaluates the divergent of F
   * @param dsol [in] vector of solution values derived with respect to the spatial variables
   * @param phi [in] matrix containig the shapefuntions
   * @param dphi [in] matrix containig the derivatives of shapefuntions
   * @param Ai [in] vector of matrices tangent to the flux components.
   * T must be a FAD class type, so the consistent Divergent jacobian is
   * evaluated
   * @param Div [out] value of the divergent
   * @param dDiv [out] an apposimation to the matrix tangent to the divergent
   *
   */

template <class T>
void Divergent(TPZFMatrix &dsol,
			   TPZFMatrix & phi,
			   TPZFMatrix & dphi,
			   TPZVec<TPZDiffMatrix<T> > & Ai,
			   TPZVec<REAL> & Div,
			   TPZDiffMatrix<REAL> * dDiv);



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
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param dim [in] Spatial dimension of the problem
   * @param Ai [in] vector of tensors tangent to the directional fluxes
   * @param Tau [out] Diffusive tensors
   */
  template <class T>
  void ComputeTau(int dim,
		 TPZFMatrix &jacinv,
		 TPZVec<T> & Sol,
		 TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau);

//----------------------Intermediate Matrices



template <class T>
inline static void RotMatrix(TPZVec<T> & sol, T & us, TPZDiffMatrix<T> &Rot, TPZDiffMatrix<T> &RotT);

template <class T>
inline static void MMatrix(TPZVec<T> & sol, T & us, REAL gamma, TPZDiffMatrix<T> &M, TPZDiffMatrix<T> &Mi);

template <class T>
inline static void RMMatrix(TPZVec<T> & sol, T & us, REAL gamma, TPZDiffMatrix<T> &RTM, TPZDiffMatrix<T> &RMi);

template <class T>
static void EigenSystemSUPG(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZDiffMatrix<T> &X, TPZDiffMatrix<T> &Xi, TPZDiffMatrix<T> &Lambda);

template <class T>
static void EigenSystemBornhaus(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZVec<REAL> & aaS, TPZDiffMatrix<T> &Y, TPZDiffMatrix<T> &Yi, TPZDiffMatrix<T> &Lambda);

template <class T>
static void ContributeBornhaus(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZVec<REAL> & aaS, TPZDiffMatrix<T> &Mat);

private:

  template <class T>
  void SUPG(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  void LS(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  void LST(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  template <class T>
  void Bornhaus(int dim, TPZFMatrix &jacinv, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

public:

//------------------ Diff setup

  /**
   * Computes the common values A B C and Tau vector of matrixes
   * for contributions
   *
   * @param dim [in]
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param U [in] vector of solutions at the point
   * @param Ai [out] matrixes A B C
   * @param Tau [out] diffusive vector of matrixes
   */
  template <class T>
  void PrepareDiff(int dim,
		 TPZFMatrix &jacinv, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);

  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div
   * @param pTaudDiv [out] pointer of a vector to matrices to store the
   * approximate derivatives of Tau_i*Div. If Null, then no approximation
   * derivative is evaluated.
   *
   */

void PrepareFastDiff(int dim,
		 TPZFMatrix &jacinv, TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv,
		 TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv = NULL);

#ifdef _AUTODIFF
  /**
   * Prepares the data to compute the diffusive term
   * as fast as possible, sparing operations
   *
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions setup with shape functions
   * @param dsol [in] derivatives of sol with respect to the dim dimensions
   * *aligned as a vector)
   * @param TauDiv [out] Vector of Vectors to store the values of Tau_i*Div
   *
   */
void PrepareFastDiff(int dim,
		 TPZFMatrix &jacinv, TPZVec<FADREAL> &sol,
                 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv);

template <int dim>
void TPZArtDiff::PrepareFastestDiff(TPZFMatrix &jacinv,
		TPZVec<REAL> &sol,
                TPZFMatrix &dsol,
		TPZFMatrix &phi,
                TPZFMatrix &dphi,
		TPZVec<TPZVec<REAL> > & TauDiv,
		TPZVec<TPZDiffMatrix<REAL> > & dTauDiv);

#endif

//-----------------Contribute

  /**
   * Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * @param deltaX [in] diameter of element (used only if fDelta > 0);
   */
void ContributeApproxImplDiff(int dim,
			   TPZFMatrix &jacinv,
                           TPZVec<REAL> &sol, TPZFMatrix &dsol,
			   TPZFMatrix &dphix,
			   TPZFMatrix &ek, TPZFMatrix &ef,
			   REAL weight, REAL timeStep,
			   double deltaX);

  /**
   * Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions
   * @param dsol [in] derivatives of U with respect to the dim dimensions
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * @param deltaX [in] diameter of element (used only if fDelta > 0);
   */
void ContributeExplDiff(int dim,
			   TPZFMatrix &jacinv,
                           TPZVec<REAL> &sol, TPZFMatrix &dsol,
			   TPZFMatrix &dphix,
			   TPZFMatrix &ef,
			   REAL weight, REAL timeStep,
			   double deltaX);

#ifdef _AUTODIFF

  /**
   * Contributes the diffusion term to the tangent matrix (ek) and residual vector (ef)
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions setup with derivatives
   * @param dsol [in] derivatives of U with respect to the dim dimensions setup with derivatives (xi components aligned as a vector)
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * @param deltaX [in] diameter of element (used only if fDelta > 0);
   */
void ContributeImplDiff(int dim,
			   TPZFMatrix &jacinv,
                           TPZVec<FADREAL> &sol, TPZVec<FADREAL> &dsol,
			   TPZFMatrix &ek, TPZFMatrix &ef,
			   REAL weight, REAL timeStep,
			   double deltaX);

  /**
   * Contributes the diffusion term to the tangent matrix (ek) and residual vector (ef) with FAD template classes
   * @param dim [in] dimension
   * @param jacinv [in] Inverse jacobian of the mapping.
   * @param sol [in] solution of the dim+2 state functions setup with derivatives
   * @param dsol [in] derivatives of U with respect to the dim dimensions setup with derivatives (xi components aligned as a vector)
   * @param phi [in] list of shape functions.
   * @param dphi [in] derivatives of shape functions.
   * @param ek [out] Tangent matrix to contribute to
   * @param ef [out] Residual vector to contribute to
   * @param weight [in] Gaussian quadrature integration weight
   * @param timeStep [in]
   * @param deltaX [in] diameter of element (used only if fDelta > 0);
   */

template <int dim>
void ContributeFastestImplDiff_dim(
                          TPZFMatrix &jacinv,
			  TPZVec<REAL> &sol, TPZFMatrix &dsol,
			  TPZFMatrix &phi, TPZFMatrix &dphi,
			  TPZFMatrix &ek, TPZFMatrix &ef,
			  REAL weight, REAL timeStep,
			  double deltaX);

void ContributeFastestImplDiff(int dim, TPZFMatrix &jacinv,
			  TPZVec<REAL> &sol, TPZFMatrix &dsol,
			  TPZFMatrix &phi, TPZFMatrix &dphi,
			  TPZFMatrix &ek, TPZFMatrix &ef,
			  REAL weight, REAL timeStep,
			  double deltaX);

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


template <class T>
void TPZArtDiff::RotMatrix(TPZVec<T> & sol, T & us, TPZDiffMatrix<T> &Rot, TPZDiffMatrix<T> &RotT)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T u, v, w, u2, v2, w2, uspuInv, usInv;

   switch (dim)
   {
      case (2):
      u = sol[1]/sol[0];
      v = sol[2]/sol[0];
      u2 = u * u;
      v2 = v * v;
      usInv= T(1.)/us;

      Rot.Redim(nstate,nstate);
      Rot(0,0) = 1.;
      Rot(3,3) = 1.;
      Rot(1,1) = u * usInv;
      Rot(1,2) = v * usInv;
      Rot(2,1) = - Rot(1,2);
      Rot(2,2) = Rot(1,1);

      Rot.Transpose(RotT);
      break;
      case (3):
      u = sol[1]/sol[0];
      v = sol[2]/sol[0];
      w = sol[3]/sol[0];
      u2 = u * u;
      v2 = v * v;
      w2 = w * w;
      uspuInv = T(1.)/(us + u);
      usInv = T(1.)/us;
      uspuInv *= usInv;

      Rot.Redim(nstate,nstate);
      Rot(0,0) = 1.;
      Rot(4,4) = 1.;
      Rot(1,1) = u * usInv;
      Rot(1,2) = v * usInv;
      Rot(1,3) = w * usInv;
      Rot(2,1) = -Rot(1,2);
      Rot(2,2) = Rot(1,1) + w2 * uspuInv;
      Rot(2,3) = - w*v * uspuInv;
      Rot(3,1) = -Rot(1,3);
      Rot(3,2) = Rot(2,3);
      Rot(3,3) = Rot(1,1) + v2 * uspuInv;

      Rot.Transpose(RotT);
      break;

      default:
      PZError << "TPZArtDiff::RotMatrix Error: Invalid Dimension\n";
   }

}


template <class T>
void TPZArtDiff::MMatrix(TPZVec<T> & sol, T & us, REAL gamma, TPZDiffMatrix<T> &M, TPZDiffMatrix<T> &Mi)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T rhoInv, usInv;

   switch (dim)
   {
      case (2):
      usInv= T(1.)/us;
      rhoInv = T(1.)/sol[0];

      M.Redim(nstate,nstate);
      M(0,0) = 1.;
      M(1,0) = us;
      M(1,1) = sol[0];
      M(2,2) = sol[0];
      M(3,0) = us * us /2.;
      M(3,1) = sol[0] * us;
      M(3,3) = T(1.)/(gamma-1.);

      Mi.Redim(nstate,nstate);
      Mi(0,0) = 1.;
      Mi(1,0) = - rhoInv * us;
      Mi(1,1) = rhoInv;
      Mi(2,2) = rhoInv;
      Mi(3,0) = T((gamma-1.)/2.) * us * us;
      Mi(3,1) = T(1.-gamma)*us;
      Mi(3,3) = gamma-1.;

      break;
      case(3):
      usInv= T(1.)/us;
      rhoInv = T(1.)/sol[0];

      M.Redim(nstate,nstate);
      M(0,0) = 1.;
      M(1,0) = us;
      M(1,1) = sol[0];
      M(2,2) = sol[0];
      M(3,3) = sol[0];
      M(4,0) = us * us /2.;
      M(4,1) = sol[0] * us;
      M(4,4) = T(1.)/(gamma-1.);

      Mi.Redim(nstate,nstate);
      Mi(0,0) = 1.;
      Mi(1,0) = - rhoInv * us;
      Mi(1,1) = rhoInv;
      Mi(2,2) = rhoInv;
      Mi(3,3) = rhoInv;
      Mi(4,0) = T((gamma-1.)/2.) * us * us;
      Mi(4,1) = T(1.-gamma)*us;
      Mi(4,4) = gamma-1.;
      break;

      default:
      PZError << "TPZArtDiff::MMatrix Error: Invalid Dimension\n";
   }
}


template <class T>
void TPZArtDiff::RMMatrix(TPZVec<T> & sol, T & us, REAL gamma, TPZDiffMatrix<T> &RTM, TPZDiffMatrix<T> &RMi)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T u, v, w, u2, v2, w2,
     rhoInv, usInv, uspu, uspuInv;
   T kt, st;

   switch (dim)
   {
      case (2):
      u = sol[1]/sol[0];
      usInv= T(1.)/us;
      rhoInv = T(1.)/sol[0];
      // kt = u/us, st = v/us
      kt = sol[1]/sol[0]/us;
      st = sol[2]/sol[0]/us;

      RTM.Redim(nstate,nstate);
      RTM(0,0) =  1.;
      RTM(1,0) =  us * kt;
      RTM(1,1) =  sol[0] * kt;
      RTM(1,2) = -sol[0] * st;
      RTM(2,0) =  us * st;
      RTM(2,1) = -RTM(1,2);
      RTM(2,2) =  RTM(1,1);
      RTM(3,0) = us * us /2.;
      RTM(3,1) = sol[0] * us;
      RTM(3,3) = 1./(gamma-1.);

      RMi.Redim(nstate,nstate);
      RMi(0,0) = 1.;
      RMi(1,0) = - rhoInv * us;
      RMi(1,1) = rhoInv * kt;
      RMi(1,2) = rhoInv * st;
      RMi(2,1) = - RMi(1,2);
      RMi(2,2) = RMi(1,1);
      RMi(3,0) = T((gamma-1.)/2.) * us * us;
      RMi(3,1) = T(1.-gamma)*us*kt;
      RMi(3,2) = T(1.-gamma)*us*st;
      RMi(3,3) = gamma-1.;

      break;
      case(3):
      u = sol[1]/sol[0];
      v = sol[2]/sol[0];
      w = sol[3]/sol[0];
      u2 = u * u;
      v2 = v * v;
      w2 = w * w;
      uspu = us + u;
      uspuInv = T(1.)/uspu;
      usInv = T(1.)/us;
      kt = T(1.)/(sol[0] * us * (uspu * u + v2 + w2));
#ifdef NOTDEFINED
      RTM.Redim(nstate,nstate);
      RTM(0,0) = 1.;
      RTM(1,0) = u2 * usInv;
      RTM(1,1) = sol[1]/*rhou*/*usInv;
      RTM(1,2) = sol[2]/*rhov*/*usInv;
      RTM(1,3) = sol[3]/*rhow*/*usInv;
      RTM(2,0) = - u * v * usInv;
      RTM(2,1) = - RTM(1,2);
      RTM(2,2) = RTM(1,1) + sol[3] * w * uspuInv * usInv;
      RTM(2,3) = - RTM(1,2) * w * uspuInv;
      RTM(3,0) = - u * w * usInv;
      RTM(3,1) = -RTM(1,3);
      RTM(3,2) = RTM(2,3);
      RTM(3,3) = RTM(1,1) + sol[2] * v * uspuInv * usInv;
      RTM(4,0) = u2/2.;
      RTM(4,1) = sol[1];
      RTM(4,4) = 1./(gamma-1.);

      RMi.Redim(nstate,nstate);
      RMi(0,0) = 1.;
      RMi(1,0) = - u / sol[0];
      RMi(1,1) = - RMi(1,0) * usInv;
      RMi(1,2) = - v * usInv / sol[0];
      RMi(1,3) = - w * usInv / sol[0];
      RMi(2,1) = - RMi(1,2);
      RMi(2,2) = (u * v2 + uspu * (u2+w2)) * kt;
      RMi(2,3) = (u - uspu) * v * w * kt;
      RMi(3,1) = - RMi(1,3);
      RMi(3,2) = RMi(2,3);
      RMi(3,3) = (uspu * (u2+v2) + u* w2) * kt;
      RMi(4,0) = u2 * T((gamma-1.)/2.);
      RMi(4,1) = - u2 * usInv * T(gamma - 1.);
      RMi(4,2) =   u * v * usInv * T(gamma -1.);
      RMi(4,3) =   u * w * usInv * T(gamma -1.);
      RMi(4,4) = gamma -1.;
#endif

      RTM.Redim(nstate,nstate);
      RTM(0,0) = 1.;
      RTM(1,0) = u;
      RTM(1,1) = sol[1]/*rhou*/*usInv;
      RTM(1,2) = - sol[2]/*rhov*/*usInv;
      RTM(1,3) = - sol[3]/*rhow*/*usInv;
      RTM(2,0) = v;
      RTM(2,1) = - RTM(1,2);
      RTM(2,2) = RTM(1,1) + sol[3] * w * uspuInv * usInv;
      RTM(2,3) = RTM(1,2) * w * uspuInv;
      RTM(3,0) = w;
      RTM(3,1) = -RTM(1,3);
      RTM(3,2) = RTM(2,3);
      RTM(3,3) = RTM(1,1) + sol[2] * v * uspuInv * usInv;
      RTM(4,0) = us * us /2.;
      RTM(4,1) = sol[0] * us;
      RTM(4,4) = 1./(gamma-1.);

      RMi.Redim(nstate,nstate);
      RMi(0,0) = 1.;
      RMi(1,0) = - us / sol[0];
      RMi(1,1) = u * usInv / sol[0];
      RMi(1,2) = v * usInv / sol[0];
      RMi(1,3) = w * usInv / sol[0];
      RMi(2,1) = - RMi(1,2);
      RMi(2,2) = (u * v2 + uspu * (u2+w2)) * kt;
      RMi(2,3) = (u - uspu) * v * w * kt;
      RMi(3,1) = - RMi(1,3);
      RMi(3,2) = RMi(2,3);
      RMi(3,3) = (uspu * (u2+v2) + u* w2) * kt;
      RMi(4,0) = us * us * T((gamma-1.)/2.);
      RMi(4,1) = u * T(1.-gamma);
      RMi(4,2) = v * T(1.-gamma);
      RMi(4,3) = w * T(1.-gamma);
      RMi(4,4) = gamma -1.;

      break;

      default:
      PZError << "TPZArtDiff::MMatrix Error: Invalid Dimension\n";
   }
}



template <class T>
void TPZArtDiff::EigenSystemSUPG(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZDiffMatrix<T> &X, TPZDiffMatrix<T> &Xi, TPZDiffMatrix<T> &Lambda)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T us2, c2, s, rho_c_us;

   c2 = c * c;
   us2 = us * us;
   rho_c_us = sol[0] * c * us;

   switch (dim)
   {
      case (2):
      s = sqrt(c2 + us2 * 16.);

      X.Redim(nstate,nstate);
      X(0,0) = 1.;
      X(0,2) = 4. * sol[0] * us / c2;
      X(0,3) = X(0,2);
      X(1,2) = - s/c -1.;
      X(1,3) =   s/c -1.;
      X(2,1) = 1.;
      X(3,2) = 4. * sol[0] * us;
      X(3,3) = X(3,2);

      Xi.Redim(nstate,nstate);
      Xi(0,0) = 1.;
      Xi(0,3) = - 1. / c2;
      Xi(1,2) = 1.;
      Xi(2,1) = - c / (2. * s);
      Xi(2,3) = (s-c)/(8. * sol[0] * s * us);
      Xi(3,1) = -Xi(2,1);
      Xi(3,3) = (s+c)/(8. * sol[0] * s * us);

      Lambda.Redim(nstate,nstate);
      Lambda(0,0) = 1./us;
      Lambda(1,1) = 1./sqrt(us2 + c2);
      Lambda(2,2) = 1./sqrt(us2 + 1.5 * c2 - .5 * c * s);
      Lambda(3,3) = 1./sqrt(us2 + 1.5 * c2 + .5 * c * s);

      break;
      case(3):

      s = sqrt(c2 + us2 * 4.);

      X.Redim(nstate,nstate);
      X(0,0) = 1.;
      X(0,3) = 1. / c2;
      X(0,4) = X(0,3);
      X(1,3) = (-s-c)/(2. * rho_c_us);
      X(1,4) = ( s-c)/(2. * rho_c_us);
      X(2,2) = 1.;
      X(3,1) = 1.;
      X(4,3) = 1.;
      X(4,4) = 1.;


      Xi.Redim(nstate,nstate);
      Xi(0,0) = 1.;
      Xi(0,4) = - 1. / c2;
      Xi(1,3) = 1.;
      Xi(2,2) = 1.;
      Xi(3,1) = - rho_c_us / s;
      Xi(3,4) = .5 - c / (2. * s);
      Xi(4,1) = -Xi(3,1);
      Xi(4,4) = .5 + c / (2. * s);

      Lambda.Redim(nstate,nstate);
      Lambda(0,0) = 1./us;
      Lambda(1,1) = 1./sqrt(us2 + c2);
      Lambda(2,2) = Lambda(1,1);
      Lambda(3,3) = 1./sqrt(us2 + 2. * c2 - c * s);
      Lambda(4,4) = 1./sqrt(us2 + 2. * c2 + c * s);
      break;

      default:
      PZError << "TPZArtDiff::EigenSystemSUPG Error: Invalid Dimension\n";
   }
}




template <class T>
void TPZArtDiff::EigenSystemBornhaus(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZVec<REAL> & aaS, TPZDiffMatrix<T> &Y, TPZDiffMatrix<T> &Yi, TPZDiffMatrix<T> &Lambda)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T k, us2, c2, rho_c, temp1, temp2, temp3, k2;

   c2 = c * c;
   us2 = us * us;
   rho_c = sol[0] * c;

   switch (dim)
   {
      case (2):
      k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1];
      k = sqrt(k2);


      Y.Redim(nstate,nstate);
      Y(0,1) = 1.;
      Y(0,2) = 1. / c2;
      Y(0,3) = Y(0,2);
      Y(1,0) = -aaS[1] ;/// aaS[0];
      Y(1,2) = -aaS[0] / (k * rho_c);
      Y(1,3) = -Y(1,2);
      Y(2,0) = aaS[0];//1.;
      Y(2,2) = -aaS[1] / (k * rho_c);
      Y(2,3) = -Y(2,2);
      Y(3,2) = 1.;
      Y(3,3) = 1.;

      Yi.Redim(nstate,nstate);
      Yi(0,1) = - /*aaS[0] * */aaS[1] / k2;
      Yi(0,2) =   /*aaS[0] * */aaS[0] / k2;
      Yi(1,0) = 1.;
      Yi(1,3) = -1. / c2;
      Yi(2,1) = - aaS[0] * rho_c / (2. * k);
      Yi(2,2) = - aaS[1] * rho_c / (2. * k);
      Yi(2,3) = .5;
      Yi(3,1) = -Yi(2,1);
      Yi(3,2) = -Yi(2,2);
      Yi(3,3) = .5;

      temp1 = aaS[0] * us;
      if(val(temp1) < 0)temp1 = -temp1;
      temp2 = aaS[0] * us - k * c;
      if(val(temp2) < 0)temp2 = -temp2;
      temp3 = aaS[0] * us + k * c;
      if(val(temp3) < 0)temp3 = -temp3;

      Lambda.Redim(nstate,nstate);
      Lambda(0,0) = temp1;
      Lambda(1,1) = temp1;
      Lambda(2,2) = temp2;
      Lambda(3,3) = temp3;

      break;
      case (3):
      k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1] + aaS[2]*aaS[2];
      k = sqrt(k2);

      Y.Redim(nstate,nstate);
      Y(0,2) = 1.;
      Y(0,3) = 1. / c2;
      Y(0,4) = Y(0,3);
      Y(1,0) = -aaS[2]/aaS[0];
      Y(1,1) = -aaS[1]/aaS[0];
      Y(1,3) = -aaS[0]/(k * rho_c);
      Y(1,4) = -Y(1,3);
      Y(2,1) = 1.;
      Y(2,3) = -aaS[1]/(k * rho_c);
      Y(2,4) = -Y(2,3);
      Y(3,0) = 1.;
      Y(3,3) = -aaS[2]/(k * rho_c);
      Y(3,4) = -Y(3,3);
      Y(4,3) = 1.;
      Y(4,4) = 1.;

      Yi.Redim(nstate,nstate);
      Yi(0,1) = - aaS[0] * aaS[2] / k2;
      Yi(0,2) = - aaS[1] * aaS[2] / k2;
      Yi(0,3) =  (aaS[0] * aaS[0] + aaS[1] * aaS[1]) / k2;
      Yi(1,1) = - aaS[0] * aaS[1] / k2;
      Yi(1,2) =  (aaS[0] * aaS[0] + aaS[2] * aaS[2]) / k2;
      Yi(1,3) =  Yi(0,2);
      Yi(2,0) = 1.;
      Yi(2,4) = -1. / c2;
      Yi(3,1) = - aaS[0] * rho_c / (2. * k);
      Yi(3,2) = - aaS[1] * rho_c / (2. * k);
      Yi(3,3) = - aaS[2] * rho_c / (2. * k);
      Yi(3,4) = .5;
      Yi(4,1) = - Yi(3,1);
      Yi(4,2) = - Yi(3,2);
      Yi(4,3) = - Yi(3,3);
      Yi(4,4) = .5;

      temp1 = aaS[0] * us;
      if(val(temp1) < 0)temp1 = -temp1;
      temp2 = aaS[0] * us - k * c;
      if(val(temp2) < 0)temp2 = -temp2;
      temp3 = aaS[0] * us + k * c;
      if(val(temp3) < 0)temp3 = -temp3;

      Lambda.Redim(nstate,nstate);
      Lambda(0,0) = temp1;
      Lambda(1,1) = temp1;
      Lambda(2,2) = temp1;
      Lambda(3,3) = temp2;
      Lambda(4,4) = temp3;

      break;
      default:
      PZError << "TPZArtDiff::EigenSystemBornhaus Error: Invalid Dimension\n";
   }
}



template <class T>
void TPZArtDiff::ContributeBornhaus(TPZVec<T> & sol, T & us, T & c, REAL gamma, TPZVec<REAL> & aaS, TPZDiffMatrix<T> &Mat)
{
   int nstate = sol.NElements();
   int dim = nstate - 2;

   T k, us2, l1, l3, l4, twoCK, c2, temp1, temp2/*, temp3*/, k2;


   c2 = c * c;
   us2 = us * us;
   //rho_c = sol[0] * c;

   switch (dim)
   {
      case (2):
      k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1];
      k = sqrt(k2);
      twoCK = c * T(k * 2.);

      l1 = aaS[0] * us;
      if(val(l1) < 0.)l1 = -l1;
      l3 = aaS[0] * us - k * c;
      if(val(l3) < 0.)l3 = -l3;
      l4 = aaS[0] * us + k * c;
      if(val(l4) < 0.)l4 = -l4;

      temp1 = (l4 - l3) * sol[0] / twoCK;
      temp2 = l3 + l4 - 2.*l1;
      Mat(0,0) += l1;
      Mat(0,1) += aaS[0] * temp1;
      Mat(0,2) += aaS[1] * temp1;
      Mat(0,3) += temp2 / (2. * c2);
      Mat(1,1) += (aaS[1]*aaS[1] * l1 + aaS[0]*aaS[0] * (l3+l4) *.5)/k2;
      Mat(1,2) += aaS[0]*aaS[1] * temp2/(2. * k2);
      Mat(1,3) += aaS[0]*(l4-l3)/(c * k * sol[0] * 2.);
      Mat(2,1) += aaS[0]*aaS[1] * temp2/(2. * k2);
      Mat(2,2) += (aaS[0]*aaS[0] * l1 + aaS[1]*aaS[1] * (l3+l4) *.5)/k2;
      Mat(2,3) += aaS[1]*(l4-l3)/(c * k * sol[0] * 2.);
      Mat(3,1) += aaS[0] * c2 * temp1;
      Mat(3,2) += aaS[1] * c2 * temp1;
      Mat(3,3) += (l3 + l4)/2.;


      break;
/*      case (3):
      k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1] + aaS[2]*aaS[2];
      k = sqrt(k2);

      Y.Redim(nstate,nstate);
      Y(0,2) = 1.;
      Y(0,3) = 1. / c2;
      Y(0,4) = Y(0,3);
      Y(1,0) = -aaS[2]/aaS[0];
      Y(1,1) = -aaS[1]/aaS[0];
      Y(1,3) = -aaS[0]/(k * rho_c);
      Y(1,4) = -Y(1,3);
      Y(2,1) = 1.;
      Y(2,3) = -aaS[1]/(k * rho_c);
      Y(2,4) = -Y(2,3);
      Y(3,0) = 1.;
      Y(3,3) = -aaS[2]/(k * rho_c);
      Y(3,4) = -Y(3,3);
      Y(4,3) = 1.;
      Y(4,4) = 1.;

      Yi.Redim(nstate,nstate);
      Yi(0,1) = - aaS[0] * aaS[2] / k2;
      Yi(0,2) = - aaS[1] * aaS[2] / k2;
      Yi(0,3) =  (aaS[0] * aaS[0] + aaS[1] * aaS[1]) / k2;
      Yi(1,1) = - aaS[0] * aaS[1] / k2;
      Yi(1,2) =  (aaS[0] * aaS[0] + aaS[2] * aaS[2]) / k2;
      Yi(1,3) =  Yi(0,2);
      Yi(2,0) = 1.;
      Yi(2,4) = -1. / c2;
      Yi(3,1) = - aaS[0] * rho_c / (2. * k);
      Yi(3,2) = - aaS[1] * rho_c / (2. * k);
      Yi(3,3) = - aaS[2] * rho_c / (2. * k);
      Yi(3,4) = .5;
      Yi(4,1) = - Yi(3,1);
      Yi(4,2) = - Yi(3,2);
      Yi(4,3) = - Yi(3,3);
      Yi(4,4) = .5;

      temp1 = aaS[0] * us;
      if(val(temp1) < 0)temp1 = -temp1;
      temp2 = aaS[0] * us - k * c;
      if(val(temp2) < 0)temp2 = -temp2;
      temp3 = aaS[0] * us + k * c;
      if(val(temp3) < 0)temp3 = -temp3;

      Lambda.Redim(nstate,nstate);
      Lambda(0,0) = temp1;
      Lambda(1,1) = temp1;
      Lambda(2,2) = temp1;
      Lambda(3,3) = temp2;
      Lambda(4,4) = temp3;

      break;*/
      default:
      PZError << "TPZArtDiff::EigenSystemBornhaus Error: Invalid Dimension\n";
   }
}

#endif

