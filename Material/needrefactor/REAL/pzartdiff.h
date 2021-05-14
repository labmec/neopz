/**
 * \file
 * @brief Contains the TPZArtDiff class which implements a numerical diffusivity coefficient.
 */

#ifndef TPZARTDIFF_HH
#define TPZARTDIFF_HH


#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzdiffmatrix.h"
#include "pzstring.h"
#include "TPZCompElDisc.h"
#include "TPZSavable.h"

#include "fadType.h"


/** @brief Returns value of the variable */ 
inline REAL val(STATE & number)
{
	return number;
}

/**
 * @ingroup material 
 * @enum TPZArtDiffType
 * @brief Enumerate to define the possible types of artificial diffusion term to stabilize the numerical scheme
 * @var None_AD
 * @brief No Artificial diffusion term is considered
 * @var LeastSquares_AD
 * @brief Use Least squares method to applied artificial diffusion term
 * @var Bornhaus_AD
 * @brief Use Bornhaus method to applied artificial diffusion term
 * @var SUPG_AD 
 * @brief Use Supg method to applied artificial diffusion term
 * @var TrnLeastSquares_AD
 * @brief Use Transpose of the Least Squares method to applied artificial diffusion term
 */
enum TPZArtDiffType
{
	None_AD = -1,
	LeastSquares_AD = 1,
	Bornhaus_AD = 2,
	SUPG_AD = 3,
	TrnLeastSquares_AD = 4
};

/**
 * @ingroup material
 * @author Erick Slis
 * @author Cedric Marcelo Augusto Ayala Bravo
 * @brief This class adds to the term of diffusion to the variacional formulation \n
 * of the differential equation partial compressible of Euler (hyperbolic).
 */
/**
 * This term is introduced to stabilize the numerical method of approach.
 */
class TPZArtDiff : public TPZSavable
{
public:
	/** @brief Simple constructor */
	TPZArtDiff();
	/** @brief Specific constructor */
	TPZArtDiff(TPZArtDiffType type, REAL gamma, REAL CFL = 0., REAL delta = 0.);
	
	virtual ~TPZArtDiff();
	
	/** @name Attributes and parameters
	 * @{
	 */
	
	/** @brief Returns the type of artifical diffusion */
	TPZArtDiffType ArtDiffType();
	
	/**
	 * @brief Configures the type of artificial diffusion
	 * @param[in] type Type of diffusive term
	 */
	void SetArtDiffType(TPZArtDiffType type);
	
	/** @brief Returns the name of diffusive term */
	TPZString DiffusionName();
	
	/** @brief Returns the best value for delta based on the interpolation degree */
	REAL OptimalDelta();
	
	/**
	 * @brief Returns the stored value for delta
	 * @param[in] sol Solution
	 * @param[in] deltaX Delta X 
	 */
	/**
	 * If fDelta is positive, then it is used as an overall
	 * scale on the deltaX. \n
	 * If fDelta is negative, then the absolute value of it
	 * is used as the delta for all elements. \n
	 * If fDelta is zero, then the routine OptimalDelta
	 * is called, but this method is useful only for
	 * the reflected shock problem. \n
	 * Prefer to use fDelta = 1.
	 */
	REAL Delta(REAL deltaX, TPZVec<STATE> & sol);
	
	/** @brief Sets the value for delta */
	void SetDelta(REAL delta);
	
	/** @brief Returns the best value for \f$CFL\f$ based on the interpolation degree */
	REAL OptimalCFL(int degree = TPZCompEl::GetgOrder());
	
	/** @brief Save the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Read the element data from a stream */
	void Read(TPZStream &buf, void *context) override;
	
	/** @brief Class identificator */
public:
int ClassId() const override;

	
	/** @brief pressure */
	template< class T >
	static void Pressure(REAL gamma, int dim, T& press, TPZVec<T> &U);
	
	/** @} */
	
	/** @name A B C matrices and operations
	 * @{ 
	 */
	
	/**
	 * @brief Operation product point in the diffusion term
	 * @param[in] dphi vector of scalars that make part of the ODot operation
	 * @param[in] M vector of matrices to operate
	 * @param[out] Result Result of the operation
	 */
	template <class T>
	void ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZDiffMatrix<T> > &M, TPZDiffMatrix<T> &Result);
	
	/**
	 * @brief Operation product point in the diffusion term
	 * @param[in] dphi vector of scalars that make part of the ODot operation
	 * @param[in] TauDiv vector of vectors to operate
	 * @param[out] Result Result of the operation
	 */
	template <class T>
	void ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZVec<T> > &TauDiv, TPZVec<T> &Result);
	
	
	/**
	 * @brief Evaluates the divergent of F
	 * @param[in] dsol vector of solution values derived with respect to the spatial variables
	 * @param[in] dphi matrix containig the derivatives of shapefuntions
	 * @param[in] Ai vector of matrices tangent to the flux components.
	 * @param[out] Div value of the divergent
	 * @param[out] dDiv an apposimation to the matrix tangent to the divergent
	 */
	
	void Divergent(TPZFMatrix<STATE> &dsol,
				   TPZFMatrix<REAL> & dphi,
				   TPZVec<TPZDiffMatrix<STATE> > & Ai,
				   TPZVec<STATE> & Div,
				   TPZDiffMatrix<STATE> * dDiv);
	

	/**
	 * @brief Evaluates the divergent of F
	 * @param[in] dsol Vector of solution values derived with respect to the spatial variables
	 * @param[in] phi Matrix containig the shapefuntions
	 * @param[in] dphi Matrix containig the derivatives of shapefuntions
	 * @param[in] Ai Vector of matrices tangent to the flux components.
	 * @param[out] Div Value of the divergent
	 * @param[out] dDiv An apposimation to the matrix tangent to the divergent
	 */
	/**
	 * @note T must be a FAD class type, so the consistent Divergent jacobian is evaluated
	 */
	template <class T>
	void Divergent(TPZFMatrix<STATE> &dsol,
				   TPZFMatrix<REAL> & phi,
				   TPZFMatrix<REAL> & dphi,
				   TPZVec<TPZDiffMatrix<T> > & Ai,
				   TPZVec<STATE> & Div,
				   TPZDiffMatrix<STATE> * dDiv);
	
	
	
	/**
	 * @brief Evaluates the divergent of F with FAD type
	 * @param[in] dsol Derivatives of the solution aligned in one vector
	 * @param[in] Ai Vector of matrices tangent to the flux components
	 * @param[out] Div Divergent with derivatives.
	 */
	void Divergent(TPZVec<FADREAL> &dsol,
				   TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
				   TPZVec<FADREAL> & Div);
	/** @} */
	
	/** Tau tensor */
	/**
	 * @brief Computes the diffusive term according to the name
	 * @param[in] dim Spatial dimension of the problem
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] Sol Solution vector
	 * @param[in] Ai Vector of tensors tangent to the directional fluxes
	 * @param[out] Tau Diffusive tensors
	 */
	template <class T>
	void ComputeTau(int dim,
					TPZFMatrix<REAL> &jacinv,
					TPZVec<T> & Sol,
					TPZVec<TPZDiffMatrix<T> > &Ai,
					TPZVec<TPZDiffMatrix<T> > &Tau);
	
	/** @name Intermediate Matrices
	 * @{
	 */
	
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
	
	/** @} */
	
private:
	
	template <class T>
	void SUPG(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);
	
	template <class T>
	void LS(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);
	
	template <class T>
	void LST(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);
	
	template <class T>
	void Bornhaus(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);
	
public:
	
	/** @name Diff setup
	 * @{
	 */
	
	/**
	 * @brief Computes the common values A B C and Tau vector of matrixes for contributions
	 * @param[in] dim Spatial dimension of the problem
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] U Vector of solutions at the point
	 * @param[out] Ai Matrixes A B C
	 * @param[out] Tau Diffusive vector of matrixes
	 */
	template <class T>
	void PrepareDiff(int dim,
					 TPZFMatrix<REAL> &jacinv, TPZVec<T> &U,
					 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau);
	
	/**
	 * @brief Prepares the data to compute the diffusive term as fast as possible, sparing operations
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions
	 * @param[in] dsol Derivatives of U with respect to the dim dimensions
	 * @param[in] dphi Derivatives of shape functions.
	 * @param[out] TauDiv Vector of Vectors to store the values of Tau_i*Div
	 * @param[out] pTaudDiv Pointer of a vector to matrices to store the
	 * approximate derivatives of \f$Tau_i*Div\f$. \n If Null, then no approximation
	 * derivative is evaluated.
	 */
	void PrepareFastDiff(int dim,
						 TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol,
						 TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> & dphi,
						 TPZVec<TPZVec<STATE> > & TauDiv,
						 TPZVec<TPZDiffMatrix<STATE> > * pTaudDiv = NULL);
	
	/**
	 * @brief Prepares the data to compute the diffusive term as fast as possible, sparing operations
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions setup with shape functions
	 * @param[in] dsol Derivatives of sol with respect to the dim (dimensions aligned as a vector)
	 * @param[out] TauDiv Vector of Vectors to store the values of \f$Tau_i*Div\f$
	 */
	void PrepareFastDiff(int dim,
						 TPZFMatrix<REAL> &jacinv, TPZVec<FADREAL> &sol,
						 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv);
	
	template <int dim>
	void PrepareFastestDiff(TPZFMatrix<REAL> &jacinv,
							TPZVec<STATE> &sol,
							TPZFMatrix<STATE> &dsol,
							TPZFMatrix<REAL> &phi,
							TPZFMatrix<REAL> &dphi,
							TPZVec<TPZVec<STATE> > & TauDiv,
							TPZVec<TPZDiffMatrix<STATE> > & dTauDiv);
	/** @} */
	

	/** @name Contribute methods
	 * @{
	 */
	
	/**
	 * @brief Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions
	 * @param[in] dsol Derivatives of U with respect to the dim dimensions
	 * @param[in] dphix Derivatives of shape functions.
	 * @param[out] ek Tangent matrix to contribute to
	 * @param[out] ef Residual vector to contribute to
	 * @param[in] weight Gaussian quadrature integration weight
	 * @param[in] timeStep Time step
	 * @param[in] deltaX Diameter of element (used only if \f$fDelta > 0\f$);
	 */
	void ContributeApproxImplDiff(int dim,
								  TPZFMatrix<REAL> &jacinv,
								  TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
								  TPZFMatrix<REAL> &dphix,
								  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
								  REAL weight, REAL timeStep,
								  REAL deltaX);
	
	/**
	 * @brief Contributes the diffusion term to the tangent matrix (ek-approximated) and residual vector (ef)
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions
	 * @param[in] dsol Derivatives of U with respect to the dim dimensions
	 * @param[in] dphix Derivatives of shape functions.
	 * @param[out] ef Residual vector to contribute to
	 * @param[in] weight Gaussian quadrature integration weight
	 * @param[in] timeStep Time step
	 * @param[in] deltaX Diameter of element (used only if fDelta > 0);
	 */
	void ContributeExplDiff(int dim,
							TPZFMatrix<REAL> &jacinv,
							TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
							TPZFMatrix<REAL> &dphix,
							TPZFMatrix<STATE> &ef,
							REAL weight, REAL timeStep,
							REAL deltaX);
	

	/**
	 * @brief Contributes the diffusion term to the tangent matrix (ek) and residual vector (ef)
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions setup with derivatives
	 * @param[in] dsol Derivatives of U with respect to the dim dimensions setup with derivatives (xi components aligned as a vector)
	 * @param[out] ek Tangent matrix to contribute to
	 * @param[out] ef Residual vector to contribute to
	 * @param[in] weight [in] Gaussian quadrature integration weight
	 * @param[in] timeStep 
	 * @param[in] deltaX Diameter of element (used only if \f$ fDelta > 0\f$);
	 */
	void ContributeImplDiff(int dim,
							TPZFMatrix<REAL> &jacinv,
							TPZVec<FADREAL> &sol, TPZVec<FADREAL> &dsol,
							TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
							REAL weight, REAL timeStep,
							REAL deltaX);
	
	/**
	 * @brief Contributes the diffusion term to the tangent matrix (ek) and residual vector (ef) with FAD template classes
	 * @param[in] dim Spatial dimension
	 * @param[in] jacinv Inverse jacobian of the mapping.
	 * @param[in] sol Solution of the dim+2 state functions setup with derivatives
	 * @param[in] dsol Derivatives of U with respect to the dim dimensions setup with derivatives (xi components aligned as a vector)
	 * @param[in] phi List of shape functions.
	 * @param[in] dphi Derivatives of shape functions.
	 * @param[out] ek Tangent matrix to contribute to
	 * @param[out] ef Residual vector to contribute to
	 * @param[in] weight Gaussian quadrature integration weight
	 * @param[in] timeStep Time step
	 * @param[in] deltaX Diameter of element (used only if \f$fDelta > 0\f$);
	 */
	void ContributeFastestImplDiff(int dim, TPZFMatrix<REAL> &jacinv,
								   TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
								   TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
								   TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
								   REAL weight, REAL timeStep,
								   REAL deltaX);

	template <int dim>
	void ContributeFastestImplDiff_dim(
									   TPZFMatrix<REAL> &jacinv,
									   TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,
									   TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
									   TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,
									   REAL weight, REAL timeStep,
									   REAL deltaX);

	/** @} */
	
private:
	
	/** @brief Kind of artificial diffusion term to apply */
	TPZArtDiffType fArtDiffType;
	
	/** @brief Ratio between specific heat is constant and the specific heat the constant volume of a polytropic gas */
	REAL fGamma;
	
	/** @brief Scalar coefficient of the element in the diffusion term */
	REAL fDelta;
	
	/** @brief \f$CFL\f$ number */
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
			s = sqrt(c2 + us2 * REAL(16.));
			
			X.Redim(nstate,nstate);
			X(0,0) = 1.;
			X(0,2) = REAL(4.) * sol[0] * us / c2;
			X(0,3) = X(0,2);
			X(1,2) = - s/c -1.;
			X(1,3) =   s/c -1.;
			X(2,1) = 1.;
			X(3,2) = REAL(4.) * sol[0] * us;
			X(3,3) = X(3,2);
			
			Xi.Redim(nstate,nstate);
			Xi(0,0) = 1.;
			Xi(0,3) = REAL(- 1.) / c2;
			Xi(1,2) = 1.;
			Xi(2,1) = - c / (REAL(2.) * s);
			Xi(2,3) = (s-c)/(REAL(8.) * sol[0] * s * us);
			Xi(3,1) = -Xi(2,1);
			Xi(3,3) = (s+c)/(REAL(8.) * sol[0] * s * us);
			
			Lambda.Redim(nstate,nstate);
			Lambda(0,0) = REAL(1.)/us;
			Lambda(1,1) = REAL(1.)/sqrt(us2 + c2);
			Lambda(2,2) = REAL(1.)/sqrt(us2 + REAL(1.5) * c2 - REAL(.5) * c * s);
			Lambda(3,3) = 1./sqrt(us2 + REAL(1.5) * c2 + REAL(.5) * c * s);
			
			break;
		case(3):
			
			s = sqrt(c2 + us2 * REAL(4.));
			
			X.Redim(nstate,nstate);
			X(0,0) = 1.;
			X(0,3) = REAL(1.) / c2;
			X(0,4) = X(0,3);
			X(1,3) = (-s-c)/(REAL(2.) * rho_c_us);
			X(1,4) = ( s-c)/(REAL(2.) * rho_c_us);
			X(2,2) = 1.;
			X(3,1) = 1.;
			X(4,3) = 1.;
			X(4,4) = 1.;
			
			
			Xi.Redim(nstate,nstate);
			Xi(0,0) = 1.;
			Xi(0,4) = REAL(- 1.) / c2;
			Xi(1,3) = 1.;
			Xi(2,2) = 1.;
			Xi(3,1) = - rho_c_us / s;
			Xi(3,4) = REAL(.5) - c / (REAL(2.) * s);
			Xi(4,1) = -Xi(3,1);
			Xi(4,4) = REAL(.5) + c / (REAL(2.) * s);
			
			Lambda.Redim(nstate,nstate);
			Lambda(0,0) = REAL(1.)/us;
			Lambda(1,1) = REAL(1.)/sqrt(us2 + c2);
			Lambda(2,2) = Lambda(1,1);
			Lambda(3,3) = REAL(1.)/sqrt(us2 + REAL(2.) * c2 - c * s);
			Lambda(4,4) = REAL(1.)/sqrt(us2 + REAL(2.) * c2 + c * s);
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
			Y(0,2) = REAL(1.) / c2;
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
			Yi(1,3) = REAL(-1.) / c2;
			Yi(2,1) = - aaS[0] * rho_c / (REAL(2.) * k);
			Yi(2,2) = - aaS[1] * rho_c / (REAL(2.) * k);
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
			Y(0,3) = REAL(1.) / c2;
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
			Yi(2,4) = REAL(-1.) / c2;
			Yi(3,1) = - aaS[0] * rho_c / (REAL(2.) * k);
			Yi(3,2) = - aaS[1] * rho_c / (REAL(2.) * k);
			Yi(3,3) = - aaS[2] * rho_c / (REAL(2.) * k);
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
	
	T k, l1, l3, l4, l5, lstar,
	lsum, ldiff, twoCK, c2, temp1, temp2/*, temp3*/, k2;
	//T us2;
	
	c2 = c * c;
	//us2 = us * us;
	//rho_c = sol[0] * c;
	
	switch (dim)
	{
		case (2):
			k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1];
			k = sqrt(k2);
			twoCK = c * T(k * REAL(2.));
			
			l1 = aaS[0] * us;
			if(val(l1) < 0.)l1 = -l1;
			l3 = aaS[0] * us - k * c;
			if(val(l3) < 0.)l3 = -l3;
			l4 = aaS[0] * us + k * c;
			if(val(l4) < 0.)l4 = -l4;
			
			temp1 = (l4 - l3) * sol[0] / twoCK;
			temp2 = l3 + l4 - REAL(2.)*l1;
			Mat(0,0) += l1;
			Mat(0,1) += aaS[0] * temp1;
			Mat(0,2) += aaS[1] * temp1;
			Mat(0,3) += temp2 / (REAL(2.) * c2);
			Mat(1,1) += (aaS[1]*aaS[1] * l1 + aaS[0]*aaS[0] * (l3+l4) *REAL(.5))/k2;
			Mat(1,2) += aaS[0]*aaS[1] * temp2/(REAL(2.) * k2);
			Mat(1,3) += aaS[0]*(l4-l3)/(c * k * sol[0] * REAL(2.));
			Mat(2,1) += aaS[0]*aaS[1] * temp2/(REAL(2.) * k2);
			Mat(2,2) += (aaS[0]*aaS[0] * l1 + aaS[1]*aaS[1] * (l3+l4) *REAL(.5))/k2;
			Mat(2,3) += aaS[1]*(l4-l3)/(c * k * sol[0] * REAL(2.));
			Mat(3,1) += aaS[0] * c2 * temp1;
			Mat(3,2) += aaS[1] * c2 * temp1;
			Mat(3,3) += (l3 + l4)/REAL(2.);
			
			
			break;
		case (3):
			k2 = aaS[0]*aaS[0] + aaS[1]*aaS[1] + aaS[2]*aaS[2];
			k = sqrt(k2);
			
			l1 = aaS[0] * us;
			if(val(l1) < 0.)l1 = -l1;
			l4 = aaS[0] * us - k * c;
			if(val(l3) < 0.)l3 = -l3;
			l5 = aaS[0] * us + k * c;
			if(val(l4) < 0.)l4 = -l4;
			
			lsum = l5 + l4;
			ldiff= l5 - l4;
			lstar = lsum - REAL(2.)*l1;
			
			twoCK = c * k * REAL(2.);
			
			Mat.Redim(nstate,nstate);
			Mat(0,0) += l1;
			Mat(0,1) += aaS[0] * sol[0] * ldiff / twoCK;
			Mat(0,2) += aaS[1] * sol[0] * ldiff / twoCK;
			Mat(0,3) += aaS[2] * sol[0] * ldiff / twoCK;
			Mat(0,4) += lstar / ( c2 * REAL(2.) );
			Mat(1,1) += aaS[0] * aaS[0] * lstar / (k2 * REAL(2.)) + l1;
			Mat(1,2) += aaS[0] * aaS[1] * lstar / (k2 * REAL(2.));
			Mat(1,3) += aaS[0] * aaS[2] * lstar / (k2 * REAL(2.));
			Mat(1,4) += aaS[0] / twoCK / sol[0] * ldiff;
			Mat(2,1) += aaS[0] * aaS[1] * lstar / (k2 * REAL(2.));
			Mat(2,2) += aaS[1] * aaS[1] * lstar / (k2 * REAL(2.)) + l1;
			Mat(2,3) += aaS[1] * aaS[2] * lstar / (k2 * REAL(2.));
			Mat(2,4) += aaS[1] / twoCK / sol[0] * ldiff;
			Mat(3,1) += aaS[0] * aaS[2] * lstar / (k2 * REAL(2.));
			Mat(3,2) += aaS[1] * aaS[2] * lstar / (k2 * REAL(2.));
			Mat(3,3) += aaS[2] * aaS[2] * lstar / (k2 * REAL(2.)) + l1;
			Mat(3,4) += aaS[2] / twoCK / sol[0] * ldiff;
			Mat(4,1) += aaS[0] * c * sol[0] * ldiff / (k * REAL(2.));
			Mat(4,2) += aaS[1] * c * sol[0] * ldiff / (k * REAL(2.));
			Mat(4,3) += aaS[2] * c * sol[0] * ldiff / (k * REAL(2.));
			Mat(4,4) += lstar/REAL(2.) + l1;
			
			break;
		default:
			PZError << "TPZArtDiff::EigenSystemBornhaus Error: Invalid Dimension\n";
	}
}

#endif

