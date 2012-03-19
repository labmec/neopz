/**
 * \file
 * @brief Contains the TPZEulerConsLaw class which implements the weak statement of the compressible euler equations.
 */
//$Id: pzeulerconslaw.h,v 1.41 2009-09-01 22:07:16 phil Exp $

#ifndef EULERCONSLAW_H
#define EULERCONSLAW_H

#include <iostream>
#include "pzmaterial.h"
#include "tpzoutofrange.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzconslaw.h"
#include "pzartdiff.h"

#include "pzlog.h"

#ifdef LOG4CXX
extern LoggerPtr fluxroe;
extern LoggerPtr fluxappr;
#endif

#ifdef _AUTODIFF
#include "fadType.h"

#define _TFAD
//#define _FAD
//#define _TINYFAD // Doesn't work well -> apparently, derivatives do not accumulate well.
#endif

/**
 * @ingroup material
 * @brief This material implements the weak statement of the compressible euler equations
 */
class TPZEulerConsLaw  : public TPZConservationLaw
{
	public :
	
	TPZEulerConsLaw(int nummat,REAL timeStep,
					 REAL gamma,int dim,
					 TPZArtDiffType artdiff);
	
	~TPZEulerConsLaw();
	
	TPZEulerConsLaw();
	
	TPZEulerConsLaw(const TPZEulerConsLaw &cp) : TPZConservationLaw(cp), fArtDiff(cp.fArtDiff),fDiff(cp.fDiff),
	fConvVol(cp.fConvVol),fConvFace(cp.fConvFace)
	{
	}
	
	TPZAutoPointer<TPZMaterial> NewMaterial()
	{
		return new TPZEulerConsLaw(*this);
	}

	/**
	 * @brief Configures the time discretization of some contributions
	 * @param[in] Diff Refers to the diffusion term
	 * @param[in] ConvVol Refers to the volume contribution to the convective term
	 * @param[in] ConvFace Refers to the numerical flux
	 */
	void SetTimeDiscr(TPZTimeDiscr Diff, TPZTimeDiscr ConvVol, TPZTimeDiscr ConvFace);
	
	/** @brief Returns a reference to the artificial diffusion term */
	TPZArtDiff & ArtDiff(){ return fArtDiff; };
	
	/**
	 * @brief returns the best value for the CFL number based on the interpolation degree.
	 */
	REAL OptimalCFL(int degree);
	
	/** @brief See declaration in base class */
	virtual REAL SetTimeStep(REAL maxveloc,REAL deltax,int degree);
	
	/** @brief See declaration in base class */
	static int NStateVariables(int dim);
	
	/** @brief Object-based overload */
	int NStateVariables();
	
	/** @brief Estimates the deltax (element diameter) based on the inverse of the jacobian. */
	/**
	 * <!> Works only for quadratic and hexahedral elements
	 * (it assumes that each element is parametrized tensorially
	 * with functions varying from -1 to 1)
	 */
	REAL DeltaX(REAL detJac);
	
	/**
	 * @brief Computes the determinant of a 2d or 3d matrix. \n
	 * Used by recompute the element size
	 */
	REAL Det(TPZFMatrix<REAL> & Mat);
	
	/**
	 * @brief Thermodynamic pressure determined by the law of an ideal gas
	 * @param[in] gamma Gamma value
	 * @param[in] dim Spatial dimension
	 * @param[out] press Computed pressure
	 * @param[in] U Vector of state variables (sol)
	 */
	template< class T >
	static void Pressure(REAL gamma, int dim, T& press, TPZVec<T> &U);
	
	/** @brief Evaluates the speed of sound in the fluid */
	template <class T>
	static void cSpeed(TPZVec<T> & sol, REAL gamma, T & c);
	
	/** \f$ u = Sqrt(u2 + v2 + w2) \f$ */
	template <class T>
	static void uRes(TPZVec<T> & sol, T & us);
	
	virtual REAL Pressure(TPZVec<REAL> &U);
	
	virtual void Print(std::ostream & out);
	
	virtual std::string Name();
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes();
	
	
	/** @name Solutions methods */
	/** @{ */
	
	/** @brief Computes the ghost state variables bsed on the BC type */
	template <class T>
	void ComputeGhostState(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> &normal, TPZBndCond &bc, int & entropyFix);

protected:
	virtual void Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
public:
	
	/** @brief returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZConservationLaw::Solution(data,var,Solout);
	}
	/** @} */

	/** @name Fluxes methods */
	/** @{ */
	
	/** @brief tensor of the three-dimensional flux of Euler */
	template < class T >
	void Flux(TPZVec<T> &U,TPZVec<T> &Fx,TPZVec<T> &Fy,TPZVec<T> &Fz);
	
	/**
	 * @brief Jacobian of the tensor flux of Euler
	 * @param[in] gamma Gamma value to use
	 * @param[in] dim Spatial dimension
	 * @param[in] U \f$dim+2\f$ solutions at given point
	 * @param[out] Ai Vector of dim tensors of \f$ (dim+2)*(dim*2)\f$,
	 * representing the derivatives of F with respect to the 
	 * dim spatial dimensions
	 */
	template <class T>
	static void JacobFlux(REAL gamma, int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai);

	/** @brief Test flux -> returns the averaged state variables across an interface */
	template <class T>
	void Test_Flux(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> & normal, REAL gamma, TPZVec<T> & flux);

	/**
	 * @brief This flux encapsulates the two and three dimensional fluxes
	 * acquired from the Mouse program
	 * @param[in] solL Vector with solution from left element 
	 * @param[in] solR Vector with solution from right element
	 * @param[in] normal Normal vector of the interface
	 * @param[in] gamma Gamma value to use
	 * @param[out] flux Vector to return flux values
	 * @param entropyFix
	 */
	template <class T>
	static void Roe_Flux(TPZVec<T> &solL, TPZVec<T> &solR,
						 TPZVec<REAL> & normal, REAL gamma,
						 TPZVec<T> & flux, int entropyFix = 1);

	/**
	 * @brief This flux encapsulates the two and three dimensional fluxes acquired from the Mouse program
	 * @note This function is called Approx because it evaluates the derivative of the Roe Flux without using automatic differentiation
	 * @param[in] solL Vector with solution from left element 
	 * @param[in] solR Vector with solution from right element
	 * @param[in] normal Normal vector of the interface
	 * @param[in] gamma Gamma value to use
	 * @param[out] flux Vector to return flux values
	 * @param entropyFix
	 */
	template <class T>
	static void ApproxRoe_Flux(TPZVec<T> &solL, TPZVec<T> &solR,
							   TPZVec<REAL> & normal, REAL gamma,
							   TPZVec<T> & flux, int entropyFix = 1);
	
	
	/** @brief Sets the delta parameter inside the artifficial diffusion term. */
	void SetDelta(REAL delta);
	
public:
	/** @brief Flux of Roe (MOUSE program) */
	template <class T>
	static void Roe_Flux(const T & rho_f,
						 const T & rhou_f,
						 const T & rhov_f,
						 const T & rhow_f,
						 const T & rhoE_f,
						 const T & rho_t,
						 const T & rhou_t,
						 const T & rhov_t,
						 const T & rhow_t,
						 const T & rhoE_t,
						 const REAL nx,
						 const REAL ny,
						 const REAL nz,
						 const REAL gam,
						 T & flux_rho,
						 T & flux_rhou,
						 T & flux_rhov,
						 T & flux_rhow,
						 T & flux_rhoE, int entropyFix = 1);
	
	template <class T>
	static void Roe_Flux(const T & rho_f,
						 const T & rhou_f,
						 const T & rhov_f,
						 const T & rhoE_f,
						 const T & rho_t,
						 const T & rhou_t,
						 const T & rhov_t,
						 const T & rhoE_t,
						 const REAL nx,
						 const REAL ny,
						 const REAL gam,
						 T &flux_rho,
						 T &flux_rhou,
						 T &flux_rhov,
						 T &flux_rhoE, int entropyFix = 1);
	/** @brief Flux of Roe (MOUSE program) */
	template <class T>
	static void ApproxRoe_Flux(const T & rho_f,
							   const T & rhou_f,
							   const T & rhov_f,
							   const T & rhow_f,
							   const T & rhoE_f,
							   const T & rho_t,
							   const T & rhou_t,
							   const T & rhov_t,
							   const T & rhow_t,
							   const T & rhoE_t,
							   const REAL nx,
							   const REAL ny,
							   const REAL nz,
							   const REAL gam,
							   T & flux_rho,
							   T & flux_rhou,
							   T & flux_rhov,
							   T & flux_rhow,
							   T & flux_rhoE, int entropyFix = 1);
	
	template <class T>
	static void ApproxRoe_Flux(const T & rho_f,
							   const T & rhou_f,
							   const T & rhov_f,
							   const T & rhoE_f,
							   const T & rho_t,
							   const T & rhou_t,
							   const T & rhov_t,
							   const T & rhoE_t,
							   const REAL nx,
							   const REAL ny,
							   const REAL gam,
							   T &flux_rho,
							   T &flux_rhou,
							   T &flux_rhov,
							   T &flux_rhoE, int entropyFix = 1);
	/** @} */
		
#ifdef _AUTODIFF
	/** @name Differentiable variables setup */
	/** @{ */

	/**
	 * @brief Converts the REAL values into FAD differentiable objects
	 * @param[in] sol Values of flattened solution (Sum(ui*phi))
	 * @param[in] dsol Values of dsol(columns) with respect to the dim dimensions (row-aligned)
	 * @param[in] phi Value of shape functions at a given point
	 * @param[in] dphi Derivatives of shape functions with respect to the spatial dimensions at a given point.
	 * @param[out] FADsol Vector of solutions and coefficient derivatives (phi) at the point
	 * @param[out] FADdsol Vector of spatial derivatives of solution and coefficient derivatives (dphi) at the point
	 */
	void PrepareFAD(TPZVec<REAL> & sol, TPZFMatrix<REAL> & dsol,
					TPZFMatrix<REAL> & phi, TPZFMatrix<REAL> & dphi,
					TPZVec<FADREAL> & FADsol,
					TPZVec<FADREAL> & FADdsol);
	
	/**
	 * @brief Converts the REAL values into FAD differentiable objects
	 * with respect to the left and right interface volumes
	 * @param[in] solL Values of flattened solution (Sum(ui*phi))
	 * @param[in] solR Values of flattened solution (Sum(ui*phi))
	 * @param[in] phiL Value of shape functions at a given point
	 * @param[in] phiR Value of shape functions at a given point
	 * @param[out] FADsolL Vector of solutions and coefficient derivatives (phi) at the point
	 * @param[out] FADsolR Vector of solutions and coefficient derivatives (phi) at the point
	 */
	void PrepareInterfaceFAD(
							 TPZVec<REAL> &solL,TPZVec<REAL> &solR,
							 TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
							 TPZVec<FADREAL> & FADsolL,
							 TPZVec<FADREAL> & FADsolR);
	
	/**
	 * @brief Converts the REAL values into FAD differentiable objects
	 * with respect to the left and right interface volumes.
	 */
	/** The derivatives are only with respect to the left and
	 * right state variables, not the shape functions.
	 * It returns a TFADET class type
	 */
	template <class T>
	void PrepareFastestInterfaceFAD(
									TPZVec<REAL> &solL,TPZVec<REAL> &solR,
									TPZVec<T> & FADsolL,
									TPZVec<T> & FADsolR);
	/** @} */
	
#endif
	
	/** @name Contributions methods */
	/** @{ */
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef);
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ef);

	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc)
	{
    	TPZDiscontinuousGalerkin::ContributeBC(data,weight,ef,bc);
	}
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
									   TPZBndCond &bc);
	
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
									   REAL weight,
									   TPZFMatrix<REAL> &ef,
									   TPZBndCond &bc);
	
	/** @name Internal contributions */
	/** @{ */
	
	void ContributeFastestBCInterface(int dim,
									  TPZVec<REAL> &x,
									  TPZVec<REAL> &solL, TPZFMatrix<REAL> &dsolL,
									  REAL weight, TPZVec<REAL> &normal,
									  TPZFMatrix<REAL> &phiL, TPZFMatrix<REAL> &dphiL,
									  TPZFMatrix<REAL> &axesleft,
									  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, TPZBndCond &bc);
	
	template <int dim>
	void ContributeFastestBCInterface_dim(
										  TPZVec<REAL> &x,
										  TPZVec<REAL> &solL, TPZFMatrix<REAL> &dsolL,
										  REAL weight, TPZVec<REAL> &normal,
										  TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &dphiL,
										  TPZFMatrix<REAL> &axesleft,
										  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc);

	virtual void ContributeLast(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
								TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
								REAL weight,
								TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
								TPZFMatrix<REAL> &ef);
	
	virtual void ContributeAdv(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
							   TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
							   REAL weight,
							   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
							   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	virtual void ContributeAdv(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
							   TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
							   REAL weight,
							   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
							   TPZFMatrix<REAL> &ef);
	
	void ContributeApproxImplDiff(TPZVec<REAL> &x,
								  TPZFMatrix<REAL> &jacinv,
								  TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
								  REAL weight,
								  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
								  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	void ContributeExplDiff(TPZVec<REAL> &x,
							TPZFMatrix<REAL> &jacinv,
							TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
							REAL weight,
							TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
							TPZFMatrix<REAL> &ef);
	
#ifdef _AUTODIFF
	void ContributeImplDiff(TPZVec<REAL> &x,
							TPZFMatrix<REAL> &jacinv,
							TPZVec<FADREAL> &sol,TPZVec<FADREAL> &dsol,
							REAL weight,
							TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	void ContributeFastestImplDiff(int dim,
								   TPZVec<REAL> &x,
								   TPZFMatrix<REAL> &jacinv,
								   TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
								   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
								   REAL weight,
								   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
#endif
	
	void ContributeExplConvFace(TPZVec<REAL> &x,
								TPZVec<REAL> &solL,TPZVec<REAL> &solR,
								REAL weight,TPZVec<REAL> &normal,
								TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
								TPZFMatrix<REAL> &ef, int entropyFix = 1);
	
	void ContributeApproxImplConvFace(TPZVec<REAL> &x, REAL faceSize,
									  TPZVec<REAL> &solL,TPZVec<REAL> &solR,
									  REAL weight,TPZVec<REAL> &normal,
									  TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
									  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix = 1
									  );
	
#ifdef _AUTODIFF
	
	void ContributeApproxImplConvFace(TPZVec<REAL> &x, REAL faceSize,
									  TPZVec<FADREAL> &solL,TPZVec<FADREAL> &solR,
									  REAL weight,TPZVec<REAL> &normal,
									  TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
									  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix = 1
									  );
	
	void ContributeImplConvFace(TPZVec<REAL> &x,
								TPZVec<FADREAL> &solL,TPZVec<FADREAL> &solR,
								REAL weight,TPZVec<REAL> &normal,
								TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
								TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix = 1);

	void ContributeFastestImplConvFace(int dim,
									   TPZVec<REAL> &x,
									   TPZVec<REAL> &solL,TPZVec<REAL> &solR,
									   REAL weight,TPZVec<REAL> &normal,
									   TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
									   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix = 1);

	template <int dim>
	void ContributeFastestImplConvFace_dim(
										   TPZVec<REAL> &x,
										   TPZVec<REAL> &solL,TPZVec<REAL> &solR,
										   REAL weight,TPZVec<REAL> &normal,
										   TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
										   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix = 1);

	template <class T>
	void ContributeFastestImplConvFace_T(TPZVec<REAL> &x,
										 TPZVec<T> &FADsolL,TPZVec<T> &FADsolR,
										 REAL weight,TPZVec<REAL> &normal,
										 TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
										 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,int entropyFix = 1);

#endif
	
	void ContributeImplConvVol(TPZVec<REAL> &x,
							   TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
							   REAL weight,
							   TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
							   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	void ContributeExplConvVol(TPZVec<REAL> &x,
							   TPZVec<REAL> &sol,
							   REAL weight,
							   TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
							   TPZFMatrix<REAL> &ef);
	
#ifdef _AUTODIFF
	void ContributeImplConvVol(TPZVec<REAL> &x,
							   TPZVec<FADREAL> &sol,TPZVec<FADREAL> &dsol,
							   REAL weight,
							   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
#endif
	
	void ContributeExplT1(TPZVec<REAL> &x,
						  TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
						  REAL weight,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
						  TPZFMatrix<REAL> &ef);
	
	void ContributeImplT1(TPZVec<REAL> &x,
						  TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
						  REAL weight,
						  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
						  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef);
	
	void ContributeExplT2(TPZVec<REAL> &x,
						  TPZVec<REAL> &sol,
						  REAL weight,
						  TPZFMatrix<REAL> &phi,
						  TPZFMatrix<REAL> &ef);
	/** @} */
	
	/** @} */
	
	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	void Read(TPZStream &buf, void *context);
	
	/** @brief Class identificator */
	int ClassId() const;

	/** @name Attributes */
	/** @{ */
protected:
	
	/** @brief diffusive term */
	TPZArtDiff fArtDiff;
	
	/** @brief variables indication whether the following terms are implicit */
	TPZTimeDiscr fDiff, fConvVol, fConvFace;
	
	/** @} */
};

inline std::string TPZEulerConsLaw::Name()
{
	return "TPZEulerConsLaw";
}

template < class T >
inline void TPZEulerConsLaw::Flux(TPZVec<T> &U,TPZVec<T> &Fx,TPZVec<T> &Fy,TPZVec<T> &Fz) {
	
	T press;
	Pressure(fGamma, fDim, press, U);
	int nstate = NStateVariables();
	if(nstate < 3 && nstate > 5){
		PZError << "TPZEulerConsLaw::Flux case not implemented\n";
		Fx.Resize(0);
		Fy.Resize(0);
		Fz.Resize(0);
		exit(-1);
	}
	
	Fx.Resize(5,0.0);//v�ido
	Fy.Resize(5,0.0);//para
	Fz.Resize(5,0.0);//R , R , R
	
	if(nstate == 5){
		Fx[0] = U[1];//ro u
		Fx[1] = (U[1]/U[0])*U[1] + press;//ro u2 + p
		Fx[2] = U[1]*(U[2]/U[0]);//ro u v
		Fx[3] = U[1]*(U[3]/U[0]);//ro u w
		Fx[4] = (U[4]+press)*(U[1]/U[0]);//(ro e + p) u
		
		Fy[0] = U[2];//ro v
		Fy[1] = U[2]*(U[1]/U[0]);//ro v u
		Fy[2] = (U[2]/U[0])*U[2] + press;//ro v2 + p
		Fy[3] = U[2]*(U[3]/U[0]);//ro v w
		Fy[4] = (U[4] + press)*(U[2]/U[0]);//(ro e + p) v
		
		Fz[0] = U[3];//ro w
		Fz[1] = U[3]*(U[1]/U[0]);//ro w u
		Fz[2] = U[3]*(U[2]/U[0]);//ro w v
		Fz[3] = (U[3]/U[0])*U[3] + press;//ro w2 + p
		Fz[4] = (U[4] + press)*(U[3]/U[0]);//(ro e + p) w
		return;
	}
	
	if(nstate == 4){
		Fx[0] = U[1];//ro u
		Fx[1] = U[1]*U[1] / U[0] + press;//ro u2 + p
		Fx[2] = U[1]*U[2] / U[0];//ro u v
		Fx[3] = (U[3]+press)*U[1] / U[0];//(E + p) u
		
		Fy[0] = U[2];//ro v
		Fy[1] = U[1]*U[2] / U[0];//ro u v
		Fy[2] = U[2]*U[2] / U[0] + press;//ro v2 + p
		Fy[3] = (U[3] + press)*U[2] / U[0];//(E + p) v
		return;
	}
	
	if(nstate == 3){
		Fx[0] = U[1];//ro u
		Fx[1] = (U[1]/U[0])*U[1] + press;//ro u2 + p
		Fx[2] = (U[2]+press)*(U[1]/U[0]);//(ro e + p) u
	}
}


template <class T>
inline void TPZEulerConsLaw::JacobFlux(REAL gamma, int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai)
{
	
	Ai.Resize(dim);
	int i;
	for(i=0;i<dim;i++)Ai[i].Redim(TPZEulerConsLaw::NStateVariables(dim), TPZEulerConsLaw::NStateVariables(dim));
	
	if(U[0] < REAL(1.e-6)) {
		PZError << "TPZEulerConsLaw::JacobFlux: Density negative "   << U[0] << std::endl;
		TPZOutofRange obj;
		throw(obj);
	}
	
	T    u,v,w,e;
	REAL gamma1 = gamma-1.;
	REAL gamma2 = gamma1/2.;
	REAL gamma3 = gamma-3;
	
	if(dim == 3){
		
		u = U[1]/U[0];
		v = U[2]/U[0];
		w = U[3]/U[0];
		e = U[4]/U[0];
		
		T    u2 = u*u;
		T    v2 = v*v;
		T    w2 = w*w;
		T    vel = u2+v2+w2;
		
		Ai[0](0,0) = 0.;
		Ai[0](0,1) = 1.;
		Ai[0](0,2) = 0.;
		Ai[0](0,3) = 0.;
		Ai[0](0,4) = 0.;
		
		Ai[0](1,0) =  gamma2*vel-u2;
		Ai[0](1,1) = -gamma3*u;
		Ai[0](1,2) = -gamma1*v;
		Ai[0](1,3) = -gamma1*w;
		Ai[0](1,4) =  gamma1;
		
		Ai[0](2,0) = -u*v;
		Ai[0](2,1) =  v;
		Ai[0](2,2) =  u;
		Ai[0](2,3) =  0.;
		Ai[0](2,4) =  0.;
		
		Ai[0](3,0) = -u*w;
		Ai[0](3,1) =  w;
		Ai[0](3,2) =  0.;
		Ai[0](3,3) =  u;
		Ai[0](3,4) =  0.;
		
		Ai[0](4,0) = -gamma*e*u + gamma1*u*vel;
		Ai[0](4,1) =  gamma*e - gamma1*u2 - gamma2*vel;
		Ai[0](4,2) = -gamma1*u*v;
		Ai[0](4,3) = -gamma1*u*w;
		Ai[0](4,4) =  gamma*u;
		
		Ai[1](0,0) = 0.;
		Ai[1](0,1) = 0.;
		Ai[1](0,2) = 1.;
		Ai[1](0,3) = 0.;
		Ai[1](0,4) = 0.;
		
		Ai[1](1,0) = -u*v;
		Ai[1](1,1) =  v;
		Ai[1](1,2) =  u;
		Ai[1](1,3) =  0.;
		Ai[1](1,4) =  0.;
		
		Ai[1](2,0) =  gamma2*vel-v2;
		Ai[1](2,1) = -gamma1*u;
		Ai[1](2,2) = -gamma3*v;
		Ai[1](2,3) = -gamma1*w;
		Ai[1](2,4) =  gamma1;
		
		Ai[1](3,0) = -v*w;
		Ai[1](3,1) =  0.;
		Ai[1](3,2) =  w;
		Ai[1](3,3) =  v;
		Ai[1](3,4) =  0.;
		
		Ai[1](4,0) = -gamma*e*v + gamma1*v*vel;
		Ai[1](4,1) = -gamma1*u*v;
		Ai[1](4,2) =  gamma*e - gamma1*v2 - gamma2*vel;
		Ai[1](4,3) = -gamma1*v*w;
		Ai[1](4,4) =  gamma*v;
		
		Ai[2](0,0) = 0.;
		Ai[2](0,1) = 0.;
		Ai[2](0,2) = 0.;
		Ai[2](0,3) = 1.;
		Ai[2](0,4) = 0.;
		
		Ai[2](1,0) = -u*w;
		Ai[2](1,1) =  w;
		Ai[2](1,2) =  0.;
		Ai[2](1,3) =  u;
		Ai[2](1,4) =  0.;
		
		Ai[2](2,0) = -v*w;
		Ai[2](2,1) =  0.;
		Ai[2](2,2) =  w;
		Ai[2](2,3) =  v;
		Ai[2](2,4) =  0.;
		
		Ai[2](3,0) =  gamma2*vel-w2;
		Ai[2](3,1) = -gamma1*u;
		Ai[2](3,2) = -gamma1*v;
		Ai[2](3,3) = -gamma3*w;
		Ai[2](3,4) =  gamma1;
		
		Ai[2](4,0) = -gamma*e*w + gamma1*w*vel;
		Ai[2](4,1) = -gamma1*u*w;
		Ai[2](4,2) = -gamma1*v*w;
		Ai[2](4,3) =  gamma*e - gamma1*w2 - gamma2*vel;
		Ai[2](4,4) =  gamma*w;
		
	} else if(dim == 2){
		
		u = U[1]/U[0];
		v = U[2]/U[0];
		e = U[3]/U[0];
		
		T    u2 = u*u;
		T    v2 = v*v;
		T    vel = u2+v2;
		
		Ai[0](0,0) = 0.;
		Ai[0](0,1) = 1.;
		Ai[0](0,2) = 0.;
		Ai[0](0,3) = 0.;
		
		Ai[0](1,0) =  gamma2*vel-u2;
		Ai[0](1,1) = -gamma3*u;
		Ai[0](1,2) = -gamma1*v;
		Ai[0](1,3) =  gamma1;
		
		Ai[0](2,0) = -u*v;
		Ai[0](2,1) =  v;
		Ai[0](2,2) =  u;
		Ai[0](2,3) =  0.;
		
		Ai[0](3,0) = -gamma*e*u + gamma1*u*vel;
		Ai[0](3,1) =  gamma*e - gamma1*u2 - gamma2*vel;
		Ai[0](3,2) = -gamma1*u*v;
		Ai[0](3,3) =  gamma*u;
		
		Ai[1](0,0) = 0.;
		Ai[1](0,1) = 0.;
		Ai[1](0,2) = 1.;
		Ai[1](0,3) = 0.;
		
		Ai[1](1,0) = -u*v;
		Ai[1](1,1) =  v;
		Ai[1](1,2) =  u;
		Ai[1](1,3) =  0.;
		
		Ai[1](2,0) =  gamma2*vel-v2;
		Ai[1](2,1) = -gamma1*u;
		Ai[1](2,2) = -gamma3*v;
		Ai[1](2,3) =  gamma1;
		
		Ai[1](3,0) = -gamma*e*v + gamma1*v*vel;
		Ai[1](3,1) = -gamma1*u*v;
		Ai[1](3,2) =  gamma*e - gamma1*v2 - gamma2*vel;
		Ai[1](3,3) =  gamma*v;
		
	} else if(dim == 1){
		
		u = U[1]/U[0];
		e = U[2]/U[0];
		
		T    u2 = u*u;
		T    vel = u2;
		
		Ai[0](0,0) = 0.;
		Ai[0](0,1) = 1.;
		Ai[0](0,2) = 0.;
		
		Ai[0](1,0) =  gamma2*vel-u2;
		Ai[0](1,1) = -gamma3*u;
		Ai[0](1,2) =  gamma1;
		
		Ai[0](2,0) = -gamma*e*u + gamma1*u*vel;
		Ai[0](2,1) =  gamma*e - gamma1*u2 - gamma2*vel;
		Ai[0](2,2) =  gamma*u;
	}
}



#ifdef _AUTODIFF
template <class T>
inline REAL val(T & number)
{
	return number.val();
}
#endif

template< class T >
inline void TPZEulerConsLaw::Pressure(REAL gamma, int dim, T & press, TPZVec<T> &U)
{
	if(fabs(val(U[0])) < 1.e-6) {
		PZError << "\nTPZEulerConsLaw::Pressure> Density negative "
		<< U[0] << std::endl;
		TPZOutofRange obj;
		throw(obj);
		//    exit(-1);
	}
	// Press� = (gam-1)*(E - ro*||(u,v,w)||/2)
	// onde aqui ro_e = E (nota�)
	
	int nstate = NStateVariables(dim);
	press = 0.0;
	if(U.NElements() != nstate) U.Resize(nstate);
	if(nstate == 5){
		//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro w , ro e)
		T rho_velocity = ( U[1]*U[1] + U[2]*U[2] + U[3]*U[3] )/U[0];
		press = ((gamma-1.)*( U[4] - REAL(0.5) * rho_velocity ));
	} else
		if(nstate == 4){
			//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro e)
			T rho_velocity = ( U[1]*U[1] + U[2]*U[2] )/U[0];
			press = ((gamma-1.)*( U[3] - REAL(0.5) * rho_velocity ));
		} else
			if(nstate == 3){
				//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro e)
				T rho_velocity = ( U[1]*U[1] )/U[0];
				press = ((gamma-1.)*( U[2] - REAL(0.5) * rho_velocity ));
			} else {
				std::cout << "\nTPZEulerConsLaw::Pressure> Unknown case - returning zero\n";
				press = 0.0;
				return;
			}
	if(val(press) < 0){
		T temp = ((T)(gamma-1.))*U[nstate-1];
		PZError << "TPZEulerConsLaw::Pressure> Negative pressure: " << press << " (gama-1)*E = " << temp << std::endl;
		TPZOutofRange obj;
		throw(obj);
	}
}

//----------------Test Flux

template <class T>
void TPZEulerConsLaw::Test_Flux(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> & normal, REAL gamma, TPZVec<T> & flux)
{
	int nState = flux.NElements();
	for(int i = 0; i < nState; i++)
	{
		flux[i] = (solR[i]-solL[i]);
	}
}


//----------------Roe Flux

template <class T>
void TPZEulerConsLaw::Roe_Flux(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> & normal, REAL gamma, TPZVec<T> & flux, int entropyFix)
{
	// Normals outgoing from the BC elements into the
	// mesh elements -> all the normals are opposited to
	// the common convention -> changing the left/right
	// elements and normals.
	int nState = solL.NElements();
	if(nState == 5)
	{
		Roe_Flux(solL[0], solL[1], solL[2], solL[3], solL[4],
				 solR[0], solR[1], solR[2], solR[3], solR[4],
				 normal[0], normal[1], normal[2],
				 gamma,
				 flux[0], flux[1], flux[2], flux[3], flux[4], entropyFix);
		
	}else if(nState == 4)
	{
		Roe_Flux(solL[0], solL[1], solL[2], solL[3],
				 solR[0], solR[1], solR[2], solR[3],
				 normal[0], normal[1],
				 gamma,
				 flux[0], flux[1], flux[2], flux[3], entropyFix);
	}else if(nState == 3)
	{
		//using the 2D expression for 1d problem
		T auxL = REAL(0.),
		auxR = REAL(0.),
		fluxaux = REAL(0.);
		auxL = flux[0];
		Roe_Flux(solL[0], solL[1], auxL, solL[2],
				 solR[0], solR[1], auxR, solR[2],
				 normal[0], 0,
				 gamma,
				 flux[0], flux[1], fluxaux, flux[2], entropyFix);
	}else
	{
		PZError << "No flux on " << nState << " state variables.\n";
	}
	
}

//left = **_f    right = **_t
template <class T>
inline void TPZEulerConsLaw::Roe_Flux(
									   const T & rho_f, const T & rhou_f, const T & rhov_f, const T & rhow_f,
									   const T & rhoE_f, const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhow_t,
									   const T & rhoE_t, const REAL nx, const REAL ny, const REAL nz, const REAL gam,
									   T & flux_rho, T &flux_rhou, T &flux_rhov,
									   T & flux_rhow, T &flux_rhoE, int entropyFix){
	
	T    alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
	T    a1,a2,a3,a4,a5,b1,b2,b3,b4,b5;
	T    ep_t, ep_f, p_t, p_f;
	T    rhouv_t, rhouv_f, rhouw_t, rhouw_f, rhovw_t, rhovw_f;
	T    lambda_f, lambda_t;
	T    delta_rho, delta_rhou, delta_rhov, delta_rhow, delta_rhoE;
	T    hnx, hny, hnz;
	T    tempo11, usc;
	
	flux_rho = 0;
	flux_rhou = 0;
	flux_rhov = 0;
	flux_rhow = 0;
	flux_rhoE = 0;
	
	REAL gam1 = gam - 1.0;
	T    irho_f = REAL(1.0)/rho_f;
	T    irho_t = REAL(1.0)/rho_t;
	
	//
	//.. Compute the ROE Averages
	//
	//.... some useful quantities
	T    coef1 = sqrt(rho_f);
	T    coef2 = sqrt(rho_t);
	T    somme_coef = coef1 + coef2;
	T    isomme_coef = REAL(1.0)/somme_coef;
	T    u_f = rhou_f*irho_f;
	T    v_f = rhov_f*irho_f;
	T    w_f = rhow_f*irho_f;
	T    h_f = (gam * rhoE_f*irho_f) - (.5*gam1) * (u_f * u_f + v_f * v_f + w_f * w_f);
	T    u_t = rhou_t*irho_t;
	T    v_t = rhov_t*irho_t;
	T    w_t = rhow_t*irho_t;
	T    h_t = (gam * rhoE_t*irho_t) - (.5*gam1) * (u_t * u_t + v_t * v_t + w_t * w_t);
	
	//.... averages
	//REAL rho_ave = coef1 * coef2;
	T    u_ave = (coef1 * u_f + coef2 * u_t) * isomme_coef;
	T    v_ave = (coef1 * v_f + coef2 * v_t) * isomme_coef;
	T    w_ave = (coef1 * w_f + coef2 * w_t) * isomme_coef;
	T    h_ave = (coef1 * h_f + coef2 * h_t) * isomme_coef;
	//
	//.. Compute Speed of sound
	T    scal = u_ave * nx + v_ave * ny + w_ave * nz;
	T    norme = sqrt(nx * nx + ny * ny + nz * nz);
	T    inorme = REAL(1.0)/norme;
	T    u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
	T    c_speed = gam1 * (h_ave - REAL(0.5) * u2pv2pw2);
	if(c_speed < REAL(1e-6)) c_speed = 1e-6;// <!> zeroes the derivatives?   // avoid division by 0 if critical
	c_speed = sqrt(c_speed);
	T    c_speed2 = c_speed * norme;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	T    eig_val1 = scal - c_speed2;
	T    eig_val2 = scal;
	T    eig_val3 = scal + c_speed2;
	//
	//.. Compute the ROE flux
	//.... In this part many tests upon the eigenvalues
	//.... are done to simplify calculations
	//.... Here we use the two formes of the ROE flux :
	//.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
	//.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
	//
	if(eig_val2 <= REAL(0.0)) {
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t +
											rhov_t * rhov_t + rhow_t * rhow_t) * irho_t);
		ep_t = rhoE_t + p_t;
		rhouv_t = rhou_t * v_t;
		rhouw_t = rhou_t * w_t;
		rhovw_t = rhov_t * w_t;
		flux_rho  = rhou_t * nx + rhov_t * ny + rhow_t * nz;
		flux_rhou = (rhou_t * u_t + p_t) * nx + rhouv_t * ny + rhouw_t * nz;
		flux_rhov = rhouv_t * nx + (rhov_t * v_t + p_t) * ny + rhovw_t * nz;
		flux_rhow = rhouw_t * nx + rhovw_t * ny + (rhow_t * w_t + p_t) * nz;
		flux_rhoE = ep_t * (u_t * nx + v_t * ny + w_t * nz);
		//
		//.... A Entropic modification
		//
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f + rhov_f * rhov_f
											+ rhow_f * rhow_f) * irho_f);
		lambda_f = u_f * nx + v_f * ny + w_f * nz + norme
		* sqrt(gam * p_f * irho_f);
		lambda_t = u_t * nx + v_t * ny + w_t * nz + norme
		* sqrt(gam * p_t * irho_t);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
		}
		//
		if (eig_val3 > REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the last
			//.. colomne of T-1 with the last row of T with eig_val3                //Cedric
			delta_rho  = rho_t - rho_f;                                             //right - left
			delta_rhou = rhou_t - rhou_f;                                           //**_t  - **_f
			delta_rhov = rhov_t - rhov_f;
			delta_rhow = rhow_t - rhow_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal * inorme;
			hnx = nx * inorme;
			hny = ny * inorme;
			hnz = nz * inorme;
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. Last columne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc + hnx;
			a3 = v_ave * usc + hny;
			a4 = w_ave * usc + hnz;
			a5 = REAL(0.5) * u2pv2pw2 * usc + REAL(2.5) * c_speed + scal;
			//.. Last row of the matrix T * eig_val3
			b1 = REAL(0.5) * (REAL(0.5) * tempo11 * u2pv2pw2 - scal);
			b2 = REAL(0.5) * (hnx - tempo11 * u_ave);
			b3 = REAL(0.5) * (hny - tempo11 * v_ave);
			b4 = REAL(0.5) * (hnz - tempo11 * w_ave);
			b5 = REAL(0.5) * tempo11;
			//
			alpha1 = b1 * delta_rho;
			alpha2 = b2 * delta_rhou;
			alpha3 = b3 * delta_rhov;
			alpha4 = b4 * delta_rhow;
			alpha5 = b5 * delta_rhoE;
			alpha  = eig_val3 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
			//
			flux_rho  -= a1 * alpha;
			flux_rhou -= a2 * alpha;
			flux_rhov -= a3 * alpha;
			flux_rhow -= a4 * alpha;
			flux_rhoE -= a5 * alpha;
		}
	}
	//
	if(eig_val2 > REAL(0.0)) {
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f +
											rhov_f * rhov_f + rhow_f * rhow_f) * irho_f);
		ep_f = rhoE_f + p_f;
		rhouv_f = rhou_f * v_f;
		rhouw_f = rhou_f * w_f;
		rhovw_f = rhov_f * w_f;
		flux_rho  = rhou_f * nx + rhov_f * ny + rhow_f * nz;
		flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny + rhouw_f * nz;
		flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny + rhovw_f * nz;
		flux_rhow = rhouw_f * nx + rhovw_f * ny + (rhow_f * w_f + p_f) * nz;
		flux_rhoE = ep_f * (u_f * nx + v_f * ny + w_f * nz);
		//
		// A Entropic modification
		//
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t +
											rhov_t * rhov_t + rhow_t * rhow_t) * irho_t);
		lambda_f = u_f * nx + v_f * ny + w_f * nz - norme
		* sqrt(gam * p_f * irho_f);
		lambda_t   = u_t * nx + v_t * ny + w_t * nz - norme
		* sqrt(gam * p_t * irho_t);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
		}
		//
		if (eig_val1 < REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the first
			//.. columne of T-1 with the first row of T with eig_val1
			delta_rho  = rho_t - rho_f;
			delta_rhou = rhou_t - rhou_f;
			delta_rhov = rhov_t - rhov_f;
			delta_rhow = rhow_t - rhow_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal * inorme;
			hnx = nx * inorme;
			hny = ny * inorme;
			hnz = nz * inorme;
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. First colomne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc - hnx;
			a3 = v_ave * usc - hny;
			a4 = w_ave * usc - hnz;
			a5 = REAL(0.5) * u2pv2pw2 * usc + REAL(2.5) * c_speed - scal;
			//.. First row of the matrix T * eig_val1
			b1 = REAL(0.5) * (REAL(0.5) * tempo11 * u2pv2pw2 + scal);
			b2 = -REAL(0.5) * (hnx + tempo11 * u_ave);
			b3 = -REAL(0.5) * (hny + tempo11 * v_ave);
			b4 = -REAL(0.5) * (hnz + tempo11 * w_ave);
			b5 = REAL(0.5) * tempo11;
			//
			alpha1 = b1 * delta_rho;
			alpha2 = b2 * delta_rhou;
			alpha3 = b3 * delta_rhov;
			alpha4 = b4 * delta_rhow;
			alpha5 = b5 * delta_rhoE;
			alpha  = eig_val1 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
			//
			flux_rho  += a1 * alpha;
			flux_rhou += a2 * alpha;
			flux_rhov += a3 * alpha;
			flux_rhow += a4 * alpha;
			flux_rhoE += a5 * alpha;
		}
	}
}

//left = **_f    right = **_t
template <class T>
inline void TPZEulerConsLaw::Roe_Flux(const T & rho_f, const T & rhou_f, const T & rhov_f, const T & rhoE_f,
									   const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhoE_t,
									   const REAL nx, const REAL ny, const REAL gam,
									   T & flux_rho, T & flux_rhou,T & flux_rhov, T & flux_rhoE, int entropyFix){
	
	T    alpha1,alpha2,alpha3,alpha4,a1,a2,a3,a4,b1,b2,b3,b4,alpha;
	T    ep_t, ep_f, p_t, p_f;
	T    rhouv_t, rhouv_f;
	T    lambda_f, lambda_t;
	T    delta_rho, delta_rhou,delta_rhov, delta_rhoE;
	REAL hnx, hny;
	T    tempo11, usc;
	
	flux_rho = 0;
	flux_rhou = 0;
	flux_rhov = 0;
	flux_rhoE = 0;
	
	REAL gam1 = gam - REAL(1.0);
	//REAL gam2 = gam * (gam - 1.0);
	//REAL igam = 1.0 / (gam - 1.0);
	
	//
	//.. Compute the ROE Averages
	//
	//.... some useful quantities
	T    coef1 = sqrt(rho_f);
	T    coef2 = sqrt(rho_t);
	T    somme_coef = coef1 + coef2;
	T    u_f = rhou_f/rho_f;
	T    v_f = rhov_f/rho_f;
	T    h_f = (gam * rhoE_f/rho_f) -  (u_f * u_f + v_f * v_f) * (gam1 / REAL(2.0));
	T    u_t = rhou_t/rho_t;
	T    v_t = rhov_t/rho_t;
	T    h_t = (gam * rhoE_t/rho_t) -  (u_t * u_t + v_t * v_t) * (gam1 / REAL(2.0));
	
	//.... averages
	//REAL rho_ave = coef1 * coef2;
	T    u_ave = (coef1 * u_f + coef2 * u_t) / somme_coef;
	T    v_ave = (coef1 * v_f + coef2 * v_t) / somme_coef;
	T    h_ave = (coef1 * h_f + coef2 * h_t) / somme_coef;
	//
	//  cout << "Correct coef1 " << coef1 << "coef2 " << coef2 << "h_f " << h_f << "h_t " << h_t << endl;
	
	//.. Compute Speed of sound
	T    scal = u_ave * nx + v_ave * ny;
	REAL norme = sqrt(nx * nx + ny * ny);
	T    u2pv2 = u_ave * u_ave + v_ave * v_ave;
	T    c_speed = gam1 * (h_ave - REAL(0.5) * u2pv2);
	//  cout << "c_speed " << c_speed << endl << "h_ave " << h_ave << endl << "u2pv2 " << u2pv2 << endl;
	
	if(c_speed < REAL(1e-6)) c_speed = REAL(1e-6);    // avoid division by 0 if critical
	c_speed = sqrt(c_speed);
	T    c_speed2 = c_speed * norme;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	T    eig_val1 = scal - c_speed2;
	T    eig_val2 = scal;
	T    eig_val3 = scal + c_speed2;
	//  cout << "Eigenvalues correct" << eig_val1 << endl << eig_val2 << endl <<
	//      eig_val3 << endl;
	
	//
	//.. Compute the ROE flux
	//.... In this part many tests upon the eigenvalues
	//.... are done to simplify calculations
	//.... Here we use the two formes of the ROE flux :
	//.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
	//.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
	//
	if(eig_val2 <= REAL(0.0)) {
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t + rhov_t * rhov_t) / rho_t);
		ep_t = rhoE_t + p_t;
		rhouv_t = rhou_t * v_t;
		flux_rho  = rhou_t * nx + rhov_t * ny;
		flux_rhou = (rhou_t * u_t + p_t) * nx + rhouv_t * ny;
		flux_rhov = rhouv_t * nx + (rhov_t * v_t + p_t) * ny;
		flux_rhoE = ep_t * (u_t * nx + v_t * ny);
		//
		//.... A Entropic modification
		//
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f + rhov_f * rhov_f) / rho_f);
		lambda_f = u_f * nx + v_f * ny + norme * sqrt(gam * p_f / rho_f);
		lambda_t   = u_t * nx + v_t * ny + norme
		* sqrt(gam * p_t / rho_t);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
		}
		//
		if (eig_val3 > REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the last
			//.. colomne of T-1 with the last row of T with eig_val3
			delta_rho  = rho_t - rho_f;
			delta_rhou = rhou_t - rhou_f;
			delta_rhov = rhov_t - rhov_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal / norme;
			hnx = nx / norme;
			hny = ny / norme;
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. Last columne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc + hnx;
			a3 = v_ave * usc + hny;
			a4 = REAL(0.5) * u2pv2 * usc + REAL(2.5) * c_speed + scal;
			//.. Last row of the matrix T * eig_val3
			b1 = REAL(0.5) * eig_val3 * (REAL(0.5) * tempo11 * u2pv2 - scal);
			b2 = REAL(0.5) * eig_val3 * (hnx - tempo11 * u_ave);
			b3 = REAL(0.5) * eig_val3 * (hny - tempo11 * v_ave);
			b4 = REAL(0.5) * eig_val3 * tempo11;
			//
			alpha1 = a1 * b1 * delta_rho;
			alpha2 = a1 * b2 * delta_rhou;
			alpha3 = a1 * b3 * delta_rhov;
			alpha4 = a1 * b4 * delta_rhoE;
			alpha = alpha1 + alpha2 + alpha3 + alpha4;
			//
			flux_rho  -= alpha;
			flux_rhou -= a2 * b1 * delta_rho + a2 * b2 * delta_rhou +
			a2 * b3 * delta_rhov + a2 * b4 * delta_rhoE;
			flux_rhov -= a3 * b1 * delta_rho + a3 * b2 * delta_rhou +
			a3 * b3 * delta_rhov + a3 * b4 * delta_rhoE;
			flux_rhoE -= a4 * b1 * delta_rho + a4 * b2 * delta_rhou +
			a4 * b3 * delta_rhov + a4 * b4 * delta_rhoE;
		}
	}
	//
	if(eig_val2 > REAL(0.0)) {
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f +
											rhov_f * rhov_f) / rho_f);
		ep_f = rhoE_f + p_f;
		rhouv_f = rhou_f * v_f;
		flux_rho  = rhou_f * nx + rhov_f * ny;
		flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny;
		flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny;
		flux_rhoE = ep_f * (u_f * nx + v_f * ny);
		//
		// A Entropic modification
		//
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t +
											rhov_t * rhov_t) / rho_t);
		lambda_f = u_f * nx + v_f * ny - norme * sqrt(gam * p_f / rho_f);
		lambda_t   = u_t * nx + v_t * ny - norme * sqrt(gam * p_t / rho_t);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
		}
		//
		if (eig_val1 < REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the first
			//.. columne of T-1 with the first row of T with eig_val1
			delta_rho  = rho_t - rho_f;
			delta_rhou = rhou_t - rhou_f;
			delta_rhov = rhov_t - rhov_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal / norme;
			hnx = nx / norme;
			hny = ny / norme;
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. First colomne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc - hnx;
			a3 = v_ave * usc - hny;
			a4 = REAL(0.5) * u2pv2 * usc + REAL(2.5) * c_speed - scal;
			//.. First row of the matrix T * eig_val1
			b1 = REAL(0.5) * eig_val1 * (REAL(0.5) * tempo11 * u2pv2 + scal);
			b2 = -REAL(0.5) * eig_val1 * (hnx + tempo11 * u_ave);
			b3 = -REAL(0.5) * eig_val1 * (hny + tempo11 * v_ave);
			b4 = REAL(0.5) * eig_val1 * tempo11;
			//
			alpha1 = a1 * b1 * delta_rho;
			alpha2 = a1 * b2 * delta_rhou;
			alpha3 = a1 * b3 * delta_rhov;
			alpha4 = a1 * b4 * delta_rhoE;
			alpha = alpha1 + alpha2 + alpha3 + alpha4;
			//
			flux_rho  += alpha;
			flux_rhou += a2 * b1 * delta_rho + a2 * b2 * delta_rhou +
			a2 * b3 * delta_rhov + a2 * b4 * delta_rhoE;
			flux_rhov += a3 * b1 * delta_rho + a3 * b2 * delta_rhou +
			a3 * b3 * delta_rhov + a3 * b4 * delta_rhoE;
			flux_rhoE += a4 * b1 * delta_rho + a4 * b2 * delta_rhou +
			a4 * b3 * delta_rhov + a4 * b4 * delta_rhoE;
		}
	}
}


template <class T>
void TPZEulerConsLaw::cSpeed(TPZVec<T> & sol, REAL gamma, T & c)
{
	if(sol[0] < REAL(1e-10))
	{
		PZError << "TPZEulerConsLaw::cSpeed Too low or negative density\n";
		TPZOutofRange obj;
		throw(obj);
	}
	
	int dim = sol.NElements() - 2;
	T press, temp;
	TPZEulerConsLaw::Pressure(gamma, dim, press, sol);
	temp = gamma * press;
	
	if(temp < REAL(1e-10)) // too low or negative
	{
		PZError << "TPZEulerConsLaw::cSpeed Too low or negative numerator\n";
	}
	c = sqrt(gamma * press/ sol[0]);
}

template <class T>
inline void TPZEulerConsLaw::uRes(TPZVec<T> & sol, T & us)
{
	if(sol[0] < REAL(1e-10))
	{
		PZError << "TPZEulerConsLaw::cSpeed Too low or negative density\n";
		TPZOutofRange obj;
		throw(obj);
		//      exit(-1);
	}
	
	T temp;
	switch(sol.NElements())
	{
		case(3):
			us = sol[1]/sol[0];
		case(4):
			temp = sol[1]*sol[1] + sol[2]*sol[2];
			if(temp < REAL(1e-40))
			{
				PZError << "TPZEulerConsLaw::uRes Zero Velocity\n";
				TPZOutofRange obj;
				throw(obj);
				//	 exit(-1);
			}
			us = sqrt(temp)/sol[0];
			break;
		case(5):
			temp = sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3];
			if(temp < REAL(1e-40))
			{
				PZError << "TPZEulerConsLaw::uRes Zero Velocity\n";
				TPZOutofRange obj;
				throw(obj);
				//	 exit(-1);
			}
			us = sqrt(temp)/sol[0];
			break;
		default:
			PZError << "TPZArtDiff::uRes Error: invalid Dimension\n";
	}
}

#endif
