/**
 * \file
 * @brief Contains the TPZDiffusionConsLaw class which implements a Euler equation where is introduced a diffusive term to stabilize.
 */

#ifndef DIFFUSIONFORCONSLAWHH
#define DIFFUSIONFORCONSLAWHH


#include "pzvec.h"
#include "pzfmatrix.h"

/**
 * @ingroup analysis
 * @brief This class adds to the term of diffusion to the variacional formularization \n
 * of the differential equation partial compressible of Euler (hyperbolic). \ref analysis Analysis
 */
/**
 * This term is introduced to stabilize the numerical method of approach.
 */
class TPZDiffusionConsLaw {
	
	
public:
	/** @brief Ratio between specific heat is constant and the specific heat the constant volume of a polytropic gas */
	static REAL fGamma;
	
	/** @brief Coefficient of size ratio of the element in the diffusion term */
	static REAL fDelta;
	
	/** @brief Parameter that it limits the condition of stability of the numerical approach */
	static REAL fCFL;
	
	/** @brief Term that adds stability to the numerical method of approach: SUPG, LS, Bornhaus */
	static std::string fArtificialDiffusion;
	
private:
	/** @brief Matrix computation derivatives of fluxes: dFx/dx, dFy/dy, dFz/dz */
	TPZFMatrix<REAL> fA,fB,fC;
	
	/** @brief Problem dimension */
	int fDimension;
	
public:
	/** @brief Default constructor */	
	TPZDiffusionConsLaw();
	/** @brief Constructor for model dimension dim */
	TPZDiffusionConsLaw(TPZVec<REAL> U,REAL gamma,int dim,const std::string &diff);
	/** @brief Destructor */
	~TPZDiffusionConsLaw();
	
	//static void SetDelta(REAL delta){TPZDiffusionConsLaw::fDelta = delta;}
	
	void GradientOfTheFlow(TPZFMatrix<REAL> &DF1,TPZFMatrix<REAL> &DF2,TPZFMatrix<REAL> &DF3);
	
	REAL DeltaOtimo();
	
	REAL Delta();
	
	REAL CFL(int degree);
	
	/** @brief Jacobiano of the tensor flux of Euler */
	static void JacobFlux(TPZVec<REAL> U,TPZFMatrix<REAL> &A,TPZFMatrix<REAL> &B,TPZFMatrix<REAL> &C);
	
	void Divergence(TPZVec<REAL> &dphi,TPZFMatrix<REAL> &diverg);
	
	/** @brief Operation product point in the diffusion term */
	void PointOperator(TPZVec<REAL> &dphi,TPZFMatrix<REAL> &diff_term);
	
	void Tau(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz);
	
	void SUPG(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz);
	
	void LS(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz);
	
	void Bornhaus(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz);
	
	/** @brief Flux of Roe (MOUSE program) */
	static void Roe_Flux(REAL rho_f, REAL rhou_f, REAL rhov_f, REAL rhow_f, REAL rhoE_f, 
						 REAL rho_t, REAL rhou_t, REAL rhov_t,REAL rhow_t, REAL rhoE_t, 
						 REAL nx, REAL ny, REAL nz, REAL gam, REAL & flux_rho, REAL & flux_rhou, 
						 REAL & flux_rhov,REAL & flux_rhow, REAL & flux_rhoE);
	
	static void Roe_Flux(REAL rho_f, REAL rhou_f, REAL rhov_f, REAL rhoE_f,
						 REAL rho_t, REAL rhou_t, REAL rhov_t, REAL rhoE_t,
						 REAL nx, REAL ny, REAL gam,
						 REAL &flux_rho, REAL &flux_rhou,
						 REAL &flux_rhov, REAL &flux_rhoE);
};

#endif
