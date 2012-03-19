/**
 * \file
 * @brief Contains the TEulerDiffusivity class which implements a numerical diffusivity coefficient for SUPG.
 */
#ifndef EULERDIFFUSIONHH
#define EULERDIFFUSIONHH

#include "pzvec.h" 
#include "pzfmatrix.h"

/** 
 * @brief Implements a numerical diffusivity coeficient for the SUPG method. \ref analysis Analysis
 * @ingroup analysis
 */
/**
 * The associated material is the bi-dimensional Euler equations for dynamic of the gases.
 * @note This class is to be used as a template argument for a different class.
 *
 * \par Variables:
 * \f$ U = (\rho, \rho * \upsilon_x, \rho * \upsilon_y, E) \f$ \n
 * \f$ \rho = \f$ density \n
 * \f$ (\upsilon_x, \upsilon_y) = \f$ velocity \n
 * \f$ \| \mbox{velocity} \| = \| v \| = \sqrt(\upsilon_x^2 + \upsilon_y^2 ) \f$ \n
 * \f$ E = \f$ energy \n
 * \f$ p = \f$ pressure \n
 * \f$ c_p = \f$ specific heat at constant pressure of the gas \n
 * \f$ c_V = \f$ specific heat at constant volume of the gas.
 *
 * \par Fluxes:
 * \f$ F_x(U) = (\rho * \upsilon_x, \rho * \upsilon_x^2 + p, \rho * \upsilon_x * \upsilon_y, \upsilon_x (\rho * E + p)) \f$ \n
 * \f$ F_y(U) = (\rho * \upsilon_y, \rho * \upsilon_x * \upsilon_y, \rho * \upsilon_y^2 + p, \upsilon_y (\rho * E + p)) \f$
 *
 * \par Equation of state:
 * \f$ p = (\gamma - 1) \star (E - \frac{1}{2} \rho \| v \| ) \f$ \n
 * \f$ \gamma = \frac{c_p}{c_V} \f$ \n
 *
 * For tests we used \f$ \gamma = 1.4 \f$ .
 */
class TEulerDiffusivity {
	/** @brief Polytropic gas constant */
	/** Ratio of the specific heat at constant pressure and the specific heat at constant volume. It is used at equation of state. */
	static  REAL  fGamma;
	
public:
	/**
	 * @brief Calculates the pressure from equation of state.\n
	 * It is expected: \f$ U = (density, density * velocityX, density * velocityY, energy) \f$
	 */
	static  REAL Pressure(TPZVec<REAL> &U);
	/** @brief Calculates the fluxes \f$ F_x \f$ and \f$ F_y \f$ */
	static  void Flux(TPZVec<REAL> &u,TPZVec<REAL> &flux);
	static  void JacobFlux(TPZVec<REAL> &u,TPZFMatrix<REAL> &Ajacob,TPZFMatrix<REAL> &Bjacob);
	static	void JacobFlux(TPZVec<REAL> &U,TPZFMatrix<REAL> &jacob,TPZVec<REAL> &normal);
	
	static  void ValJacobFlux(TPZVec<REAL> &u,TPZFMatrix<REAL> &valjacob,TPZVec<REAL> &normal);
	
	static  void MatrixDiff(TPZVec<REAL> &sol, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &jacinv,TPZFMatrix<REAL>
							&ATauA,TPZFMatrix<REAL> &ATauB,TPZFMatrix<REAL> &BTauA,TPZFMatrix<REAL> &BTauB);
	
	static  void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv);
	
	static  void InverseJacob(TPZFMatrix<REAL> &jac);
	/// Static main for test
	static int main();
	
};

#endif
