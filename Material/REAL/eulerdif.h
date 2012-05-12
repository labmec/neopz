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
	static  STATE  fGamma;
	
public:
	/**
	 * @brief Calculates the pressure from equation of state.\n
	 * It is expected: \f$ U = (density, density * velocityX, density * velocityY, energy) \f$
	 */
	static  STATE Pressure(TPZVec<STATE> &U);
	/** @brief Calculates the fluxes \f$ F_x \f$ and \f$ F_y \f$ */
	static  void Flux(TPZVec<STATE> &u,TPZVec<STATE> &flux);
	static  void JacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &Ajacob,TPZFMatrix<STATE> &Bjacob);
	static	void JacobFlux(TPZVec<STATE> &U,TPZFMatrix<STATE> &jacob,TPZVec<REAL> &normal);
	
	static  void ValJacobFlux(TPZVec<STATE> &u,TPZFMatrix<STATE> &valjacob,TPZVec<REAL> &normal);
	
	static  void MatrixDiff(TPZVec<STATE> &sol, TPZFMatrix<REAL> &axes, TPZFMatrix<REAL> &jacinv,TPZFMatrix<STATE>
							&ATauA,TPZFMatrix<STATE> &ATauB,TPZFMatrix<STATE> &BTauA,TPZFMatrix<STATE> &BTauB);
	
	static  void InvJacob2d(TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &jacinv);
	
	static  void InverseJacob(TPZFMatrix<STATE> &jac);
	/// Static main for test
	static int main();
	
};

#endif
