/**
 * @file
 * @brief Contains the TPZGaussRule class which implements the Gaussian quadrature.
 */

#ifndef TPZGAUSSRULE_H
#define TPZGAUSSRULE_H

#include "pzvec.h"

#define PZINTEGRAL_MAXITERATIONS_ALLOWED 100000

/**
 * @ingroup integral
 * @brief Implements the Gaussian quadrature. \ref integral "Numerical Integration"
 * Abstract class.
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGaussRule {
	/**
	 * @brief fType = 1 (Gauss Lobatto quadrature), 
	 * fType = 2 (Gauss Jacobi quadrature) for any alpha and beta the parameters of Jacobi polynomials,
	 * fType = 3 (Gauss Chebyshev) is a special Gauss Jacobi quadrature for alpha = -.5 and beta = -0.5
	 * fType = 0(default) (Gauss Legendre quadrature) is a special Gauss Jacobi quadrature for alpha = beta = 0.0
	 */
	int fType;

	/** @brief A instance of the cubature rule for pyramid can access computing methods to construct its points and weights. */
	friend class TPZIntRuleP3D;

	/** @brief The list can to access the constructor of the current class. */
    friend class TPZIntRuleList;
	
	/** @brief Number of integration points for this object */
    short	   fNumInt;
	/** @brief Location of the integration point */
    TPZVec<long double>	fLocation;
	/** @brief Weight of the integration point */
    TPZVec<long double>	fWeight;

	/** @brief ALPHA is the exponent of (1-X) in the quadrature rule: weight = (1-X)^ALPHA * (1+X)^BETA */
	long double fAlpha;
	/** @brief BETA is the exponent of (1+X) in the quadrature rule: weight = (1-X)^ALPHA * (1+X)^BETA */
	long double fBeta;

	/**
	 * @brief Constructor of quadrature rule 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 * @param type Identifies the type of integration points constructed
	 * @param alpha Exponent of the factor (1 - ksi) to Jacobi type
	 * @param beta Exponent of the factor (1 + ksi) to Jacobi type
	 * @note If alpha = beta = 0.0L the Jacobi type is identical to Legendre type
	 */
    TPZGaussRule(int order,int type = 0,long double alpha = 0.0L,long double beta = 0.0L);
	/** @brief Default destructor */
    ~TPZGaussRule();

public:
    enum {NRULESLOBATTO_ORDER = 645, NRULESLEGENDRE_ORDER = 1000};
	
	/** @brief Returns number of integration points */
	int NInt() const{ return fNumInt;}
	
	/** @brief Returns the gaussian quadrature type (Legendre, Lobatto, Jacobi) */
	int Type() { return fType; }
	
	/** @brief Returns location of the ith point */
	long double Loc(int i) const;
	/** @brief Returns weight for the ith point */
	long double W(int i) const;

	/** @brief Returns ALPHA: it is the exponent of (1-X) in the quadrature rule: \f$ weight = (1-X)^\alpha * (1+X)^\beta \f$ */
	long double Alpha() { return fAlpha; }
	/** @brief Returns BETA: it is the exponent of (1-X) in the quadrature rule:  \f$ weight = (1-X)^ALPHA * (1+X)^BETA \f$ */
	long double Beta() { return fBeta; }
	
	/**
	 * @brief Sets alpha and beta, the parameters of the Jacobi polynomials 
	 * @param alpha Exponent of the factor (1 - ksi) to Jacobi type
	 * @param beta Exponent of the factor (1 + ksi) to Jacobi type
	 */
	void SetParametersJacobi(long double alpha, long double beta) {
		fAlpha = alpha; fBeta = beta;
	}

	/**
	 * @brief Sets a gaussian quadrature: 0 - Gauss Legendre, 1 - Gauss Lobatto, 
	 * 2 - Gauss Jacobi with alpha = beta = 1., 3 - Gauss Chebyshev for alpha = beta = -0.5
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 * @param type Identifies the type of integration points constructed
	 */
	void SetType(int &type,int order);

	/** @brief Prints the number of integration points, all points and weights (as one dimension) */
	void Print(std::ostream & out = std::cout);

	/**
	 * @brief Checks sum of the weights is equal than measure of the master element, 
	 * and all of integration points belong to the master element.
	 * @return Returns false if one integration point is outside of the master element or the sum of weights is not one.
	 */
	bool CheckCubatureRule();
	
protected:
	/**
	 * @brief Computes the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0] 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 */
	void ComputingGaussLegendreQuadrature(int order);
	/** 
	 * @brief Computes the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0] 
	 * It is util to be called by another rule that need the Gaussian quadrature based on number of points not on order.
	 * @param npoints Number of integration points
	 * @param Location Vector of integration points ksi
	 * @param Weight Vector of corresponding weights
	 */
	void ComputingGaussLegendreQuadrature(int npoints,TPZVec<long double> &Location,TPZVec<long double> &Weight);
	/** 
	 * @brief Computes the points and weights for Gauss Lobatto Quadrature over the parametric 1D element \f$ [-1.0, 1.0] \f$
	 * It is to integrate functions \f$ f(x) \f$, but the first and last integration points are \f$ -1.0 \f$ and \f$1.0\f$ respectively
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 * @note The implementation follow Karniadakis and Sherwin: "Spectral/hp element methods for computational fluid dynamics
	 * (Oxford University Press, 2005).
	 */
	void ComputingGaussLobattoQuadrature(int order);
	/** 
	 * @brief Computes the points and weights for Gauss Lobatto Quadrature over the parametric 1D element \f$ [-1.0, 1.0] \f$
	 * It is to integrate functions \f$ f(x) \f$, but one of integration point is \f$-1.0\f$ or \f$1.0\f$ 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 */
	void ComputingGaussRadauQuadrature(int order) {
	}

	/** 
	 * @brief Compute numerical integration points and weight for Gauss Jacobi quadrature over \f$ [-1,1] \f$
	 * It is to integrate functions as \f$ (1-x)^\alpha (1+x)^\beta f(x) \f$ 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 * @param alpha Exponent of the factor (1 - ksi) to Jacobi type
	 * @param beta Exponent of the factor (1 + ksi) to Jacobi type
	 */
	void ComputingGaussJacobiQuadrature(int order,long double alpha,long double beta);
	/** 
	 * @brief Compute numerical integration points and weight for Gauss Jacobi quadrature over \f$ [-1,1] \f$
	 * It is to integrate functions as \f$ (1-x)^\alpha (1+x)^\beta f(x) \f$. It is util to be called by another rule.
	 * Do not stores the integration points and weights in the current instance.
	 * @param npoints Number of integration points
	 * @param alpha Exponent of factor (1+x)
	 * @param beta Exponent of factor (1-x)
	 * @param Location Vector of integration points ksi
	 * @param Weight Vector of corresponding weights
	 */
	void ComputingGaussJacobiQuadrature(int npoints,long double alpha, long double beta,TPZVec<long double> &Location,TPZVec<long double> &Weight);

	/**
	 * @brief Evaluates the Jacobi polinomial for real x.
	 * @param x Real number to compute J(x)
	 * @param order Order of the polinomial
	 * @param alpha Exponent of factor (1+x)
	 * @param beta Exponent of factor (1-x)
	 * @param b Recursion coefficients (Jacobi matrix)
	 * @param c Recursion coefficients 
	 * @param dp2 Returns the value of J'(x)
	 * @param p1 Returns the value of J[order-1](x), where J[order-1] represent the Jacobi polinomial to degree order-1.
	 */
	long double JacobiPolinomial(long double x,int order,long double alpha,long double beta,
								 TPZVec<long double> &b,TPZVec<long double> &c,long double *dp2,long double *p1);
	
	/**
	 * @brief Evaluates the polinomial Jacobi on x point.
	 * @param x point on the polinomial Jacobi is evaluated, J(x)
	 * @param n is the order of the polinomial
	 * @param alpha is the exponent of the factor (1-x)
	 * @param beta is the exponent of the factor (1+x)
	 * @note The Jacobi polinomial is evaluated using a recursion formulae
	 */
	long double JacobiPolinomial(long double x,int alpha,int beta,unsigned int n);

	/** 
	 * @brief Computes the points and weights for Gauss Chebyshev Quadrature over the parametric 1D element [-1.0, 1.0] 
	 * It is to integrate functions as \f$ f(x)/\sqrt(1-x^2)) \f$ 
	 * @param order Order of the polinomial will be integrated exactly with this cubature rule
	 */
	void ComputingGaussChebyshevQuadrature(int order);
};

/** @addtogroup integral */
/** @{ */

/** @name Auxiliar methods to evaluate the Jacobi polinomial
 * @{ */

/** Returns the machine precision of the computer's arithmetic as long double number */
long double machinePrecision();

/**
 * @brief Evaluate the factorial of a integer
 * @param n The integer to compute the factorial
 * @note The argument n cann't be negative.
 */
long double gamma(unsigned int n);

/** 
 * @brief Evaluates the Gamma function for a real number x, Gamma(x). The computation is based on an algorithm 
 * outlined in "An Overview of Software Development for Special Functions, in Numerical Analysis" (Dundee, 1975) by William Cody.
 * @param x real number to evaluate Gamma(x)
 * @note The implementation uses rational functions that approximate the gamma function as long double precision (19 significant digits).
 */
long double gamma(long double x);

/** @} */

/** @} */

#endif
