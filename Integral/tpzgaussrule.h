/**
 * @file
 * @brief Contains the TPZGaussRule class which implements the Gaussian quadrature.
 */
//
// C++ Interface: tpzintrule
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2008
//
// Copyright: See COPYING file that comes with this distribution
//

#ifndef TPZGAUSSRULE_H
#define TPZGAUSSRULE_H

#include "pzreal.h"
#include "pzvec.h"

/**
 * @ingroup integral
 * @brief Implements the Gaussian quadrature. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGaussRule {

	/**
	 * @brief Computes the points and weights for Gauss Legendre Quadrature over the parametric 1D element [-1.0, 1.0]
	 * @param[in] Npoints Number of points for the quadrature rule
	 * @param[out] points Pointer to vector of locations on \f$[-1.0 1.0]\f$
	 * @param[out] weights Pointer to vector of weights
	 */
	void ComputingGaussLegendreQuadrature(int Npoints,TPZVec<REAL> &points,TPZVec<REAL> &weights);

public:

	enum {NUMINT_RULES = 200};
    friend class TPZIntRuleList;
	/** @brief Number of integration points for this object */
    short	   fNumInt;
	/** @brief Location of the integration point */
    TPZVec<REAL>	fLocation;
	/** @brief Weight of the integration point */
    TPZVec<REAL>	fWeight;

    TPZGaussRule(int i);
    ~TPZGaussRule();

public:
	/** @brief Returns number of integration points */
	int NInt() const{ return fNumInt;}
	/** @brief Returns location of the ith point */
	REAL Loc(int i) const;
	/** @brief Returns weight for the ith point */
	REAL W(int i) const;

};

#endif
