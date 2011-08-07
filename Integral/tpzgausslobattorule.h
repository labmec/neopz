/**
 * @file
 * @brief Contains the TPZGaussLobattoRule class which implements similar Gaussian quadrature with two differences.
 */
// $Id: tpzgausslobattorule.h,v 1.1 2009-06-16 02:48:55 erick Exp $
#ifndef TPZGAUSSLOBATTORULE_H
#define TPZGAUSSLOBATTORULE_H

#include "pzreal.h"

/**
 * @ingroup integral
 * @author Erick Slis
 * @brief Similar to Gaussian quadrature but with two differences. \ref integral "Numerical Integration" 
 */
/**
 * Similar to Gaussian quadrature with following differences: \n
 * The integration points include the end points of the integration interval. \n
 * It is accurate for polynomials up to degree 2n–3, where n is the number of integration points.
 */
class TPZGaussLobattoRule {
public:
	
	enum {NUMINT_RULES = 11};
    friend class TPZIntRuleList;
	/** @brief Number of integration points for this object */
    short	   fNumInt;
	/** @brief Location of the integration point */
    REAL	*fLocation;
	/** @brief Weight of the integration point */
    REAL	*fWeight;
	
	/** @brief Constructor indicating the order */
    TPZGaussLobattoRule(int i);
	/** @brief Destructor */
    ~TPZGaussLobattoRule();
	
public:
	/** @brief Return number of integration points */
	short NInt(){ return fNumInt;}
	/** @brief Return location of the ith point */
	REAL Loc(int i);
	/** @brief Return weight for the ith point */
	REAL W(int i);
};

#endif
