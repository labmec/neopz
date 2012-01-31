/**
 * @file
 * @brief Contains the TPZIntRule class which is the integration rule basic class.
 */
//$Id: tpzintrule.h,v 1.3 2010-06-17 13:07:27 phil Exp $
#ifndef TPZINTRULE_H
#define TPZINTRULE_H

#include "pzreal.h"
#include "tpzgaussrule.h"
#include "tpzgausslobattorule.h"

/** @ingroup integral */
typedef TPZGaussRule INTRULE_PARENT;
//typedef TPZGaussLobattoRule INTRULE_PARENT;   // Gauss Lobatto implemented by Erick Slis

/**
 * @ingroup integral
 * @brief Integration rule basic class. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRule : public INTRULE_PARENT {
public:
	
    friend class TPZIntRuleList;
	/**
	 * @brief Numerical Integration Rule Constructor 
	 * @param precision Polinomial degree integrated exactly
	 */
    TPZIntRule(int precision);
	
	/// Default destructor
    ~TPZIntRule();
	
public:
	
	///return number of integration points
	int NInt() const;
	/**
	 * @brief Returns location of the ith point
	 * @param i Integration point for the rule of integration corresponding (i < NInt()). Zero based.
	 */
	REAL Loc(int i) const;
	/**
	 * @brief Returns weight of the ith point
	 * @param i Integration point for the rule of integration corresponding (i < NInt()). Zero based.
	 */
	REAL W(int i) const;
};

#endif
