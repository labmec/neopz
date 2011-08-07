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
//typedef TPZGaussLobattoRule INTRULE_PARENT;

/**
 * @ingroup integral
 * @brief Integration rule basic class. \ref integral "Numerical Integration"
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZIntRule : public INTRULE_PARENT{
public:
	
    friend class TPZIntRuleList;
	
    TPZIntRule(int i);
    ~TPZIntRule();
	
public:
	
	///return number of integration points
	int NInt() const;
	/// return location of the ith point
	REAL Loc(int i) const;
	///return weight for the ith point
	REAL W(int i) const;
};

#endif
