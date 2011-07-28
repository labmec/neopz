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
//
#ifndef TPZGAUSSRULE_H
#define TPZGAUSSRULE_H

#include "pzreal.h"

/**
 * @ingroup integral
 *
 * @brief Implements the Gaussian quadrature.
 * @author Philippe R. B. Devloo <phil@fec.unicamp.br>
 */
class TPZGaussRule{
public:
	
	enum {NUMINT_RULES = 38};
    friend class TPZIntRuleList;
	/// number of integration points for this object
    short	   fNumInt;
	/// location of the integration point
    REAL	*fLocation;
	/// weight of the integration point
    REAL	*fWeight;
	
    TPZGaussRule(int i);
    ~TPZGaussRule();
	
public:
	///return number of integration points
	int NInt() const{ return fNumInt;}
	///return location of the ith point
	REAL Loc(int i) const;
	///return weight for the ith point
	REAL W(int i) const;
};


#endif
