// $Id: tpzgausslobattorule.h,v 1.1 2009-06-16 02:48:55 erick Exp $
#ifndef TPZGAUSSLOBATTORULE_H
#define TPZGAUSSLOBATTORULE_H

#include "pzreal.h"

/**
 * @ingroup integral
 *
 * @brief  Similar to Gaussian quadrature with the following differences: \n
 * The integration points include the end points of the integration interval.
 * 
 * It is accurate for polynomials up to degree 2n–3, where n is the number of integration points.
 *
 * @author Erick Slis
 */
class TPZGaussLobattoRule {
public:
	
	enum {NUMINT_RULES = 11};
    friend class TPZIntRuleList;
	/// number of integration points for this object
    short	   fNumInt;
	/// location of the integration point
    REAL	*fLocation;
	/// weight of the integration point
    REAL	*fWeight;
	
	/// Constructor indicating the order
    TPZGaussLobattoRule(int i);
	/// Destructor
    ~TPZGaussLobattoRule();
	
public:
	///return number of integration points
	short NInt(){ return fNumInt;}
	///return location of the ith point
	REAL Loc(int i);
	///return weight for the ith point
	REAL W(int i);
};

#endif
