// $Id: tpzgausslobattorule.h,v 1.1 2009-06-16 02:48:55 erick Exp $
#ifndef TPZGAUSSLOBATTORULE_H
#define TPZGAUSSLOBATTORULE_H

#include "pzreal.h"

/**
	@author Erick Slis
*/
class TPZGaussLobattoRule{
public:

  enum {NUMINT_RULES = 11};
    friend class TPZIntRuleList;
    short	   fNumInt;		// number of integration points for this object
    REAL	*fLocation;	// location of the integration point
    REAL	*fWeight;	// weight of the integration point

    TPZGaussLobattoRule(int i);
    ~TPZGaussLobattoRule();

    public:

      short NInt(){ return fNumInt;}	//return number of integration points

      REAL Loc(int i);						//return location of the ith pot

      REAL W(int i);						//return weight for the ith point
  };

#endif
