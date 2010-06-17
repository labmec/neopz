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
	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZGaussRule{
public:

  enum {NUMINT_RULES = 38};
    friend class TPZIntRuleList;
    short	   fNumInt;		// number of integration points for this object
    REAL	*fLocation;	// location of the integration point
    REAL	*fWeight;	// weight of the integration point

    TPZGaussRule(int i);
    ~TPZGaussRule();

    public:

      int NInt() const{ return fNumInt;}	//return number of integration points

      REAL Loc(int i) const;						//return location of the ith pot

      REAL W(int i) const;						//return weight for the ith point
  };

#endif
