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
#ifndef TPZINTRULE_H
#define TPZINTRULE_H

#include "pzreal.h"

/**
	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntRule{
public:

  enum {NUMINT_RULES = 38};
    friend class TPZIntRuleList;
    short	   fNumInt;		// number of integration points for this object
    REAL	*fLocation;	// location of the integration point
    REAL	*fWeight;	// weight of the integration point

    TPZIntRule(int i);
    ~TPZIntRule();

    public:

      short NInt(){ return fNumInt;}	//return number of integration points

      REAL Loc(int i);						//return location of the ith pot

      REAL W(int i);						//return weight for the ith point
  };

#endif
