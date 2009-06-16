//$Id: tpzintrule.h,v 1.2 2009-06-16 02:48:55 erick Exp $
#ifndef TPZINTRULE_H
#define TPZINTRULE_H

#include "pzreal.h"
#include "tpzgaussrule.h"
#include "tpzgausslobattorule.h"


typedef TPZGaussRule INTRULE_PARENT;
//typedef TPZGaussLobattoRule INTRULE_PARENT;

/**
	@author Philippe R. B. Devloo <phil@fec.unicamp.br>
*/
class TPZIntRule : public INTRULE_PARENT{
public:

    friend class TPZIntRuleList;

    TPZIntRule(int i);
    ~TPZIntRule();

    public:

      short NInt();							//return number of integration points

      REAL Loc(int i);						//return location of the ith pot

      REAL W(int i);						//return weight for the ith point
  };

#endif
