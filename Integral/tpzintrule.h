//$Id: tpzintrule.h,v 1.3 2010-06-17 13:07:27 phil Exp $
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

      int NInt() const;							//return number of integration points

      REAL Loc(int i) const;						//return location of the ith pot

      REAL W(int i) const;						//return weight for the ith point
  };

#endif
