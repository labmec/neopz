//$Id: tpzintrule.cpp,v 1.3 2010-06-17 13:07:27 phil Exp $
#include "tpzintrule.h"
#include "tpzintrulelist.h"
#include "pzerror.h"


TPZIntRule::TPZIntRule(int precision) : INTRULE_PARENT(precision){ ; }

TPZIntRule::~TPZIntRule(){

}

REAL TPZIntRule::Loc(int i) const {
	return INTRULE_PARENT::Loc(i);
}

REAL TPZIntRule::W(int i) const{
	return INTRULE_PARENT::W(i);
}

int TPZIntRule::NInt() const{
	return INTRULE_PARENT::NInt();
}
