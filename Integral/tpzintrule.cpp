//$Id: tpzintrule.cpp,v 1.2 2009-06-16 02:48:55 erick Exp $
#include "tpzintrule.h"
#include "tpzintrulelist.h"
#include "pzerror.h"


TPZIntRule::TPZIntRule(int precision) : INTRULE_PARENT(precision){ ; }

TPZIntRule::~TPZIntRule(){

}

REAL TPZIntRule::Loc(int i) {
	return INTRULE_PARENT::Loc(i);
}

REAL TPZIntRule::W(int i) {
	return INTRULE_PARENT::W(i);
}

short TPZIntRule::NInt(){
	return INTRULE_PARENT::NInt();
}
