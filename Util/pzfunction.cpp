/** 
 * @file 
 * @brief Contains the implementation of the methods to TPZFunction class. 
 */

#include "pzfunction.h"

TPZFunction::TPZFunction()
{
}

TPZFunction::~TPZFunction()
{
}

int TPZFunction::ClassId() const{
	return TPZFUNCTIONID;
}

void TPZFunction::Write(TPZStream &buf, int withclassid){
	TPZSaveable::Write(buf, withclassid);
}

void TPZFunction::Read(TPZStream &buf, void *context){
	TPZSaveable::Read(buf, context);
}

