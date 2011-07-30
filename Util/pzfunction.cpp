/** 
 * @file 
 * @brief Contains the implementation of the methods to TPZFunction class. 
 */
//$Id: pzfunction.cpp,v 1.1 2007-09-04 12:35:22 tiago Exp $

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

