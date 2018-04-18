// $Id: TPZThermoForceA.cpp,v 1.2 2008-03-08 03:12:52 erick Exp $

#include "TPZThermoForceA.h"
#include "Hash/TPZHash.h"

int TPZThermoForceA::ClassId() const{
    return Hash("TPZThermoForceA");
}