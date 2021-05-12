#include "TPZMatGeneralisedEigenVal.h"
#include "Hash/TPZHash.h"

int TPZMatGeneralisedEigenVal::ClassId() const{
    return Hash("TPZMatGeneralisedEigenVal");
}