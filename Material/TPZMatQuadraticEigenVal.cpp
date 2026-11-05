#include "TPZMatQuadraticEigenVal.h"
#include "Hash/TPZHash.h"

int TPZMatQuadraticEigenVal::ClassId() const{
    return Hash("TPZMatQuadraticEigenVal");
}