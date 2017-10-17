
#include "TPZMohrCoulomb.h"

int TPZMohrCoulomb::ClassId() const{
    return Hash("TPZMohrCoulomb") ^ MOHRCOULOMBPARENT::ClassId() << 1;
}