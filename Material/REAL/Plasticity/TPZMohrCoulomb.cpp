
#include "TPZMohrCoulomb.h"

int TPZMohrCoulomb::ClassId(){
    return MOHRCOULOMBPARENT::ClassId() ^ Hash("TPZMohrCoulomb");
}