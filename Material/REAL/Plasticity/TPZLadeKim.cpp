
#include "TPZLadeKim.h"

int TPZLadeKim::ClassId() const{
    return Hash("TPZLadeKim") ^ LADEKIMPARENT::ClassId() << 1;
}