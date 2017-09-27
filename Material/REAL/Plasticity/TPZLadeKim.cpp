
#include "TPZLadeKim.h"

int TPZLadeKim::ClassId() {
    return LADEKIMPARENT::ClassId() ^ Hash("TPZLadeKim");
}