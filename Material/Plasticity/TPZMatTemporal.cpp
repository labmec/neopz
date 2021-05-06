
#include "TPZMatTemporal.h"

int TPZMatTemporal::ClassId() const{
    return Hash("TPZMatTemporal") ^ TPZMaterialData::ClassId() << 1;
}