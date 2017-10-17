
#include "TPZYCSandlerDimaggioL.h"

int TPZYCSandlerDimaggioL::ClassId() const{
    return Hash("TPZYCSandlerDimaggioL") ^ TPZYCSandlerDimaggio::ClassId() << 1;
}