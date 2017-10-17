
#include "TPZYCSandlerDimaggioL2.h"
#include "TPZYCSandlerDimaggioL.h"

int TPZYCSandlerDimaggioL2::ClassId() const{
    return Hash("TPZYCSandlerDimaggioL2") ^ TPZYCSandlerDimaggioL::ClassId() << 1;
}