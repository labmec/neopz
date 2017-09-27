
#include "TPZYCSandlerDimaggioL.h"

int TPZYCSandlerDimaggioL::ClassId(){
    return TPZYCSandlerDimaggio::ClassId() ^ Hash("TPZYCSandlerDimaggioL");
}