
#include "tpzyctrescaregularized.h"

int TPZYCTrescaRegularized::ClassId() const{
    return Hash("TPZYCTrescaRegularized") ^ TPZYCTresca::ClassId() << 1;
}