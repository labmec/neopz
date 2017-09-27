
#include "tpzyctrescaregularized.h"

int TPZYCTrescaRegularized::ClassId(){
    return TPZYCTresca::ClassId() ^ Hash("TPZYCTrescaRegularized");
}