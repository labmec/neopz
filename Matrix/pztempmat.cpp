
#include "pztempmat.h"

/*
TPZTempFMatrix operator+(TPZFMatrix &first, TPZTempFMatrix other) {
        other.Object() += first;
   return other;
}
*/

TPZTempFMatrix operator+(TPZTempFMatrix first, TPZTempFMatrix other) {
        other.Object() += first.Object();
   return other;
}

//TPZTempFMatrix operator*(const REAL val, TPZTempFMatrix first) {
//	first.Object()*=val;
//   return first;
//}


TPZFMatrix &operator+=(TPZFMatrix &first, TPZTempFMatrix other) {
        return first += other.Object();
}


