#ifndef TPZTEMPFMATRIX
#define TPZTEMPFMATRIX

#include <iostream>
#include "pzfmatrix.h"
#include "pzreal.h"


class TPZSFMatrix;

class TPZTempFMatrix {

        struct MIntern {
        TPZFMatrix fObject;
		int fNRef;

		MIntern() : fObject() {
                fNRef = 1;
        }
		MIntern(const TPZFMatrix &init) : fObject(init) {
                fNRef = 1;
        }

   } *fp;

public:

inline   TPZTempFMatrix() {
        fp = new MIntern;
   }


inline   TPZTempFMatrix(const TPZFMatrix &basis) {
        fp = new MIntern(basis);
   }

inline   ~TPZTempFMatrix() {
        fp->fNRef--;
      if(!fp->fNRef) delete fp;
   }

inline   TPZTempFMatrix(const TPZTempFMatrix &other) {
        fp = other.fp;
      fp->fNRef++;
   }


inline TPZFMatrix &Object() { return fp->fObject;}

inline TPZTempFMatrix operator*(const REAL val) {
	Object() *= val;
   return *this;
}

inline TPZTempFMatrix operator-(const TPZFMatrix &other) {
	Object() -= other;
   return *this;
}

};

//TPZTempFMatrix operator+(TPZFMatrix &first, TPZTempFMatrix other);
TPZTempFMatrix operator+(TPZTempFMatrix first, TPZTempFMatrix other);
//TPZTempFMatrix operator+(const TPZMatrix &A, const TPZMatrix & B );
TPZTempFMatrix operator-(const TPZMatrix & A, const TPZMatrix & B );
TPZTempFMatrix operator*(const TPZMatrix & A, const TPZFMatrix & B );
inline TPZTempFMatrix operator*(const REAL val, TPZTempFMatrix other) {
	other.Object() *= val;
   return other;
}


/* inline TPZTempFMatrix operator*(TPZFMatrix &first, TPZTempFMatrix other) { */
/*         TPZTempFMatrix result; */
/*    first.Multiply(other.Object(),result.Object()); */
/*    return result; */
/* } */


//template<class TSub>
//inline TPZTempFMatrix operator*(TPZTempFMatrix &first, TSub &other) {
//        return first.Object()*other;
//}

TPZTempFMatrix operator*(const REAL value, const TPZFMatrix &A );





//template class TPZTempFMatrix<REAL>;
//template class TPZTempFMatrix<TPZSFMatrix>;

#endif


