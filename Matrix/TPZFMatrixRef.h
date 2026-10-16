#include "pzfmatrix.h"


/** @brief TPZFMatrix wrapper for acessing storage of a given matrix*/
// it assumes that the storage of a matrix is a vector of data
// it takes the form of a full matrix with a single column and size of the data
template<class TVar>
class TPZFMatrixRef : public TPZFMatrix<TVar>{
public:
  TPZFMatrixRef(const int64_t rows, TVar * &buf) :
    TPZFMatrix<TVar>(rows,1,buf,rows), fOrigRef(buf) {}
  TPZFMatrixRef(const TPZFMatrixRef<TVar> &ref) = default;
  TPZFMatrixRef(TPZFMatrixRef<TVar> &&ref) = default;
  TPZFMatrixRef &operator=(const TPZFMatrixRef<TVar> &ref) = delete;
  TPZFMatrixRef &operator=(TPZFMatrixRef<TVar> &&rval) = delete;
  TPZFMatrixRef &operator=(const TPZFMatrix<TVar> &ref) = delete;
  TPZFMatrixRef &operator=(TPZFMatrix<TVar> &&rval)
  {
    TPZFMatrix<TVar>::operator=(rval);
    return *this;
  }
  ~TPZFMatrixRef()
  {
    fOrigRef = this->Elem();
    this->Elem() = nullptr;
  }

private:
  TVar * &fOrigRef;
};

