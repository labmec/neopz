#ifndef TPZMATRIXWINDOW_H
#define TPZMATRIXWINDOW_H

#include <pzmatrix.h>

template<class TVar>
class TPZFMatrix;
/**
   @brief Creates a window to an existing TPZFMatrix<T>, allowing to perform MultAdd on a contiguous block
   of an existing matrix.
 */
template<class TVar>
class TPZMatrixWindow : public TPZMatrix<TVar>{
public:
  /**
     @brief Constructs a window based on a ROW MAJOR memory area
     @param mem_area First position of block to be windowed
     @param nrows Number of rows of window
     @param ncols Number of cols of window
     @param leading_dim Distance in memory between two sucessive rows
   */
  TPZMatrixWindow(TVar* mem_area, const int nrows, const int ncols, const int leading_dim);
  /**
     @brief Constructs a window based on a memory area
     @param mat original TPZFmatrix
     @param i row position of window beginning
     @param j col position of window beginning
     @param nrows Number of rows of window
     @param ncols Number of cols of window
   */
  TPZMatrixWindow(TPZFMatrix<TVar> &mat, const int i, const int j, const int nrows, const int ncols);

  /**
   * @brief It computes z = beta * y + alpha * opt_a(this)*opt_x(x) but z and x can not overlap in memory.
   * @param x Is x on the above operation
   * @param y Is y on the above operation
   * @param z Is z on the above operation
   * @param alpha Is alpha on the above operation
   * @param beta Is beta on the above operation
   * @param opt_a Indicates if A(this) is Transpose or not
   * @param opt_x Indicates if x is Transpose or not
   */
  void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
               const TVar alpha=1.,const TVar beta = 0.,const int opt_a = 0, const int opt_x = 0) const;
private:
  int64_t fLeadingDim{-1};
  TVar *fStorage{nullptr};
};

#endif