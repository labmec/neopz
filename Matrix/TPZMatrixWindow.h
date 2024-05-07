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
  //!Copy constructor (what is the expected behaviour?)
  TPZMatrixWindow(const TPZMatrixWindow<TVar> &mat){
    DebugStop();
  }
  //!Move constructor (what is the expected behaviour?)
  TPZMatrixWindow(TPZMatrixWindow<TVar> &&mat){
    DebugStop();
  }
  //! Copy assignment operator. Both matrices must have the same size!
  TPZMatrixWindow<TVar> &operator=(const TPZMatrixWindow<TVar> &mat);
  //!Move assignment operator (what is the expected behaviour?)
  TPZMatrixWindow<TVar> &operator=(TPZMatrixWindow<TVar> &&mat){
    DebugStop();
    return *this;
  }
  //! Copy assignment from TPZFMatrix<T> operator. Both matrices must have the same size!
  TPZMatrixWindow<TVar> &operator=(const TPZFMatrix<TVar> &mat);
  //! Default destructor (no memory deallocation)
  ~TPZMatrixWindow() = default;
  /**
     @brief Constructs a window based on a ROW MAJOR memory area
     @param mem_area First position of block to be windowed
     @param nrows Number of rows of window
     @param ncols Number of cols of window
     @param leading_dim Distance in memory between two sucessive rows
     @param size size of memory area (must be able to accomodate nrows*ncols)
   */
  TPZMatrixWindow(TVar* mem_area, const int nrows, const int ncols, const int leading_dim, const int size_mem);
  /**
     @brief Constructs a window based on a memory area
     @param mat original TPZFmatrix
     @param i row position of window beginning
     @param j col position of window beginning
     @param nrows Number of rows of window
     @param ncols Number of cols of window
   */
  TPZMatrixWindow(TPZFMatrix<TVar> &mat, const int i, const int j, const int nrows, const int ncols);

  //! Multiply itself by a given scalar
  TPZMatrixWindow<TVar> &operator*=(const TVar val) override;
  
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
  void MultAdd(const TPZMatrixWindow<TVar> &x,const TPZMatrixWindow<TVar> &y, TPZMatrixWindow<TVar> &z,
               const TVar alpha=1.,const TVar beta = 0.,const int opt_a = 0, const int opt_x = 0) const;
  /**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 */
  inline void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
                      const TVar alpha=1., const TVar beta = 0., const int opt = 0) const override{
    //is it safe to do so?
    const TPZMatrixWindow<TVar> x_window(const_cast<TVar*>(x.Elem()),x.Rows(),x.Cols(),x.Rows(), x.Rows()*x.Cols());
    const TPZMatrixWindow<TVar> y_window(const_cast<TVar*>(y.Elem()),y.Rows(),y.Cols(),y.Rows(), y.Rows()*y.Cols());
    if(beta!=(TVar)0){
      z=y;
    }else{
      if(opt==0){
        z.Redim(this->Rows(),x.Cols());
      }else{
        z.Redim(this->Cols(),x.Cols());
      }
    }
    TPZMatrixWindow<TVar> z_window(z,0,0,z.Rows(),z.Cols());
    //forwards to general implementation (opt_x is always zero)
    MultAdd(x_window,y_window,z_window,alpha,beta,opt,0);
  }
  
  inline TVar &operator()(const int64_t row, const int64_t col) {
#ifndef PZNODEBUG
	if ( (row >= this->Rows()) || (col >= this->Cols()) || row <0 || col<0 ) {
		PZError<<__PRETTY_FUNCTION__
           <<"\nIndex out of range"<<std::endl;
    DebugStop();
	}
#endif
    return fStorage[fLeadingDim*col+row];
  }
  inline TVar operator()(const int64_t row, const int64_t col)const {
#ifndef PZNODEBUG
	if ( (row >= this->Rows()) || (col >= this->Cols()) || row <0 || col<0 ) {
		PZError<<__PRETTY_FUNCTION__
           <<"\nIndex out of range"<<std::endl;
    DebugStop();
	}
#endif
    return fStorage[fLeadingDim*col+row];
  }

private:
  TVar &g(const int64_t row, const int64_t col) const {
    return fStorage[fLeadingDim*col+row];
  }
public:

  const TVar GetVal(const int64_t row,const int64_t col ) const override{
    return fStorage[fLeadingDim*col+row];
  }

  int PutVal(const int64_t row, const int64_t col,const TVar & value ) override{
    fStorage[col*fLeadingDim+row ] = value;
    return( 1 );
  }

  //! We do not allow resizing a TPZMatrixWindow 
  int Redim(const int64_t newRows, const int64_t newCols ) override {
    if(newRows!=this->fRow || newCols!=this->fCol){
      DebugStop();
    }
    Zero();
    return 1;
  }
  //! We do not allow resizing a TPZMatrixWindow 
	int Resize(const int64_t newRows, const int64_t newCols ) override{
		if(newRows!=this->fRow || newCols!=this->fCol){
      DebugStop();
    }
		return 1;
	}
  //! Sets all elements to zero
  int Zero() override{
    for(int ic = 0; ic < this->fCol; ic++){
      for(int ir = 0; ir < this->fRow; ir++){
        this->PutVal(ir,ic,(TVar)0);
      }
    }
    return 1;
  }
  
  //! Not available for this matrix type 
  inline int64_t Size() const override{
    /*
      Theoretically these methods could be implemented,
      however, if both Size() and Elem() are implemented,
      then it could lead to the Storage() method, defined at
      TPZMatrix<T>, being called with no warnings whatsoever.
      Let us keep it safe for now.
    */
    DebugStop();
    return -1;
  }
  //! Not available for this matrix type 
  inline TVar* &Elem() override{
    /*
      Theoretically these methods could be implemented,
      however, if both Size() and Elem() are implemented,
      then it could lead to the Storage() method, defined at
      TPZMatrix<T>, being called with no warnings whatsoever.
      Let us keep it safe for now.
    */
    DebugStop();
    static TVar* myvar{nullptr};
    return myvar;
  }
  //! Not available for this matrix type 
  inline const TVar* Elem() const override{
    /*
      Theoretically these methods could be implemented,
      however, if both Size() and Elem() are implemented,
      then it could lead to the Storage() method, defined at
      TPZMatrix<T>, being called with no warnings whatsoever.
      Let us keep it safe for now.
    */
    DebugStop();
    return nullptr;
  }
  //! Not available for this matrix type 
  TPZMatrix<TVar> *Clone() const override{
    DebugStop();
    return nullptr;
  }
  //! Not available for this matrix type
  TPZMatrix<TVar> *NewMatrix() const override{
    DebugStop();
    return nullptr;
  }
  //! Not available for this matrix type
  int Decompose(const DecomposeType dt) override{
    DebugStop();
    return -1;
  }
  //! Not available for this matrix type
  void CopyFrom(const TPZMatrix<TVar> *mat) override{
    DebugStop();
  }
  //! Not available for this matrix type
  int SolveDirect(TPZFMatrix<TVar>& F , const DecomposeType dt) override{
    DebugStop();
    return -1;
  }
  //! Not available for this matrix type
  int SolveDirect(TPZFMatrix<TVar>& F , const DecomposeType dt) const override{
    DebugStop();
    return -1;
  }
private:
  //! Auxiliary function for checking input param
  void CheckConstructor(const int i, const int j, const int nr, const int nc, const int nr_orig, const int nc_orig);
  int64_t fLeadingDim{-1};
  TVar *fStorage{nullptr};
};

#endif