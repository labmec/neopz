/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a symmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 * Some of the functionalities of this class depends on the MKL library and thus needs the NeoPZ library
 * to be configured using USING_MKL=ON during the CMake process. Search on this header for MKL to see which functionalities are these.
 */

#ifndef SYSMPMATH
#define SYSMPMATH

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZPardisoSolver.h"

 /**
  * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrix : public TPZMatrix<TVar>{
	
    friend class TPZPardisoSolver<TVar>;
    
public :
    /** @brief Constructor based on number of rows and columns */
  TPZSYsmpMatrix();
	/** @brief Constructor based on number of rows and columns */
  TPZSYsmpMatrix(const int64_t rows, const int64_t cols );
	/** @brief Copy constructor */
  TPZSYsmpMatrix(const TPZSYsmpMatrix<TVar> &cp) = default;
  /** @brief Move constructor*/
  TPZSYsmpMatrix(TPZSYsmpMatrix<TVar> &&cp) = default;
  /** @brief Copy-assignment operator*/
  TPZSYsmpMatrix &operator=(const TPZSYsmpMatrix<TVar> &copy) = default;
  /** @brief Move-assignment operator*/
  TPZSYsmpMatrix &operator=(TPZSYsmpMatrix<TVar> &&copy) = default;
  inline TPZSYsmpMatrix<TVar>*NewMatrix() const override {return new TPZSYsmpMatrix<TVar>{};}
  CLONEDEF(TPZSYsmpMatrix)
	/** @brief Destructor */
	virtual ~TPZSYsmpMatrix();

  /** @brief Creates a copy from another TPZSYsmpMatrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override
  {                                                           
    auto *from = dynamic_cast<const TPZSYsmpMatrix<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
  
  /** @brief Checks if the current matrix is symmetric */
  virtual int IsSymmetric() const  override { return 1; }
  /** @brief Checks if current matrix is square */
  inline int IsSquare() const { return 1;}
    
  /** @brief Zeroes the matrix */
  virtual int Zero() override;

  /** @brief Zeroes the matrix */
  virtual int Redim(int64_t rows, int64_t cols) override
  {
    if(rows == this->fRow && cols == this->fCol)
      {
        Zero();
      }
    else
      {
        DebugStop();
      }
    return 0;
  }

  void SetIsDecomposed(int val) override;
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric) override;
	  
	  /** @brief Get the matrix entry at (row,col) without bound checking */
	  virtual const TVar GetVal(const int64_t row, const int64_t col ) const override;
    
    /** @brief Put values without bounds checking \n
     *  This method is faster than "Put" if DEBUG is defined.
     */
    virtual int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val ) override;
  /** @name Arithmetic*/
  /** @{ */
  TPZSYsmpMatrix<TVar> operator+(const TPZSYsmpMatrix<TVar> & A) const;
  TPZSYsmpMatrix<TVar> operator-(const TPZSYsmpMatrix<TVar> & A) const;
  TPZSYsmpMatrix<TVar> operator*(const TVar alpha) const;
  TPZSYsmpMatrix<TVar> &operator+=(const TPZSYsmpMatrix<TVar> &A );
  TPZSYsmpMatrix<TVar> &operator-=(const TPZSYsmpMatrix<TVar> &A );
  TPZSYsmpMatrix<TVar> &operator*=(const TVar val) override;
  
	/** @brief Computes z = beta * y + alpha * opt(this)*x */
	/** @note z and x cannot overlap in memory */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const override;
	/** @} */
	/** @brief Sets data to the class */
	virtual void SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A );

    virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex) override;
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex) override;
    /// Access function for the coefficients
    TPZVec<TVar> &A()
    {
        return fA;
    }
    
    TPZVec<int64_t> &IA()
    {
        return fIA;
    }
    
    TPZVec<int64_t> &JA()
    {
        return fJA;
    }

	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const override;
    
    /**
     * @name Factorization
     * @brief Those member functions perform the matrix factorization
     * @{
     */
    

    /**
     * @brief Decomposes the current matrix using LDLt. Depends on MKL. \n
     * The current matrix has to be symmetric.
     * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
     */
    virtual int Decompose_LDLt(std::list<int64_t> &singular) override;
    /** @brief Decomposes the current matrix using LDLt. Depends on MKL.*/
    virtual int Decompose_LDLt() override;

    /** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. Depends on MKL.*/
    virtual int Decompose_Cholesky() override;
    /**
     * @brief Decomposes the current matrix using Cholesky method.Depends on MKL.
     * @param singular
     */
    virtual int Decompose_Cholesky(std::list<int64_t> &singular)  override;
    

    /** @} */
    
    /**
     * @name Substitutions
     * @brief Substitutions forward and backward
     * @{
     */
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.Depends on MKL.
     * @param b right hand side and result after all
     */
    virtual int Subst_LForward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.Depends on MKL.
     * @param b right hand side and result after all
     */
    virtual int Subst_LBackward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.Depends on MKL.
     * @param b right hand side and result after all
     */
    virtual int Subst_Diag( TPZFMatrix<TVar>* b ) const override;

    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular.Depends on MKL.
     * @param b right hand side and result after all
     */
    virtual int Subst_Forward( TPZFMatrix<TVar>* b ) const override;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular.Depends on MKL.
     * @param b right hand side and result after all
     */
    virtual int Subst_Backward( TPZFMatrix<TVar>* b ) const override;
    

    /** @} */
    

    int ClassId() const override;

    void ComputeDiagonal();
    //! Gets reference to TPZPardisoSolver instance for fine-tuning
    TPZPardisoSolver<TVar> & GetPardisoControl()
    {return fPardisoControl;}
protected:
  void CheckTypeCompatibility(const TPZMatrix<TVar>*aPtr,
                              const TPZMatrix<TVar>*bPtr)const override;
  inline TVar *&Elem() override
  {
    return fA.begin();
  }
  inline const TVar *Elem() const override
  {
    return fA.begin();
  }
  inline int64_t Size() const override
  {
    return fA.size();
  }
private:
	
	
	TPZVec<int64_t>  fIA;
	TPZVec<int64_t>  fJA;
	TPZVec<TVar> fA;
	
  TPZPardisoSolver<TVar> fPardisoControl;
	
	TPZVec<TVar> fDiag;
};

template<class TVar>
inline void TPZSYsmpMatrix<TVar>::SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A )
{
	// Pass the data to the class.
	fIA = IA;
	fJA = JA;
	fA  =  A;
	ComputeDiagonal();
}

#endif
