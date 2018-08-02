/**
 * @file
 * @brief Contains TPZSYsmpMatrix class which implements a symmetric sparse matrix. \n
 * Purpose: Defines operations on symmetric sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */

#ifndef SYSMPMATH
#define SYSMPMATH

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef USING_MKL

#include "TPZPardisoControl.h"
#endif
 /**
  * @brief Implements a symmetric sparse matrix. \ref matrix "Matrix"
  * @ingroup matrix
  */
template<class TVar>
class TPZSYsmpMatrix : public TPZMatrix<TVar>{
	
#ifdef USING_MKL
    friend class TPZPardisoControl<TVar>;
#endif
    
public :
    /** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrix();
	/** @brief Constructor based on number of rows and columns */
    TPZSYsmpMatrix(const int64_t rows, const int64_t cols );
	/** @brief Copy constructor */
    TPZSYsmpMatrix(const TPZSYsmpMatrix<TVar> &cp) : 
    TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
    TPZMatrix<TVar>(cp), fIA(cp.fIA), fJA(cp.fJA), fA(cp.fA), fDiag(cp.fDiag)
#ifdef USING_MKL
    , fPardisoControl(cp.fPardisoControl)
#endif
    {
#ifdef USING_MKL
        fPardisoControl.SetMatrix(this);
#endif
    }
    
    TPZSYsmpMatrix &operator=(const TPZSYsmpMatrix<TVar> &copy);
    
    CLONEDEF(TPZSYsmpMatrix)
	/** @brief Destructor */
	virtual ~TPZSYsmpMatrix();
    
    /** @brief Checks if the current matrix is symmetric */
    virtual int IsSimetric() const    { return 1; }
    /** @brief Checks if current matrix is square */
    inline int IsSquare() const { return 1;}
    
    /** @brief Zeroes the matrix */
    virtual int Zero(){
        fA.Fill(0.);
        fDiag.Fill(0.);
#ifndef USING_MKL
        TPZMatrix<TVar>::fDecomposed = ENoDecompose;
#endif
        return 0;
        
    }

    
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric);
	
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar &GetVal(const int64_t row, const int64_t col ) const;
    
    /** @brief Put values without bounds checking \n
     *  This method is faster than "Put" if DEBUG is defined.
     */
    virtual int PutVal(const int64_t /*row*/,const int64_t /*col*/,const TVar & val );

	
	/** @brief Computes z = beta * y + alpha * opt(this)*x */
	/** @note z and x cannot overlap in memory */
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0) const ;
	
	/** @brief Sets data to the class */
	virtual void SetData(const TPZVec<int64_t> &IA,const TPZVec<int64_t> &JA, const TPZVec<TVar> &A );
	
    /// Access function for the coefficients
    TPZVec<TVar> &A()
    {
        return fA;
    }
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout ,const MatrixOutputFormat = EFormatted ) const;
    
#ifdef USING_MKL
    /**
     * @name Factorization
     * @brief Those member functions perform the matrix factorization
     * @{
     */
    

    /**
     * @brief Decomposes the current matrix using LDLt. \n
     * The current matrix has to be symmetric.
     * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
     */
    virtual int Decompose_LDLt(std::list<int64_t> &singular);
    /** @brief Decomposes the current matrix using LDLt. */
    virtual int Decompose_LDLt();

    /** @brief Decomposes the current matrix using Cholesky method. The current matrix has to be symmetric. */
    virtual int Decompose_Cholesky() ;
    /**
     * @brief Decomposes the current matrix using Cholesky method.
     * @param singular
     */
    virtual int Decompose_Cholesky(std::list<int64_t> &singular) ;
    

    /** @} */
    
    /**
     * @name Substitutions
     * @brief Substitutions forward and backward
     * @{
     */
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LForward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
     * @param b right hand side and result after all
     */
    virtual int Subst_LBackward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is diagonal matrix.
     * @param b right hand side and result after all
     */
    virtual int Subst_Diag( TPZFMatrix<TVar>* b ) const;

    /**
     * @brief Computes B = Y, where A*Y = B, A is lower triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Forward( TPZFMatrix<TVar>* b ) const;
    
    /**
     * @brief Computes B = Y, where A*Y = B, A is upper triangular.
     * @param b right hand side and result after all
     */
    virtual int Subst_Backward( TPZFMatrix<TVar>* b ) const;
    

    /** @} */
    

#endif
    public:
virtual int ClassId() const;

    void ComputeDiagonal();

private:
	
	
	TPZVec<int64_t>  fIA;
	TPZVec<int64_t>  fJA;
	TPZVec<TVar> fA;
	
#ifdef USING_MKL
    TPZPardisoControl<TVar> fPardisoControl;
#endif
	
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
