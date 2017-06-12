/**
 * @file
 * @brief Contains the TPZFYsmpMatrix class which implements a non symmetric sparse matrix.
 */

#ifndef YSMPMATH
#define YSMPMATH
#ifdef USING_BLAS
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#elif USING_MKL
#include <mkl.h>
#else
extern "C"{
     #include "cblas.h"
     };
#endif
#endif

template<class TVar>
class TPZVerySparseMatrix;

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef USING_MKL
#include "TPZPardisoControl.h"
#endif

/**
 * @brief Implements a non symmetric sparse matrix (Yale Sparse Matrix Storage). \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Defines operations on general sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
template<class TVar>
class TPZFYsmpMatrix : public TPZMatrix<TVar> {
	
	public :
	
	/** @brief An auxiliary structure to hold the data of the subset \n of equations used to multiply in a multitrheaded environment */
	/**
	 In future versions this structure should be defined in a derived class
	 */
	struct TPZMThread {
		const TPZFYsmpMatrix<TVar> *target;
		long fFirsteq;
		long fLasteq;
		const TPZFMatrix<TVar> *fX;
		TPZFMatrix<TVar> *fZ;
		TVar fAlpha;
		int fOpt;
	};
	
private:
	
	static void * ExecuteMT(void *entrydata);
	
public:
    
    TPZFYsmpMatrix();
    
    TPZFYsmpMatrix(const long rows,const long cols );
	
	TPZFYsmpMatrix(const TPZVerySparseMatrix<TVar> &cp);
	
	TPZFYsmpMatrix &operator=(const TPZFYsmpMatrix<TVar> &copy);
	
	TPZFYsmpMatrix &operator=(const TPZVerySparseMatrix<TVar> &cp);
    
	CLONEDEF(TPZFYsmpMatrix)
	
	virtual ~TPZFYsmpMatrix();	
	
    /** @brief Fill matrix storage with randomic values */
    /** This method use GetVal and PutVal which are implemented by each type matrices */
    void AutoFill(long nrow, long ncol, int symmetric);
    

    
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const TVar &GetVal(const long row,const long col ) const;
	
	long NumTerms()
	{
		return fIA[this->Rows()];
	}
	
	int PutVal(const long row, const long col, const TVar &Value);
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0., const int opt = 0) const;
	
	virtual void MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						   const TVar alpha=1.,const TVar beta = 0., const int opt = 0);
	
	virtual int GetSub(const long sRow,const long sCol,const long rowSize,
					   const long colSize, TPZFMatrix<TVar> & A ) const;
	
	void GetSub(const TPZVec<long> &indices,TPZFMatrix<TVar> &block) const;
	
	/** @brief Pass the data to the class. */
	virtual void SetData( long *IA, long *JA, TVar *A );
    
    /** @brief Pass the data to the class. */
    virtual void SetData( TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<TVar> &A );
	
	/** @brief Print the matrix along with a identification title */
	virtual void Print(const char *title, std::ostream &out = std::cout , const MatrixOutputFormat form = EFormatted) const;
	
	/**
	 * @name Solvers
	 * @brief Linear system solvers. \n
	 */
	 /** For symmetric decompositions lower triangular matrix is used. \n
	 * Solves a system A*X = B returning X in B
	 */  
	//@{
	/**
	 * @brief Solves the linear system using Jacobi method. \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveJacobi(long & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result,
							 TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent = 0) ;
	
	void SolveSOR(long &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
				  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 ) ;    
	// @}
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * putting it on destination indexes position
	 */
	virtual void AddKelOld(
						   TPZFMatrix<TVar> & elmat //! Member stiffness matrix beeing added
						   , TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
						   );    
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<long> & destinationindex);
	
	virtual void AddKel(TPZFMatrix<TVar> & elmat, TPZVec<long> & sourceindex, TPZVec<long> & destinationindex);
	
	void MultiplyDummy(TPZFYsmpMatrix<TVar> & B, TPZFYsmpMatrix<TVar> & Res);
	
	virtual int Zero();
	
	/**
	 * @name Factorization
	 * @brief Those member functions are related to matrices factorization
	 */
	//@{
	/**
	 * @brief Decomposes the current matrix using LU decomposition.
	 */
	virtual int Decompose_LU(std::list<long> &singular);
	virtual int Decompose_LU();
	
	//@}
	
	/**
	 * @name Substitutions
	 * @brief Substitutions forward and backward
	 */
	//@{  
	/**
	 * @brief Computes Forward and Backward substitution for a "LU" decomposed matrix.
	 * @param B right hand side and result after all
	 */
	virtual int Substitution( TPZFMatrix<TVar> * B ) const;
	
	//@}
	
	
private:
	
	void ComputeDiagonal();
	
	/*
	 * @brief Perform row update of the sparse matrix
	 */
	void RowLUUpdate(long sourcerow, long destrow);
	
protected:
	TPZVec<long>  fIA;
	TPZVec<long>  fJA;
	TPZVec<TVar> fA;
	
	TPZVec<TVar> fDiag;
	
	int   fSymmetric;
	
#ifdef USING_MKL
    friend class TPZPardisoControl<TVar>;
    
    TPZPardisoControl<TVar> fPardisoControl;
#endif
protected:
	
	/**
	 * @brief Implements a initialization method for the sparse structure. It sets the initial value for the fIA and fJA.
	 */ 
	/**
	 * -fIA will contain the initial positions for all the equations
	 * -fJA will contain (-1) on all its positions
	 * -fA will contain 0 on all its value 
	 */
	void InitializeData();
};


template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( long *IA, long *JA, TVar *A ) {
	// Pass the data to the class.
    int nel = this->Rows()+1;
    fIA.resize(nel);
    memccpy(&fIA[0], IA, nel, sizeof(long));
    long nval = fIA[nel-1];
    fJA.resize(nval);
    memccpy(&fJA[0], JA, nval, sizeof(long));
    fA.resize(nval);
    memccpy(&fA[0], A, nval, sizeof(TVar));
	ComputeDiagonal();
}

/** @brief Pass the data to the class. */
template<class TVar>
inline void TPZFYsmpMatrix<TVar>::SetData( TPZVec<long> &IA, TPZVec<long> &JA, TPZVec<TVar> &A ){
    
    if (IA.size() != this->Rows() + 1 ) {
        DebugStop();
    }
    
    if (JA.size() != IA[this->Rows()]) {
        DebugStop();
    }
    
    fIA = IA;
    fJA = JA;
    fA = A;
}

#endif
