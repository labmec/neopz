/**
 * @file
 * @brief Contains the TPZFYsmpMatrix class which implements a non symmetric sparse matrix.
 */
/******************************************************************************
 *
 * Class definition:    TPZFYsmpMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on general sparse matrices stored
 *                      in the (old) Yale Sparse Matrix Package format.
 *
 * Solvers:             SOR
 *
 *****************************************************************************/

#ifndef YSMPMATH
#define YSMPMATH
#ifdef USING_BLAS
extern "C"{
     #include "cblas.h"
     };
#endif

#include "pzmatrix.h"
#include "tpzverysparsematrix.h" 


class TPZFMatrix;

/**
 * @brief Implements a non symmetric sparse matrix (Yale Sparse Matrix Storage). \ref matrix "Matrix"
 * @ingroup matrix
 */
/**
 * Defines operations on general sparse matrices stored in the (old) Yale Sparse Matrix Package format.
 */
class TPZFYsmpMatrix : public TPZMatrix {
	
	public :
	
	/** @brief An auxiliary structure to hold the data of the subset \n of equations used to multiply in a multitrheaded environment */
	/**
	 In future versions this structure should be defined in a derived class
	 */
	struct TPZMThread {
		const TPZFYsmpMatrix *target;
		int fFirsteq;
		int fLasteq;
		const TPZFMatrix *fX;
		TPZFMatrix *fZ;
		REAL fAlpha;
		int fOpt;
		int fStride;
	};
	
private:
	
	static void * ExecuteMT(void *entrydata);
	
public: 
    TPZFYsmpMatrix(const int rows,const int cols );
	
    //TPZFYsmpMatrix(const TPZFYsmpMatrix &cp) : TPZMatrix(cp), fSymmetric(cp.fSymmetric)
	//    {
	//      int fjasize = fIA[Rows()];
	//      fIA = new int[Rows()+1];
	//      fDiag = new REAL[Rows()];
	//      fJA = new int[fjasize];
	//			fA = new REAL[fjasize];
	//			memcpy(fIA,cp.fIA,(Rows()+1)*sizeof(int));
	//			memcpy(fJA,cp.fJA,fjasize*sizeof(int));
	//			memcpy(fDiag,cp.fDiag,Rows()*sizeof(REAL));
	//			memcpy(fA,cp.fA,fjasize*sizeof(REAL));
	//		}
	
	// Replace the above destructor
	
	TPZFYsmpMatrix(const TPZVerySparseMatrix &cp);
	
	TPZFYsmpMatrix &operator=(const TPZFYsmpMatrix &copy);
	
	TPZFYsmpMatrix &operator=(const TPZVerySparseMatrix &cp);
    
	CLONEDEF(TPZFYsmpMatrix)
	
	virtual ~TPZFYsmpMatrix();
	
	
	/** @brief Get the matrix entry at (row,col) without bound checking */
	virtual const REAL &GetVal(const int row,const int col ) const;
	
	int NumTerms() 
	{
		return fIA[Rows()];
	}
	
	int PutVal(const int row, const int col, const REAL &Value);
	
	virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						 const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 ) const;
	
	virtual void MultAddMT(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						   const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 );
	
	virtual int GetSub(const int sRow,const int sCol,const int rowSize,
					   const int colSize, TPZFMatrix & A ) const;
	
	void GetSub(const TPZVec<int> &indices,TPZFMatrix &block) const;
	
	/** @brief Pass the data to the class. */
	virtual void SetData( int *IA, int *JA, REAL *A );
	
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
	virtual void SolveJacobi(int & numiterations, const TPZFMatrix & F, TPZFMatrix & result,
							 TPZFMatrix * residual, TPZFMatrix & scratch, REAL & tol, const int FromCurrent = 0) ;
	
	void SolveSOR(int &numiterations, const TPZFMatrix &rhs, TPZFMatrix &x,
				  TPZFMatrix *residual, TPZFMatrix &scratch,
				  const REAL overrelax, REAL &tol,
				  const int FromCurrent = 0,const int direction = 1 ) ;    
	// @}
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * putting it on destination indexes position
	 */
	virtual void AddKelOld(
						   TPZFMatrix & elmat //! Member stiffness matrix beeing added
						   , TPZVec < int > & destinationindex //! Positioning of such members on global stiffness matrix
						   );    
	
	void AddKel(TPZFMatrix & elmat, TPZVec<int> & destinationindex);    
	
	void AddKel(TPZFMatrix & elmat, TPZVec<int> & sourceindex, TPZVec<int> & destinationindex);    
	
	void MultiplyDummy(TPZFYsmpMatrix & B, TPZFYsmpMatrix & Res);
	
	virtual int Zero();
	
	/**
	 * @name Factorization
	 * @brief Those member functions are related to matrices factorization
	 */
	//@{
	/**
	 * @brief Decomposes the current matrix using LU decomposition.
	 */
	virtual int Decompose_LU(std::list<int> &singular);
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
	virtual int Substitution( TPZFMatrix * B ) const;
	
	//@}
	
	
private:
	
	void ComputeDiagonal();
	
	/*
	 * @brief Perform row update of the sparse matrix
	 */
	void RowLUUpdate(int sourcerow, int destrow);
	
protected:
	int  *fIA;
	int  *fJA;
	REAL *fA;
	
	REAL *fDiag;
	
	int   fSymmetric;
	
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



inline void TPZFYsmpMatrix::SetData( int *IA, int *JA, REAL *A ) {
	// Pass the data to the class.
	fIA = IA;
	fJA = JA;
	fA  =  A;
	ComputeDiagonal();
}

#endif
