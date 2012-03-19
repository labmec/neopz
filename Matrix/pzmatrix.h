/**
 * @file
 * @brief Contains TPZMatrix<TVar>class, root matrix class.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatrix.hh
//
// Class:  TPZMatrix
//
// Obs.:   Implements abstract matrix class.
//
// Versao: 04 / 1996.
//


#ifndef _TMATRIXHH_
#define _TMATRIXHH_


#include "pzstream.h"
#include "pzreal.h"
#include "pzsave.h"

#include <list>
#include <sstream>

/** @brief To create clone matrix */
#define CLONEDEF(A) virtual TPZMatrix<TVar>*Clone() const { return new A(*this); }

/** @brief Gets maxime value between a and b */
#define MAX( a, b )   ( (a) > (b) ? (a) : (b) )
/** @brief Gets minime value between a and b */
#define MIN( a, b )   ( (a) < (b) ? (a) : (b) )
template<class TVar>
class TPZMatrix;

template<class TVar>
class TPZFMatrix;

template<class TVar>
class TPZSolver;

template<class T>
class TPZVec;

extern "C"{
	/// Extern BLAS FUNCTION 
	double ddot(int *N, double *X, int *INCX, double *Y, int *INCY);
}

/** \addtogroup matrix
 * @{
 */
/**
 * @enum DecomposeType
 * @brief Defines decomposition type for any matrix classes
 * @param ENoDecompose Not decomposed
 * @param ELU Decomposed using LU method
 * @param ECholesky Decomposed using Cholesky method
 * @param ELDLt Decomposed using LDLt method
 */
enum DecomposeType {ENoDecompose, ELU, ELUPivot, ECholesky, ELDLt};

/** @brief Defines output format */
enum MatrixOutputFormat {EFormatted, EInputFormat, EMathematicaInput, EMatlabNonZeros};

/** @brief Root matrix class (abstract). \ref matrix "Matrix" */
/** Abstract class TPZMatrix<TVar>which defines interface of derived matrix classes. */
template<class TVar=REAL>
class TPZMatrix: public TPZSaveable

{
public:
	/** @brief Simple constructor */
	TPZMatrix() {
		fDecomposed = 0;
		fDefPositive = 0;
		fRow = 0;
		fCol = 0;
		gZero = 0.;
	}
	
	TPZMatrix<TVar>(const TPZMatrix<TVar>&cp) : fRow(cp.fRow), fCol(cp.fCol), fDecomposed(cp.fDecomposed),fDefPositive(cp.fDefPositive)
	{
	}
	/** @brief Simple destructor */
	virtual ~TPZMatrix();
	
	virtual TPZMatrix<TVar>*Clone() const = 0;
	
	/** @brief Fill matrix storage with randomic values */
	/** This method use GetVal and PutVal which are implemented by each type matrices */
	void AutoFill();
	
	/** @brief Checks if current matrix value is symmetric */
	virtual int VerifySymmetry(REAL tol = 1.e-13) const;
	
	/**
     @brief Put values with bounds checking if DEBUG variable is defined.
     @param row Row number.
     @param col Column number.
     @param value Value being put.
	 */
	virtual int Put(const int row,const int col,const TVar & value );
	/**
	 * @brief Get value with bound checking
     * @param row Row number.
     * @param col Column number.
	 */
	virtual const TVar &Get(const int row,const int col ) const;
	
	/**
	 * @brief Substitution for the () operator when const arguments are needed
     * @param row Row number.
     * @param col Column number.
	 */
	const TVar &g(const int row, const int col) const {return Get(row,col);}
	
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
	 * @param col Column number.
	 */
	TVar &operator() (const int row,const int col );
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
     * @param col Column number.
	 */
	virtual TVar &s(const int row, const int col);
	/**
	 * @brief The operators check on the bounds if the DEBUG variable is defined
     * @param row Row number.
	 */
	TVar &operator()(const int row);
	
	/** @brief Put values without bounds checking \n
	 *  This method is faster than "Put" if DEBUG is defined.
	 */
	virtual int PutVal(const int /*row*/,const int /*col*/,const TVar & /*val*/ ) { return 0; }
	/** @brief Get values without bounds checking \n
	 *  This method is faster than "Get" if DEBUG is defined.
	 */
	virtual const TVar &GetVal(const int /*row*/, const int /*col*/ ) const  { return gZero; }
	
	/** @name Algebraic
	 *  @brief Implements algebraic operations with matrices
	 */
	//@{
	
	/**
	 * @brief It mutiplies itself by TPZMatrix<TVar>A putting the result in res
	 * @param A TPZMatrix<TVar>object to multiplied to
	 * @param res TPZFMatrix<TVar>containing the result
	 * @param opt Indicates if is Transpose or not
	 * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
	 */
	virtual void Multiply(const TPZFMatrix<TVar>& A,TPZFMatrix<TVar>& res, int opt = 0, int stride = 1) const;
	/**
	 * @brief It adds itself to TPZMatrix<TVar>A putting the result in res
	 * @param A TPZMatrix<TVar>to added to current matrix
	 * @param res Contains the result
	 */
	virtual void Add(const TPZMatrix<TVar>& A,TPZMatrix<TVar>& res) const;
	/**
	 * @brief It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
	 * @param x Is x on the above operation
	 * @param y Is y on the above operation
	 * @param z Is z on the above operation
	 * @param alpha Is alpha on the above operation
	 * @param beta Is beta on the above operation
	 * @param opt Indicates if is Transpose or not
	 * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
	 */
	virtual void MultAdd(const TPZFMatrix<TVar> & x,const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,
						 const TVar alpha=1., const TVar beta = 0., const int opt = 0, const int stride = 1 ) const;
	
	/**
	 * @brief Computes res = rhs - this * x
	 */
	virtual void Residual(const TPZFMatrix<TVar>& x,const TPZFMatrix<TVar>& rhs, TPZFMatrix<TVar>& res ) ;
	/**
	 * @brief It substracts A from storing the result in result
	 */
	virtual void Substract(const TPZMatrix<TVar>& A,TPZMatrix<TVar>& result) const;
	
	
	/** @brief Converts the matrix in an identity matrix*/
	virtual void Identity();
	/**
	 * @brief It makes *T the transpose of current matrix.
	 */
	virtual void Transpose(TPZMatrix<TVar>*const T) const;
	
	/**
	 * @brief It makes Inv =[this].
	 * IMPORTANT OBSERVATION --> The original matrix (calling object) no is more equal. 
	 * It containts the some decomposition (LU or Cholesky or ...)
	 */
	int Inverse(TPZFMatrix<TVar>&Inv);
	
	/**
	 * @brief Computes the matrix norm of this
	 * @param p interpolation order
	 * @param numiter is used by 2-norm calculation in the SolveEigenvaluesJacobi method required to compute the maximum eigenvalue
	 * @param tol - same of numiter
	 * @note It is available p-norm = 1, 2 and infinity.\n
	 * p=1 is the maximum absolute column sum norm \n
	 * p=2 is the spectral norm wich is the square root of the maximum eigenvalue of Tranpose[this].this \n
	 * p=infinity is the maximum absolute row sum norm - p infinity is implemented with p = 0 \n
	 * These operations are defined on the website of the Mathematica software: \n
	 * http://mathworld.wolfram.com/MatrixNorm.html \n
	 * Be careful when choosing 2-norm. It has a high computational cost.
	 
	 */
	TVar MatrixNorm(int p, int numiter = 2000000, REAL tol = 1.e-10) const; // -<
	
	/**
	 * @brief Computes the matrix condition number of this
	 *
	 * It is available p-norm = 1, 2 and infinity.
	 * p=1 is the maximum absolute column sum norm
	 * p=2 is the spectral norm wich is the square root of the maximum eigenvalue of Tranpose[this].this
	 * p=infinity is the maximum absolute row sum norm - p infinity is implemented with p = 0
	 * These operations are defined on the website of the Mathematica software:
	 * http://mathworld.wolfram.com/MatrixNorm.html
	 * All norms require the computation of the inverse matrix.
	 * It has a high computational cost and a high memory requirement.
	 */
	TVar ConditionNumber(int p, int numiter = 2000000, REAL tol = 1.e-10);
	
	//@}
	
	/**
	 * @name Generic
	 * @brief Performs generic operations on matrices
	 */
	
	//@{
	/** @brief Input operation*/
	virtual void Input(std::istream & in = std::cin );
	
	/** @brief Input operation*/
	template <class TT> friend std::istream & operator >> (std::istream& in,TPZMatrix<TT>& A) ;
	
	/** @brief It prints the matrix data in a MatrixFormat Rows X Cols */
	virtual void Print(const char *name = 0, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	
	void PrintMath(const char *name = 0, std::ostream &out = std::cout);
	
	/** @brief Returns number of rows */
	int Rows() const;
	/** @brief Returns number of cols */
	int Cols() const;
	
	/**
	 * @brief Returns the dimension of the matrix if the matrix is square. \n
	 * If the matrix is not square, returns an error
	 */
	inline virtual int Dim() const;
	
	
	/**
	 * @brief Redimensions a matriz keeping the previous values
	 * @param newRows Specifies the new number of rows in matrix
	 * @param newCols Specifies the new number of Columns in matrix
	 */
	virtual int Resize(const int newRows, const int newCols ) {
		fRow = newRows;
		fCol = newCols;
		return 0;
	}
	
	/**
	 * @brief Redimensions the matrix reinitializing it with zero
	 * @param newRows Specifies the new number of rows in matrix.
	 * @param newCols Specifies the new number of Columns in matrix.
	 */
	virtual int Redim(const int newRows, const int newCols ) {
		fRow = newRows;
		fCol = newCols;
		return 0;
	}
	
	/** @brief Zeroes the matrix */
	virtual int Zero(){
		std::cout << "WARNING! TPZMatrix<TVar>::Zero is called\n";
		return 0; }
	
	//@}
	/**
	 * @name SubMatrices
	 * Operations with SUB MATRICES.
	 */
	
	//@{
	/**
	 * @brief It puts submatrix Source on actual matrix structure.
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param Source The matrix to be inserted
	 */
	virtual int PutSub( const int sRow, const int sCol, const TPZFMatrix<TVar>& Source );
	
	/**
	 * @brief Gets submatrix storing it on Target.
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param Target The matrix to be aquired.
	 */
	virtual int GetSub( const int sRow, const int sCol, const int rowSize,
					   const int colSize, TPZFMatrix<TVar>& Target ) const;
	
	/**
	 * @brief It adds Source matrix on current matrix from position (sRow, sCol)
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param Source The matrix to be added
	 */
	virtual int AddSub(const int sRow, const int sCol, const TPZFMatrix<TVar>& Source );
	
	//insere uma submatriz do objeto corrente para *Target sem redimencionar-la
	/**
	 * @brief Inserts a submatrix from current object on matrix *Target with no \n
	 * redimentioning
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param pRow Specifies final row on current object
	 * @param pCol Specifies final column on current object.
	 * @param Target The matrix to be inserted.
	 */
	virtual int InsertSub(const int sRow,const int sCol,const int rowSize,const int colSize,
						  const int pRow,const int pCol, TPZMatrix<TVar>* Target ) const;
	
	//adiciona uma submatriz do objeto corrente para *target
	//10-10-95
	/**
	 * @brief Adds a submatrix from current object in *Target
	 * @param sRow Specifies starting row on current object
	 * @param sCol Specifies starting column on current object.
	 * @param rowSize Specifies the amount of rows from sRow
	 * @param colSize Specifies the amount of columns from sCol
	 * @param pRow Specifies starting row on pA
	 * @param pCol Specifies starting column on pA.
	 * @param pA The matrix to be added.
	 */
	virtual int AddSub(const int sRow,const int sCol,const int rowSize,
					   const int colSize,const int pRow,const int pCol, TPZMatrix<TVar>* pA ) const;
	
	
	//@}
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * @param elmat Element matrix to be contributed
	 * @param destinationindex Contains destine indexes on current matrix
	 */
	virtual  void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int> &destinationindex);
	
	/**
	 * @brief Add a contribution of a stiffness matrix
	 * @param elmat Element matrix to be contributed
	 * @param sourceindex Contains source indexes on current matrix
	 * @param destinationindex Contains destine indexes on current matrix
	 */
	virtual  void AddKel(TPZFMatrix<TVar>&elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex);
	
	/**
	 * @name Inquire
	 * @brief Returns information of the current object
	 */
	
	/**
	 * @brief Updates the values of the matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > /* mat*/)    
	{
		std::cout << "TPZMatrix<TVar>::UdateFrom is not implemented\n";
        DebugStop();
	}
	//@{
	
	/**
	 * @brief Checks if the current matrix is symmetric
	 */
	
	virtual int IsSimetric() const    { return 0; }
	/**
	 * @brief Checks if current matrix is square
	 */
	inline int IsSquare() const { return fRow == fCol;}
	
	/**
	 * @brief Checks if current matrix is definite positive
	 */
	virtual int IsDefPositive() const{ return 0; }
	/**
	 * @brief Checks if current matrix is already decomposed
	 */
	int IsDecomposed() const         { return fDecomposed; }
	//@}
	/**
	 * @brief Sets current matrix to decomposed state
	 */
	void SetIsDecomposed(int val) {fDecomposed = (char) val; }
	
	
	/**
	 * @name Solvers
	 * @brief Linear system solvers. \n
	 * 
	 * For symmetric decompositions lower triangular matrix is used. \n
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
	virtual void SolveJacobi(int & numiterations, const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							 TPZFMatrix<TVar>* residual, TPZFMatrix<TVar>& scratch, REAL & tol, const int FromCurrent = 0);
	
	/**
	 * @brief Solves the linear system using Successive Over Relaxation method (Gauss Seidel). \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param overrelax The over relaxation parameter
	 * @param tol The tolerance value..
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 * @param direction Indicates interaction direction, from first to last (default 1) or from last to first (-1)
	 */
	virtual void SolveSOR(int & numiterations, const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						  TPZFMatrix<TVar>* residual,TPZFMatrix<TVar>& scratch,const TVar overrelax, TVar & tol,
						  const int FromCurrent = 0,const int direction = 1) ;
	/**
	 * @brief Solves the linear system using Symmetric Successive Over Relaxation method (Gauss Seidel). \n
	 * @param numiterations The number of interations for the process.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param scratch Available manipulation area on memory.
	 * @param overrelax The over relaxation parameter
	 * @param tol The tolerance value..
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveSSOR(int & numiterations,const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						   TPZFMatrix<TVar>* residual, TPZFMatrix<TVar>& scratch, const TVar overrelax, TVar & tol,
						   const int FromCurrent = 0) ;
	
	/**
	 * @brief Solves the linear system using Conjugate Gradient method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent.
	 */
	virtual void SolveCG(int & numiterations, TPZSolver<TVar> & preconditioner,
						 const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						 TPZFMatrix<TVar>* residual, TVar & tol,
						 const int FromCurrent = 0) ;
	/**
	 * @brief Solves the linear system using Bi-Conjugate Gradient method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param tol The tolerance value.
	 */
	virtual void SolveBICG(int & numiterations, TPZSolver<TVar> & preconditioner,
						   const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						   TVar & tol) ;
	
	/**
	 * @brief Solves the linear system using Bi-Conjugate Gradient stabilized method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent.
	 */
	virtual void SolveBICGStab(int & numiterations, TPZSolver<TVar> & preconditioner,
							   const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							   TPZFMatrix<TVar>* residual, TVar & tol,
							   const int FromCurrent = 0) ;
	
	/**
	 * @brief Solves the linear system using Generalized Minimal Residual (GMRES) method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param H The right hand side of the system
	 * @param numvectors The number of vectors involved
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveGMRES(int & numiterations, TPZSolver<TVar> & preconditioner,
							TPZFMatrix<TVar>& H, int & numvectors,
							const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
							TPZFMatrix<TVar>* residual, TVar & tol,const int FromCurrent) ;
	
	/**
	 * @brief Solves the linear system using IR method. \n
	 * @param numiterations The number of interations for the process.
	 * @param preconditioner The preconditioner attribute used.
	 * @param F The right hand side of the system.
	 * @param result The solution.
	 * @param residual Returns F - A*U which is the solution residual.
	 * @param tol The tolerance value.
	 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
	 */
	virtual void SolveIR(int & numiterations, TPZSolver<TVar> & preconditioner,
						 const TPZFMatrix<TVar>& F, TPZFMatrix<TVar>& result,
						 TPZFMatrix<TVar>* residual, TVar & tol,
						 const int FromCurrent = 0);
	
	/** @brief Transforms this matrix in a diagonal matrix, where the diagonal values are its eigenvalues.
	 * This method is efficient only for small matrices.
	 * @param numiterations The number of interations for the process.
	 * @param tol The tolerance value.
	 * @param Sort diagonal values from big to small
	 * @return Returns true if tolerance is achieved or false otherwise.
	 */
	virtual bool SolveEigenvaluesJacobi(int &numiterations, REAL &tol, TPZVec<REAL> * Sort = 0);
	
	/** @brief Compute Eigenvalues and Eigenvectors of this matrix. \n
	 * This method is efficient only for small matrices.
	 * @param numiterations The number of interations for the process.
	 * @param tol The tolerance value.
	 * @param Eigenvalues ordered from big to small
	 * @param Eigenvectors: each row represent one eigenvector. It is in same order of eigenvalues.
	 * @return Returns true if tolerance is achieved or false otherwise.
	 */
	virtual bool SolveEigensystemJacobi(int &numiterations, REAL & tol, TPZVec<REAL> & Eigenvalues, TPZFMatrix<TVar>& Eigenvectors) const;
	
	/**
	 * @brief Solves the linear system using Direct methods
	 * @param F The right hand side of the system and where the solution is stored.
	 * @param dt Indicates type of decomposition
	 * @param singular
	 */
	virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt, std::list<int> &singular);
	/**
	 * @brief Solves the linear system using Direct methods
	 * @param F The right hand side of the system and where the solution is stored.
	 * @param dt Indicates type of decomposition
	 */
	virtual int SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt);
	
    /** @brief decompose the system of equations acording to the decomposition scheme */
    virtual int Decompose(const DecomposeType dt, std::list<int> &singular)
    {
        switch (dt) {
            case ELU:
                return Decompose_LU(singular);
                break;
            case ELDLt:
                return Decompose_LDLt(singular);
                break;
            case ECholesky:
                return Decompose_Cholesky(singular);
                break;                
            default:
                DebugStop();
                break;
        }
        return -1;
    }
	/**
	 * @brief Retorna o valor mais proximo a "val" (exceto valores no intervalo -tol <= val <= +tol) contido no vetor Vec
	 */
	static TVar ReturnNearestValue(TVar val, TPZVec<REAL> &Vec, TVar tol);
	
	/**
	 * @brief Solves the linear system using LU method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 * @param singular
	 */
	int Solve_LU ( TPZFMatrix<TVar>* B, std::list<int> &singular );
	/**
	 * @brief Solves the linear system using LU method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 */
	int Solve_LU ( TPZFMatrix<TVar>* B );
	
	//para usar estos solves e' responsabilidade do usuario
	//que a matriz corrente seja simetrica
	/**
	 * @brief Solves the linear system using Cholesky method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 */
	int Solve_Cholesky( TPZFMatrix<TVar>* B);
	/**
	 * @brief Solves the linear system using Cholesky method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 * @param singular
	 */
	int Solve_Cholesky( TPZFMatrix<TVar>* B, std::list<int> &singular );
	/**
	 * @brief Solves the linear system using LDLt method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 * @param singular
	 */
	int Solve_LDLt    ( TPZFMatrix<TVar>* B, std::list<int> &singular );
	/**
	 * @brief Solves the linear system using LDLt method\n
	 * @param B The right hand side of the system and where the solution is stored.
	 */
	int Solve_LDLt    ( TPZFMatrix<TVar>* B);
	
	//@}
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
	
	// Decompoe a matriz em GGt (onde G e' triangular inferior).
	//a matriz corrente tem que ser simetrica.
	/**
	 * @brief Decomposes the current matrix using Cholesky method. \n
	 * The current matrix has to be symmetric.
	 */
	virtual int Decompose_Cholesky() ;
	/**
	 * @brief Decomposes the current matrix using Cholesky method.
	 * @param singular
	 */
	virtual int Decompose_Cholesky(std::list<int> &singular) ;
	
	// Decompoe a matriz em LDLt (onde L e' triangular inferior com
	//  1.0 na diagonal, e D e' uma matriz diagonal).
	//a matriz a ser decomposta tem que ser simetrica
	/**
	 * @brief Decomposes the current matrix using LDLt. \n
	 * The current matrix has to be symmetric.
	 * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
	 */
	virtual int Decompose_LDLt(std::list<int> &singular);
	/** @brief Decomposes the current matrix using LDLt. */
	virtual int Decompose_LDLt();
	
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
	virtual int Substitution( TPZFMatrix<TVar>* B ) const;
	
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
	//@}
	/**
	 * @name TPZSaveable
	 * Methods which would make TPZMatrix<TVar>compliant with TPZSaveable
	 */
	//@{
	/**
	 * @brief Unpacks the object structure from a stream of bytes
	 * @param buf The buffer containing the object in a packed form
	 * @param context 
	 */
	virtual void  Read(TPZStream &buf, void *context );
	/**
	 * @brief Packs the object structure in a stream of bytes
	 * @param buf Buffer which will receive the bytes
	 * @param withclassid
	 */
	virtual void Write( TPZStream &buf, int withclassid );
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	
	/**
	 * @brief Extract the block indicated by the indices from the matrix
	 */
	virtual void GetSub(const TPZVec<int> &indices,TPZFMatrix<TVar>&block) const;
	
	/**
	 * @brief Compare values of this to B, with a precision tolerance tol.
	 */
    bool CompareValues(TPZMatrix<TVar>&M, TVar tol);
	
protected:
	
	/**
	 * @brief Is an auxiliar method used by MultiplyAdd
	 * @see MultAdd
	 */
	void PrepareZ(const TPZFMatrix<TVar>& y, TPZFMatrix<TVar>& z,const TVar beta,const int opt,const int stride) const;
	
	/**
	 * @brief Constructor
	 * @param row Number of rows
	 * @param col Number of cols
	 */
	inline  TPZMatrix<TVar>(const int row,const int col )
	{ fRow = row; fCol = col;fDefPositive=0; fDecomposed = 0;}
	
public:
	/**
	 * @brief Returns error messages
	 * @param msg First message.
	 * @param msg2 Second message.
	 */
	static int Error(const char *msg ,const char *msg2 = 0);
	
protected:
	/**
	 * @brief It clears data structure.
	 */
	virtual int Clear() { return 0; }
	
	/**
	 * @brief Swaps contents of a in b and b in a
	 */
	static void Swap(int *a, int *b);
	/* inline void
	 TPZMatrix<TVar>::Swap( int *a, int *b ) const
	 {
	 int aux = *a;
	 *a = *b;
	 *b = aux;
	 }
	 */
	/** @brief Number of rows in matrix */
	int fRow;
	/** @brief Number of cols in matrix */
	int fCol;
	/** @brief Decomposition type used to decompose the current matrix */
	char  fDecomposed;
	/** @brief Definite Posistiveness of current matrix */
	char fDefPositive;
	static TVar gZero;
};

/** @} */

/** @brief Overload << operator to print entries of the matrix ***/
template<class TVar>
std::ostream & operator<<(std::ostream& out, const TPZMatrix<TVar> & A);

/******** Inline ********/


template<class TVar>
inline int TPZMatrix<TVar>::Rows() const {
	return fRow;
}


template<class TVar>
inline int TPZMatrix<TVar>::Cols() const {
	return fCol;
}

template<class TVar>
inline void TPZMatrix<TVar>::Residual(const TPZFMatrix<TVar>& x,const TPZFMatrix<TVar>& rhs, TPZFMatrix<TVar>& res )  {
	MultAdd( x, rhs, res, -1.0, 1.0 );
}

/***********/
/*** Put ***/

template<class TVar>
inline int TPZMatrix<TVar>::Put(const int row,const int col,const TVar & value ) {     
	// verificando se o elemento a inserir esta dentro da matriz
#ifdef DEBUG
	if ( row >= Rows() || col >= Cols() || row <0 || col < 0 ) {
		std::stringstream sout;
		sout << "TPZMatrix<TVar>::Put" << " Index out of range row = " << row << " col = " << col << " Rows() " << Rows() << " Cols() " << Cols() << std::endl;
		Error(sout.str().c_str());
	}
#endif
	
	return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

template<class TVar>
inline const TVar &TPZMatrix<TVar>::Get(const int row, const int col ) const {
	// verificando se o elemento pedido esta dentro da matriz
#ifdef DEBUG
	if ( (row >= Rows()) || (col >= Cols()) || row <0 || col <0 ) {
		Error("TPZMatrix::Get", "Index out of range");
		gZero=0.;
		return gZero;
	}
#endif
	return( GetVal( row, col ) );
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::operator()(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
#ifndef NODEBUG
	if ( (row >= Rows()) || (col >= Cols()) || row <0 || col<0 ) {
		gZero = 0.;
		Error("TPZMatrix<TVar>::Operator()","Index out of range");
		return gZero;
	}
#endif
	return s(row,col);
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::s(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
	std::cout << "TPZMatrix<TVar>::s not implemented\n";
	return gZero;
}

template<class TVar>
inline TVar &TPZMatrix<TVar>::operator()(const int row) {
	return operator()(row,0);
}

template<class TVar>
inline int TPZMatrix<TVar>::Dim() const{
	if ( IsSquare() ) return Rows();
	Error( "matrix is not square" );
	return ( 0 );
}
//***Solve LU ***/

template<class TVar>
inline int TPZMatrix<TVar>::Solve_LU( TPZFMatrix<TVar>*B, std::list<int> &singular) {
	if ( IsSimetric() ) Error( "LU decomposition is a not symetric decomposition" );
	return ( ( !Decompose_LU(singular) )?  0 : Substitution( B )  );
}

template<class TVar>
inline int TPZMatrix<TVar>::Solve_LU( TPZFMatrix<TVar>*B ) {
	if ( IsSimetric() ) Error( "LU decomposition is a not symetric decomposition" );
	return ( ( !Decompose_LU() )?  0 : Substitution( B )  );
}
/**********************/
/*** Solve Cholesky ***/
//
//  Se nao conseguir resolver por Cholesky retorna 0 e a matriz
//   sera' modificada (seu valor perdera' o sentido).
//
template<class TVar>
inline int TPZMatrix<TVar>::Solve_Cholesky( TPZFMatrix<TVar>* B )
{
	return(
		   ( !Decompose_Cholesky() )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
		   );
}

template<class TVar>
inline int TPZMatrix<TVar>::Solve_Cholesky( TPZFMatrix<TVar>* B, std::list<int> &singular ) {
	return(
		   ( !Decompose_Cholesky(singular) )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
		   );
}



/******************/
/*** Solve LDLt ***/

template<class TVar>
inline int TPZMatrix<TVar>::Solve_LDLt( TPZFMatrix<TVar>* B ) {
	
	return(
		   ( !Decompose_LDLt() )? 0 :
		   ( Subst_LForward( B ) && Subst_Diag( B ) && Subst_LBackward( B ) )
		   );
}



template<class TVar>
inline int TPZMatrix<TVar>::Solve_LDLt( TPZFMatrix<TVar>* B, std::list<int> &singular ) {
	
	return(
		   ( !Decompose_LDLt(singular) )? 0 :
		   ( Subst_LForward( B ) && Subst_Diag( B ) && Subst_LBackward( B ) )
		   );
}

template<class TVar>
inline void
TPZMatrix<TVar>::Swap( int *a, int *b )
{
	int aux = *a;
	*a = *b;
	*b = aux;
}
/******************/
/***  Multiply ***/


#endif
// _TMATRIXH_

