
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatrix.hh
//
// Class:  TPZMatrix
//
// Obs.:   Implementa matrizes cheias (normais).
//
// Versao: 04 / 1996.
//


#ifndef _TMATRIXHH_
#define _TMATRIXHH_


#include "pzstream.h"
#include "pzreal.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

class TPZFMatrix;
class TPZSolver;

template<class T>
class TPZVec;

extern "C"{
//BLAS FUNTIONS
double ddot(int *N, double *X, int *INCX, double *Y, int *INCY);
}

/**
   @enum DecomposeType
 * Defines decomposition type for any matrix classes
 * @param ENoDecompose Not decomposed
 * @param ELU Decomposed using LU method
 * @param ECholesky Decomposed using Cholesky method
 * @param ELDLt Decomposed using LDLt method
 */
enum DecomposeType {ENoDecompose, ELU, ECholesky, ELDLt};

enum MatrixOutputFormat {EFormatted, EInputFormat};


/**
 * Abstract class TPZMatrix which defines interface of derived matrix classes. 
 * @ingroup matrix
 */
class TPZMatrix
#ifdef OOPARLIB
  : public TSaveable
#endif

{
public:
  /**
   *Simple constructor
   */
  TPZMatrix() {
    fDecomposed = 0;
    fDefPositive = 0;
    fRow = 0;
    fCol = 0;
    gZero = 0.;
  }  
  /**
   *Simple destructor
   */

  virtual ~TPZMatrix() {
    fDecomposed = 0;
    fDefPositive = 0;
    fRow = 0;
    fCol = 0; 
  }

  /**
     Put values with bounds checking if DEBUG variable is defined.
     @param row Row number.
     @param col Column number.
     @param value Value being put.
   */
  virtual int    Put(const int row,const int col,const REAL & value );
  /** 
   * Get value with bound checking 
     @param row Row number.
     @param col Column number.
   */
  virtual const REAL &Get(const int row,const int col ) const;

  /** 
   * Substitution for the () operator when const arguments are needed 
     @param row Row number.
     @param col Column number.
   */
  const REAL &g(const int row, const int col) const {return Get(row,col);}

  /**
   * The operators check on the bounds if the DEBUG variable is defined
     @param row Row number.
     @param col Column number.
   */
  REAL &operator() (const int row,const int col );
  /** 
   * The operators check on the bounds if the DEBUG variable is defined
     @param row Row number.
     @param col Column number.
   */
  virtual REAL &s(const int row, const int col);
  /**
   * The operators check on the bounds if the DEBUG variable is defined
     @param row Row number.
   */
  REAL &operator()(const int row);

  /** Put values without bounds checking \n
   *  this method is faster than "Put" if DEBUG is defined.
   */
  virtual int PutVal(const int /*row*/,const int /*col*/,const REAL & /*val*/ ) { return 0; }
  /** Get values without bounds checking \n
   *  this method is faster than "Get" if DEBUG is defined.
   */
  virtual const REAL &GetVal(const int /*row*/, const int /*col*/ ) const  { return gZero; }

  /** @name Algebraic
   *  Implements algebraic operations with matrices
   */
  //@{

  /**
   * It mutiplies itself by TPZMatrix A putting the result in res
   * @param A TPZMatrix object to multiplied to
   * @param res TPZFMatrix containing the result
   * @param opt Indicates if is Transpose or not
   * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
   */
  virtual void Multiply(const TPZFMatrix & A,TPZFMatrix & res,const int opt = 0,const int stride = 1) const;
  /**
   * It adds itself to TPZMatrix A putting the result in res
   * @param A TPZMatrix to added to current matrix
   * @param res Contains the result
   */
  virtual void Add(const TPZMatrix & A,TPZMatrix & res) const;
  /**
   * It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
   * @param x Is x on the above operation
   * @param y Is y on the above operation
   * @param z Is z on the above operation
   * @param alpha Is alpha on the above operation
   * @param beta Is beta on the above operation
   * @param opt Indicates if is Transpose or not
   * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
   */
  virtual void MultAdd(const TPZFMatrix & x,const TPZFMatrix & y, TPZFMatrix & z,
		       const REAL alpha=1., const REAL beta = 0., const int opt = 0, const int stride = 1 ) const;

  /**
   * Computes res = rhs - this * x 
   */
  virtual void Residual(const TPZFMatrix & x,const TPZFMatrix & rhs, TPZFMatrix & res ) const;
  /**
   * It substracts A from storing the result in result
   */
  virtual void Substract(const TPZMatrix & A,TPZMatrix & result) const;


  /**Converts the matrix in an identity matrix*/
  virtual void Identity();
  /** 
   * It makes *T the transpose of current matrix.
   */ 
  virtual void Transpose(TPZMatrix *const T) const;
  //@}

  /**
   * @name Generic
   * Performs generic operations on matrices
   */

  //@{
  /**Input operation*/
  virtual void Input(istream & in = cin );

  /**Input operation*/
  friend istream & operator>>(istream& in,TPZMatrix & A);

  /** It prints the matrix data in a MatrixFormat Rows X Cols */
  virtual void Print(const char *name = 0, ostream &out = cout ,const MatrixOutputFormat form = EFormatted) const;

  /**Returns number of rows */
  int Rows() const;
  /** Returns number of cols */
  int Cols() const;

  /**
   * Returns the dimension of the matrix if the matrix is square. \n
   * If the matrix is not square, returns an error
   */
  inline virtual int Dim() const;


  /** 
   * Redimensions a matriz keeping the previous values 
   * @param newRow Specifies the new number of rows in matrix
   * @param newCol Specifies the new number of Columns in matrix
   */
  virtual int Resize(const int newRows, const int newCols ) { 
    fRow = newRows;
    fCol = newCols;
    return 0; 
  }

  /** 
   * Redimensions the matrix reinitializing it with zero 
   * @param newRow Specifies the new number of rows in matrix.
   * @param newCol Specifies the new number of Columns in matrix.
   */
  virtual int Redim(const int newRows, const int newCols ) { 
    fRow = newRows;
    fCol = newCols;
    return 0; 
  }

  /** Zeroes the matrix */
  virtual int Zero(){
    cout << "WARNING! TPZMatrix::Zero is called\n";
    return 0; }

  //@}
  /**
   * @name SubMatrices
   * Operations with SUB MATRICES.
   */

  //@{
  /**
   * It puts submatrix Source on actual matrix structure.
   * @param sRow Specifies starting row on current object
   * @param sCol Specifies starting column on current object.
   * @param Source The matrix to be inserted
   */
  virtual int PutSub( const int sRow, const int sCol, const TPZFMatrix & Source );

  /**
   * Gets submatrix storing it on Target.
   * @param sRow Specifies starting row on current object
   * @param sCol Specifies starting column on current object.
   * @param rowSize Specifies the amount of rows from sRow
   * @param colSize Specifies the amount of columns from sCol
   * @param Target The matrix to be aquired.
   */
  virtual int GetSub( const int sRow, const int sCol, const int rowSize,
		      const int colSize, TPZFMatrix & Target ) const;

  /**
   * It adds Source matrix on current matrix from position (sRow, sCol)
   * @param sRow Specifies starting row on current object
   * @param sCol Specifies starting column on current object.
   * @param Source The matrix to be added
   */
  virtual int AddSub(const int sRow, const int sCol, const TPZFMatrix & Source );

  //insere uma submatriz do objeto corrente para *Target sem redimencionar-la
  /**
   * Inserts a submatrix from current object on matrix *Target with no \n
   * redimentioning
   * @param sRow Specifies starting row on current object
   * @param sCol Specifies starting column on current object.
   * @param rowSize Specifies the amount of rows from sRow
   * @param colSize Specifies the amount of columns from sCol
   * @param Target The matrix to be inserted.
   */
  virtual int InsertSub(const int sRow,const int sCol,const int rowSize,const int colSize,
			const int pRow,const int pCol, TPZMatrix * Target ) const;

  //adiciona uma submatriz do objeto corrente para *target
  //10-10-95
  /**
   * Adds a submatrix from current object in *Target
   * @param sRow Specifies starting row on current object
   * @param sCol Specifies starting column on current object.
   * @param rowSize Specifies the amount of rows from sRow
   * @param colSize Specifies the amount of columns from sCol
   * @param pRow Specifies starting row on pA
   * @param pCol Specifies starting column on pA.
   * @param pA The matrix to be added.
   */
  virtual int AddSub(const int sRow,const int sCol,const int rowSize,
		     const int colSize,const int pRow,const int pCol, TPZMatrix * pA ) const;


  //@}

  /**
   * Add a contribution of a stiffness matrix
   * @param elmat Element matrix to be contributed
   * @param destinationindex Contains destine indexes on current matrix
   */
  virtual  void AddKel(TPZFMatrix &elmat, TPZVec<int> &destinationindex);

  /**
   * Add a contribution of a stiffness matrix
   * @param elmat Element matrix to be contributed
   * @param sourceindex Contains source indexes on current matrix
   * @param destinationindex Contains destine indexes on current matrix
   */
  virtual  void AddKel(TPZFMatrix &elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex);

  /**
   * @name Inquire
   * Returns information of the current object
   */
  //@{

  /**
   * Checks if the current matrix is symmetric
   */
  
  virtual int IsSimetric() const    { return 0; }
  /**
   * Checks if current matrix is square
   */
  inline int IsSquare() const { return fRow == fCol;}

  /**
   * Checks if current matrix is definite positive
   */
  virtual int IsDefPositive() const{ return 0; }
  /**
   * Checks if current matrix is already decomposed
   */
  int IsDecomposed() const         { return fDecomposed; }
  //@}
  /**
   * Sets current matrix to decomposed state
   */
  void SetIsDecomposed(int val) {fDecomposed = (char) val; }

  
  /**
   * @name Solvers
   * Linear system solvers. \n
   * For symmetric decompositions lower triangular matrix is used. \n
   * Solves a system A*X = B returning X in B
   */  
  //@{
  /**
   * Solves the linear system using Jacobi method. \n
   * @param numinterations The number of interations for the process.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param scratch Available manipulation area on memory.
   * @param tol The tolerance value.
   * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
   */
  virtual void SolveJacobi(int & numiterations, const TPZFMatrix & F, TPZFMatrix & result,
			   TPZFMatrix * residual, TPZFMatrix & scratch, REAL & tol, const int FromCurrent = 0) const;

  /**
   * Solves the linear system using Successive Over Relaxation method (Gauss Seidel). \n
   * @param numinterations The number of interations for the process.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param scratch Available manipulation area on memory.
   * @param overrelax The over relaxation parameter  
   * @param tol The tolerance value..
   * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
   * @param direction Indicates interaction direction, from first to last (default 1) or from last to first (-1)
   */
  virtual void SolveSOR(int & numiterations, const TPZFMatrix & F, TPZFMatrix & result,
			TPZFMatrix * residual,TPZFMatrix & scratch,const REAL overrelax, REAL & tol,
			const int FromCurrent = 0,const int direction = 1) const;
  /**
   * Solves the linear system using Symmetric Successive Over Relaxation method (Gauss Seidel). \n
   * @param numinterations The number of interations for the process.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param scratch Available manipulation area on memory.
   * @param overrelax The over relaxation parameter   
   * @param tol The tolerance value..
   * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
   */
  virtual void SolveSSOR(int & numiterations,const TPZFMatrix & F, TPZFMatrix & result,
			 TPZFMatrix * residual, TPZFMatrix & scratch, const REAL overrelax, REAL & tol,
			 const int FromCurrent = 0) const;

  /**
   * Solves the linear system using Conjugate Gradient method. \n
   * @param numinterations The number of interations for the process.
   * @param preconditioner The preconditioner attribute used.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param tol The tolerance value.
   */
  virtual void SolveCG(int & numiterations, TPZSolver & preconditioner,
		       const TPZFMatrix & F, TPZFMatrix & result,
		       TPZFMatrix * residual, REAL & tol,
		       const int FromCurrent = 0) const;
  /**
   * Solves the linear system using Bi-Conjugate Gradient method. \n
   * @param numinterations The number of interations for the process.
   * @param preconditioner The preconditioner attribute used.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param tol The tolerance value.
   */
  virtual void SolveBICG(int & numiterations, TPZSolver & preconditioner,
		       const TPZFMatrix & F, TPZFMatrix & result,
		       REAL & tol) const;

  /**
   * Solves the linear system using Generalized Minimal Residual (GMRES) method. \n
   * @param numinterations The number of interations for the process.
   * @param preconditioner The preconditioner attribute used.
   * @param H The right hand side of the system
   * @param numvectors The number of vectors involved
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param tol The tolerance value.
   * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
   */
  virtual void SolveGMRES(int & numiterations, TPZSolver & preconditioner,
			  TPZFMatrix & H, int & numvectors,
			  const TPZFMatrix & F, TPZFMatrix & result,
			  TPZFMatrix * residual, REAL & tol,const int FromCurrent) const;

  /**
   * Solves the linear system using IR method. \n
   * @param numinterations The number of interations for the process.
   * @param preconditioner The preconditioner attribute used.
   * @param F The right hand side of the system.
   * @param result The solution.
   * @param residual Returns F - A*U which is the solution residual.
   * @param tol The tolerance value.
   * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
   */
  virtual void SolveIR(int & numiterations, TPZSolver & preconditioner,
		       const TPZFMatrix & F, TPZFMatrix & result,
		       TPZFMatrix * residual, REAL & tol,
		       const int FromCurrent = 0) const;

  /**
   * Solves the linear system using Direct methods\n
   * @param F The right hand side of the system and where the solution is stored.
   * @param DecomposeType Indicates type of decomposition
   */
  virtual int SolveDirect ( TPZFMatrix & F , const DecomposeType dt);


  /**
   * Solves the linear system using LU method\n
   * @param B The right hand side of the system and where the solution is stored.
   */
  int Solve_LU ( TPZFMatrix * B );

  //para usar estos solves e' responsabilidade do usuario
  //que a mtriz corrente seja simetrica
  /**
   * Solves the linear system using Cholesky method\n
   * @param B The right hand side of the system and where the solution is stored.
   */  
  int Solve_Cholesky( TPZFMatrix * B );
  /**
   * Solves the linear system using LDLt method\n
   * @param B The right hand side of the system and where the solution is stored.
   */  
  int Solve_LDLt    ( TPZFMatrix * B );

  //@}
  /**
   * @name Factorization
   * Those member functions are related to matrices factorization
   */
  //@{
  /**
   * Decomposes the current matrix using LU decomposition.
   */
  virtual int Decompose_LU();

  // Decompoe a matriz em GGt (onde G e' triangular inferior).
  //a matriz corrente tem que ser simetrica.
  /**
   * Decomposes the current matrix using Cholesky method. \n
   * The current matrix has to be symmetric.
   */
  virtual int Decompose_Cholesky() ;

  // Decompoe a matriz em LDLt (onde L e' triangular inferior com
  //  1.0 na diagonal, e D e' uma matriz diagonal).
  //a matriz a ser decomposta tem que ser simetrica
  /**
   * Decomposes the current matrix using LDLt. \n
   * The current matrix has to be symmetric.
   * "L" is lower triangular with 1.0 in its diagonal and "D" is a Diagonal matrix.
   */
  virtual int Decompose_LDLt();
  //@}
 
  /**
   * @name Substitutions
   * Substitutions forward and backward
   */
  //@{  
  /**
   * Computes Forward and Backward substitution for a "LU" decomposed matrix.
   * @param B right hand side and result after all
   */
  virtual int Substitution( TPZFMatrix * B ) const;

  /**
   * Computes B = Y, where A*Y = B, A is lower triangular.
   * @param b right hand side and result after all
   */
  virtual int Subst_Forward( TPZFMatrix * b ) const;

  /**
   * Computes B = Y, where A*Y = B, A is upper triangular.
   * @param b right hand side and result after all
   */
  virtual int Subst_Backward( TPZFMatrix * b ) const;

  /** 
   * Computes B = Y, where A*Y = B, A is lower triangular with A(i,i)=1.
   * @param b right hand side and result after all
   */
  virtual int Subst_LForward( TPZFMatrix * b ) const;

  /**
   * Computes B = Y, where A*Y = B, A is upper triangular with A(i,i)=1.
   * @param b right hand side and result after all
   */
  virtual int Subst_LBackward( TPZFMatrix * b ) const;

  /**
   * Computes B = Y, where A*Y = B, A is diagonal matrix.
   * @param b right hand side and result after all
   */
  virtual int Subst_Diag( TPZFMatrix * b ) const;
  //@}
#ifdef OOPARLIB
  /**
   * @name Saveable
   * Methods which would make TPZMatrix compliant with TSaveable
   */
  //@{
  /**
   * Returns Class ID
   */
  virtual long GetClassID() const{ return TMATRIX_ID; }

  /**
   * Unpacks the object structure from a stream of bytes
   * @param buf The buffer containing the object in a packed form
   */
  virtual int Unpack( TReceiveStorage * buf );


  /**
   * Static method which will recompose an object of type TPZMatrix based on buf
   * @param buf Stream of bytes where the object is stored
   */
  static TSaveable *Restore(TReceiveStorage * buf);
  /**
   * Packs the object structure in a stream of bytes
   * @param buf Buffer which will receive the bytes
   */
  virtual int Pack( TSendStorage * buf ) const;

  /**
   * Returns class name
   */
  virtual char *ClassName() const{ return( "TPZMatrix" ); }

  /**
   * Returns true if the object belongs to a class which is derived from Classid
   * @param Classid Inquired base class ID
   */
  virtual int DerivedFrom(const long Classid) const;

  /**
   * Returns true if the object belongs to a class which is derived from classname
   * @param classname Inquired base class name
   */
  virtual int DerivedFrom(const char * classname) const; // a class with name classname
  //@}
#endif

protected:

  /**
   * Is an auxiliar method used by MultiplyAdd
   * @see MultAdd
   */
  void PrepareZ(const TPZFMatrix & y, TPZFMatrix & z,const REAL beta,const int opt,const int stride) const;
  
  /**
   * Constructor
   * @param row Number of rows
   * @param col Number of cols
   */
  TPZMatrix (const int row,const int col )
  { fRow = row; fCol = col;fDefPositive=0; fDecomposed = 0;}

  /**
   * Returns error messages
   * @param msg First message.
   * @param msg2 Second message.
   */
  virtual int Error(const char *msg ,const char *msg2 = 0) const;
  /**
   * It clears data structure.
   */
  virtual int Clear() { return 0; }

  /**
   * Swaps contents of a in b and b in a
   */
  static void Swap(int *a, int *b);
  /* inline void
    TPZMatrix::Swap( int *a, int *b ) const
    {
      int aux = *a;
      *a = *b;
      *b = aux;
    }
  */
  /** Number of rows in matrix */
  int fRow;
  /** Number of cols in matrix */
  int fCol;
  /** Decomposition type used to decompose the current matrix */
  char  fDecomposed;
  /** Definite Posistiveness of current matrix */
  char fDefPositive;
  static REAL gZero;
};


ostream & operator<<(ostream& out, const TPZMatrix & A);

/******** Inline ********/

inline int TPZMatrix::Rows() const {
         return fRow;
}

inline int TPZMatrix::Cols() const {
         return fCol;
}

inline void TPZMatrix::Residual(const TPZFMatrix & x,const TPZFMatrix & rhs, TPZFMatrix & res ) const {
         MultAdd( x, rhs, res, -1.0, 1.0 );
}

/***********/
/*** Put ***/

inline int TPZMatrix::Put(const int row,const int col,const REAL & value ) {
  // verificando se o elemento a inserir esta dentro da matriz
#ifdef DEBUG
         if ( row >= Rows() || col >= Cols() || row <0 || col < 0 ) {
         	Error("TPZMatrix::Put","Index out of range");
         }
#endif

         return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

inline const REAL &TPZMatrix::Get(const int row, const int col ) const {
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


inline REAL &TPZMatrix::operator()(const int row, const int col) {
         // verificando se o elemento a inserir esta dentro da matriz
#ifndef NODEBUG
        if ( (row >= Rows()) || (col >= Cols()) || row <0 || col<0 ) {
                gZero = 0.;
                Error("TPZMatrix::Operator()","Index out of range");
                return gZero;
        }
#endif
	return s(row,col);
}

inline REAL &TPZMatrix::s(const int row, const int col) {
         // verificando se o elemento a inserir esta dentro da matriz
	cout << "TPZMatrix::s not implemented\n";
	return gZero;
}

inline REAL &TPZMatrix::operator()(const int row) {
        return operator()(row,0);
}


inline int TPZMatrix::Dim() const{
        if ( IsSquare() ) return Rows();
        Error( "matrix is not square" );
        return ( 0 );
}
//***Solve LU ***/

inline int TPZMatrix::Solve_LU( TPZFMatrix *B) {
        if ( IsSimetric() ) Error( "LU decomposition is a not symetric decomposition" );
        return ( ( !Decompose_LU() )?  0 : Substitution( B )  );
}
/**********************/
/*** Solve Cholesky ***/
//
//  Se nao conseguir resolver por Cholesky retorna 0 e a matriz
//   sera' modificada (seu valor perdera' o sentido).
//
inline int TPZMatrix::Solve_Cholesky( TPZFMatrix * B ) {
        return(
            ( !Decompose_Cholesky() )?  0 :( Subst_Forward( B ) && Subst_Backward( B ) )
                          );
}



/******************/
/*** Solve LDLt ***/
inline int TPZMatrix::Solve_LDLt( TPZFMatrix * B ) {

        return(
              ( !Decompose_LDLt() )? 0 :
                       ( Subst_LForward( B ) && Subst_Diag( B ) && Subst_LBackward( B ) )
              );
}

inline void
TPZMatrix::Swap( int *a, int *b )
{
  int aux = *a;
  *a = *b;
  *b = aux;
}


#endif
// _TMATRIXH_


