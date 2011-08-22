/**
 * @file
 * @brief Contains TPZFMatrix class which implements full matrix.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tfullmat.hh
//
// Class:  TPZFMatrix
//
// Obs.:   Implements matrix classes (normais).
//
// Versao: 04 / 1996.
//
#ifndef _TMATRIXHH_
#include "pzmatrix.h"
#endif


#ifndef _TFULLMATRIXH_
#define _TFULLMATRIXH_


#include <iostream>
#include <memory.h>

#include <math.h>


#include "pzsave.h"
#include "pzmatrixid.h"


template <class T>
class TPZVec;

/** \addtogroup matrix
 * @{ */

/** @brief MACRO to get MAT(row,col) entry */
#define GETVAL(MAT,rows,row,col) MAT->fElem[((unsigned)col)*rows+row]
/** @brief MACRO to put value val into MAT(row,col) entry */
#define PUTVAL(MAT,rows,row,col,val) MAT->fElem[((unsigned)col)*rows+row]=val
/** @brief MACRO to get the entry of the vector (ptr[col*rows+row]) as matrix ( ptr(row,col) )*/
#define SELECTEL(ptr,rows,row,col) ptr[col*rows+row]

/**
 * @brief Full matrix class. \ref matrix "Matrix"
 * @note The full matrix class is special in that the data is stored column wise
 */
class TPZFMatrix : public TPZMatrix {
	
public:
	/**
	 * @brief Simple constructor
	 */
	TPZFMatrix () : TPZMatrix( 0, 0 ), fElem(0),fGiven(0),fSize(0) {}
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param buf Preallocated memory area which can be used by the matrix object
     @param size Size of the area pointed to by buf
	 */
	TPZFMatrix (const int rows ,const int columns, REAL * buf,const int size);
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param val Inital value fill all elements
	 */
	TPZFMatrix (const int rows ,const int columns,const REAL & val );
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
	 */
	inline  TPZFMatrix(const int rows ,const int columns = 1) : TPZMatrix(rows,columns), fElem(0),fGiven(0),fSize(0) {
		if(rows*columns) fElem = new REAL[rows*columns];
	}
	//@{
	/**
     * @brief Copy constructor
     * @param refmat Used as a model for current object
	 */
	TPZFMatrix (const TPZFMatrix & refmat);
	
	CLONEDEF(TPZFMatrix)
	
	TPZFMatrix(const TPZMatrix &refmat); // copy the elements one by one
	//@}
	/**
     @brief Constructor that uses a temporary matrix
	 */
	//  TPZFMatrix(TPZTempFMatrix);
	/**
     @brief Simple destructor
	 */
	virtual  ~TPZFMatrix();
	
	int PutVal(const int row,const int col,const REAL & value );
	const REAL &GetVal(const int row,const int col ) const;
	
	virtual REAL &s(const int row, const int col);
	
	REAL &g(const int row, const int col) const;
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix &rhs,TPZVec<int> &destination);
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param source Source index on rhs
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix &rhs,TPZVec<int> &source, TPZVec<int> &destination);
	
	/**
	 * B = this * X
	 * If opt = 1 then B = Transpose[this] * X
	 */
	void ConstMultiply(const TPZFMatrix & x,TPZFMatrix & B,const int opt = 0) const;
	
	virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						 const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 ) const ;
	
	static void MultAdd(const REAL *ptr, int rows, int cols, const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						const REAL alpha=1.,const REAL beta = 0.,const int opt = 0,const int stride = 1 );
	
	//@{
	/**
     @brief Generic operator with REAL type
	 */
	REAL &operator()(const int row,const int col);
	REAL &operator()(const int row);
	//@}
	/**
     @name FULL
     @brief Operations with FULL matrices
	 */
	//@{
	/**
     @brief Generic operator with FULL matrices
	 */
	virtual TPZFMatrix &operator= (const TPZFMatrix &A );
	//TPZFMatrix &operator= ( TPZTempFMatrix A);
	TPZFMatrix operator+  (const TPZFMatrix &A ) const;
	//TPZTempFMatrix operator+ (TPZTempFMatrix A);
	TPZFMatrix operator-  (const TPZFMatrix &A ) const;
	TPZFMatrix operator*  ( TPZFMatrix A ) const ;
	TPZFMatrix &operator+=(const TPZFMatrix &A );
	//    TPZFMatrix &operator+=(TPZTempFMatrix A );
	TPZFMatrix &operator-=(const TPZFMatrix &A );
	//@}
	/**
	 * @brief Performs an ZAXPY operation being *this += alpha * p
	 * @param alpha Being alpha on above opereation
	 * @param p Being p on above operation
	 */
	void ZAXPY(const REAL alpha, const TPZFMatrix &p);
	/**
	 * @brief Performs an operation *this = this * beta + z
	 * @param beta Being beta on above opereation
	 * @param z Being z on above operation
	 */
	void TimesBetaPlusZ(const REAL beta, const TPZFMatrix &z);
	
	/// Operations with matrices GENERICAS.
	/**
     @name Generics
     @brief Generic operators with matrices
	 */
	//@{
	/**
     @brief Generic operator with matrices
	 */
	TPZFMatrix &operator= (const TPZMatrix &A );
	//  TPZFMatrix operator+  (const TPZMatrix &A ) const;
	//  TPZFMatrix operator-  (const TPZMatrix &A ) const;
	//  TPZFMatrix operator*  (const TPZMatrix &A ) const;
	//  TPZFMatrix &operator+=(const TPZMatrix &A );
	//  TPZFMatrix &operator-=(const TPZMatrix &A );
	//@}
	// Operations with values NUMERICOS.
	
	/**
     @name Numerics
     @brief Numeric operations with matrices
	 */
	//@{
	/**
     Numeric operator with matrices
	 */
	TPZFMatrix &operator= (const REAL val );
	TPZFMatrix operator+  (const REAL val ) const;
	TPZFMatrix operator-  (const REAL val ) const;// { return operator+( -val ); }
	TPZFMatrix operator*  (const REAL val ) const;
	TPZFMatrix &operator+=(const REAL val );
	TPZFMatrix &operator-=(const REAL val )  { return operator+=( -val ); }
	TPZFMatrix &operator*=(const REAL val );
	
	TPZFMatrix operator-() const;// { return operator*( -1.0 ); }
	//@}
	//  void Input( istream & in = cin );
	
	/// Redimension a matrix, but maintain your elements.
	int Resize(const int newRows,const int wCols );
	
	/// Redimension the matrix doing nothing with the elements
	int SetSize(int newRows, int newCols);
	
	/// Remodel the shape of the  matrix, but keeping the same dimension.
	int Remodel(const int newRows,const int wCols );
	
	
	/// Redimension a matrix and ZERO your elements.
	int Redim(const int newRows,const int newCols );
	
	/// Makes Zero all the elements
	int Zero();
	
	/** 
	 * @brief This method implements a Gram Schimidt method. \n this = Orthog.TransfToOrthog
	 * @param Orthog [out] each column represents a vector orthogonalized with respect to the first vector (first column of *this). Vectors are normalized
	 * @param TransfToOrthog [out] is the basis change from *this to Orthog 
	 * @author Caju
	 * @since 2007
	 */
	void GramSchmidt(TPZFMatrix &Orthog, TPZFMatrix &TransfToOrthog);
	
	void DeterminantInverse(REAL &determinant, TPZFMatrix &inverse);
	
	void Transpose(TPZMatrix *const T) const;
	/**
     @see TPZMatrix::Transpose
	 */
	void Transpose();
	
	/*** @name Solve some systems ***/
	/** @{ */
	
	/** @brief LU Decomposition. Stores L and U matrices at the storage of the same matrix */
	virtual int Decompose_LU(std::list<int> &singular);
	virtual int Decompose_LU();
	
	static int Substitution(const REAL *ptr, int rows, TPZFMatrix *B);
	
	
	virtual int Substitution( TPZFMatrix *B ) const;
	
	/** @brief LU Decomposition using pivot
	 * @author Edimar Cesar Rylo
	 */
	virtual int Decompose_LU(TPZVec<int> &index);
	
	/** @brief LU substitution using pivot.
	 * @author Edimar Cesar Rylo
	 */
	virtual int Substitution( TPZFMatrix *B, TPZVec<int> &index ) const;
	
	/** @brief LU substitution using pivot. Static version.
	 * @author Edimar Cesar Rylo - Philippe Devloo
	 */
	static int Substitution(const REAL *ptr, int rows,  TPZFMatrix *B, TPZVec<int> &index );
	
	/**
	 * @}
	 */
	
	/// routines to send and receive messages
	virtual int ClassId() const;
	
	virtual void Read( TPZStream &buf, void *context );
	virtual void Write(TPZStream &buf, int withclassid );
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/// Compare the object for identity with the object pointed to, eventually copy the object
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * generate an interupt if the override flag is true and the objects are different
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	
	operator const REAL*() const { return fElem; }
	
	static void PrintStatic(const REAL *ptr, int rows, int cols, const char *name, std::ostream& out,const MatrixOutputFormat form);
	
	
private:
	
	static int Error(const char *msg1,const char *msg2=0 );
	int Clear();
	
	REAL *fElem;
	REAL *fGiven;
	int fSize;
};

/** @} */

inline TPZFMatrix::TPZFMatrix(const int rows,const int cols,REAL * buf,const int sz)
: TPZMatrix( rows, cols ), fElem(buf),fGiven(buf),fSize(sz) {
    int size = rows * cols;
	if(size == 0)
	{
		fElem = NULL;
	}
	else if(size > sz) 
	{
		fElem=new REAL[size];
#ifndef NODEBUG
		if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
	}
}

inline TPZFMatrix::TPZFMatrix(const int rows,const int cols,const REAL & val )
: TPZMatrix( rows, cols ), fElem(0), fGiven(0), fSize(0) {
	int size = rows * cols;
	if(!size) return;
	fElem=new REAL[size];
#ifdef DEBUG
	if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
	for(int i=0;i<size;i++) fElem[i] = val;
}

/** @brief Implements a scalar product val*A */
inline TPZFMatrix operator*(REAL val, const TPZFMatrix &A)
{
	return A*val;
}


/*******************************/
/*** Operator*( TPZFMatrix & ) ***/

inline TPZFMatrix TPZFMatrix::operator*( TPZFMatrix A ) const {
	if ( Cols() != A.Rows() )
		Error( "Operator* <matrixs with incompatible dimensions>" );
	TPZFMatrix res;
	res.Redim( Rows(), A.Cols() );
	MultAdd(A,A,res,1.,0.,0);
	return( res );
}

/**
 * @brief Non abstract class which implements full matrices with preallocated storage. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<int N>
class TPZFNMatrix : public TPZFMatrix {
	REAL fBuf[N+1];
	
public:
	/*
	 * Constructor which does not initialize the data
	 * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
	 * @param row Number of rows
	 * @param col Number of cols
	 */
	inline TPZFNMatrix(int row, int col) : TPZFMatrix(row,col,fBuf,N) {}
	
	inline TPZFNMatrix() : TPZFMatrix(0,0,fBuf,N)
	{
	}
	
	inline TPZFNMatrix(const TPZFMatrix &copy) : TPZFMatrix(0,0,fBuf,N)
	{
		*this = copy;
	}
	
	virtual ~TPZFNMatrix()
	{
	}
	
	CLONEDEF(TPZFNMatrix)
	/*
	 * Constructor which initializes the data
	 * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
	 * @param row Number of rows
	 * @param col Number of cols
	 * @param val initial value of the matrix elements
	 */
	inline  TPZFNMatrix(int row, int col, const REAL &val) : TPZFMatrix(row,col,fBuf,N) {
		TPZFMatrix::operator=(val);
	}
	
	inline  TPZFMatrix &operator=(const TPZFMatrix &copy) {
		return TPZFMatrix::operator=(copy);
	}
	inline  TPZFNMatrix<N> &operator=(const TPZFNMatrix<N> &copy) {
		TPZFMatrix::operator=(copy);
		return *this;
	}
	
};

/** @brief Returns a dot product to matrices */
REAL Dot(const TPZFMatrix &A,const TPZFMatrix &B);

/** @brief Returns the norm of the matrix A */
REAL Norm(const TPZFMatrix &A);





/**************/
/*** PutVal ***/
inline int TPZFMatrix::PutVal(const int row, const int col,const REAL & value ) {
	fElem[ ((unsigned)col) * Rows() + row ] = value;
	return( 1 );
}

/******************/
/*** Destructor ***/

inline TPZFMatrix::~TPZFMatrix () {
	if(fElem && fElem != fGiven) delete[]( fElem );
	fElem = 0;
	fSize = 0;
}


/**************/
/*** GetVal ***/
inline const REAL &TPZFMatrix::GetVal( const int row, const int col ) const {
#ifdef DEBUG
	if(row >= Rows() || row<0 || col >= Cols() || col<0) {
		Error("TPZFMatrix::operator() "," Index out of bounds");
		DebugStop();
		return gZero;
	}
#endif
	return( fElem[ ((unsigned)col) * Rows() + row ] );
}


inline REAL &TPZFMatrix::operator()( const int row, const int col) {
#ifndef NODEBUG
	if(row >= Rows() || row<0 || col >= Cols() || col<0) {
		Error("TPZFMatrix::operator() "," Index out of bounds");
		DebugStop();
		return gZero;
	}
#endif
	return *(fElem+col*fRow+row);
}

inline REAL &TPZFMatrix::s(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}


inline REAL &TPZFMatrix::g( const int row, const int col) const {
#ifdef DEBUG
	if(row >= Rows() || row<0 || col >= Cols() || col<0) {
		Error("TPZFMatrix::operator() "," Index out of bounds");
		DebugStop();
		return gZero;
	}
#endif
	return *(fElem+col*fRow+row);
}


inline REAL &TPZFMatrix::operator()(const int row) {
#ifdef DEBUG
	if(row >= Rows() || row<0) {
		Error("TPZFMatrix::operator() "," Index out of bounds");
		DebugStop();
		return gZero;
	}
#endif
	return *(fElem+row);
}

inline int TPZFMatrix::Redim(const int newRows,const int newCols) {
	int newsize = newRows*newCols;
	int size = fRow*fCol;
	if ( newsize == size) {
		fRow = newRows;
		fCol = newCols;
		Zero();
		return( 1 );
	}
	if(fElem && fElem != fGiven) delete []fElem;
	
	if(fGiven && newsize <= fSize) {
		fElem = fGiven;
	} else if(newsize == 0) {
		fElem = NULL;
	} else {
		fElem = new REAL[ newsize ] ;
	}
#ifndef NODEBUG
	if (newsize && fElem == NULL )
		Error( "Resize <memory allocation error>." );
#endif
	
	fRow  = newRows;
	fCol  = newCols;
	
	Zero();
	
	return( 1 );
}

/***************/
/****Zero*******/

inline int TPZFMatrix::Zero() {
	int size = fRow * fCol * sizeof(REAL);
	memset(fElem,'\0',size);
	fDecomposed = 0;
	return( 1 );
}


/**************************/
/*** Operations Global ***/



//inline TPZFMatrix &TPZFMatrix::operator+=(TPZTempFMatrix A) {
//	return (*this) += A.Object();
//}

inline REAL Norm(const TPZFMatrix &A) {
	return sqrt(Dot(A,A));
}



#endif


