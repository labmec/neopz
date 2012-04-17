/**
 * @file
 * @brief Contains TPZMatrixclass which implements full matrix.
 */

#ifndef _TMATRIXHH_
#include "pzmatrix.h"
#endif


#ifndef _TFULLMATRIXH_
#define _TFULLMATRIXH_


#include <iostream>
#include <memory.h>

#include <math.h>

#include "pzsave.h"

template <class T>
class TPZVec;

/** 
 * @addtogroup matrix
 * @{ 
 */

/** @brief Id of full matrix */
const int TPZFMATRIXID = 100;

/** @brief MACRO to get MAT(row,col) entry */
#define GETVAL(MAT,rows,row,col) MAT->fElem[((unsigned)col)*rows+row]
/** @brief MACRO to put value val into MAT(row,col) entry */
#define PUTVAL(MAT,rows,row,col,val) MAT->fElem[((unsigned)col)*rows+row]=val
/** @brief MACRO to get the entry of the vector (ptr[col*rows+row]) as matrix ( ptr(row,col) )*/
#define SELECTEL(ptr,rows,row,col) ptr[col*rows+row]

/**
 * @brief Full matrix class. \ref matrix "Matrix"
 * @note The full matrix class is special in that the data is stored column wise
 * @author Misael Mandujano
 * @since 04/1996
 */

template<class TVar=REAL>
class TPZFMatrix: public TPZMatrix<TVar> {
	
public:
	/** @brief Simple constructor */
	TPZFMatrix() : TPZMatrix<TVar>( 0, 0 ), fElem(0),fGiven(0),fSize(0) {}
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param buf Preallocated memory area which can be used by the matrix object
     @param size Size of the area pointed to by buf
	 */
	TPZFMatrix(const int rows ,const int columns, TVar* buf,const int size);
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param val Inital value fill all elements
	 */
	TPZFMatrix(const int rows ,const int columns,const TVar & val );
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
	 */
	inline  TPZFMatrix(const int rows ,const int columns = 1) : TPZMatrix<TVar>(rows,columns), fElem(0),fGiven(0),fSize(0) {
		if(rows*columns) fElem = new TVar[rows*columns];
	}

	/**
     * @brief Copy constructor
     * @param refmat Used as a model for current object
	 */
	TPZFMatrix(const TPZFMatrix<TVar> & refmat);
	
	CLONEDEF(TPZFMatrix<TVar>)
	TPZFMatrix(const TPZMatrix<TVar> & refmat);

	/** @brief Simple destructor */
	virtual  ~TPZFMatrix();
	
	int PutVal(const int row,const int col,const TVar & value );
	const TVar &GetVal(const int row,const int col ) const;
	
	virtual TVar &s(const int row, const int col);
	
	TVar &g(const int row, const int col) const;
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int> &destination);
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param source Source index on rhs
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int> &source, TPZVec<int> &destination);
	
	/**
	 * B = this * X
	 * If opt = 1 then B = Transpose[this] * X
	 */
	void ConstMultiply(const TPZFMatrix<TVar> & x,TPZFMatrix<TVar> & B,const int opt = 0) const;
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const ;
	
	static void MultAdd(const TVar *ptr, int rows, int cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 );  
	
	/**
     * @name Generic operator with TVar type
	 * @{
	 */
	TVar &operator()(const int row,const int col);
	TVar &operator()(const int row);
	/** @} */
	
	/**
     * @name Operations with FULL matrices
	 * @{
	 */
	
	/** @brief Generic operator with FULL matrices */
	virtual TPZFMatrix&operator= (const TPZFMatrix<TVar> &A );
	TPZFMatrix<TVar> operator+  (const TPZFMatrix<TVar> &A ) const;
	TPZFMatrix<TVar> operator-  (const TPZFMatrix<TVar> &A ) const;
	TPZFMatrix<TVar> operator*  ( TPZFMatrix<TVar> A ) const ;
	TPZFMatrix<TVar> &operator+=(const TPZFMatrix<TVar> &A );
	TPZFMatrix<TVar> &operator-=(const TPZFMatrix<TVar> &A );
	
	/** @} */
	
	/**
	 * @brief Performs an ZAXPY operation being *this += alpha * p
	 * @param alpha Being alpha on above opereation
	 * @param p Being p on above operation
	 */
	void ZAXPY(const TVar alpha, const TPZFMatrix<TVar> &p);
	/**
	 * @brief Performs an operation *this = this * beta + z
	 * @param beta Being beta on above opereation
	 * @param z Being z on above operation
	 */
	void TimesBetaPlusZ(const TVar beta, const TPZFMatrix<TVar> &z);
	
	/** @brief Generic operator with matrices */
	TPZFMatrix<TVar> &operator= (const TPZMatrix<TVar> &A );
	
	/**
     * @name Numerics
     * @brief Numeric operations with matrices
	 * @{
	 */
	
	/** @brief Numeric operator with matrices */
	TPZFMatrix<TVar> &operator= (const TVar val );
	TPZFMatrix<TVar> operator+  (const TVar val ) const;
	TPZFMatrix<TVar> operator-  (const TVar val ) const;
	TPZFMatrix<TVar> operator*  (const TVar val ) const;
	TPZFMatrix<TVar> &operator+=(const TVar val );
	TPZFMatrix<TVar> &operator-=(const TVar val )  { return operator+=( -val ); }
	TPZFMatrix<TVar> &operator*=(const TVar val );
	
	TPZFMatrix<TVar> operator-() const;// { return operator*( -1.0 ); }

	/** @} */
	
	/** @brief Redimension a matrix, but maintain your elements. */
	int Resize(const int newRows,const int wCols );
	
	/** @brief Redimension the matrix doing nothing with the elements */
	int SetSize(int newRows, int newCols);
	
	/** @brief Remodel the shape of the  matrix, but keeping the same dimension. */
	int Remodel(const int newRows,const int wCols );

	/** @brief Redimension a matrix and ZERO your elements. */
	int Redim(const int newRows,const int newCols );
	
	/** @brief Makes Zero all the elements */
	int Zero();
	
	/** 
	 * @brief This method implements a Gram Schimidt method. \n this = Orthog.TransfToOrthog
	 * @param Orthog [out] each column represents a vector orthogonalized with respect to the first vector (first column of *this). Vectors are normalized
	 * @param TransfToOrthog [out] is the basis change from *this to Orthog 
	 * @author Caju
	 * @since 2007
	 */
	void GramSchmidt(TPZFMatrix<TVar> &Orthog, TPZFMatrix<TVar> &TransfToOrthog);
	
	void DeterminantInverse(TVar &determinant, TPZFMatrix<TVar> &inverse);
	
	void Transpose(TPZMatrix<TVar> *const T) const;
	
	/** @see TPZMatrix<TVar>::Transpose */
	
	void Transpose();
	
	/*** @name Solve some systems ***/
	/** @{ */
	
	/** @brief LU Decomposition. Stores L and U matrices at the storage of the same matrix */
	virtual int Decompose_LU(std::list<int> &singular);
	virtual int Decompose_LU();
	
	static int Substitution(const TVar *ptr, int rows, TPZFMatrix<TVar> *B);
	
	virtual int Substitution( TPZFMatrix<TVar> *B ) const;
	
	/** @brief LU Decomposition using pivot */
	virtual int Decompose_LU(TPZVec<int> &index);
	
	/** @brief LU substitution using pivot. */
	virtual int Substitution( TPZFMatrix<TVar> *B, TPZVec<int> &index ) const;
	
	/** @brief LU substitution using pivot. Static version. */
	static int Substitution(const TVar *ptr, int rows,  TPZFMatrix<TVar> *B, TPZVec<int> &index );
	
	/** @} */
	
	/** @brief Routines to send and receive messages */
	virtual int ClassId() const;
	
	virtual void Read( TPZStream &buf, void *context );
	virtual void Write(TPZStream &buf, int withclassid );
	
	/** @brief Compare the object for identity with the object pointed to, eventually copy the object */
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * overwrite the calling object if the override flag is true
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false);
	
	/** @brief Compare the object for identity with the object pointed to, eventually copy the object */
	/**
	 * compare both objects bitwise for identity. Put an entry in the log file if different
	 * generate an interupt if the override flag is true and the objects are different
	 */
	virtual bool Compare(TPZSaveable *copy, bool override = false) const;
	
	operator const TVar*() const { return fElem; }
	
	static void PrintStatic(const TVar *ptr, int rows, int cols, const char *name, std::ostream& out,const MatrixOutputFormat form);

private:
	
	static int Error(const char *msg1,const char *msg2=0 );
	int Clear();
	
	TVar *fElem;
	TVar *fGiven;
	int fSize;
};

/** @} */

template<class TVar>
inline TPZFMatrix<TVar>::TPZFMatrix(const int rows,const int cols,TVar * buf,const int sz)
: TPZMatrix<TVar>( rows, cols ), fElem(buf),fGiven(buf),fSize(sz) {
    int size = rows * cols;
	if(size == 0)
	{
		fElem = NULL;
	}
	else if(size > sz) 
	{
		fElem=new TVar[size];
#ifndef NODEBUG
		if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
	}
}
template<class TVar>
inline TPZFMatrix<TVar>::TPZFMatrix(const int rows,const int cols,const TVar & val )
: TPZMatrix<TVar>( rows, cols ), fElem(0), fGiven(0), fSize(0) {
	int size = rows * cols;
	if(!size) return;
	fElem=new TVar[size];
#ifdef DEBUG
	if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
	for(int i=0;i<size;i++) fElem[i] = val;
}

template<class TVar>
/** @brief Implements a scalar product val*A */
inline TPZFMatrix<TVar> operator*(TVar val, const TPZFMatrix<TVar> &A)
{
	return A*val;
}


/*******************************/
/*** Operator*( TPZMatrix<TVar> & ) ***/
template<class TVar>
inline TPZFMatrix<TVar> TPZFMatrix<TVar>::operator*( TPZFMatrix<TVar> A ) const {
	if ( this->Cols() != A.Rows() )
		Error( "Operator* <matrixs with incompatible dimensions>" );
	TPZFMatrix<TVar> res;
	res.Redim( this->Rows(), A.Cols() );
	MultAdd(A,A,res,1.,0.,0);
	return( res );
}

/**
 * @brief Non abstract class which implements full matrices with preallocated storage. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<int N, class TVar=REAL>
class TPZFNMatrix : public TPZFMatrix<TVar> {
	TVar fBuf[N+1];
	
public:
	/*
	 * Constructor which does not initialize the data
	 * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
	 * @param row Number of rows
	 * @param col Number of cols
	 */
	inline TPZFNMatrix(int row, int col) : TPZFMatrix<TVar>(row,col,fBuf,N) {}
	
	inline TPZFNMatrix() : TPZFMatrix<TVar>(0,0,fBuf,N)
	{
	}
	
	inline TPZFNMatrix(const TPZFMatrix<TVar> &copy) : TPZFMatrix<TVar>(0,0,fBuf,N)
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
	inline  TPZFNMatrix(int row, int col, const TVar &val) : TPZFMatrix<TVar>(row,col,fBuf,N) {
		TPZFMatrix<TVar>::operator=(val);
	}
	
	inline  TPZFMatrix<TVar> &operator=(const TPZFMatrix<TVar> &copy) {
		return TPZFMatrix<TVar>::operator=(copy);
	}
	inline  TPZFNMatrix<N, TVar> &operator=(const TPZFNMatrix<N, TVar> &copy) {
		TPZFMatrix<TVar>::operator=(copy);
		return *this;
	}
	
};

/** @brief Returns a dot product to matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A,const TPZFMatrix<TVar> &B);


/** @brief Returns the norm of the matrix A */
template<class TVar>
TVar Norm(const TPZFMatrix<TVar> &A);





/**************/
/*** PutVal ***/
template<class TVar>
inline int TPZFMatrix<TVar>::PutVal(const int row, const int col,const TVar & value ) {
	fElem[ ((unsigned)col) * this->Rows() + row ] = value;
	return( 1 );
}

/******************/
/*** Destructor ***/
template<class TVar>
inline TPZFMatrix<TVar>::~TPZFMatrix<TVar>() {
	if(fElem && fElem != fGiven) delete[]( fElem );
	fElem = 0;
	fSize = 0;
}


/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &TPZFMatrix<TVar>::GetVal( const int row, const int col ) const {
#ifdef DEBUG
	if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
		Error("TPZFMatrix::operator() "," Index out of bounds");
		DebugStop();
		return this->gZero;
	}
#endif
	return( fElem[ ((unsigned)col) *  this->Rows() + row ] );
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::operator()( const int row, const int col) {
#ifndef NODEBUG
	if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
		Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
		DebugStop();
		return this->gZero;
	}
#endif
	return *(this->fElem+col*this->fRow+row);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::s(const int row, const int col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::g( const int row, const int col) const {
#ifdef DEBUG
	if(row >=  this->Rows() || row<0 || col >=  this->Cols() || col<0) {
		Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
		DebugStop();
		return this->gZero;
	}
#endif
	return *(this->fElem+col*this->fRow+row);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::operator()(const int row) {
#ifdef DEBUG
	if(row >=  this->Rows() || row<0) {
		Error("TPZFMatrix<TVar>::operator() "," Index out of bounds");
		DebugStop();
		return this->gZero;
	}
#endif
	return *(this->fElem+row);
}

template<class TVar>
inline int TPZFMatrix<TVar>::Redim(const int newRows,const int newCols) {
	int newsize = newRows*newCols;
	int size = this->fRow*this->fCol;
	if ( newsize == size) {
		this->fRow = newRows;
		this->fCol = newCols;
		Zero();
		return( 1 );
	}
	if(this->fElem && this->fElem != this->fGiven) delete []this->fElem;
	
	if(this->fGiven && newsize <= this->fSize) {
		this->fElem = this->fGiven;
	} else if(newsize == 0) {
		this->fElem = NULL;
	} else {
		this->fElem = new TVar[ newsize ] ;
	}
#ifndef NODEBUG
	if (newsize && this->fElem == NULL )
		Error( "Resize <memory allocation error>." );
#endif
	
	this->fRow  = newRows;
	this->fCol  = newCols;
	
	Zero();
	
	return( 1 );
}

/***************/
/****Zero*******/

template<class TVar>
inline int TPZFMatrix<TVar>::Zero() {
	int size = this->fRow * this->fCol * sizeof(TVar);
	memset(this->fElem,'\0',size);
	this->fDecomposed = 0;
	return( 1 );
}

/**************************/
/*** Operations Global ***/

inline int Norm(const TPZFMatrix<int> &A) {
	return sqrt(Dot(A,A));
}

inline float Norm(const TPZFMatrix<float> &A) {
	return sqrt(Dot(A,A));
}

inline double Norm(const TPZFMatrix<double> &A) {
	return sqrt(Dot(A,A));
}

inline long double Norm(const TPZFMatrix<long double> &A) {
	return sqrt(Dot(A,A));
}

inline std::complex <float> Norm(const TPZFMatrix< std::complex <float> > &A) {
	return sqrt(Dot(A,A));
}

inline std::complex <double> Norm(const TPZFMatrix< std::complex <double> > &A) {
	return sqrt(Dot(A,A));
}

inline std::complex <long double> Norm(const TPZFMatrix< std::complex <long double> > &A) {
	return sqrt(Dot(A,A));
}

#endif


