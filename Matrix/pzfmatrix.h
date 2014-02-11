/**
 * @file
 * @brief Contains TPZMatrixclass which implements full matrix (using column major representation).
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
#include "pzmatrixid.h"


template <class T>
class TPZVec;

template <class TVar>
class TPZVerySparseMatrix;

/**
 * @addtogroup matrix
 * @{
 */

/** @brief MACRO to get MAT(row,col) entry */
#define GETVAL(MAT,rows,row,col) MAT->fElem[((unsigned)col)*rows+row]
/** @brief MACRO to put value val into MAT(row,col) entry */
#define PUTVAL(MAT,rows,row,col,val) MAT->fElem[((unsigned)col)*rows+row]=val
/** @brief MACRO to get the entry of the vector (ptr[col*rows+row]) as matrix ( ptr(row,col) )*/
#define SELECTEL(ptr,rows,row,col) ptr[col*rows+row]


/** @brief Returns a dot product to matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A,const TPZFMatrix<TVar> &B);

/** @brief Returns the norm of the matrix A */
template<class TVar>
TVar Norm(const TPZFMatrix<TVar> &A);


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
	TPZFMatrix(const long rows ,const long columns, TVar* buf,const long size);
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
     @param val Inital value fill all elements
	 */
	TPZFMatrix(const long rows ,const long columns,const TVar & val );
	/**
     @brief Constructor with initialization parameters
     @param rows Initial number of rows
     @param columns Number of columns
	 */
	inline  TPZFMatrix(const long rows ,const long columns = 1) : TPZMatrix<TVar>(rows,columns), fElem(0),fGiven(0),fSize(0) {
		if(rows*columns) fElem = new TVar[rows*columns];
	}
    
    /**
     * @brief Copy constructor specialized form TPZVerySparseMatrix
     * @param refmat Used as a model for current object
	 */
    TPZFMatrix(TPZVerySparseMatrix<TVar> const & A);
	
	/**
     * @brief Copy constructor
     * @param refmat Used as a model for current object
	 */
	TPZFMatrix(const TPZFMatrix<TVar> & refmat);
	
    
	
	CLONEDEF(TPZFMatrix<TVar>)
	TPZFMatrix(const TPZMatrix<TVar> & refmat);
    
	/** @brief Simple destructor */
	virtual  ~TPZFMatrix();
    
    long MemoryFootprint() const
    {
        return (sizeof(TVar)*this->Rows()*this->Cols());
	}
    
    TVar *Adress()
    {
        return fElem;
	}
	int PutVal(const long row,const long col,const TVar & value );
	const TVar &GetVal(const long row,const long col ) const;
	
	virtual TVar &s(const long row, const long col);
	
	TVar &g(const long row, const long col) const;
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<long> &destination);
	/**
	 * @brief Performs a right hand side assemblage
	 * @param rhs Load vector
	 * @param source Source index on rhs
	 * @param destination Destine index on current matrix
	 */
	void AddFel(TPZFMatrix<TVar> &rhs,TPZVec<long> &source, TPZVec<long> &destination);
	
	/**
	 * B = this * X
	 * If opt = 1 then B = Transpose[this] * X
	 */
	void ConstMultiply(const TPZFMatrix<TVar> & x,TPZFMatrix<TVar> & B,const int opt = 0) const;
	
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
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const ;
    
	
	static void MultAdd(const TVar *ptr, long rows, long cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 );
	
	/**
     * @name Generic operator with TVar type
	 * @{
	 */
	TVar &operator()(const long row,const long col);
	TVar &operator()(const long row);
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
	
    //	TPZFMatrix<TVar> operator-() const;// { return operator*( -1.0 ); }
    
	/** @} */
	
	/** @brief Redimension a matrix, but maintain your elements. */
	int Resize(const long newRows,const long wCols );
	
	/** @brief Redimension the matrix doing nothing with the elements */
	int SetSize(long newRows, long newCols);
	
	/** @brief Remodel the shape of the  matrix, but keeping the same dimension. */
	int Remodel(const long newRows,const long wCols );
    
	/** @brief Redimension a matrix and ZERO your elements. */
	int Redim(const long newRows,const long newCols );
	
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
	
	/** @brief Cholesky Decomposition Optmized. for walks in the direction of the vector that composes the matrix */
	virtual int Decompose_Cholesky();
	virtual int Decompose_Cholesky(std::list<long> &singular);
	
	/** @brief LU Decomposition. Stores L and U matrices at the storage of the same matrix */
	virtual int Decompose_LU(std::list<long> &singular);
	virtual int Decompose_LU();
	
	static int Substitution(const TVar *ptr, long rows, TPZFMatrix<TVar> *B);
	
	virtual int Substitution( TPZFMatrix<TVar> *B ) const;
	
	/** @brief LU Decomposition using pivot */
	virtual int Decompose_LU(TPZVec<long> &index);
	
	/** @brief LU substitution using pivot. */
	virtual int Substitution( TPZFMatrix<TVar> *B, TPZVec<long> &index ) const;
	
	/** @brief LU substitution using pivot. Static version. */
	static int Substitution(const TVar *ptr, long rows,  TPZFMatrix<TVar> *B, TPZVec<long> &index );
	
	/** @} */
	
	/** @brief Routines to send and receive messages */
	virtual int ClassId() const;
	
	virtual void Read( TPZStream &buf, void *context );
	virtual void Write(TPZStream &buf, int withclassid );
	virtual void Write(TPZStream &buf, int withclassid ) const;
	
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
	
	static void PrintStatic(const TVar *ptr, long rows, long cols, const char *name, std::ostream& out,const MatrixOutputFormat form);
    
private:
	
	static int Error(const char *msg1,const char *msg2=0 );
	int Clear();
	
	TVar *fElem;
	TVar *fGiven;
	long fSize;
};

/** @} */

template<class TVar>
inline TPZFMatrix<TVar>::TPZFMatrix(const long rows,const long cols,TVar * buf,const long sz)
: TPZMatrix<TVar>( rows, cols ), fElem(buf),fGiven(buf),fSize(sz) {
    long size = rows * cols;
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
inline TPZFMatrix<TVar>::TPZFMatrix(const long rows,const long cols,const TVar & val )
: TPZMatrix<TVar>( rows, cols ), fElem(0), fGiven(0), fSize(0) {
	long size = rows * cols;
	if(!size) return;
	fElem=new TVar[size];
#ifdef DEBUG
	if ( fElem == NULL && size) Error( "Constructor <memory allocation error>." );
#endif
	for(long i=0;i<size;i++) fElem[i] = val;
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

/**************/
/*** PutVal ***/
template<class TVar>
inline int TPZFMatrix<TVar>::PutVal(const long row, const long col,const TVar & value ) {
 	fElem[ ((unsigned)col) * this->Rows() + row ] = value;
	return( 1 );
}

/******************/
/*** Destructor ***/
template<class TVar>
inline TPZFMatrix<TVar>::~TPZFMatrix() {
	if(fElem && fElem != fGiven) delete[]( fElem );
	fElem = 0;
	fSize = 0;
}


/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &TPZFMatrix<TVar>::GetVal( const long row, const long col ) const {
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
inline TVar &TPZFMatrix<TVar>::operator()( const long row, const long col) {
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
inline TVar &TPZFMatrix<TVar>::s(const long row, const long col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

template<class TVar>
inline TVar &TPZFMatrix<TVar>::g( const long row, const long col) const {
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
inline TVar &TPZFMatrix<TVar>::operator()(const long row) {
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
inline int TPZFMatrix<TVar>::Redim(const long newRows,const long newCols) {
	long newsize = newRows*newCols;
	long size = this->fRow*this->fCol;
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
	long size = this->fRow * this->fCol * sizeof(TVar);
	memset(this->fElem,'\0',size);
	this->fDecomposed = 0;
	return( 1 );
}

/**************************/
/*** Operations Global ****/

inline long Norm(const TPZFMatrix<long> &A) {
	return (long)sqrt((REAL)Dot(A,A));
}

inline int Norm(const TPZFMatrix<int> &A) {
	return (int)sqrt((REAL)Dot(A,A));
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

inline float Norm(const TPZFMatrix< std::complex <float> > &A) {
	return sqrt(Dot(A,A).real());
}

inline double Norm(const TPZFMatrix< std::complex <double> > &A) {
	return sqrt(Dot(A,A).real());
}

inline long double Norm(const TPZFMatrix< std::complex <long double> > &A) {
	return sqrt(Dot(A,A).real());
}

inline TPZFlopCounter Norm(const TPZFMatrix<TPZFlopCounter> &A)
{
    return sqrt(Dot(A, A));
}
/**
 * @brief Non abstract class which implements full matrices with preallocated storage with (N+1) entries. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<int N, class TVar=REAL>
class TPZFNMatrix : public TPZFMatrix<TVar> {
    
	TVar fBuf[N+1];
	
public:
	/*
	 * @brief Constructor which does not initialize the data. \n
	 * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
	 * @param row Number of rows
	 * @param col Number of cols
	 */
	inline TPZFNMatrix(long row, long col) : TPZFMatrix<TVar>(row,col,fBuf,N) {}
	
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
	 * @brief Constructor which initializes the data. \n
	 * WARNING : this class will dynamically allocate memory if the template parameter N is smaller than row*col
	 * @param row Number of rows
	 * @param col Number of cols
	 * @param val initial value of the matrix elements
	 */
	inline  TPZFNMatrix(long row, long col, const TVar &val) : TPZFMatrix<TVar>(row,col,fBuf,N) {
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
    
    
#endif
