/**
 * @file
 * @brief Contains TPZSpMatrix class which defines sparse matrix class.
 */

#ifndef TESPMATH
#define TESPMATH

#include "pzlink.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef WORKPOOL

#include "pzworkpool.h"

#endif

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

template<class TVar>
class TPZSSpMatrix;

/**
 * @brief Defines sparce matrix class. \ref matrix "Matrix"
 * @ingroup matrix
 * @author Misael Santana Mandujano
 * @since 04/1996
 */
/** Stores data as a linked list of nonzero elements */
template<class TVar>
class TPZSpMatrix : public TPZMatrix<TVar>
{
	friend class TPZSSpMatrix<TVar>;
	
public:
	
	/**
     @struct TPZNode
     @brief Defines a node
     @param elem Element on node
     @param col Column number
	 */
	struct TPZNode
	{
		REAL elem;
		int64_t    col;
	};
	
public:
	/** @brief Simple constructor */
	TPZSpMatrix ()
    : TPZRegisterClassId(&TPZSpMatrix::ClassId),
    TPZMatrix<TVar>( 0, 0 )  { fElem = NULL; }
	/**
     * @brief Constructor with initialization parameters
     * @param rows Number of rows
     * @param cols Number of columns
	 */
	TPZSpMatrix ( int64_t rows, int64_t cols );
	/**
     * @brief Copy constructor
     * @param A Model object
	 */
	TPZSpMatrix (const TPZSpMatrix<TVar> &A )
    : TPZRegisterClassId(&TPZSpMatrix::ClassId),
    TPZMatrix<TVar>( A )  { fCopy( &A ); }
	
	CLONEDEF(TPZSpMatrix)
	/** @brief Simple destructor */
	~TPZSpMatrix();
	
	int    Put(const int64_t row,const int64_t col,const TVar & value ) override;
	const TVar &Get(const int64_t row,const int64_t col ) const override;
	
	// Nao verifica limites da matriz (e' mais rapido).
	int    PutVal(const int64_t row,const int64_t col,const TVar & element ) override;
	const TVar &GetVal(const int64_t row,const int64_t col ) const override;
	
	/**
	 * @name Operadores com matrizes ESPARSAS NAO simetricas.
	 * @{
	 */
	
	/** @brief Generic overloaded operator */ 
	TPZSpMatrix &operator= (const TPZSpMatrix<TVar> &A );
	TPZSpMatrix operator+  (const TPZSpMatrix<TVar> &A ) const;
	TPZSpMatrix operator-  (const TPZSpMatrix<TVar> &A ) const;
	TPZSpMatrix &operator+=(const TPZSpMatrix<TVar> &A );
	TPZSpMatrix &operator-=(const TPZSpMatrix<TVar> &A );
	
	/** @} */
	
	/**
	 * @name Operadores com matrizes GENERICAS.
	 * @{
	 */
	
	//TPZSpMatrix & operator= (const TPZMatrix<TVar> & A );
	
	/** @} */
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const REAL alpha ,const REAL beta = 0.,const int opt = 0) const override;
	
	
	/**
	 * @name Operadores com valores NUMERICOS.
	 * @{
	 */
	
	/** @brief Numerical values operator */
	TPZSpMatrix operator*  (const TVar v ) const;
	TPZSpMatrix &operator*=(const TVar v );
	
	TPZSpMatrix operator-() const { return operator*(-1.0); }
	
	/** @} */
	
	/// Desallocate all the elements of the matrix
	TPZSpMatrix &Reset();
	
	/// Redimension the matrix to newRows x newCols arrange.
	int Resize(const int64_t newRows,const int64_t newCols ) override;
	/**
     * @brief Redimensions current matrix keeping its elements
     * @param newDim New matrix dimensio
	 */
	int Resize(const int64_t newDim )   { return Resize( newDim, newDim ); }
	
	/// Redimension the matrix to newRows x newCols arrange and zeroes its elements
	int Redim(const int64_t newRows,const int64_t newCols ) override;
	/**
     * @brief Redimensions matrix deleting all its elements
     * @param newDim New matrix dimension
	 */
	int Redim(const int64_t newDim )    { return Redim( newDim, newDim ); }
	
	/// Zeroes all elements of the matrix
	int Zero() override;
	
	/*** Resolucao de sistemas ***/
	
#ifdef OOPARLIB
	
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *CreateInstance(TReceiveStorage *buf);
	inline virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSpMatrix" ); }
	virtual int DerivedFrom(const int64_t Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif
	public:
int ClassId() const override;

protected:
	
	/**  @see TPZSpMatrix::operator+= */
	
	int fAdd(const TPZSpMatrix<TVar> *const A );               // operator +=
//	int fAdd(const TPZSpMatrix<TVar> *A, TPZSpMatrix<TVar> *B ); // operator +
	
	int fSub(const TPZSpMatrix<TVar> *const A );               // operator -=
//	int fSub( TPZSpMatrix<TVar> *A, TPZSpMatrix<TVar> *B ); // operator -
	
	int fCopy(const TPZSpMatrix<TVar> *const A );              // operator =
	
	int fMult( TVar val );                // operator *
	
	/** Swap (troca) the values of the variables */
	inline void Swap( int64_t *a, int64_t *b );
	
	int  Clear() override { delete [] fElem; return ( 1 );}
	
	/**
     * @brief Computes dot product with respect to lines row_i e row_j using only elements that belongs to columns less than k.
	 * @note Afterwards lines row_i and row_j would be positioned on the element of column k. (where it should be).
     * @param row_i Ith row to be used
     * @param row_j Jth row to be used
	 * @param k number of column limit
	 */
	REAL ProdEsc( TPZLink<TPZNode> *row_i, TPZLink<TPZNode> *row_j,
				 int64_t k );
	
	
	TPZLink<TPZNode> *fElem;
	
#ifdef WORKPOOL
	TPZWorkPool fWp;
#endif
};

/*** Swap ***/
template<class TVar>
inline void
TPZSpMatrix<TVar>::Swap( int64_t *a, int64_t *b )
{
	int64_t c = *a;
	*a = *b;
	*b = c;
}

#include "pzsespmat.h"

#endif
