/**
 * @file
 * @brief Contains the implementation of the TPZSBMatrixLapack methods.
 */

#include <math.h>
#include <stdlib.h>

#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "TPZSBMatrixLapack.h"

#ifdef USING_LAPACK
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#else
#include <clapack.h>
#endif
#endif

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzsbmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSBMatrixLapack ***/

/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

//final_ok
template<class TVar> TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack() : TPZMatrix<TVar>(0,0) , fDiag() {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  fBand = 0;
}
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack( long dim, long band )
: TPZMatrix<TVar>( dim, dim )
{
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
	fBand = ( band > (dim - 1) ? (dim - 1) : band );
	fDiag.Resize( Size() );
	if ( fDiag.begin() == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrixLapack( dim ) <Error creating Matrix>" );
	
	Zero();
}

//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack(const TPZSBMatrixLapack<TVar> &A ) : TPZMatrix<TVar>(A)  {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  Copy(A);
}


/**************/
/*** PutVal ***/
//final_ok
template <class TVar>
int
TPZSBMatrixLapack<TVar>::PutVal(const long r,const long c,const TVar& value )
{
  // initialises row and col in order to work with upper triangular matrix
  long row(r),col(c);
  if ( row > col )
    this->Swap( &row, &col );
  
  long index;
  if ( (index = col-row) > fBand )
  {
#ifdef PZDEBUG
    if (value != TVar(0.) ) {
      DebugStop();
    }
#endif
    return( 0 );        // Element out of bounds
  }
  fDiag[ fBand * ( col + 1) + row ] = value;
  this->fDecomposed = 0;
  return( 1 );
}

/**************/
/*** GetVal ***/
//final_ok
template<class TVar>
const TVar
&TPZSBMatrixLapack<TVar>::GetVal(const long r,const long c ) const
{
	
	// initialises row and col in order to work with upper triangular matrix
  long row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	long index;
	if ( (index = col-row) > fBand )
		return( this->gZero );        // Element out of bounds
	const long dim = this->fCol;//square matrix
	return( fDiag[ fBand * ( col + 1) + row ] );
}

/*************/
/*** Print ***/
//final_ok
template<class TVar>
void
TPZSBMatrixLapack<TVar> ::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	out.width( 8 );
	out.precision( 4 );
	
	out << "Writing matrix '" << name;
	out << "' (" << this->Rows() << " x " << this->Cols() << ")  Bandwith = "<<fBand<<"\n";
	TPZMatrix<TVar>::Print(0,out,form);
}

/** @brief Overload << operator to output entries of TPZSBMatrixLapack matrix ***/
//final_ok
template<class TVar>
std::ostream&
operator<<(std::ostream& out,TPZSBMatrixLapack<TVar>  &A)
{
	out.width( 8 );
	out.precision( 4 );
	
	out <<"\n(" << A.Rows() << " x " << A.Cols()
	<< ")  Bandwith = "<< A.GetBand()<<"\n";
	
	for ( long row = 0; row < A.Rows(); row++)
    {
		out << "\t";
		for ( long col = 0; col < A.Cols(); col++ )
			out << A.GetVal( row, col) << "  ";
		out << "\n";
    }
	
	return  out << "\n";
}

/******** Operacoes com matrizes BANDA SIMETRICA  ********/

/******************/
/*** Operator = ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator=(const TPZSBMatrixLapack<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}

/******************/
/*** Operator + ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator+(const TPZSBMatrixLapack<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-( TPZSBMatrixLapack ) <incompatible dimensions>" );
    
    long nel = Size();
    TPZSBMatrixLapack<TVar> res( this->Dim(), fBand );
    
    for (long el=0; el<nel; el++) {
        res.fDiag[el] = fDiag[el]+A.fDiag[el];
    }
    return res;
}

/******************/
/*** Operator - ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator-(const TPZSBMatrixLapack<TVar> &A ) const
{
  if ( this->Dim() != A.Dim() || fBand != A.fBand )
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-( TPZSBMatrixLapack ) <incompatible dimensions>" );
  
    long nel = Size();
    TPZSBMatrixLapack<TVar> res( this->Dim(), fBand );
    
    for (long el=0; el<nel; el++) {
        res.fDiag[el] = fDiag[el]-A.fDiag[el];
    }
    return res;
}

/*******************/
/*** Operator += ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator+=(const TPZSBMatrixLapack<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+=( TPZSBMatrixLapack ) <incompatible dimensions>" );
    
    long nel = Size();
    
    for (long el=0; el<nel; el++) {
        fDiag[el] = fDiag[el]+A.fDiag[el];
    }
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return *this;

}

/*******************/
/*** Operator -= ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator-=(const TPZSBMatrixLapack<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator-=( TPZSBMatrixLapack ) <incompatible dimensions>" );
    
    long nel = Size();
    
    for (long el=0; el<nel; el++) {
        fDiag[el] = fDiag[el]-A.fDiag[el];
    }
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return *this;
}
//NOK
template<class TVar>
void TPZSBMatrixLapack<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
								const TVar alpha,const TVar beta ,const int opt,const int stride ) const {
  DebugStop();
}

/******** Operacoes com valores NUMERICOS ********/
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//

/*****************************/
/*** Operator * ( REAL ) ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator*(const TVar value ) const
{
	TPZSBMatrixLapack<TVar> res( this->Dim(), fBand );
	
	TVar *pr  = res.fDiag.begin();
	TVar *pm  = fDiag.begin();
	TVar *end = fDiag.begin() + Size() * sizeof( TVar );
	while ( pm < end )
		*pr++ = (*pm++) * value;
	return( res );
}

/******************************/
/*** Operator += ( REAL ) ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator+=(const TVar value )
{
	TVar *pm  = fDiag.begin();
	TVar *end = fDiag.begin() + Size() * sizeof( TVar );
	for(long iCol = 0 ; iCol < this->Dim() ; iCol ++){
		long nElemCol = iCol >= fBand ? fBand + 1 : 1 + iCol;
		for (int i = 0 ; i < fBand + 1 - nElemCol; i++) {
			pm++;//skip zeros,if any
		}
		for(int i = 0 ; i < nElemCol ; i++){
			*pm++ += value;
		}
	}
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( *this );
}

/******************************/
/*** Operator *= ( REAL ) ***/
//final_ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator*=(const TVar value )
{
	TVar *pm  = fDiag.begin();
	TVar *end = fDiag.begin() + Size() * sizeof( TVar );
	while ( pm < end )
		*pm++ *= value;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( *this );
}

/**************/
/*** Resize ***/
//
// Changes matrix dimensions, but keeps its current values.
// New positions are created with zeros.
//
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Resize(const long newDim ,const long)
{
	if ( newDim == this->Dim() )
		return( 1 );
	
	// creates new matrix
	long newBand = fBand > newDim - 1 ? newDim - 1 : fBand;
	long newSize  = newDim * (newBand + 1);
	
	TVar *newDiag = new TVar[newSize] ;
	TVar *oldDiag = fDiag.begin();
	for (long iCol = 0; iCol <= newDim; iCol++) {
		for (int iSkip = 0; iSkip < fBand - newBand; iSkip++) {
			oldDiag++;
		}
		if( iCol < this->Dim() ){
			for (int iCopy = 0; iCopy < newBand + 1 - (fBand - newBand); iCopy++) {
				*newDiag++ = *oldDiag++;
			}
		}
		else{
			for(int iZero = 0 ; iZero < newBand + 1 ; iZero ++){
				*newDiag++ = this->gZero;
			}
		}
	}
	// Deletes old matrix and takes new one
	if ( fDiag.begin() != NULL )
		fDiag.clear();
	fDiag =  *newDiag;
	delete newDiag;
	this->fCol = this->fRow = newDim;
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}

/*************/
/*** Redim ***/
//
// Changes matrix dimwensions and erases all its elements
//
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Redim(const long newDim ,const long)
{
	if (newDim != this->Dim()){
		this->fRow = this->fCol = newDim;
		fDiag.Resize( newDim );
	}
	fDiag.Fill( this->gZero , 0 , Size() );
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Zero()
{
	TVar *dst = this->fDiag.begin();
	TVar *end = this->fDiag.begin() + Size() * sizeof( TVar );
	while ( dst < end )
		*dst++ = this->gZero;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}


/****************/
/*** Set Band ***/
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::SetBand(const long newBand )
{
	if ( this->fBand == newBand )
		return( 1 );
	
	if ( this->fBand > (this->Dim() - 1) )
		return( 0 );
	
	TVar *newDiag = new TVar[this->Dim() * (newBand + 1)];
	
	// Copia os elementos antigos para a nova alocacao.
	TVar * const pNew = newDiag;
	TVar * const pOld = fDiag.begin();
	long diffSize = abs( newBand - fBand );
	TVar *src = pOld;
	TVar *dst = pNew;
	
	if( newBand > fBand ){
		for (int iCol = 0 ; iCol < this->Dim() ; iCol ++) {
			for( int iZero = 0 ; iZero < diffSize ; iZero++) {
				*dst++ = this->gZero;
			}
			for(int iCopy = 0 ; iCopy < fBand + 1 ; iCopy++){
				*dst++ = *src++;
			}
		}
	}
	else{
		for (int iCol = 0 ; iCol < this->Dim() ; iCol ++) {
			for( int iSkip = 0 ; iSkip < diffSize ; iSkip++) {
				dst++;
			}
			for(int iCopy = 0 ; iCopy < newBand + 1 ; iCopy++){
				*dst++ = *src++;
			}
		}
	}
	fDiag.clear();
	fDiag = * newDiag;
	delete newDiag;
	fBand = newBand;
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}

/********************* Resolucao de sistemas *********************/

/**************************/
/*** Decompose Cholesky ***/
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky(std::list<long> &singular)
{
	return Decompose_Cholesky();
}

//final_ok
template<>
int
TPZSBMatrixLapack<std::complex< float > >::Decompose_Cholesky()
{
    if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<std::complex<float> >::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
    if (this->fDecomposed) {
        return 1;
    }
    
#ifdef USING_LAPACK
    char uplo = 'u';
    int n = this->Dim();
    int lda = this->Dim();
    int info = -666;
    cpotrf_(&uplo, &n, (__CLPK_complex*) fDiag.begin(), &lda, &info);
#endif
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}

//final_ok
template<>
int
TPZSBMatrixLapack<std::complex< double > >::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<std::complex<double> >::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;
	zpotrf_(&uplo, &n, (__CLPK_doublecomplex *) fDiag.begin(), &lda, &info);
#endif
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}
//final_ok
template<>
int
TPZSBMatrixLapack<double>::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;

	dpotrf_(&uplo, &n, fDiag.begin(), &lda, &info);

#endif
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}

//final_ok
template<>
int
TPZSBMatrixLapack<float>::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;
	
	spotrf_(&uplo, &n, fDiag.begin(), &lda, &info);

#endif
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky()
{
	TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <LAPACK does not support this specific data type>" );
	DebugStop();
	return 0;
}

/**********************/
/*** Decompose LDLt ***/
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_LDLt(std::list<long> &singular)
{
	return Decompose_LDLt();
}
//final_ok
template< >
int
TPZSBMatrixLapack< complex<double> >::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix< complex<double> >::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
#ifdef USING_LAPACK
  char uplo = 'u';
  int n = this->Dim();
  int lda = this->Dim();
  int info = -666;
  int *ipiv = new int [n];
  int lwork = this->Dim();
  complex<double> *work = new complex<double> [n];
	
  zsytrf_(&uplo, &n, ( __CLPK_doublecomplex * )fDiag.begin(), &lda, ipiv, ( __CLPK_doublecomplex * )work, &lwork, &info);
  
  if ( work != NULL ) delete[] work;
  if( ipiv != NULL ) delete[] ipiv;
  
#endif
  
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}

//final_ok
template< >
int
TPZSBMatrixLapack< complex<float> >::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix< complex<float> >::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;
	int *ipiv = new int [n];
	int lwork = this->Dim();
	complex<float> *work = new complex<float> [n];
	
	csytrf_(&uplo, &n, ( __CLPK_complex * )fDiag.begin(), &lda, ipiv, ( __CLPK_complex * )work, &lwork, &info);
	
	if ( work != NULL ) delete[] work;
	if( ipiv != NULL ) delete[] ipiv;
	
#endif
	
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}

//final_ok
template< >
int
TPZSBMatrixLapack< double >::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix< double >::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;
	int *ipiv = new int [n];
	int lwork = this->Dim();
	double *work = new double [n];
	
	dsytrf_(&uplo, &n, fDiag.begin(), &lda, ipiv, work, &lwork, &info);
	
	if ( work != NULL ) delete[] work;
	if( ipiv != NULL ) delete[] ipiv;
	
#endif
	
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}

//final_ok
template< >
int
TPZSBMatrixLapack< float >::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix< float >::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
#ifdef USING_LAPACK
	char uplo = 'u';
	int n = this->Dim();
	int lda = this->Dim();
	int info = -666;
	int *ipiv = new int [n];
	int lwork = this->Dim();
	float *work = new float [n];
	
	ssytrf_(&uplo, &n, fDiag.begin(), &lda, ipiv, work, &lwork, &info);
	
	if ( work != NULL ) delete[] work;
	if( ipiv != NULL ) delete[] ipiv;
	
#endif
	
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}
//final_ok
template< class TVar >
int
TPZSBMatrixLapack< TVar >::Decompose_LDLt()
{
	
	TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <LAPACK does not support this specific data type>" );
	DebugStop();
	return 0;
}

/*********************/
/*** Subst Forward ***/
//
//  Does Ax = b, where A is lower triangular matrix
//
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
	if ( ( B->Rows() != this->Dim() ) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag.begin() + fBand * sizeof ( TVar );//first element of main diagonal
	for ( long k = 0; k < this->Dim(); k++, row_k += k >= fBand ? fBand + 1 : fBand + 1 - k ){
		for ( long j = 0; j < B->Cols(); j++ )
		{
			// Does sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			TVar sum = this->gZero;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i - k > fBand ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current++;
			}
			// Does B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	}
	return( 1 );
}
//nok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
	if ( ( B->Rows() != this->Dim() ) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag.begin() + ( this->Size() - 1 ) * sizeof ( TVar );//last element of main diagonal
	for ( long k = this->Dim() ; k >= 0; k--, row_k -= fBand + 1 ){
		for ( long j = B->Cols(); j >= 0; j-- )
		{
			
			TVar sum = this->gZero;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i < k ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current+= fBand;
			}
			
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	}
	return( 1 );
	
}

/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//
//final_ok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
	if ( ( B->Rows() != this->Dim() ) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag.begin() + fBand * sizeof ( TVar );//first element of main diagonal
	for ( long k = 0; k < this->Dim(); k++, row_k += k >= fBand ? fBand + 1 : fBand + 1 - k ){
		for ( long j = 0; j < B->Cols(); j++ ){
			// Does sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			TVar sum = this->gZero;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i - k > fBand ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current++;
			}
			// Does B[k,j] = (B[k,j] - sum)
			//
			B->PutVal( k, j, B->GetVal(k, j) - sum  );
		}
	}
	return( 1 );
}

//final_ok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
	if ( ( B->Rows() != this->Dim() ) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag.begin() + ( this->Size() - 1 ) * sizeof ( TVar );//last element of main diagonal
	for ( long k = this->Dim() ; k >= 0; k--, row_k -= fBand + 1 ){
		for ( long j = B->Cols(); j >= 0; j-- )
		{
			
			TVar sum = this->gZero;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i < k ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current+= fBand;
			}
			
			B->PutVal( k, j, B->GetVal(k, j) - sum  );
		}
	}
	return( 1 );
	
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
//final_ok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
	
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
	
	
	TVar *row_k = fDiag.begin() + fBand * sizeof( TVar );//first element of main diagonal
	for ( long k = 0; k < this->Dim(); k++, row_k += ( fBand + 1 ) * sizeof( TVar ) )
		for ( long j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j) / *row_k );
	
	return( 1 );
}

/**************************** PRIVATE ****************************/

/*************/
/*** Clear ***/
//final_ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Clear()
{
	if ( this->fDiag.begin() != NULL ){
		fDiag.clear() ;
	}
	this->fRow = this->fCol = 0;
	this->fDecomposed = 0;
	return( 1 );
}

/************/
/*** Copy ***/
//final_ok
template<class TVar>
void
TPZSBMatrixLapack<TVar>::Copy(const TPZSBMatrixLapack<TVar> &A )
{
	this->fBand = A.fBand;
	this->fRow = this->fCol = A.Dim();
	this->fDiag.clear();
	this->fDiag = TPZVec <TVar> ( Size() , this->gZero );
	this->fDecomposed  = A.fDecomposed;
	this->fDefPositive = A.fDefPositive;
	
	if ( fDiag.begin() == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Copy( TPZSBMatrixLapack ) <memory allocation error>" );
	
	TVar *dst = fDiag.begin();
	TVar *src = A.fDiag.begin();
	TVar *end = A.fDiag.begin() + sizeof(TVar) * Size();
	while ( src < end )
		*dst++ = *src++;
}


// Inicializando os templates
template class TPZSBMatrixLapack< float >;
template class TPZSBMatrixLapack< double >;
template class TPZSBMatrixLapack< complex<double> >;
template class TPZSBMatrixLapack< complex<float> >;

