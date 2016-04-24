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

//ok
template<class TVar> TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack() : TPZMatrix<TVar>(0,0) {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  fDiag = NULL;
  fBand = 0;
}
//ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack( long dim, long band )
: TPZMatrix<TVar>( dim, dim )
{
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
	fBand = ( band > (dim - 1) ? (dim - 1) : band );
	fDiag = new TVar[Size()];
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrixLapack( dim ) <Error creating Matrix>" );
	
	Zero();
}

//ok
template<class TVar>
TPZSBMatrixLapack<TVar>::TPZSBMatrixLapack(const TPZSBMatrixLapack<TVar> &A ) : TPZMatrix<TVar>(A)  {
#ifndef USING_LAPACK //for symmetric band matrices without lapack see the TPZSBMatrix class
  DebugStop();
#endif
  Copy(A);
}


/**************/
/*** PutVal ***/
//ok
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
  const long dim = this->fCol;//square matrix
  fDiag[ dim * (fBand - index) + col ] = value;
  this->fDecomposed = 0;
  return( 1 );
}

/**************/
/*** GetVal ***/
//ok
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
	return( fDiag[ dim * (fBand - index) + col ] );
}

/*************/
/*** Print ***/
//ok
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
//ok
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
//ok
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
//ok
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
//ok
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
//ok
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
    return *this;

}

/*******************/
/*** Operator -= ***/
//ok
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
//ok
template<class TVar>
TPZSBMatrixLapack<TVar>
TPZSBMatrixLapack<TVar>::operator*(const TVar value ) const
{
	TPZSBMatrixLapack<TVar> res( this->Dim(), fBand );
	
	TVar *pr  = res.fDiag;
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	while ( pm < end )
		*pr++ = (*pm++) * value;
	return( res );
}

/******************************/
/*** Operator += ( REAL ) ***/
//ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator+=(const TVar value )
{
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	for(long iDiag = 0 ; iDiag < fBand + 1 ; iDiag ++){
		long nElemDiag = this->Dim() - fBand + iDiag;
		*pm += fBand - iDiag;//skip zeros
		for(int i = 0 ; i < nElemDiag ; i++){
			*pm++ += value;
		}
	}while ( pm < end )
	
	return( *this );
}

/******************************/
/*** Operator *= ( REAL ) ***/
//ok
template<class TVar>
TPZSBMatrixLapack<TVar> &
TPZSBMatrixLapack<TVar>::operator*=(const TVar value )
{
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	while ( pm < end )
		*pm++ *= value;
	
	this->fDecomposed = 0;
	return( *this );
}

/**************/
/*** Resize ***/
//
// Changes matrix dimensions, but keeps its current values.
// New positions are created with zeros.
//
//ok
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
  
  //if newBand < fBand, some diagonals are going to be skipped
  long nSkipDiags = fBand - newBand;
  //iDiag begins in the smaller diagonal
  for (long iDiag = 0; iDiag < newBand + 1; iDiag ++) {
    TVar *ini = newDiag + newDim * iDiag;
    long pos = 0;
    for ( ; pos < newBand - iDiag ; pos ++) {
      *ini++ = 0.; //fills the new positions with zeros
    }
    TVar *oldIni = fDiag + this->Dim() * (iDiag + nSkipDiags) + fBand - iDiag;//first valid element of corresponding diagonal
    for ( ; pos < newDim ; pos++){
      *ini++ = *oldIni++;
    }
  }
	// Deletes old matrix and takes new one
	if ( fDiag != NULL )
    DebugStop(); //it really shouldnt be null
	fDiag = newDiag;
	this->fCol = this->fRow = newDim;
	this->fDecomposed = 0;
	return( 1 );
}

/*************/
/*** Redim ***/
//
// Changes matrix dimwensions and erases all its elements
//
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Redim(const long newDim ,const long)
{
	if ( newDim != this->Dim() )
    {
		// Aloca a nova matriz.
		this->fRow = this->fCol = newDim;
		if ( this->fDiag != NULL )
			delete( this->fDiag );
		this->fDiag = new TVar[Size()];
		if ( fDiag == NULL )
			TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>" );
    }
	
	TVar *dst = fDiag;
	TVar *end = &fDiag[Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Zero()
{
	TVar *dst =this->fDiag;
	TVar *end = &fDiag[this->Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}


/****************/
/*** Set Band ***/
//ok
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
	TVar *pNew = newDiag;
	TVar *pOld = fDiag;
	long newSize  = newBand + 1;
	long oldSize  = fBand + 1;
	long diffSize = abs( newSize - oldSize );
	long minSize = MIN(newSize , oldSize);
	TVar *src = newSize > oldSize ? pOld : pOld + diffSize * this->Dim();
	TVar *dst = newSize > oldSize ? pNew + diffSize * this->Dim() : pNew;
	
	for ( int i = 0; i < minSize * this->Dim() ; i++ ){
		*dst++ = *src++;
	}
	if ( newBand > fBand ){
		for(int i = 0 ; i < diffSize * this->Dim() ; i++){
			*pNew++ = 0.;
		}
	}
	else{
		this->fDecomposed = 0;
	}
	
	if ( fDiag == NULL )
		DebugStop();//it really wasnt supposed to be null
	fDiag = newDiag;
	fBand = newBand;
	
	return( 1 );
}

/********************* Resolucao de sistemas *********************/

/**************************/
/*** Decompose Cholesky ***/
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky(std::list<long> &singular)
{
	return Decompose_Cholesky();
}
//NOK

template<>
int
TPZSBMatrixLapack<std::complex<double> >::Decompose_Cholesky()
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
    cpotrf_(&uplo, &n, fDiag, &lda, &info); //nao compila
#endif
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}


template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
  
#ifdef USING_LAPACK
  char uplo = 'u';
  int n = this->Dim();
  int lda = this->Dim();
  int info = -666;
#ifdef STATE_COMPLEX
  //cpotrf_(&uplo, &n, fDiag, &lda, &info); //nao compila
#else
  spotrf_(&uplo, &n, fDiag, &lda, &info);
#endif
#endif
  
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}

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
#ifdef STATE_COMPLEX
    //cpotrf_(&uplo, &n, fDiag, &lda, &info); //nao compila
#else
    spotrf_(&uplo, &n, fDiag, &lda, &info);
#endif
#endif
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}

/**********************/
/*** Decompose LDLt ***/
//OK
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_LDLt(std::list<long> &singular)
{
	return Decompose_LDLt();
}
//NOK
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
#ifdef USING_LAPACK
  char uplo = 'u';
  int n = this->Dim();
  int lda = this->Dim();
  int info = -666;
  int *ipiv = new int [n];
  int lwork = this->Dim();
  TVar *work = new TVar [n];
  
#ifdef STATE_COMPLEX
  
  //csytrf_(&uplo, &n, fDiag, &lda, ipiv, work, &lwork, &info);
  
#else
  
  dsytrf_(&uplo, &n, fDiag, &lda, ipiv, work, &lwork, &info);
  
#endif
  
  if ( work != NULL ) delete[] work;
  if( ipiv != NULL ) delete[] ipiv;
  
#endif
  
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}

/*********************/
/*** Subst Forward ***/
//
//  Does Ax = b, where A is lower triangular matrix
//
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag + fBand * this->Dim();//first element of main diagonal
	for ( long k = 0; k < this->Dim(); k++, row_k += this->Dim() - k )
		for ( long j = 0; j < B->Cols(); j++ )
		{
			// Does sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			TVar sum = 0.0;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i - k > fBand ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current += this->Dim();
			}
			// Does B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	
	return( 1 );
}
//nok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	long k,j,i,jmax,stepcol=fBand+2;
	for(k=0; k<B->Cols() ; k++)
    {
		for(i=this->Rows()-1; i>=0; i--)
		{
			TVar *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>this->Rows()-1)? this->Rows()-1 : i+fBand;
			for(j=i+1;j<=jmax;j++)
			{
				diagk+=stepcol;
				B->operator()(i,k)-=B->GetVal(j,k)*(*diagk);
			}
			B->operator()(i,k)/=GetVal(i,i);
			
		}
		
    }
	
	return ( 1 ) ;
	
}

/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//
//ok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
	
	TVar *row_k = fDiag + fBand * this->Dim();//first element of main diagonal
	for ( long k = 0; k < this->Dim(); k++, row_k += this->Dim() - k )
		for ( long j = 0; j < B->Cols(); j++ )
		{
			// Does sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			TVar sum = 0.0;
			TVar *current = row_k;
			for ( long i = 0; i < k ; i++ ){
				if ( i - k > fBand ) continue; //empty positions
				sum += (*current) * B->GetVal( i, j );
				current += this->Dim();
			}
			// Does B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, B->GetVal(k, j) - sum );
		}
	
	return( 1 );
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
//ok
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
	
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
	
	
	TVar *row_k = fDiag + fBand * this->Dim();
	for ( long k = 0; k < this->Dim(); k++, row_k ++ )
		for ( long j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j) / *row_k );
	
	return( 1 );
}
//NOK
template<class TVar>
int TPZSBMatrixLapack<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
	
	long k,j,i,jmax,stepcol=fBand+2;
	
	for(k=0; k<B->Cols() ; k++)
    {
		for(i=this->Rows()-1; i>=0; i--)
		{
			TVar *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>this->Rows()-1)? this->Rows()-1 : i+fBand;
			for(j=i+1;j<=jmax;j++)
			{
				diagk+=stepcol;
				B->operator()(i,k)-=B->GetVal(j,k)*(*diagk);
			}
		}
		
    }
	
	return 1;
	
}

/**************************** PRIVATE ****************************/

/*************/
/*** CLear ***/
//ok
template<class TVar>
int
TPZSBMatrixLapack<TVar>::Clear()
{
	if ( this->fDiag != NULL )
		delete []fDiag ;
	this->fRow = this->fCol = 0;
	fDiag = NULL;
	this->fDecomposed = 0;
	return( 1 );
}

/************/
/*** Copy ***/
//ok
template<class TVar>
void
TPZSBMatrixLapack<TVar>::Copy(const TPZSBMatrixLapack<TVar> &A )
{
	this->fBand = A.fBand;
	this->fRow = this->fCol = A.Dim();
	this->fDiag = new TVar[Size()];
	this->fDecomposed  = A.fDecomposed;
	this->fDefPositive = A.fDefPositive;
	
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Copy( TPZSBMatrixLapack ) <memory allocation error>" );
	
	TVar *dst = fDiag;
	TVar *src = A.fDiag;
	TVar *end = A.fDiag + Size();
	while ( src < end )
		*dst++ = *src++;
}


// Inicializando os templates
template class TPZSBMatrixLapack< float >;
template class TPZSBMatrixLapack< double >;
template class TPZSBMatrixLapack< long double >;
template class TPZSBMatrixLapack< complex<double> >;
template class TPZSBMatrixLapack< complex<float> >;

