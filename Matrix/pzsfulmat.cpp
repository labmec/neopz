/**
 * @file
 * @brief Contains the implementation of the TPZSFMatrix methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tsfulmat.cc
//
// Class:  TPZSFMatrix
//
// Obs.:   Implementa matrizes cheias (normais).
//
// Versao: 04 / 1996.
//

#include <math.h>
#include <stdlib.h>

#include "pzfmatrix.h"
#include "pzsfulmat.h"
//#include "pzerror.h"

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzsfmatrix"));
#endif

using namespace std;

/*******************/
/*** Constructor ***/
template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix(const int dim )
: TPZMatrix<TVar>( dim, dim )
{
	fElem = new( TVar[ Size() ] );
	
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Zera a Matriz.
	Zero();
}

/*********************************/
/*** Constructor( TPZSFMatrix& ) ***/
template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix (const TPZSFMatrix<TVar>  & A)
: TPZMatrix<TVar> ( A.Dim(), A.Dim() )
{
	fElem = new( REAL[Size()] );
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	TVar *end = &fElem[Size()];
	while ( dst < end )
		*dst++ = *src++;
}



/*******************************/
/*** Constructor( TPZMatrix<TVar> & ) ***/

template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix(const  TPZMatrix<TVar>  &A )
: TPZMatrix<TVar> ( A )
{
	fElem = new( TVar[Size()] );
	
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *dst = fElem;
	for ( int r = 0; r < this->Dim(); r++ )
		for ( int c = 0; c <= r; c++ )
			*dst++ = A.GetVal( r, c );
}



/******************/
/*** Destructor ***/

template<class TVar>
TPZSFMatrix<TVar> ::~TPZSFMatrix ()
{
	if ( fElem != NULL )
		delete []fElem;
}



/******** Operacoes com matrizes FULL simetricas ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSFMatrix<TVar> &
TPZSFMatrix<TVar> ::operator=(const TPZSFMatrix<TVar>  &A )
{
	if ( this->Dim() != A.Dim() )
    {
		if ( fElem != NULL )
			delete( fElem );
		fElem = new( TVar[ A.Size() ] );
		if ( fElem == NULL )
			TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Operator= <memory allocation error>." );
    }
	
	this->fRow = this->fCol =  A.Dim();
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = *src++;
	
	this->fDecomposed = A.fDecomposed;
	
	return *this;
}



/*******************************/
/*** Operator+( TPZSFMatrix & ) ***/

template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator+(const TPZSFMatrix<TVar>  &A ) const
{
	if ( A.Dim() != this->Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
	TPZSFMatrix<TVar>  res( this->Dim() );
	TVar *pm  = fElem;
	TVar *pa  = A.fElem;
	TVar *pr  = res.fElem;
	TVar *end = &(res.fElem[Size()]);
	
	while ( pr < end )
		*pr++ = (*pm++) + (*pa++);
	
	return( res );
}



/*******************************/
/*** Operator-( TPZSFMatrix & ) ***/
template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator-(const TPZSFMatrix<TVar>  &A ) const
{
	if ( A.Dim() != this->Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
	TPZSFMatrix<TVar>  res( this->Dim() );
	TVar *pm  = fElem;
	TVar *pa  = A.fElem;
	TVar *pr  = res.fElem;
	TVar *end = &(res.fElem[Size()]);
	
	while ( pr < end )
		*pr++ = (*pm++) - (*pa++);
	
	return( res );
}



/********************************/
/*** Operator+=( TPZSFMatrix & ) ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator+=(const TPZSFMatrix<TVar>  &A )
{
	if ( A.Dim() != this->Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Operator+= <matrixs with different dimensions>" );
	
	TVar *src = A.fElem;
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ += *src++;
	
	this->fDecomposed = 0;
	return( *this );
}



/*******************************/
/*** Operator-=( TPZSFMatrix & ) ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator-=(const TPZSFMatrix<TVar>  &A )
{
	if ( A.Dim() != this->Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator-= <matrixs with different dimensions>" );
	
	TVar *src = A.fElem;
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ -= *src++;
	
	this->fDecomposed = 0;
	return( *this );
}



/******** Operacoes com matrizes GENERICAS ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator=(const TPZMatrix<TVar>  &A )
{
	int newDim = Min( A.Rows(), A.Cols() );
	int size   = newDim * (newDim + 1) / 2;
	
	if ( newDim != this->Dim() )
    {
		if ( fElem != NULL )
			delete( fElem );
		fElem = new REAL[ size ] ;
    }
	
	// Copia a matriz.
	TVar *dst = fElem;
	for ( int c = 0; c < newDim; c++ )
		for ( int r = 0; r <= c; r++ )
			*dst++ = A.Get( r, c );
	
	this->fRow = this->fCol = newDim;
	
	// Ajusta Status de decomposicao.
	
	SetIsDecomposed(A.IsDecomposed());
	this->fDefPositive = (char)A.IsDefPositive();
	
	return( *this );
}



/******************/
/*** Operator + ***/

template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator+(const TPZMatrix<TVar>  &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator+ <matrixs with incoompatible dimensions" );
	
	TPZSFMatrix<TVar>  res( this->Dim() );
	TVar *pm = fElem;
	TVar *pr = res.fElem;
	for ( int c = 0; c < this->Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pr++ = (*pm++) + A.Get( r, c );
	
	return( *this );
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator-( const TPZMatrix<TVar>  &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator- <matrixs with incoompatible dimensions" );
	
	TPZSFMatrix<TVar>  res( this->Dim() );
	TVar *pm = fElem;
	TVar *pr = res.fElem;
	for ( int c = 0; c < this->Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pr++ = (*pm++) - A.Get( r, c );
	
	return( *this );
}



/*******************/
/*** Operator += ***/
template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator+=(const TPZMatrix<TVar>  &A )
{
	if (this->Dim() != A.Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator+= (TPZMatrix<>&) <different dimensions>" );
	
	TVar *pm = fElem;
	for ( int c = 0; c < this->Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pm++ += A.Get( r, c );
	
	this->fDecomposed = 0;
	return( *this );
}



/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator-=(const TPZMatrix<TVar>  &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator-= (TPZMatrix<>&) <different dimensions>" );
	
	TVar *pm = fElem;
	for ( int c = 0; c < this->Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pm++ -= A.Get( r, c );
	
	this->fDecomposed = 0;
	return( *this );
}



/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSFMatrix<TVar> &
TPZSFMatrix<TVar> ::operator=(const TVar value )
{
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = value;
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( *this );
}



/**************************/
/*** Operator+( value ) ***/

template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator+(const TVar value ) const
{
	TPZSFMatrix<TVar>  res( this->Dim() );
	
	TVar *dst = res.fElem;
	TVar *src = fElem;
	TVar *end = &fElem[ Size() ];
	while ( src < end )
		*dst++ = (*src++) + value;
	
	return( res );
}



/**************************/
/*** Operator*( value ) ***/

template<class TVar>
TPZSFMatrix<TVar> 
TPZSFMatrix<TVar> ::operator*(const TVar value ) const
{
	TPZSFMatrix<TVar>  res( this->Dim() );
	
	TVar *dst = res.fElem;
	TVar *src = fElem;
	TVar *end = &fElem[ Size() ];
	while ( src < end )
		*dst++ = (*src++) * value;
	
	return( res );
}



/***************************/
/*** Operator+=( value ) ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator+=( TVar value )
{
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ += value;
	this->fDecomposed = 0;
	return( *this );
}



/***************************/
/*** Operator*=( value ) ***/

template<class TVar>
TPZSFMatrix<TVar>  &
TPZSFMatrix<TVar> ::operator*=( TVar value )
{
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ *= value;
	this->fDecomposed = 0;
	return( *this );
}



/**************/
/*** Resize ***/
template<class TVar>
int
TPZSFMatrix<TVar> ::Resize( int newDim , int )
{
	if ( newDim == this->Dim() )
		return( 1 );
	
	int newSize = newDim * (newDim + 1) / 2;
	int oldSize = Size();
	TVar *newElem = new TVar[newSize] ;
	if ( newElem == NULL )
		return TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	int minSize  = Min( newSize, oldSize );
	TVar *src = fElem;
	TVar *dst = newElem;
	TVar *end = &fElem[ minSize ];
	
	// Copia os elementos antigos para a nova matriz.
	while ( src < end )
		*dst++ = *src++;
	
	// Preenche os elementos que sobrarem (se sobrarem) com ZEROS.
	end = &newElem[ newSize ];
	while ( dst < end )
		*dst++ = 0.0;
	
	if ( fElem != NULL )
		delete( fElem );
	fElem = newElem;
	this->fRow = this->fCol = newDim;
	this->fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/

template<class TVar>
int
TPZSFMatrix<TVar> ::Redim( int newDim , int)
{
	// Se for preciso, desaloca a matriz antiga e aloca uma
	//  nova com o novo tamanho.
	if ( newDim != this->Dim() )
    {
		this->fRow = this->fCol = newDim;
		if ( fElem != NULL )
			delete( fElem );
		fElem = new( REAL[Size()] );
    }
	
	// Zera a matriz.
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = 0.;
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	
	
	return( 1 );
}


template<class TVar>
int
TPZSFMatrix<TVar> ::Zero()
{
	TVar *dst = fElem;
	TVar *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = 0.;
	this->fDecomposed = 0;
	this->fDefPositive = 0;
	return( 1 );
}
/******** Resolucao de Sistemas ********/

/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_Cholesky(std::list<int> &singular)
{
	return Decompose_Cholesky();
}

template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_Cholesky()
{
	if (  this->fDecomposed )  TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	TVar *ptr_k = fElem;
	for ( int k = 0; k < this->Dim(); k++, ptr_k += k  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum  = 0.0;
		TVar *pk  = ptr_k;
		TVar *pkk = ptr_k + k;
		for ( ; pk < pkk;  pk++ )
			sum += (*pk) * (*pk);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		if ( (*pk -= sum) < 1.e-10 )
			return( 0 );
		*pk = sqrt( *pk );
		
		// Loop para i = k+1 ... Dim().
		//
		TVar *ptr_i = ptr_k;
		TVar *pi;
		for ( int i = k+1; i <this->Dim(); i++ )
		{
			// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1, ..., k-1.
			//
			sum   =  0.0;
			ptr_i += i;
			pk    =  ptr_k;
			pi    =  ptr_i;
			for ( int p = 0; p < k; p++ )
				sum += (*pk++) * (*pi++);
			
			// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
			//
			*pi = (*pi - sum) / *pk;
		}
    }
	
	this->fDecomposed  = ECholesky;
	this->fDefPositive = 1;
	return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_LDLt(std::list<int> &singular)
{
	return Decompose_LDLt();
}

template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_LDLt()
{
	/* REAL *aux;
     aux=new REAL [ Dim() - 1 ];
     REAL *ptr_k = fElem;
	 
     for ( int k = 0; k < Dim(); k++, ptr_k += k  )
     {
     // Faz aux(p) = A(p,p) * A(k,p), p = 1, ..., k-1.
     //
     REAL *pk    = ptr_k;
     REAL *pDiag = fElem;
     REAL *pAux  = aux;
     for ( int p = 0; p < k; p++, pDiag += p )
     *pAux++ = (*pDiag++) * (*pk++);
	 
	 // Faz sum = SOMA( A(k,p) * aux(p) ), p = 1, ..., k-1.
	 //
	 REAL sum  = 0.0;
	 REAL *end = ptr_k + k;
	 pk   = ptr_k;
	 pAux = aux;
	 while ( pk < end )
	 sum += (*pk++) * (*pAux++);
	 
	 // Faz A(k,k) = A(k,k) - sum;
	 //
	 *pk -= sum;
	 if ( IsZero(*pk) )
	 return( 0 );
	 
	 // Loop para i = k+1 ... Dim().
	 //
	 REAL *ptr_i = ptr_k;
	 REAL *pi;
	 for ( int i = k+1; i < Dim(); i++ )
	 {
	 // Faz sum = SOMA( A(i,p) * aux(p) ), p = 1, ..., k-1.
	 //
	 sum   =  0.0;
	 ptr_i += i;
	 pAux  =  aux;
	 pi    =  ptr_i;
	 end   =  ptr_i + k;
	 while ( pi < end )
	 sum += (*pi++) * (*pAux++);
	 
	 // Faz A(i,k) = (A(i,k) - sum) / A(k,k).
	 //
	 *pi = (*pi - sum) / *pk;
	 }
	 }
	 */
	if (  this->fDecomposed && this->fDecomposed != ELDLt)
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Matrix already Decomposed with a different scheme>" );
	else if ( this->fDecomposed ) return 0;
	
	int j,k,l;
	//	REAL sum;
	
	
	for ( j = 0; j < this->Dim(); j++ )
    {
		//	sum=0.;
		for ( k=0; k<j; k++)
			//			sum=sum-GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
			PutVal( j,j,GetVal(j,j)-GetVal(k,k)*GetVal(k,j)*GetVal(k,j) );
		//		PutVal(j,j,GetVal(j,j)+sum);
		for ( k=0; k<j; k++)
		{
			
			for( l=j+1; l<this->Dim();l++)
				PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
		}
		
		if ( IsZero(GetVal(j,j)) )TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Zero on diagonal>" );
		
		for( l=j+1; l<this->Dim();l++) PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
		
    }
	this->fDecomposed  = ELDLt;
	this->fDefPositive = 0;
	
	return( 1 );
}



/*********************/
/*** Subst Forward ***/

template<class TVar>
int
TPZSFMatrix<TVar> ::Subst_Forward( TPZFMatrix<TVar>  *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_Forward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = fElem;
	for ( int k = 0; k < this->Dim(); k++ )
    {
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1,.., k-1.
			//
			REAL *pk = ptr_k;
			REAL sum = 0.0;
			for ( int i = 0; i < k; i++ )
				sum += (*pk++) * B->GetVal( i, j );
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal( k, j ) - sum) / *pk );
		}
		ptr_k += k + 1;
    }
	
	return( 1 );
}



/**********************/
/*** Subst Backward ***/
template<class TVar>
int
TPZSFMatrix<TVar> ::Subst_Backward( TPZFMatrix<TVar>  *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_Backward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = &fElem[ Size()-1 ];
	for ( int k = this->Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			TVar sum = 0.0;
			TVar *pk = ptr_k;
			for ( int i = this->Dim()-1; i > k; i-- )
			{
				sum += *pk * B->GetVal( i, j );
				pk -= i;
			}
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *pk );
		}
	
	return( 1 );
}



/***********************/
/*** Subst L Forward ***/
template<class TVar>
int
TPZSFMatrix<TVar> ::Subst_LForward( TPZFMatrix<TVar>  *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_LForward <the matrix result can not be simetric>" );
	
	REAL *ptr_k = fElem;
	for ( int k = 0; k < this->Dim(); k++, ptr_k += k )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			TVar *pk = ptr_k;
			TVar sum = 0.0;
			for ( int i = 0; i < k; i++ )
				sum += (*pk++) * B->GetVal( i, j );
			
			// Faz b[k] = (b[k] - sum).
			//
			B->PutVal( k, j, B->GetVal( k, j ) - sum );
		}
	
	return( 1 );
}



/************************/
/*** Subst L Backward ***/

template<class TVar>
int
TPZSFMatrix<TVar> ::Subst_LBackward( TPZFMatrix<TVar>  *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_LBackward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = &fElem[ Size()-1 ];
	for ( int k = this->Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			TVar sum = 0.0;
			TVar *pk = ptr_k;
			for ( int i = this->Dim()-1; i > k; i-- )
			{
				sum += *pk * B->GetVal( i, j );
				pk -= i;
			}
			
			// Faz B[k,j] = B[k,j] - sum.
			//
			B->PutVal( k, j, B->GetVal(k, j) - sum );
		}
	
	return( 1 );
}



/******************/
/*** Subst Diag ***/

template<class TVar>
int
TPZSFMatrix<TVar> ::Subst_Diag( TPZFMatrix<TVar>  *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		return( 0 );
	
	TVar *pDiag = fElem;
	for ( int k = 0; k < this->Dim(); k++ )
    {
		for ( int j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j ) / *pDiag );
		pDiag += (k + 2);
    }
	
	return( 1 );
}




/************************** Private **************************/

/*************/
/*** Error ***/

/*int
 TPZSFMatrix::Error(const char *msg1,const char *msg2 )
 {
 ostringstream out;  
 out << "TPZSFMatrix::" << msg1 << msg2 << ".\n";
 //pzerror.Show();
 LOGPZ_ERROR (logger, out.str().c_str());
 DebugStop();
 return 0;
 }*/



/*************/
/*** Clear ***/

template<class TVar>
int
TPZSFMatrix<TVar> ::Clear()
{
	if ( fElem != NULL )
		delete( fElem );
	
	fElem = NULL;
	this->fRow = this->fCol = 0;
	return( 1 );
}

#ifdef OOPARLIB

template<class TVar>
int TPZSFMatrix<TVar> ::Unpack( TReceiveStorage *buf ){
	TSaveable::Unpack(buf);
	buf->UpkDouble(fElem);
	return 1;
}


template<class TVar>
TSaveable *TPZSFMatrix<TVar>::Restore(TReceiveStorage *buf) {
	TPZSFMatrix<TVar> *m = new TPZSFMatrix<TVar>();
	m->Unpack(buf);
	return m;
}

template<class TVar>
int TPZSFMatrix<TVar>::Pack( TSendStorage *buf ){
	TSaveable::Pack(buf);
	buf->PkDouble(fElem);
	return 1;
}

template<class TVar>
int TPZSFMatrix<TVar>::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSaveable::DerivedFrom(Classid);
}
template<class TVar>
int TPZSFMatrix<TVar>::DerivedFrom(char *classname){
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TSaveable::DerivedFrom(classname);
}

#endif

