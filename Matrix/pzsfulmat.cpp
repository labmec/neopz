/**
 * @file
 * @brief Contains the implementation of the TPZSFMatrix methods.
 */

#include <math.h>
#include <stdlib.h>

#include "pzfmatrix.h"
#include "pzsfulmat.h"
//#include "pzerror.h"

#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzsfmatrix");
#endif

using namespace std;

/*******************/
/*** Constructor ***/
template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix(const int64_t dim )
: TPZRegisterClassId(&TPZSFMatrix::ClassId),
TPZMatrix<TVar>( dim, dim )
{
    int64_t size = Size();
	fElem = new TVar[ size ] ;
	
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Zera a Matriz.
	Zero();
}

/*********************************/
/*** Constructor( TPZSFMatrix& ) ***/
template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix (const TPZSFMatrix<TVar>  & A)
: TPZRegisterClassId(&TPZSFMatrix::ClassId),
TPZMatrix<TVar> ( A.Dim(), A.Dim() )
{
    int64_t size = Size();
	fElem = new TVar[size] ;
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	TVar *end = &fElem[Size()];
	while ( dst < end )
		*dst++ = *src++;
}

/*** Constructor( TPZSFMatrix& ) ***/
template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix (TPZSFMatrix<TVar>  && A)
: TPZRegisterClassId(&TPZSFMatrix::ClassId),
	TPZMatrix<TVar> (A),fElem(A.fElem)
{
	A.fElem=nullptr;
}


/*******************************/
/*** Constructor( TPZMatrix<TVar> & ) ***/

template<class TVar>
TPZSFMatrix<TVar> ::TPZSFMatrix(const  TPZMatrix<TVar>  &A )
: TPZRegisterClassId(&TPZSFMatrix::ClassId),
TPZMatrix<TVar> ( A )
{
    int64_t size = Size();
	fElem = new TVar[size] ;
	
	if ( fElem == NULL )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *dst = fElem;
	for ( int64_t r = 0; r < this->Dim(); r++ )
		for ( int64_t c = 0; c <= r; c++ )
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


template<class TVar>
TPZSFMatrix<TVar> &
TPZSFMatrix<TVar> ::operator=(const TPZSFMatrix<TVar>  &A )
{
	if ( this->Dim() != A.Dim() )
    {
		if ( fElem != NULL )
        {
			delete( fElem );
        }
        int64_t size = A.Size();
		fElem = new TVar[ size ] ;
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

template<class TVar>
TPZSFMatrix<TVar> &
TPZSFMatrix<TVar> ::operator=(TPZSFMatrix<TVar>  &&A )
{
	fElem=A.fElem;
	A.fElem=nullptr;
	return *this;
}
/******** Operacoes com matrizes FULL simetricas ********/

/******************/




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
	int64_t newDim = Min( A.Rows(), A.Cols() );
	int64_t size   = newDim * (newDim + 1) / 2;
	
	if ( newDim != this->Dim() )
    {
		if ( fElem != NULL )
			delete( fElem );
		fElem = new TVar[ size ] ;
    }
	
	// Copia a matriz.
	TVar *dst = fElem;
	for ( int64_t c = 0; c < newDim; c++ )
		for ( int64_t r = 0; r <= c; r++ )
			*dst++ = A.Get( r, c );
	
	this->fRow = this->fCol = newDim;
	
	// Ajusta Status de decomposicao.
	
	this->SetIsDecomposed(A.IsDecomposed());
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
	
	auto  res(*this);
	TVar *pr = res.fElem;
	for ( int64_t c = 0; c < this->Dim(); c++ )
		for ( int64_t r = 0; r <= c; r++ )
			*pr++ += A.Get( r, c );
	
	return res;
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSFMatrix<TVar>
TPZSFMatrix<TVar> ::operator-( const TPZMatrix<TVar>  &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Operator- <matrixs with incoompatible dimensions" );
	
	auto res(*this);
	TVar *pr = res.fElem;
	for ( int64_t c = 0; c < this->Dim(); c++ )
		for ( int64_t r = 0; r <= c; r++ )
			*pr++ -= A.Get( r, c );
	
	return res;
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
	for ( int64_t c = 0; c < this->Dim(); c++ )
		for ( int64_t r = 0; r <= c; r++ )
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
	for ( int64_t c = 0; c < this->Dim(); c++ )
		for ( int64_t r = 0; r <= c; r++ )
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
	auto  res(*this);
	
	TVar *dst = res.fElem;
	TVar *end = dst+Size();
	while ( dst < end )
		*dst++ *= value;
	
	return res;
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
TPZSFMatrix<TVar> ::Resize( int64_t newDim , int64_t )
{
	if ( newDim == this->Dim() )
		return( 1 );
	
	int64_t newSize = newDim * (newDim + 1) / 2;
	int64_t oldSize = Size();
	TVar *newElem = new TVar[newSize] ;
	if ( newElem == NULL )
		return TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	int64_t minSize  = Min( newSize, oldSize );
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
TPZSFMatrix<TVar> ::Redim( int64_t newDim , int64_t)
{
	// Se for preciso, desaloca a matriz antiga e aloca uma
	//  nova com o novo tamanho.
	if ( newDim != this->Dim() )
    {
		this->fRow = this->fCol = newDim;
		if ( fElem != NULL )
			delete( fElem );
		fElem = new TVar[Size()] ;
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
TPZSFMatrix<TVar> ::Decompose_Cholesky(std::list<int64_t> &singular)
{
	return Decompose_Cholesky();
}

template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_Cholesky()
{
	if (  this->fDecomposed )  TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	TVar *ptr_k = fElem;
	for ( int64_t k = 0; k < this->Dim(); k++, ptr_k += k  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum  = 0.0;
		TVar *pk  = ptr_k;
		TVar *pkk = ptr_k + k;
		for ( ; pk < pkk;  pk++ ){
			if constexpr(is_complex<TVar>::value){
				sum += std::conj(*pk) * (*pk);
			}else{
				sum += (*pk) * (*pk);
			}
		}
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		if (IsZero(*pk -= sum))
			return( 0 );
		*pk = sqrt( *pk );
		
		// Loop para i = k+1 ... Dim().
		//
		TVar *ptr_i = ptr_k;
		TVar *pi;
		for ( int64_t i = k+1; i <this->Dim(); i++ )
		{
			// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1, ..., k-1.
			//
			sum   =  0.0;
			ptr_i += i;
			pk    =  ptr_k;
			pi    =  ptr_i;
			for ( int64_t p = 0; p < k; p++ ){
				if constexpr(is_complex<TVar>::value){
					sum += std::conj(*pk++) * (*pi++);
				}else{
					sum += (*pk++) * (*pi++);
				}
			}
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
TPZSFMatrix<TVar> ::Decompose_LDLt(std::list<int64_t> &singular)
{
	return Decompose_LDLt();
}

template<class TVar>
int
TPZSFMatrix<TVar> ::Decompose_LDLt()
{
	if (  this->fDecomposed && this->fDecomposed != ELDLt)
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Matrix already Decomposed with a different scheme>" );
	else if ( this->fDecomposed ) return 0;
	
	//	TVar sum;
	
	
	for ( auto j = 0; j < this->Dim(); j++ )
    {
		TVar sum=0.;
		for ( auto k=0; k<j; k++){
			if constexpr (is_complex<TVar>::value){
				sum+=GetVal(k,k)*std::conj(GetVal(k,j))*GetVal(k,j);
			}else{
				sum+=GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
			}
		}
		PutVal(j,j,GetVal(j,j)-sum);
		if ( IsZero(GetVal(j,j)) )TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Zero on diagonal>" );
		for( auto l=j+1; l<this->Dim();l++)
		{
			sum = 0;
			for ( auto k=0; k<j; k++){
				if constexpr(is_complex<TVar>::value){
					sum+=GetVal(k,k)*std::conj(GetVal(j,k))*GetVal(l,k);
				}else{
					sum+=GetVal(k,k)*GetVal(j,k)*GetVal(l,k);
				}
			}
			PutVal(l,j, (GetVal(l,j)-sum)/GetVal(j,j) );
		}
		
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
	
	if ( B->IsSymmetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_Forward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = fElem;
	for ( int64_t k = 0; k < this->Dim(); k++ )
    {
		for ( int64_t j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1,.., k-1.
			//
			TVar *pk = ptr_k;
			TVar sum = 0.0;
			for ( int64_t i = 0; i < k; i++ )
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
	
	if ( B->IsSymmetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_Backward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = &fElem[ Size()-1 ];
	for ( int64_t k = this->Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int64_t j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			TVar sum = 0.0;
			TVar *pk = ptr_k;
			for ( int64_t i = this->Dim()-1; i > k; i-- )
			{
				if constexpr(is_complex<TVar>::value){
					sum += std::conj(*pk) * B->GetVal( i, j );
				}else{
					sum += (*pk) * B->GetVal( i, j );
				}
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
	
	if ( B->IsSymmetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_LForward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = fElem;
	for ( int64_t k = 0; k < this->Dim(); k++, ptr_k += k )
		for ( int64_t j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			TVar *pk = ptr_k;
			TVar sum = 0.0;
			for ( int64_t i = 0; i < k; i++ )
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
	
	if ( B->IsSymmetric() )
		TPZMatrix<TVar> ::Error(__PRETTY_FUNCTION__, "Subst_LBackward <the matrix result can not be simetric>" );
	
	TVar *ptr_k = &fElem[ Size()-1 ];
	for ( int64_t k = this->Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int64_t j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			TVar sum = 0.0;
			TVar *pk = ptr_k;
			for ( int64_t i = this->Dim()-1; i > k; i-- )
			{
				if constexpr (is_complex<TVar>::value){
					sum += std::conj(*pk) * B->GetVal( i, j );
				}else{
					sum += *pk * B->GetVal( i, j );
				}
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
	for ( int64_t k = 0; k < this->Dim(); k++ )
    {
		for ( int64_t j = 0; j < B->Cols(); j++ )
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
TSaveable *TPZSFMatrix<TVar>::CreateInstance(TReceiveStorage *buf) {
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
int TPZSFMatrix<TVar>::DerivedFrom(int64_t Classid){
	return TSaveable::DerivedFrom(Classid);
}
template<class TVar>
int TPZSFMatrix<TVar>::DerivedFrom(char *classname){
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TSaveable::DerivedFrom(classname);
}

#endif

template<class TVar>
int TPZSFMatrix<TVar>::ClassId() const{
    return Hash("TPZSFMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

template class TPZSFMatrix<float>;
template class TPZSFMatrix<double>;
template class TPZSFMatrix<long double>;
template class TPZSFMatrix<std::complex<float>>;
template class TPZSFMatrix<std::complex<double>>;
template class TPZSFMatrix<std::complex<long double>>;