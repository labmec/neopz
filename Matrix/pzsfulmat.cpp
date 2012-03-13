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

TPZSFMatrix::TPZSFMatrix(const int dim )
: TPZMatrix( dim, dim )
{
	fElem = new REAL[ Size() ] ;
	
	if ( fElem == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Zera a Matriz.
	Zero();
}

/*********************************/
/*** Constructor( TPZSFMatrix& ) ***/

TPZSFMatrix::TPZSFMatrix (const TPZSFMatrix & A)
: TPZMatrix( A.Dim(), A.Dim() )
{
	fElem = new REAL[Size()] ;
	if ( fElem == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	REAL *src = A.fElem;
	REAL *dst = fElem;
	REAL *end = &fElem[Size()];
	while ( dst < end )
		*dst++ = *src++;
}



/*******************************/
/*** Constructor( TPZMatrix& ) ***/

TPZSFMatrix::TPZSFMatrix(const  TPZMatrix &A )
: TPZMatrix( A )
{
	fElem = new REAL[Size()] ;
	
	if ( fElem == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	REAL *dst = fElem;
	for ( int r = 0; r < Dim(); r++ )
		for ( int c = 0; c <= r; c++ )
			*dst++ = A.GetVal( r, c );
}



/******************/
/*** Destructor ***/

TPZSFMatrix::~TPZSFMatrix ()
{
	if ( fElem != NULL )
		delete []fElem;
}



/******** Operacoes com matrizes FULL simetricas ********/

/******************/
/*** Operator = ***/

TPZSFMatrix&
TPZSFMatrix::operator=(const TPZSFMatrix &A )
{
	if ( Dim() != A.Dim() )
    {
		if ( fElem != NULL )
			delete( fElem );
		fElem = new REAL[ A.Size() ] ;
		if ( fElem == NULL )
			TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator= <memory allocation error>." );
    }
	
	fRow = fCol =  A.Dim();
	
	// Copia a matriz
	REAL *src = A.fElem;
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = *src++;
	
	fDecomposed = A.fDecomposed;
	
	return *this;
}



/*******************************/
/*** Operator+( TPZSFMatrix & ) ***/

TPZSFMatrix
TPZSFMatrix::operator+(const TPZSFMatrix &A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
	TPZSFMatrix res( Dim() );
	REAL *pm  = fElem;
	REAL *pa  = A.fElem;
	REAL *pr  = res.fElem;
	REAL *end = &(res.fElem[Size()]);
	
	while ( pr < end )
		*pr++ = (*pm++) + (*pa++);
	
	return( res );
}



/*******************************/
/*** Operator-( TPZSFMatrix & ) ***/

TPZSFMatrix
TPZSFMatrix::operator-(const TPZSFMatrix &A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
	TPZSFMatrix res( Dim() );
	REAL *pm  = fElem;
	REAL *pa  = A.fElem;
	REAL *pr  = res.fElem;
	REAL *end = &(res.fElem[Size()]);
	
	while ( pr < end )
		*pr++ = (*pm++) - (*pa++);
	
	return( res );
}



/********************************/
/*** Operator+=( TPZSFMatrix & ) ***/

TPZSFMatrix &
TPZSFMatrix::operator+=(const TPZSFMatrix &A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator+= <matrixs with different dimensions>" );
	
	REAL *src = A.fElem;
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ += *src++;
	
	fDecomposed = 0;
	return( *this );
}



/*******************************/
/*** Operator-=( TPZSFMatrix & ) ***/

TPZSFMatrix &
TPZSFMatrix::operator-=(const TPZSFMatrix &A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Operator-= <matrixs with different dimensions>" );
	
	REAL *src = A.fElem;
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ -= *src++;
	
	fDecomposed = 0;
	return( *this );
}



/******** Operacoes com matrizes GENERICAS ********/

/******************/
/*** Operator = ***/

TPZSFMatrix &
TPZSFMatrix::operator=(const TPZMatrix &A )
{
	int newDim = Min( A.Rows(), A.Cols() );
	int size   = newDim * (newDim + 1) / 2;
	
	if ( newDim != Dim() )
    {
		if ( fElem != NULL )
			delete( fElem );
		fElem = new REAL[ size ] ;
    }
	
	// Copia a matriz.
	REAL *dst = fElem;
	for ( int c = 0; c < newDim; c++ )
		for ( int r = 0; r <= c; r++ )
			*dst++ = A.Get( r, c );
	
	fRow = fCol = newDim;
	
	// Ajusta Status de decomposicao.
	
	SetIsDecomposed(A.IsDecomposed());
	fDefPositive = (char)A.IsDefPositive();
	
	return( *this );
}



/******************/
/*** Operator + ***/

TPZSFMatrix
TPZSFMatrix::operator+(const TPZMatrix &A ) const
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Operator+ <matrixs with incoompatible dimensions" );
	
	TPZSFMatrix res( Dim() );
	REAL *pm = fElem;
	REAL *pr = res.fElem;
	for ( int c = 0; c < Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pr++ = (*pm++) + A.Get( r, c );
	
	return( *this );
}



/******************/
/*** Operator - ***/

TPZSFMatrix
TPZSFMatrix::operator-( const TPZMatrix &A ) const
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Operator- <matrixs with incoompatible dimensions" );
	
	TPZSFMatrix res( Dim() );
	REAL *pm = fElem;
	REAL *pr = res.fElem;
	for ( int c = 0; c < Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pr++ = (*pm++) - A.Get( r, c );
	
	return( *this );
}



/*******************/
/*** Operator += ***/

TPZSFMatrix &
TPZSFMatrix::operator+=(const TPZMatrix &A )
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Operator+= (TPZMatrix &) <different dimensions>" );
	
	REAL *pm = fElem;
	for ( int c = 0; c < Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pm++ += A.Get( r, c );
	
	fDecomposed = 0;
	return( *this );
}



/*******************/
/*** Operator -= ***/

TPZSFMatrix &
TPZSFMatrix::operator-=(const TPZMatrix &A )
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Operator-= (TPZMatrix &) <different dimensions>" );
	
	REAL *pm = fElem;
	for ( int c = 0; c < Dim(); c++ )
		for ( int r = 0; r <= c; r++ )
			*pm++ -= A.Get( r, c );
	
	fDecomposed = 0;
	return( *this );
}



/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/

TPZSFMatrix&
TPZSFMatrix::operator=(const REAL value )
{
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = value;
	fDecomposed = 0;
	fDefPositive = 0;
	return( *this );
}



/**************************/
/*** Operator+( value ) ***/

TPZSFMatrix
TPZSFMatrix::operator+(const REAL value ) const
{
	TPZSFMatrix res( Dim() );
	
	REAL *dst = res.fElem;
	REAL *src = fElem;
	REAL *end = &fElem[ Size() ];
	while ( src < end )
		*dst++ = (*src++) + value;
	
	return( res );
}



/**************************/
/*** Operator*( value ) ***/

TPZSFMatrix
TPZSFMatrix::operator*(const REAL value ) const
{
	TPZSFMatrix res( Dim() );
	
	REAL *dst = res.fElem;
	REAL *src = fElem;
	REAL *end = &fElem[ Size() ];
	while ( src < end )
		*dst++ = (*src++) * value;
	
	return( res );
}



/***************************/
/*** Operator+=( value ) ***/

TPZSFMatrix &
TPZSFMatrix::operator+=( REAL value )
{
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ += value;
	fDecomposed = 0;
	return( *this );
}



/***************************/
/*** Operator*=( value ) ***/

TPZSFMatrix &
TPZSFMatrix::operator*=( REAL value )
{
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ *= value;
	fDecomposed = 0;
	return( *this );
}



/**************/
/*** Resize ***/

int
TPZSFMatrix::Resize( int newDim , int )
{
	if ( newDim == Dim() )
		return( 1 );
	
	int newSize = newDim * (newDim + 1) / 2;
	int oldSize = Size();
	REAL *newElem = new REAL[newSize] ;
	if ( newElem == NULL )
		return TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	int minSize  = Min( newSize, oldSize );
	REAL *src = fElem;
	REAL *dst = newElem;
	REAL *end = &fElem[ minSize ];
	
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
	fRow = fCol = newDim;
	fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/

int
TPZSFMatrix::Redim( int newDim , int)
{
	// Se for preciso, desaloca a matriz antiga e aloca uma
	//  nova com o novo tamanho.
	if ( newDim != Dim() )
    {
		fRow = fCol = newDim;
		if ( fElem != NULL )
			delete( fElem );
		fElem = new REAL[Size()] ;
    }
	
	// Zera a matriz.
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = 0.;
	fDecomposed = 0;
	fDefPositive = 0;
	
	
	return( 1 );
}

int
TPZSFMatrix::Zero()
{
	REAL *dst = fElem;
	REAL *end = &fElem[ Size() ];
	while ( dst < end )
		*dst++ = 0.;
	fDecomposed = 0;
	fDefPositive = 0;
	return( 1 );
}
/******** Resolucao de Sistemas ********/

/**************************/
/*** Decompose Cholesky ***/
int
TPZSFMatrix::Decompose_Cholesky(std::list<int> &singular)
{
	return Decompose_Cholesky();
}

int
TPZSFMatrix::Decompose_Cholesky()
{
	if (  fDecomposed )  TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	REAL *ptr_k = fElem;
	for ( int k = 0; k < Dim(); k++, ptr_k += k  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		REAL sum  = 0.0;
		REAL *pk  = ptr_k;
		REAL *pkk = ptr_k + k;
		for ( ; pk < pkk;  pk++ )
			sum += (*pk) * (*pk);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		if ( (*pk -= sum) < 1.e-10 )
			return( 0 );
		*pk = sqrt( *pk );
		
		// Loop para i = k+1 ... Dim().
		//
		REAL *ptr_i = ptr_k;
		REAL *pi;
		for ( int i = k+1; i < Dim(); i++ )
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
	
	fDecomposed  = ECholesky;
	fDefPositive = 1;
	return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
int
TPZSFMatrix::Decompose_LDLt(std::list<int> &singular)
{
	return Decompose_LDLt();
}


int
TPZSFMatrix::Decompose_LDLt()
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
	if (  fDecomposed && fDecomposed != ELDLt)
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Matrix already Decomposed with a different scheme>" );
	else if ( fDecomposed ) return 0;
	
	int j,k,l;
	//	REAL sum;
	
	
	for ( j = 0; j < Dim(); j++ )
    {
		//	sum=0.;
		for ( k=0; k<j; k++)
			//			sum=sum-GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
			PutVal( j,j,GetVal(j,j)-GetVal(k,k)*GetVal(k,j)*GetVal(k,j) );
		//		PutVal(j,j,GetVal(j,j)+sum);
		for ( k=0; k<j; k++)
		{
			
			for( l=j+1; l<Dim();l++)
				PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
		}
		
		if ( IsZero(GetVal(j,j)) )TPZMatrix::Error(__PRETTY_FUNCTION__,"Decompose_LDLt <Zero on diagonal>" );
		
		for( l=j+1; l<Dim();l++) PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
		
    }
	fDecomposed  = ELDLt;
	fDefPositive = 0;
	
	return( 1 );
}



/*********************/
/*** Subst Forward ***/

int
TPZSFMatrix::Subst_Forward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Subst_Forward <the matrix result can not be simetric>" );
	
	REAL *ptr_k = fElem;
	for ( int k = 0; k < Dim(); k++ )
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

int
TPZSFMatrix::Subst_Backward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Subst_Backward <the matrix result can not be simetric>" );
	
	REAL *ptr_k = &fElem[ Size()-1 ];
	for ( int k = Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			REAL sum = 0.0;
			REAL *pk = ptr_k;
			for ( int i = Dim()-1; i > k; i-- )
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

int
TPZSFMatrix::Subst_LForward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Subst_LForward <the matrix result can not be simetric>" );
	
	REAL *ptr_k = fElem;
	for ( int k = 0; k < Dim(); k++, ptr_k += k )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			//
			REAL *pk = ptr_k;
			REAL sum = 0.0;
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

int
TPZSFMatrix::Subst_LBackward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		return( 0 );
	
	if ( B->IsSimetric() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Subst_LBackward <the matrix result can not be simetric>" );
	
	REAL *ptr_k = &fElem[ Size()-1 ];
	for ( int k = Dim()-1; k >= 0; k--, ptr_k-- )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ); i = N, ..., k+1.
			//
			REAL sum = 0.0;
			REAL *pk = ptr_k;
			for ( int i = Dim()-1; i > k; i-- )
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

int
TPZSFMatrix::Subst_Diag( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		return( 0 );
	
	REAL *pDiag = fElem;
	for ( int k = 0; k < Dim(); k++ )
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

int
TPZSFMatrix::Clear()
{
	if ( fElem != NULL )
		delete( fElem );
	
	fElem = NULL;
	fRow = fCol = 0;
	return( 1 );
}

#ifdef OOPARLIB

int TPZSFMatrix::Unpack( TReceiveStorage *buf ){
	TSaveable::Unpack(buf);
	buf->UpkDouble(fElem);
	return 1;
}



TSaveable *TPZSFMatrix::Restore(TReceiveStorage *buf) {
	TPZSFMatrix *m = new TPZSFMatrix();
	m->Unpack(buf);
	return m;
}

int TPZSFMatrix::Pack( TSendStorage *buf ){
	TSaveable::Pack(buf);
	buf->PkDouble(fElem);
	return 1;
}


int TPZSFMatrix::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSaveable::DerivedFrom(Classid);
}

int TPZSFMatrix::DerivedFrom(char *classname){
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TSaveable::DerivedFrom(classname);
}

#endif

