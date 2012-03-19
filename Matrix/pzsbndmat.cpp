/**
 * @file
 * @brief Contains the implementation of the TPZSBMatrix methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tsbndmat.h
//
// Class:  TPZSBMatrix
//
// Obs.:   Esta classe gerencia matrizes do tipo Band simetrica.
//
// Versao: 10 / 1996.
//

#include <math.h>
#include <stdlib.h>

#include "pzfmatrix.h"
#include "pzsbndmat.h"
//#include "pzerror.h" 

#include <sstream>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzsbmatrix"));
#endif

using namespace std;

/*******************/
/*** TPZSBMatrix ***/

/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

template<class TVar>
TPZSBMatrix<TVar>::TPZSBMatrix( int dim, int band )
: TPZMatrix<TVar>( dim, dim )
{
	fBand = ( band > (dim - 1) ? (dim - 1) : band );
	fDiag = new( TVar[Size()] );
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrix( dim ) <Error creating Matrix>" );
	
	Zero();
}

/**************/
/*** PutVal ***/

template <class TVar>
int
TPZSBMatrix<TVar>::PutVal(const int r,const int c,const TVar& value )
{
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	int index;
	if ( (index = col-row) > fBand )
	{
#ifdef DEBUG
		if (value) {
			DebugStop();
		}
#endif
		return( 0 );        // O elemento esta fora da banda.
	}
	fDiag[ col * (fBand+1) + index ] = value;
	this->fDecomposed = 0;
	return( 1 );
}



/**************/
/*** GetVal ***/
template<class TVar>
const TVar
&TPZSBMatrix<TVar>::GetVal(const int r,const int c ) const
{
	
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		this->Swap( &row, &col );
	
	int index;
	this->gZero=0.0;
	if ( (index = col-row) > fBand )
		return( this->gZero );        // O elemento esta fora da banda.
	
	return( fDiag[ col * (fBand+1) + index ] );
}


/*************/
/*** Print ***/

template<class TVar>
void
TPZSBMatrix<TVar> ::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	out.width( 8 );
	out.precision( 4 );
	
	out << "Writing matrix '" << name;
	out << "' (" << this->Rows() << " x " << this->Cols() << ")  Bandwith = "<<fBand<<"\n";
	TPZMatrix<TVar>::Print(0,out,form);
	/*
	 for ( int row = 0; row < Rows(); row++)
	 {
	 out << "\t";
	 for ( int col = 0; col < Cols(); col++ )
	 out << Get( row, col) << "  ";
	 out << "\n";
	 }
	 
	 out << "\n";
	 */
}


/** @brief Overload << operator to output entries of TPZSBMatrix matrix ***/
template<class TVar>
std::ostream&
operator<<(std::ostream& out,TPZSBMatrix<TVar>  &A)
{
	out.width( 8 );
	out.precision( 4 );
	
	out <<"\n(" << A.Rows() << " x " << A.Cols()
	<< ")  Bandwith = "<< A.GetBand()<<"\n";
	
	for ( int row = 0; row < A.Rows(); row++)
    {
		out << "\t";
		for ( int col = 0; col < A.Cols(); col++ )
			out << A.Get( row, col) << "  ";
		out << "\n";
    }
	
	return  out << "\n";
}



/******** Operacoes com matrizes BANDA SIMETRICA  ********/

/******************/
/*** Operator = ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator=(const TPZSBMatrix<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}



/******************/
/*** Operator + ***/

template<class TVar>
TPZSBMatrix<TVar>
TPZSBMatrix<TVar>::operator+(const TPZSBMatrix<TVar> &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
	
	// Define os ponteiros e tamanhos para os elementos da maior e da
	//  menor banda.
	TVar *elemMax;
	TVar *elemMin;
	int  sizeMax;
	int  sizeMin;
	if ( fBand > A.fBand )
    {
		sizeMax = fBand + 1;
		elemMax = fDiag;
		sizeMin = A.fBand + 1;
		elemMin = A.fDiag;
    }
	else
    {
		sizeMax = A.fBand + 1;
		elemMax = A.fDiag;
		sizeMin = fBand + 1;
		elemMin = fDiag;
    }
	
	TPZSBMatrix<TVar> res( this->Dim(), sizeMax-1 );
	
	// Efetua a SOMA.
	TVar *dest = res.fDiag;
	int col,i;
	for ( col = 0; col < this->Dim(); col++ )
    {
		for (  i = 0; i < sizeMin; i++ )
			*dest++ = (*elemMax++) + (*elemMin++);
		for ( ; i < sizeMax; i++ )
			*dest++ = *elemMax++;
    }
	
	return( res );
}



/******************/
/*** Operator - ***/

template<class TVar>
TPZSBMatrix<TVar>
TPZSBMatrix<TVar>::operator-(const TPZSBMatrix<TVar> &A ) const
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSBMatrix ) <incompatible dimensions>" );
	
	// Define o tamanho e os elementos das colunas das 2 matrizes.
	int  sizeThis  = fBand + 1;
	int  sizeA     = A.fBand + 1;
	TVar *elemThis = fDiag;
	TVar *elemA    = A.fDiag;
	
	TPZSBMatrix<TVar> res( this->Dim(), MAX(sizeThis, sizeA) );
	
	// Efetua a SUBTRACAO.
	TVar *dest = res.fDiag;
	int col,i;
	for (col = 0; col < this->fRow; col++ )
    {
		for (i = 0; (i < sizeThis) && (i < sizeA); i++ )
			*dest++ = (*elemThis++) - (*elemA++);
		if ( i == sizeA )
			for ( ; i < sizeThis; i++ )
				*dest++ = *elemThis++;
		else
			for ( ; i < sizeA; i++ )
				*dest++ = -(*elemA++);
    }
	
	return( res );
}



/*******************/
/*** Operator += ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator+=(const TPZSBMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+=( TPZSBMatrix ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais).
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		TVar *pm = fDiag;
		TVar *pa = A.fDiag;
		TVar *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ += *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		TVar *pThis  = fDiag;
		TVar *pA     = A.fDiag;
		int sizeThis = fBand + 1;
		int sizeA    = A.fBand + 1;
		int inc = sizeThis - sizeA;
		for ( int col = 0; col < this->Dim(); col++, pThis += inc )
			for ( int i = 0; i < sizeA; i++ )
				*pThis++ += *pA++;
    }
	
	this->fDecomposed = 0;
	return( *this );
}



/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator-=(const TPZSBMatrix<TVar> &A )
{
	if ( this->Dim() != A.Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-=( TPZSBMatrix ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais)
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		TVar *pm = fDiag;
		TVar *pa = A.fDiag;
		TVar *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ -= *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		TVar *pThis  = fDiag;
		TVar *pA     = A.fDiag;
		int sizeThis = fBand + 1;
		int sizeA    = A.fBand + 1;
		int inc = sizeThis - sizeA;
		for ( int col = 0; col < this->Dim(); col++, pThis += inc )
			for ( int i = 0; i < sizeA; i++ )
				*pThis++ -= *pA++;
    }
	
	this->fDecomposed = 0;
	return( *this );
}

template<class TVar>
void TPZSBMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						  const REAL alpha,const REAL beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && this->Cols()*stride != x.Rows()) || this->Rows()*stride != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSBMatrix::MultAdd incompatible dimensions\n");
	}
	PrepareZ(y,z,beta,opt,stride);
	int rows = this->Rows();
	int xcols = x.Cols();
	int ic, r;
	for (ic = 0; ic < xcols; ic++) {
		int begin, end;
		for ( r = 0; r < rows; r++ ) {
			begin = MAX( r - fBand, 0 );
			end   = MIN( r + fBand + 1, rows );
			TVar val = 0.;
			// Calcula um elemento da resposta.
			for ( int i = begin ; i < end; i++ ) val += GetVal( r, i ) * x.GetVal(stride * i, ic );
			val *= alpha;
			val += z.GetVal(r*stride,ic);
			z.PutVal( r*stride , ic, val );
		}
	}
}


/******** Operacoes com MATRIZES GENERICAS ********/


// Estas operacoes com matrizes genericas, usam a parte triangular
// inferior da maior matriz quadrada de A. Ex.:
//
//  Se A = 01 02 03 04   A matriz usada sera':  01 05 09
//         05 06 07 08                          05 06 10
//         09 10 11 12                          09 10 11
//


/******** Operacoes com valores NUMERICOS ********/
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//



/*****************************/
/*** Operator * ( REAL ) ***/

template<class TVar>
TPZSBMatrix<TVar>
TPZSBMatrix<TVar>::operator*(const TVar value ) const
{
	TPZSBMatrix<TVar> res( this->Dim(), fBand );
	
	TVar *pr  = res.fDiag;
	TVar *pm  = fDiag;
	TVar *end = &fDiag[Size()];
	while ( pm < end )
		*pr++ = (*pm++) * value;
	return( res );
}



/******************************/
/*** Operator += ( REAL ) ***/


/******************************/
/*** Operator *= ( REAL ) ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator*=(const TVar value )
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
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
template<class TVar>
int
TPZSBMatrix<TVar>::Resize(const int newDim ,const int)
{
	if ( newDim == this->Dim() )
		return( 1 );
	
	if (fBand>this->Dim()-1) fBand=this->Dim()-1;//misael:19/10/96
	// Cria nova matrix.
	int  newSize  = newDim * (fBand + 1);
	TVar *newDiag = new TVar[newSize] ;
	
	// Copia os elementos para a nova matriz.
	TVar *src = fDiag;
	TVar *dst = newDiag;
	TVar *end = &newDiag[ MIN(Size(), newSize) ];
	while ( dst < end )
		*dst++ = *src++;
	
	// Zera as posicoes que sobrarem (se sobrarem).
	end = &newDiag[newSize];
	while ( dst < end )
		*dst++ = 0.0;
	
	// Descarta a matriz antiga e valida a nova matriz.
	if ( fDiag != NULL )
		delete( fDiag );
	fDiag = newDiag;
	this->fCol = this->fRow = newDim;
	this->fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSBMatrix<TVar>::Redim(const int newDim ,const int)
{
	if ( newDim != this->Dim() )
    {
		// Aloca a nova matriz.
		this->fRow = this->fCol = newDim;
		if ( this->fDiag != NULL )
			delete( this->fDiag );
		this->fDiag = new( TVar[Size()] );
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

template<class TVar>
int
TPZSBMatrix<TVar>::Zero()
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
template<class TVar>
int
TPZSBMatrix<TVar>::SetBand(const int newBand )
{
	if ( this->fBand == newBand )
		return( 1 );
	
	if ( this->fBand > (this->Dim() - 1) )
		return( 0 );
	
	TVar *newDiag = new( TVar[this->Dim() * (newBand + 1)] );
	
	// Copia os elementos antigos para a nova alocacao.
	REAL *pNew = newDiag;
	REAL *pOld = fDiag;
	int newSize  = newBand + 1;
	int oldSize  = fBand + 1;
	int minSize  = MIN( newSize, oldSize );
	int i;
	for ( i = 0; i < minSize; i++ )
    {
		TVar *end = pNew + (this->Dim() * newSize);
		TVar *dst = pNew++;
		TVar *src = pOld++;
		for ( ; dst < end; src += oldSize, dst += newSize )
			*dst = *src;
    }
	
	// Preenche com zero os elementos que sobrarem (se sobrarem).
	for ( ; i < newSize; i++ )
    {
		TVar *end = pNew + (this->Dim() * newSize);
		TVar *dst = pNew++;
		for ( ; dst < end; dst += newSize )
			*dst = 0.0;
    }
	
	if ( fDiag == NULL )
		delete( fDiag );
	fDiag = newDiag;
	
	if ( newBand < fBand )
		this->fDecomposed = 0;
	fBand = newBand;
	
	return( 1 );
}



/********************* Resolucao de sistemas *********************/

/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky(std::list<int> &singular)
{
	return Decompose_Cholesky();
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky()
{
	if (  this->fDecomposed && this->fDecomposed != ECholesky )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (this->fDecomposed) {
		return 1;
	}
	int size = fBand + 1; // Tamanho da banda.
	int next = fBand + 1; // O quanto se soma para "andar" na diagonal.
	TVar *row_k   = fDiag;
	TVar *row_end = fDiag + Size();
	
	for ( int k = 0; k <this-> Dim(); k++, row_k += next  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		TVar sum = 0.0;
		TVar *elem_k = row_k + 1;
		TVar *end_k  = row_k + next;
		for ( ; elem_k < end_k; elem_k++ )
			sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		if ( (*row_k -= sum) < 1.e-10 )
			return( 0 );
		*row_k = sqrt( *row_k );
		
		// Loop para i = k+1 ... Dim().
		//
		TVar *row_i = row_k + next;
		for ( int j = 2; row_i < row_end; j++, row_i += next )
		{
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( j < size )
			{
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				TVar *elem_i = row_i + j;
				TVar *end_i  = row_i + size;
				elem_k = row_k + 1;
				while ( elem_i < end_i )
					sum += (*elem_i++) * (*elem_k++);
				
				// Faz A(i,k) = (A(i,k) - sum) / A(k,k)
				row_i[j-1] = (row_i[j-1] -sum) / (*row_k);
			}
			
			// Se nao tiverem estes elementos, sum = 0.0.
			else if ( j == size )
				row_i[j-1] /= (*row_k);
			
			// Se nao existir nem o elemento A(i,k), nao faz nada.
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
TPZSBMatrix<TVar>::Decompose_LDLt(std::list<int> &singular)
{
	return Decompose_LDLt();
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_LDLt()
{
	
	if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
	int j,k,l, begin,end;
	TVar sum;
	
	
	for ( j = 0; j < this->Dim(); j++ )
    {
		//Print("curernt");
		sum=0.;
		
		begin = MAX( int(j - fBand), 0 );
		//cout<<"begin="<<begin<<"\n";
		for ( k=begin; k<j; k++)
		{
			sum=sum-GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
			//cout<<"(k,j)"<<k<<" "<<j<<"\n";
		}
		
		
		//	 operator()(j,j)=GetVal(j,j)+sum;
		PutVal(j,j,GetVal(j,j)+sum);
		//cout<<"\n(j,j)"<<j<<" "<<j<<"\n\n";
		for ( k=0; k<j; k++)
		{
			end   = MIN( int(k + fBand )+1, this->Dim() );
			for( l=j+1; l<end;l++)
			{
				PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
				/*cout<<"end="<<end<<"\n";
				 cout<<"(l,j)"<<l<<" "<<j<<"\n";
				 cout<<"(j,k)"<<j<<" "<<k<<"\n";
				 cout<<"(l,k)"<<l<<" "<<k<<"\n\n";
				 */
			}
		}
		
		if ( IsZero(GetVal(j,j)) ) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Zero on diagonal>" );
		end  = MIN( int(j + fBand )+1, this->Dim() );
		//cout<<"end="<<end<<"\n";
		for( l=j+1; l<end;l++)
		{
			//cout<<"(l,j)"<<l<<" "<<j<<"\n";
			PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
		}
    }
	this->fDecomposed  = 1;
	this->fDefPositive = 0;
	
	return( 1 );
	
}



/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
template<class TVar>
int
TPZSBMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	
	int size = fBand + 1;
	TVar *row_k = fDiag;
	for ( int k = 0; k < this->Dim(); k++, row_k += size )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			int end=(k-fBand>0)? fBand:k;//misael
			TVar sum = 0.0;
			TVar *elem_ki = row_k + 1;
			TVar *end_ki  = row_k + end;  //misael
			for ( int i = k-1; elem_ki < end_ki ; i-- )
				sum += (*elem_ki++) * B->GetVal( i, j );
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	
	return( 1 );
}
template<class TVar>
int
TPZSBMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	/*  int i,j,k;
	 
	 for(k=0; k<B->Cols() ; k++){
	 for(i=Rows()-1; i>=0; i--){
	 short jmax;
	 jmax=( i+fBand>Rows()-1 )? Rows()-1 : i+fBand;
	 for(j=jmax ; j>i ; j--) B->operator()(i,k)-=B->GetVal(j,k)*GetVal(i,j);
	 if (  IsZero(GetVal(i,i)) ) TPZMatrix::Error(__PRETTY_FUNCTION__,"Suubst_Backward->diagonal element = zero");
	 B->operator()(i,k)/=GetVal(i,i);
	 
	 }
	 }
	 
	 return 1;
	 */
	int k,j,i,jmax,stepcol=fBand+2;
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
template<class TVar>
int
TPZSBMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
	
	int size = fBand + 1;
	TVar *row_k = fDiag;
	int i,j,k;
	for ( k = 0; k < this->Dim(); k++, row_k += size )
		for ( j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			
			int end=(k-fBand>0)? fBand:k;  //misael
			TVar sum = 0.0;
			TVar *elem_ki = row_k + 1;
			TVar *end_ki  = row_k + end;
			for ( i = k-1; elem_ki <= end_ki ; i-- )//misael
				sum += (*elem_ki++) * B->GetVal( i, j );
			
			// Faz b[k] = (b[k] - sum).
			//
			B->PutVal( k, j, (B->GetVal( k, j ) - sum) );
			
		}
	
	return( 1 );
}



/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int
TPZSBMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
	
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
	
	
	int size = fBand + 1;
	TVar *row_k = fDiag;
	for ( int k = 0; k < this->Dim(); k++, row_k += size )
		for ( int j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j) / *row_k );
	
	return( 1 );
}

template<class TVar>
int //misael 19/10/96
TPZSBMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
	if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
	
	int k,j,i,jmax,stepcol=fBand+2;
	
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
/*** Error ***/
/*int
 TPZSBMatrix::Error(const char *msg1,const char* msg2 ) 
 {
 ostringstream out;
 out << "TPZSBMatrix::" << msg1 << msg2 << ".\n";
 //pzerror.show();
 LOGPZ_ERROR (logger, out.str().c_str());
 DebugStop();
 return 1;
 }*/



/*************/
/*** CLear ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Clear()
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

template<class TVar>
void
TPZSBMatrix<TVar>::Copy(const TPZSBMatrix<TVar> &A )
{
	this->fBand = A.fBand;
	this->fRow = this->fCol = A.Dim();
	this->fDiag = new( TVar[Size()] );
	this->fDecomposed  = A.fDecomposed;
	this->fDefPositive = A.fDefPositive;
	
	if ( fDiag == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Copy( TPZSBMatrix ) <memory allocation error>" );
	
	TVar *dst = fDiag;
	TVar *src = A.fDiag;
	TVar *end = A.fDiag + Size();
	while ( src < end )
		*dst++ = *src++;
}

#ifdef OOPARLIB

template<class TVar>
int TPZSBMatrix<TVar>::Unpack (TReceiveStorage *buf ){
	TSimMatrix::Unpack(buf);
	buf->UpkInt(&fBand);
	buf->UpkDouble(fDiag,fBand);
	return 1;
}


template<class TVar>
TSaveable *TPZSBMatrix<TVar>::Restore(TReceiveStorage *buf) {
	TPZSBMatrix<TVar> *m = new TPZSBMatrix<TVar>();
	m->Unpack(buf);
	return m;
}

template<class TVar>
int TPZSBMatrix<TVar>::Pack( TSendStorage *buf ){
	TSimMatrix::Pack(buf);
	buf->PkInt(&fBand);
	buf->PkDouble(fDiag,fBand);
	return 1;
}

template<class TVar>
int TPZSBMatrix<TVar>::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSimMatrix::DerivedFrom(Classid);
}

template<class TVar>
int TPZSBMatrix<TVar>::DerivedFrom(char *classname){
	if(!strcmp(ClassName(),classname)) return 1;
	return TSimMatrix::DerivedFrom(classname);
}

#endif

