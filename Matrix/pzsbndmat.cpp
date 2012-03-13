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

TPZSBMatrix::TPZSBMatrix( int dim, int band )
: TPZMatrix( dim, dim )
{
	fBand = ( band > (dim - 1) ? (dim - 1) : band );
	fDiag = new REAL[Size()] ;
	if ( fDiag == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "TPZSBMatrix( dim ) <Error creating Matrix>" );
	
	Zero();
}

/**************/
/*** PutVal ***/

int
TPZSBMatrix::PutVal(const int r,const int c,const REAL& value )
{
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		Swap( &row, &col );
	
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
	fDecomposed = 0;
	return( 1 );
}



/**************/
/*** GetVal ***/

const REAL
&TPZSBMatrix::GetVal(const int r,const int c ) const
{
	
	// inicializando row e col para trabalhar com a triangular superior
	int row(r),col(c);
	if ( row > col )
		Swap( &row, &col );
	
	int index;
	gZero=0.0;
	if ( (index = col-row) > fBand )
		return( gZero );        // O elemento esta fora da banda.
	
	return( fDiag[ col * (fBand+1) + index ] );
}


/*************/
/*** Print ***/
void
TPZSBMatrix ::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	out.width( 8 );
	out.precision( 4 );
	
	out << "Writing matrix '" << name;
	out << "' (" << Rows() << " x " << Cols() << ")  Bandwith = "<<fBand<<"\n";
	TPZMatrix::Print(0,out,form);
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
std::ostream&
operator<<(std::ostream& out,TPZSBMatrix  &A)
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

TPZSBMatrix &
TPZSBMatrix::operator=(const TPZSBMatrix &A )
{
	Clear();
	Copy( A );
	return( *this );
}



/******************/
/*** Operator + ***/

TPZSBMatrix
TPZSBMatrix::operator+(const TPZSBMatrix &A ) const
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
	
	// Define os ponteiros e tamanhos para os elementos da maior e da
	//  menor banda.
	REAL *elemMax;
	REAL *elemMin;
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
	
	TPZSBMatrix res( Dim(), sizeMax-1 );
	
	// Efetua a SOMA.
	REAL *dest = res.fDiag;
	int col,i;
	for ( col = 0; col < Dim(); col++ )
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

TPZSBMatrix
TPZSBMatrix::operator-(const TPZSBMatrix &A ) const
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "operator-( TPZSBMatrix ) <incompatible dimensions>" );
	
	// Define o tamanho e os elementos das colunas das 2 matrizes.
	int  sizeThis  = fBand + 1;
	int  sizeA     = A.fBand + 1;
	REAL *elemThis = fDiag;
	REAL *elemA    = A.fDiag;
	
	TPZSBMatrix res( Dim(), MAX(sizeThis, sizeA) );
	
	// Efetua a SUBTRACAO.
	REAL *dest = res.fDiag;
	int col,i;
	for (col = 0; col < fRow; col++ )
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

TPZSBMatrix &
TPZSBMatrix::operator+=(const TPZSBMatrix &A )
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "operator+=( TPZSBMatrix ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais).
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		REAL *pm = fDiag;
		REAL *pa = A.fDiag;
		REAL *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ += *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		REAL *pThis  = fDiag;
		REAL *pA     = A.fDiag;
		int sizeThis = fBand + 1;
		int sizeA    = A.fBand + 1;
		int inc = sizeThis - sizeA;
		for ( int col = 0; col < Dim(); col++, pThis += inc )
			for ( int i = 0; i < sizeA; i++ )
				*pThis++ += *pA++;
    }
	
	fDecomposed = 0;
	return( *this );
}



/*******************/
/*** Operator -= ***/

TPZSBMatrix &
TPZSBMatrix::operator-=(const TPZSBMatrix &A )
{
	if ( Dim() != A.Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "operator-=( TPZSBMatrix ) <incompatible dimensions>" );
	
	// No caso de as bandas serem iguais (ou se "tornarem" iguais)
	if ( fBand <= A.fBand )
    {
		if ( fBand < A.fBand )
			SetBand( A.fBand );
		REAL *pm = fDiag;
		REAL *pa = A.fDiag;
		REAL *end = &fDiag[ Size() ];
		while ( pm < end )
			*pm++ -= *pa++;
    }
	else
    {
		// Se a banda desta matriz for maior...
		REAL *pThis  = fDiag;
		REAL *pA     = A.fDiag;
		int sizeThis = fBand + 1;
		int sizeA    = A.fBand + 1;
		int inc = sizeThis - sizeA;
		for ( int col = 0; col < Dim(); col++, pThis += inc )
			for ( int i = 0; i < sizeA; i++ )
				*pThis++ -= *pA++;
    }
	
	fDecomposed = 0;
	return( *this );
}

void TPZSBMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						  const REAL alpha,const REAL beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && Cols()*stride != x.Rows()) || Rows()*stride != x.Rows())
		TPZMatrix::Error(__PRETTY_FUNCTION__, "TPZSBMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix::Error(__PRETTY_FUNCTION__,"TPZSBMatrix::MultAdd incompatible dimensions\n");
	}
	PrepareZ(y,z,beta,opt,stride);
	int rows = Rows();
	int xcols = x.Cols();
	int ic, r;
	for (ic = 0; ic < xcols; ic++) {
		int begin, end;
		for ( r = 0; r < rows; r++ ) {
			begin = MAX( r - fBand, 0 );
			end   = MIN( r + fBand + 1, rows );
			REAL val = 0.;
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

TPZSBMatrix
TPZSBMatrix::operator*(const REAL value ) const
{
	TPZSBMatrix res( Dim(), fBand );
	
	REAL *pr  = res.fDiag;
	REAL *pm  = fDiag;
	REAL *end = &fDiag[Size()];
	while ( pm < end )
		*pr++ = (*pm++) * value;
	return( res );
}



/******************************/
/*** Operator += ( REAL ) ***/


/******************************/
/*** Operator *= ( REAL ) ***/

TPZSBMatrix &
TPZSBMatrix::operator*=(const REAL value )
{
	REAL *pm  = fDiag;
	REAL *end = &fDiag[Size()];
	while ( pm < end )
		*pm++ *= value;
	
	fDecomposed = 0;
	return( *this );
}



/**************/
/*** Resize ***/
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//
int
TPZSBMatrix::Resize(const int newDim ,const int)
{
	if ( newDim == Dim() )
		return( 1 );
	
	if (fBand>Dim()-1) fBand=Dim()-1;//misael:19/10/96
	// Cria nova matrix.
	int  newSize  = newDim * (fBand + 1);
	REAL *newDiag = new REAL[newSize] ;
	
	// Copia os elementos para a nova matriz.
	REAL *src = fDiag;
	REAL *dst = newDiag;
	REAL *end = &newDiag[ MIN(Size(), newSize) ];
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
	fCol = fRow = newDim;
	fDecomposed = 0;
	return( 1 );
}



/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
int
TPZSBMatrix::Redim(const int newDim ,const int)
{
	if ( newDim != Dim() )
    {
		// Aloca a nova matriz.
		fRow = fCol = newDim;
		if ( fDiag != NULL )
			delete( fDiag );
		fDiag = new REAL[Size()] ;
		if ( fDiag == NULL )
			TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>" );
    }
	
	REAL *dst = fDiag;
	REAL *end = &fDiag[Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	fDecomposed = 0;
	fDefPositive = 0;
	return( 1 );
}

int
TPZSBMatrix::Zero()
{
	REAL *dst = fDiag;
	REAL *end = &fDiag[Size()];
	while ( dst < end )
		*dst++ = 0.;
	
	fDecomposed = 0;
	fDefPositive = 0;
	return( 1 );
}


/****************/
/*** Set Band ***/
int
TPZSBMatrix::SetBand(const int newBand )
{
	if ( fBand == newBand )
		return( 1 );
	
	if ( fBand > (Dim() - 1) )
		return( 0 );
	
	REAL *newDiag = new REAL[Dim() * (newBand + 1)] ;
	
	// Copia os elementos antigos para a nova alocacao.
	REAL *pNew = newDiag;
	REAL *pOld = fDiag;
	int newSize  = newBand + 1;
	int oldSize  = fBand + 1;
	int minSize  = MIN( newSize, oldSize );
	int i;
	for ( i = 0; i < minSize; i++ )
    {
		REAL *end = pNew + (Dim() * newSize);
		REAL *dst = pNew++;
		REAL *src = pOld++;
		for ( ; dst < end; src += oldSize, dst += newSize )
			*dst = *src;
    }
	
	// Preenche com zero os elementos que sobrarem (se sobrarem).
	for ( ; i < newSize; i++ )
    {
		REAL *end = pNew + (Dim() * newSize);
		REAL *dst = pNew++;
		for ( ; dst < end; dst += newSize )
			*dst = 0.0;
    }
	
	if ( fDiag == NULL )
		delete( fDiag );
	fDiag = newDiag;
	
	if ( newBand < fBand )
		fDecomposed = 0;
	fBand = newBand;
	
	return( 1 );
}



/********************* Resolucao de sistemas *********************/

/**************************/
/*** Decompose Cholesky ***/
int
TPZSBMatrix::Decompose_Cholesky(std::list<int> &singular)
{
	return Decompose_Cholesky();
}

int
TPZSBMatrix::Decompose_Cholesky()
{
	if (  fDecomposed && fDecomposed != ECholesky )  TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_Cholesky <Matrix already Decomposed>" );
	if (fDecomposed) {
		return 1;
	}
	int size = fBand + 1; // Tamanho da banda.
	int next = fBand + 1; // O quanto se soma para "andar" na diagonal.
	REAL *row_k   = fDiag;
	REAL *row_end = fDiag + Size();
	
	for ( int k = 0; k < Dim(); k++, row_k += next  )
    {
		// Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
		//
		REAL sum = 0.0;
		REAL *elem_k = row_k + 1;
		REAL *end_k  = row_k + next;
		for ( ; elem_k < end_k; elem_k++ )
			sum += (*elem_k) * (*elem_k);
		
		// Faz A(k,k) = sqrt( A(k,k) - sum ).
		//
		if ( (*row_k -= sum) < 1.e-10 )
			return( 0 );
		*row_k = sqrt( *row_k );
		
		// Loop para i = k+1 ... Dim().
		//
		REAL *row_i = row_k + next;
		for ( int j = 2; row_i < row_end; j++, row_i += next )
		{
			// Se tiverem elementos na linha 'i' cuja coluna e'
			//  menor do que 'K'...
			if ( j < size )
			{
				// Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
				sum = 0.0;
				REAL *elem_i = row_i + j;
				REAL *end_i  = row_i + size;
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
	
	fDecomposed  = ECholesky;
	fDefPositive = 1;
	return( 1 );
}



/**********************/
/*** Decompose LDLt ***/
int
TPZSBMatrix::Decompose_LDLt(std::list<int> &singular)
{
	return Decompose_LDLt();
}

int
TPZSBMatrix::Decompose_LDLt()
{
	
	if (  fDecomposed )  TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
	
	int j,k,l, begin,end;
	REAL sum;
	
	
	for ( j = 0; j < Dim(); j++ )
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
			end   = MIN( int(k + fBand )+1, Dim() );
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
		
		if ( IsZero(GetVal(j,j)) ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Zero on diagonal>" );
		end  = MIN( int(j + fBand )+1, Dim() );
		//cout<<"end="<<end<<"\n";
		for( l=j+1; l<end;l++)
		{
			//cout<<"(l,j)"<<l<<" "<<j<<"\n";
			PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
		}
    }
	fDecomposed  = 1;
	fDefPositive = 0;
	
	return( 1 );
	
}



/*********************/
/*** Subst Forward ***/
//
//  Faz Ax = b, onde A e' triangular inferior.
//
int
TPZSBMatrix::Subst_Forward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
	
	int size = fBand + 1;
	REAL *row_k = fDiag;
	for ( int k = 0; k < Dim(); k++, row_k += size )
		for ( int j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			int end=(k-fBand>0)? fBand:k;//misael
			REAL sum = 0.0;
			REAL *elem_ki = row_k + 1;
			REAL *end_ki  = row_k + end;  //misael
			for ( int i = k-1; elem_ki < end_ki ; i-- )
				sum += (*elem_ki++) * B->GetVal( i, j );
			
			// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
			//
			B->PutVal( k, j, (B->GetVal(k, j) - sum) / *row_k );
		}
	
	return( 1 );
}
int
TPZSBMatrix::Subst_Backward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
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
		for(i=Rows()-1; i>=0; i--)
		{
			REAL *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>Rows()-1)? Rows()-1 : i+fBand;
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
int
TPZSBMatrix::Subst_LForward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
	
	int size = fBand + 1;
	REAL *row_k = fDiag;
	int i,j,k;
	for ( k = 0; k < Dim(); k++, row_k += size )
		for ( j = 0; j < B->Cols(); j++ )
		{
			// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
			
			int end=(k-fBand>0)? fBand:k;  //misael
			REAL sum = 0.0;
			REAL *elem_ki = row_k + 1;
			REAL *end_ki  = row_k + end;
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
int
TPZSBMatrix::Subst_Diag( TPZFMatrix *B ) const
{
	
	if ( (B->Rows() != Dim()) || !fDecomposed )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
	
	
	int size = fBand + 1;
	REAL *row_k = fDiag;
	for ( int k = 0; k < Dim(); k++, row_k += size )
		for ( int j = 0; j < B->Cols(); j++ )
			B->PutVal( k, j, B->GetVal( k, j) / *row_k );
	
	return( 1 );
}

int //misael 19/10/96
TPZSBMatrix::Subst_LBackward( TPZFMatrix *B ) const
{
	if ( (B->Rows() != Dim()) || !fDecomposed )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
	
	int k,j,i,jmax,stepcol=fBand+2;
	
	for(k=0; k<B->Cols() ; k++)
    {
		for(i=Rows()-1; i>=0; i--)
		{
			REAL *diagk=fDiag+(fBand+1)*i;
			jmax=( (i+fBand)>Rows()-1)? Rows()-1 : i+fBand;
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

int
TPZSBMatrix::Clear()
{
	if ( fDiag != NULL )
		delete []fDiag ;
	fRow = fCol = 0;
	fDiag = NULL;
	fDecomposed = 0;
	return( 1 );
}



/************/
/*** Copy ***/

void
TPZSBMatrix::Copy(const TPZSBMatrix &A )
{
	fBand = A.fBand;
	fRow = fCol = A.Dim();
	fDiag = new REAL[Size()] ;
	fDecomposed  = A.fDecomposed;
	fDefPositive = A.fDefPositive;
	
	if ( fDiag == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__,"Copy( TPZSBMatrix ) <memory allocation error>" );
	
	REAL *dst = fDiag;
	REAL *src = A.fDiag;
	REAL *end = A.fDiag + Size();
	while ( src < end )
		*dst++ = *src++;
}

#ifdef OOPARLIB

int TPZSBMatrix::Unpack (TReceiveStorage *buf ){
	TSimMatrix::Unpack(buf);
	buf->UpkInt(&fBand);
	buf->UpkDouble(fDiag,fBand);
	return 1;
}



TSaveable *TPZSBMatrix::Restore(TReceiveStorage *buf) {
	TPZSBMatrix *m = new TPZSBMatrix();
	m->Unpack(buf);
	return m;
}

int TPZSBMatrix::Pack( TSendStorage *buf ){
	TSimMatrix::Pack(buf);
	buf->PkInt(&fBand);
	buf->PkDouble(fDiag,fBand);
	return 1;
}


int TPZSBMatrix::DerivedFrom(long Classid){
	if(Classid == GetClassID()) return 1;
	return TSimMatrix::DerivedFrom(Classid);
}

int TPZSBMatrix::DerivedFrom(char *classname){
	if(!strcmp(ClassName(),classname)) return 1;
	return TSimMatrix::DerivedFrom(classname);
}

#endif

