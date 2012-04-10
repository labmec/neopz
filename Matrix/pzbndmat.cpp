/**
 * @file
 * @brief Contains the implementation of the TPZFBMatrix methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tbndmat.cc
//
// Class:  TPZFBMatrix
//
// Obs.:   Implementa matrizes Banda.
//
// Versao: 04 / 1996.
//

#include <math.h>
#include "pzfmatrix.h"
#include "pzbndmat.h"

//#include "pzerror.h"
#include <stdlib.h>

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzfbmatrix"));
#endif

using namespace std;

/*******************/
/*** Constructor ***/
template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix()
: TPZMatrix<TVar>( 0, 0 )
{
	fElem = NULL;
	fBand = 0;
}

/********************/
/*** Constructors ***/

template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix( int dim, int band_width )
: TPZMatrix<TVar>( dim, dim )
{
	int size;
	fBand = band_width;
	size  = dim * (2*fBand + 1);
	if(size) {
		fElem = new TVar[ size ] ;
		if ( fElem == NULL ) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	} else {
		fElem = NULL;
	}
	
	// Zera a Matriz.
	TVar *p   = fElem;
	TVar *end = fElem + size;
	while ( p < end )
		*p++ = 0.0;
}



/*********************************/
/*** Constructor( TPZFBMatrix& ) ***/

template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix (const TPZFBMatrix<TVar> & A)
: TPZMatrix<TVar>( A.Dim(), A.Dim() )
{
	// Philippe 20/10/97
	int size = A.Dim()*(2 * A.fBand + 1);
	fBand = A.fBand;
	fElem = new TVar[ size ] ;
	
	if ( fElem == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	for ( int i = 0; i < size; i++ )
		*dst++ = *src++;
}



/******************/
/*** Destructor ***/
template<class TVar>
TPZFBMatrix<TVar>::~TPZFBMatrix ()
{
	if ( fElem != NULL )
		delete []fElem ;
}



/***********/
/*** Put ***/

template<class TVar>
int
TPZFBMatrix<TVar>::Put(const int row,const int col,const TVar& value )
{
	if ( (row >= Dim()) || (col >= Dim()) || row<0 || col<0)
    {
		cout << "TPZFBMatrix::Put: " << row << "," << col << "," << Dim();
		cout << "\n";
		return( 0 );
    }
	
	return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

template<class TVar>
const TVar&
TPZFBMatrix<TVar>::Get(const int row,const int col ) const
{
	if ( (row >= Dim()) || (col >= Dim()) )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}



/******** Operacoes com matrizes FULL BAND  ********/

/******************/
/*** Operator = ***/


template<class TVar>
TPZFBMatrix<TVar>&
TPZFBMatrix<TVar>::operator=(const TPZFBMatrix<TVar> & A )
{
	if(this == &A) return *this;
	//  int size = A.Dim() * A.fBand;
	//Philippe 23/10/97
	int size = A.Dim() * (1+2*A.fBand);
	int oldsize = Dim()+2*fBand;
	if(oldsize != size && fElem) delete [] fElem;
	
	if(oldsize != size && size) fElem = new TVar[size] ;
	else if(size == 0) fElem = 0;
	if ( size && fElem == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator= <memory allocation error>." );
	
	this->fRow  = A.fRow;
	this->fCol  = A.fCol;
	fBand = A.fBand;
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	memcpy(dst,src,size*sizeof(TVar));
	//  for ( int i = 0; i < size; i++ )
	//	 *dst++ = *src++;
	
	return *this;
}



/*******************************/
/*** Operator+( TPZFBMatrix & ) ***/
// DEPENDS ON THE STORAGE FORMAT
template<class TVar>
TPZFBMatrix<TVar>
TPZFBMatrix<TVar>::operator+(const TPZFBMatrix<TVar> & A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
	int newBand, minBand;
	int inc_pa, inc_pm, inc_pr;
	if ( fBand < A.fBand )
    {
		newBand = A.fBand;
		minBand = 2 * fBand + 1;
		inc_pm  = 0;
		inc_pa  = inc_pr = 2 * (A.fBand - fBand);
    }
	else
    {
		newBand = fBand;
		minBand = 2 * A.fBand + 1;
		inc_pa  = 0;
		inc_pm  = inc_pr = 2 * (fBand - A.fBand);
    }
	
	TPZFBMatrix<TVar> res( Dim(), newBand );
	TVar *pm = fElem + (inc_pm >> 1);
	TVar *pa = A.fElem + (inc_pa >> 1);
	TVar *pr = res.fElem + (inc_pr >> 1);
	
	for ( int row = 0; row < Dim(); row++ )
    {
		for ( int col = 0; col < minBand; col++ )
			*pr++ = (*pm++) + (*pa++);
		pr += inc_pr;
		pm += inc_pm;
		pa += inc_pa;
    }
	
	return( res );
}



/*******************************/
/*** Operator-( TPZFBMatrix & ) ***/
// DEPENDS ON THE STORAGE FORMAT

template<class TVar>
TPZFBMatrix<TVar>
TPZFBMatrix<TVar>::operator-(const TPZFBMatrix<TVar> & A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
	int newBand, minBand;
	int inc_pa, inc_pm, inc_pr;
	if ( fBand < A.fBand )
    {
		newBand = A.fBand;
		minBand = 2 * fBand + 1;
		inc_pm  = 0;
		inc_pa  = inc_pr = 2 * (A.fBand - fBand);
    }
	else
    {
		newBand = fBand;
		minBand = 2 * A.fBand + 1;
		inc_pa  = 0;
		inc_pm  = inc_pr = 2 * (fBand - A.fBand);
    }
	
	TPZFBMatrix<TVar> res( Dim(), newBand );
	TVar *pm = fElem + (inc_pm >> 1);
	TVar *pa = A.fElem + (inc_pa >> 1);
	TVar *pr = res.fElem + (inc_pr >> 1);
	
	for ( int row = 0; row < Dim(); row++ )
    {
		for ( int col = 0; col < minBand; col++ )
			*pr++ = (*pm++) - (*pa++);
		pr += inc_pr;
		pm += inc_pm;
		pa += inc_pa;
    }
	
	return( res );
}



/*******************************/
/*** Operator*( TPZFBMatrix & ) ***/



/********************************/
/*** Operator+=( TPZFBMatrix & ) ***/
// DEPENDS ON THE STORAGE FORMAT

template <class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator+=(const TPZFBMatrix<TVar> & A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
	int newBand, minBand;
	int inc_pa, inc_pm;
	if ( fBand < A.fBand )
    {
		newBand = A.fBand;
		minBand = 2 * fBand + 1;
		inc_pm  = 0;
		inc_pa  = 2 * (A.fBand - fBand);
    }
	else
    {
		newBand = fBand;
		minBand = 2 * A.fBand + 1;
		inc_pa  = 0;
		inc_pm  = 2 * (fBand - A.fBand);
    }
	
	SetBand( newBand );
	
	TVar *pm = fElem + (inc_pm >> 1);
	TVar *pa = A.fElem + (inc_pa >> 1);
	
	for ( int row = 0; row < Dim(); row++ )
    {
		for ( int col = 0; col < minBand; col++ )
			*pm++ += *pa++;
		pm += inc_pm;
		pa += inc_pa;
    }
	
	return( *this );
}



/*******************************/
/*** Operator-=( TPZFBMatrix & ) ***/
// DEPENDS ON THE STORAGE FORMAT

template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator-=(const TPZFBMatrix<TVar> & A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
	int newBand, minBand;
	int inc_pa, inc_pm;
	if ( fBand < A.fBand )
    {
		newBand = A.fBand;
		minBand = 2 * fBand + 1;
		inc_pm  = 0;
		inc_pa  = 2 * (A.fBand - fBand);
    }
	else
    {
		newBand = fBand;
		minBand = 2 * A.fBand + 1;
		inc_pa  = 0;
		inc_pm  = 2 * (fBand - A.fBand);
    }
	
	SetBand( newBand );
	
	TVar *pm = fElem + (inc_pm >> 1);
	TVar *pa = A.fElem + (inc_pa >> 1);
	
	for ( int row = 0; row < Dim(); row++ )
    {
		for ( int col = 0; col < minBand; col++ )
			*pm++ -= *pa++;
		pm += inc_pm;
		pa += inc_pa;
    }
	
	return( *this );
}



/******** Operacoes com MATRIZES GENERICAS ********/


/*******************/
/*** MultiplyAdd ***/
//
//  perform a multiply add operation to be used by iterative solvers
//

template<class TVar>
void TPZFBMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						  const TVar alpha,const TVar beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && this->Cols()*stride != x.Rows()) || this->Rows()*stride != x.Rows())
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZFBMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::MultAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt,stride);
	int rows = this->Rows();
	int xcols = x.Cols();
	int ic, r;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			int begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBand, 0 );
				end   = MIN( r + fBand + 1, Dim() );
				TVar val = z.GetVal(r*stride,ic);
				// Calcula um elemento da resposta.
				for ( int i = begin ; i < end; i++ ) val += alpha * GetVal( r, i ) * x.GetVal( i*stride, ic );
				z.PutVal( r*stride, ic, val );
			}
		}
	} else {
		for (ic = 0; ic < xcols; ic++) {
			int begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBand, 0 );
				end   = MIN( r + fBand + 1, Dim() );
				TVar val = z.GetVal(r*stride,ic);
				// Calcula um elemento da resposta.
				for ( int i = begin ; i < end; i++ ) val += alpha * GetVal( i, r ) * x.GetVal( i*stride, ic );
				z.PutVal( r*stride, ic, val );
			}
		}
	}
}




/*******************/
/*** Operator += ***/
//

template<class TVar>
TPZFBMatrix<TVar> TPZFBMatrix<TVar>::operator-() const {
	TPZFBMatrix<TVar> temp(*this);
	temp *= -1.;
	return temp;
}


/******** Operacoes com valores NUMERICOS ********/



/**************************/
/*** Operator*( value ) ***/

template<class TVar>
TPZFBMatrix<TVar>
TPZFBMatrix<TVar>::operator*(const TVar value ) const
{
	TPZFBMatrix<TVar> res( Dim(), fBand );
	int size = Dim() * (2*fBand + 1);
	
	TVar *dst = res.fElem;
	TVar *src = fElem;
	for ( int i = 0; i < size; i++ )
		*dst++ = (*src++) * value;
	
	return( res );
}





/***************************/
/*** Operator*=( value ) ***/

template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator*=(const TVar value )
{
	if ( value != 1.0 )
    {
		int size = Dim() * (2*fBand + 1);
		TVar *dst = fElem;
		for ( int i = 0; i < size; i++ )
			*dst++ *= value;
    }
	
	return( *this );
}



/**************/
/*** Resize ***/
// DEPENDS ON THE STORAGE FORMAT

template<class TVar>
int
TPZFBMatrix<TVar>::Resize(const int newRows,const int newCols)
{
	if ( newRows != newCols )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	if ( !fBand )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	
	if ( newRows == this->Dim() )
		return( 1 );
	
	int bandSize = 2 * fBand + 1;
	TVar *newElem = new TVar[ newRows * bandSize ] ;
	if ( !newElem )
		return TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	int minDim = ( Dim() < newRows ? Dim() : newRows );
	TVar *src = fElem;
	TVar *dst = newElem;
	int r,c;
	
	// Copia as linhas da matriz antiga para a nova.
	for ( r = 0; r < minDim; r++ )
		// Copia os elementos de uma linha.
		for ( c = 0; c < bandSize; c++ )
			*dst++ = *src++;
	
	// Preenche as linha que sobrarem (se sobrarem) com ZEROS.
	for ( ; r < newRows; r++ )
		for ( c = 0; c < bandSize; c++ )
			*dst++ = 0.0;
	
	delete( fElem );
	fElem = newElem;
	this->fRow  = newRows;
	this->fCol  = newCols;
	return( 1 );
}



/*************/
/*** Redim ***/
template<class TVar>
int
TPZFBMatrix<TVar>::Redim(const int newRows,const int newCols )
{
	if ( newRows != newCols )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	//  if ( !fBand ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	if ( fElem  )  delete []fElem;
	
	int size = newRows * (2 * fBand + 1);
	if(size) {
		fElem = new TVar[ size ] ;
		if ( fElem == NULL ) TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	} else {
		fElem = NULL;
	}
	
	/*  REAL *dst = fElem;
	 for ( int i = 0; i < size; i++ )
	 *dst++ = 0.0;   */
	
	this->fRow = this->fCol = newRows;
	Zero();
	
	return( 1 );
}


/***************/
/**** Zero ****/
template<class TVar>
int
TPZFBMatrix<TVar>::Zero()
{
	int size = Dim() * (2 * fBand + 1);
	
	TVar *p = fElem,*plast=fElem+size;
	while(p < plast) *p++ = 0.0;
	
	this->fDecomposed = 0;
	
	return( 1 );
}


/***************/
/*** SetBand ***/
// DEPENDS ON THE STORAGE FORMAT
template<class TVar>
int
TPZFBMatrix<TVar>::SetBand( int newBand )
{
	if ( newBand >= Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SetBand <the band must be lower than the matrix dimension " );
	
	int newSize = Dim() * (2 * newBand + 1);
	TVar *newElem = new TVar[ newSize ] ;
	if ( newElem == NULL )
		return TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	TVar *dst = newElem;
	TVar *src = fElem;
	for ( int r = 0; r < Dim(); r++ )
		for ( int i = -newBand; i <= newBand; i++ )
			*dst++ = ( (i >= -fBand) && (i <= fBand) ) ? *src++ : 0.0;
	
	if(fElem) delete []fElem;
	fElem = newElem;
	fBand = newBand;
	
	return( 1 );
}



/********************/
/*** Transpose () ***/
template<class TVar>
void
TPZFBMatrix<TVar>::Transpose (TPZMatrix<TVar> *const T) const
{
	T->Resize( Dim(), Dim() );
	
	int end,
    begin;
	//REAL *p = fElem;
	for ( int r = 0; r < Dim(); r++ )
    {
		begin = MAX( r - fBand, 0 );
		end   = MIN( r + fBand + 1, Dim() );
		for ( int c = begin; c < end; c++ )
		{
			T->PutVal( c, r, GetVal( r, c ) );
			//			cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		}
    }
}


/*****************/
/*** Decompose_LU ***/
//fElem[ fBand * (2*row + 1) + col ]
template<class TVar>
int
TPZFBMatrix<TVar>::Decompose_LU(std::list<int> &singular)
{
    if (  this->fDecomposed && this->fDecomposed == ELU) {
        return ELU;
    } else if(this->fDecomposed) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
    }
    int rows = this->Rows();
    int  min = rows - 1;
    int imax;
    TVar nn;
    
    TVar *kFirstPtr = fElem+fBand+1;
    int k;
    for ( k = 0; k < min ; k++ )
    {
        TVar *iPtr = fElem+fBand*(2*k+1)+k;
        TVar pivot = *iPtr;
        if ( IsZero(pivot))
        {
            (*iPtr)++;
            pivot++;
            std::cout << __PRETTY_FUNCTION__ << " at row " << k << " is singular\n";
			//            TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
        }
        
        //       if(IsZero(pivot)){
        //         if(pivot < 0.) *iPtr = -1.e-10;
        //         else *iPtr = +1.e-10;
        //         pivot = *iPtr;
        //       }
        
        TVar *jFirstPtr = fElem+fBand*(2*k+3)+k+1;
        imax = k+fBand+1;
        if(imax > rows) imax = rows;
        TVar *kLastPtr = fElem+fBand*(2*k+1)+imax;
        TVar *kPtr, *jPtr;
        for ( int i = k+1; i < imax; i++ )
        {
            iPtr += 2*fBand;
            *iPtr /= pivot;
            nn = *iPtr;
            kPtr = kFirstPtr;
            jPtr = jFirstPtr;
            while(kPtr < kLastPtr) *jPtr++ -= nn * *kPtr++;
            jFirstPtr += 2*fBand;
            //			for ( int j = k+1; j < Cols(); j++ )
            //				PutVal( i, j, GetVal( i, j ) - nn * GetVal( k, j ) );
        }
        kFirstPtr += 2*fBand+1;
    }
    TVar *iPtr = fElem+fBand*(2*k+1)+k;
    TVar pivot = *iPtr;
    if ( IsZero(pivot))
    {
        (*iPtr)++;
        pivot++;
        std::cout << __PRETTY_FUNCTION__ << " at row " << k << " is singular\n";
        //            TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
    }
	
    this->fDecomposed = ELU;
    return 1;
	
}

template<class TVar>
int
TPZFBMatrix<TVar>::Decompose_LU()
{
	if (  this->fDecomposed && this->fDecomposed == ELU) {
		return ELU;
	} else if(this->fDecomposed) {
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
	}
	int rows = this->Rows();
	int  min = rows - 1;
	int imax;
	TVar nn;
	
	TVar *kFirstPtr = fElem+fBand+1;
	for ( int k = 0; k < min ; k++ )
    {
		TVar *iPtr = fElem+fBand*(2*k+1)+k;
		TVar pivot = *iPtr;
		if ( IsZero(pivot) )
		{
			TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
		}
		
		//       if(IsZero(pivot)){
		//         if(pivot < 0.) *iPtr = -1.e-10;
		//         else *iPtr = +1.e-10;
		//         pivot = *iPtr;
		//       }
		
		TVar *jFirstPtr = fElem+fBand*(2*k+3)+k+1;
		imax = k+fBand+1;
		if(imax > rows) imax = rows;
		TVar *kLastPtr = fElem+fBand*(2*k+1)+imax;
		TVar *kPtr, *jPtr;
		for ( int i = k+1; i < imax; i++ )
		{
			iPtr += 2*fBand;
			*iPtr /= pivot;
			nn = *iPtr;
			kPtr = kFirstPtr;
			jPtr = jFirstPtr;
			while(kPtr < kLastPtr) *jPtr++ -= nn * *kPtr++;
			jFirstPtr += 2*fBand;
			//			for ( int j = k+1; j < Cols(); j++ )
			//				PutVal( i, j, GetVal( i, j ) - nn * GetVal( k, j ) );
		}
		kFirstPtr += 2*fBand+1;
    }
	this->fDecomposed = ELU;
	return 1;
}

template<class TVar>
int TPZFBMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const{
	
    int rowb = B->Rows();
    int colb = B->Cols();
	TVar *BfElem = &(B->operator ()(0,0));
    if ( rowb != this->Rows() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SubstitutionLU <incompatible dimensions>" );
	int i;
    for ( i = 0; i < rowb; i++ ) {
        for ( int col = 0; col < colb; col++ ) {
			int firstj = i-fBand;
			if(firstj < 0) firstj=0;
			TVar *jptr = fElem+fBand*(2*i+1)+firstj;
			TVar *iiptr = fElem+fBand*(2*i+1)+i;
			TVar *Biptr = BfElem+i+rowb*col;
			TVar *Bjptr = BfElem+firstj+rowb*col;
			while(jptr < iiptr) {
				*Biptr -= *jptr++ * *Bjptr++;
			}
			//            for ( int j = 0; j < i; j++ ) {
			//                B->PutVal( i, col, B->GetVal(i, col)-GetVal(i, j) * B->GetVal(j, col) );
			//            }
        }
    }
    for (int col=0; col<colb; col++) {
        for ( i = rowb-1; i >= 0; i-- ) {
			int jlast = i+fBand+1;
			if(jlast > rowb) jlast = rowb;
			TVar *jlastptr = fElem + fBand*(2*i+1)+jlast;
			TVar *jptr = fElem + fBand*(2*i+1)+i+1;
			TVar *Biptr = BfElem+i+rowb*col;
			int j=i+1;
			TVar *Bjptr = BfElem+j+rowb*col;
			while(jptr < jlastptr) {
				*Biptr -= *jptr++ * *(Bjptr++);
			}
			//            for ( int j = i+1; j < rowb ; j++ ) {
			//                B->PutVal( i, col, B->GetVal(i, col) -
			//                    GetVal(i, j) * B->GetVal(j, col) );
			//            }
            if ( IsZero( GetVal(i, i) ) ) {
				//                     TPZMatrix::Error(__PRETTY_FUNCTION__, "BackSub( SubstitutionLU ) <Matrix is singular" );
            }
            B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
    }
    return( 1 );
}


/************************** Private **************************/

/*************/
/*** Error ***/
/*
 int
 TPZFBMatrix::Error(const char *msg1,const char *msg2 ) 
 {
 ostringstream out;
 out << "TPZFBMatrix::" << msg1 << msg2 << ".\n";
 LOGPZ_ERROR (logger, out.str().c_str());
 // pzerror.Show();
 DebugStop();
 return 0;
 }
 */


/*************/
/*** Clear ***/

template<class TVar>
int
TPZFBMatrix<TVar>::Clear()
{
	if ( fElem != NULL )
		delete( fElem );
	
	fElem = NULL;
	this->fRow  = this->fCol = 0;
	fBand = 0;
	return( 1 );
}

#ifdef OOPARLIB

template<class TVar>
int TPZFBMatrix<TVar>::Unpack( TReceiveStorage *buf ){
	TPZMatrix<TVar>::Unpack(buf);
	fElem= new TVar[fBand];    //dim *(2*fBand+1)??
	buf->UpkDouble(fElem,fBand);
	return 1;
}


template<class TVar>
TSaveable *TPZFBMatrix<TVar>::Restore(TReceiveStorage *buf) {
	TPZFBMatrix<TVar> *m = new TPZFBMatrix<TVar>();
	m->Unpack(buf);
	return m;
}

template<class TVar>
int TPZFBMatrix<TVar>::Pack( TSendStorage *buf ) const {
	TPZMatrix<TVar>::Pack(buf);
	buf->PkDouble(fElem,fBand);
	return 1;
}

template<class TVar>
int TPZFBMatrix<TVar>::DerivedFrom(const long Classid) const {
	if(Classid == GetClassID()) return 1;
	return TPZMatrix<TVar>::DerivedFrom(Classid);
}

template<class TVar>
int TPZFBMatrix<TVar>::DerivedFrom(const char *classname) const {
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TPZMatrix<TVar>::DerivedFrom(classname);
}

#endif

template class TPZFBMatrix<long double>;
template class TPZFBMatrix<double>;
template class TPZFBMatrix<int>;
template class TPZFBMatrix<float>;
