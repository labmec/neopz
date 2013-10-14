/**
 * @file
 * @brief Contains the implementation of the TPZFBMatrix methods.
 */

#include <math.h>
#include "pzfmatrix.h"
#include "pzbndmat.h"

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
TPZFBMatrix<TVar>::TPZFBMatrix( long dim, long band_width )
: TPZMatrix<TVar>( dim, dim )
{
	unsigned long size;
	fBand = band_width;
	size  = (2*fBand + 1);
	size *= dim;
	if(size) {
		fElem = new TVar[ size ] ;
		if ( fElem == NULL ) 
			TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
		// Zera a Matriz.
		TVar *p   = fElem;
		TVar *end = fElem + size;
		while ( p < end )
			*p++ = (TVar)(0.0);
	} else {
		fElem = NULL;
	}	
}



/*********************************/
/*** Constructor( TPZFBMatrix& ) ***/

template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix (const TPZFBMatrix<TVar> & A)
: TPZMatrix<TVar>( A.Dim(), A.Dim() )
{
	// Philippe 20/10/97
	unsigned long size = ((unsigned long)A.Dim())*(2 * A.fBand + 1);
	fBand = A.fBand;
	fElem = new TVar[ size ] ;
	
	if ( fElem == NULL )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	TVar *src = A.fElem;
	TVar *dst = fElem;
	for (unsigned long i = 0L; i < size; i++ )
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
TPZFBMatrix<TVar>::Put(const long row,const long col,const TVar& value )
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
TPZFBMatrix<TVar>::Get(const long row,const long col ) const
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
	unsigned long size = ((unsigned long)A.Dim()) * (1+2*A.fBand);
	unsigned long oldsize = ((unsigned long)Dim())+2*fBand;
	if(oldsize != size && fElem) delete [] fElem;
	
	if(oldsize != size && size) fElem = new TVar[size] ;
	else if(size == 0L) fElem = 0;
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
	
	long newBand, minBand;
	unsigned long inc_pa, inc_pm, inc_pr;
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
	
	for ( long row = 0; row < Dim(); row++ )
    {
		for ( long col = 0; col < minBand; col++ )
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
	
	long newBand, minBand;
	unsigned long inc_pa, inc_pm, inc_pr;
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
	
	for ( long row = 0; row < Dim(); row++ )
    {
		for ( long col = 0; col < minBand; col++ )
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
	
	long newBand, minBand;
	unsigned long inc_pa, inc_pm;
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
	
	for ( long row = 0; row < Dim(); row++ )
    {
		for ( long col = 0; col < minBand; col++ )
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
	
	long newBand, minBand;
	unsigned long inc_pa, inc_pm;
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
	
	for ( long row = 0; row < Dim(); row++ )
    {
		for ( long col = 0; col < minBand; col++ )
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
	long rows = this->Rows();
	long xcols = x.Cols();
	long ic, r;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			long begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBand, 0 );
				end   = MIN( r + fBand + 1, Dim() );
				TVar val = z.GetVal(r*stride,ic);
				// Calcula um elemento da resposta.
				for ( long i = begin ; i < end; i++ ) val += alpha * GetVal( r, i ) * x.GetVal( i*stride, ic );
				z.PutVal( r*stride, ic, val );
			}
		}
	} else {
		for (ic = 0; ic < xcols; ic++) {
			long begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBand, 0 );
				end   = MIN( r + fBand + 1, Dim() );
				TVar val = z.GetVal(r*stride,ic);
				// Calcula um elemento da resposta.
				for ( long i = begin ; i < end; i++ ) val += alpha * GetVal( i, r ) * x.GetVal( i*stride, ic );
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
	temp *= (TVar)(-1.);
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
	unsigned long size = ((unsigned long)Dim()) * (2*fBand + 1);
	
	TVar *dst = res.fElem;
	TVar *src = fElem;
	for (unsigned long i = 0L; i < size; i++ )
		*dst++ = (*src++) * value;
	
	return( res );
}





/***************************/
/*** Operator*=( value ) ***/

template<>
TPZFBMatrix<std::complex<float> > &
TPZFBMatrix<std::complex<float> >::operator*=(const std::complex<float> value )
{
	if ( value.real() != 1.0 || value.imag() != 0. )
    {
		unsigned long size = (2*fBand + 1);
        size *= Dim();
        std::complex<float> *dst = fElem;
		for ( unsigned long i = 0L; i < size; i++ )
			*dst++ *= value;
    }
	
	return( *this );
}

template<>
TPZFBMatrix<std::complex<double> > &
TPZFBMatrix<std::complex<double> >::operator*=(const std::complex<double> value )
{
	if ( value.real() != 1.0 || value.imag() != 0. )
    {
		unsigned long size = (2*fBand + 1);
        size *= Dim();
        std::complex<double> *dst = fElem;
		for ( unsigned long i = 0L; i < size; i++ )
			*dst++ *= value;
    }
	
	return( *this );
}



template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator*=(const TVar value )
{
	if ( value != (TVar)1.0 )
    {
		unsigned long size = (2*fBand + 1);
        size *= Dim();
		TVar *dst = fElem;
		for ( unsigned long i = 0L; i < size; i++ )
			*dst++ *= value;
    }
	
	return( *this );
}



/**************/
/*** Resize ***/
// DEPENDS ON THE STORAGE FORMAT

template<class TVar>
int
TPZFBMatrix<TVar>::Resize(const long newRows,const long newCols)
{
	if ( newRows != newCols )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	if ( !fBand )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	
	if ( newRows == this->Dim() )
		return( 1 );
	
	long bandSize = 2 * fBand + 1;
	unsigned long size = bandSize;
    size *= newRows;
	TVar *newElem = new TVar[ size ] ;
	if ( !newElem )
		return TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	long minDim = ( Dim() < newRows ? Dim() : newRows );
	TVar *src = fElem;
	TVar *dst = newElem;
	long r,c;
	
	// Copia as linhas da matriz antiga para a nova.
	for ( r = 0; r < minDim; r++ )
		// Copia os elementos de uma linha.
		for ( c = 0; c < bandSize; c++ )
			*dst++ = *src++;
	
	// Preenche as linha que sobrarem (se sobrarem) com ZEROS.
	for ( ; r < newRows; r++ )
		for ( c = 0; c < bandSize; c++ )
			*dst++ = (TVar)(0.0);
	
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
TPZFBMatrix<TVar>::Redim(const long newRows,const long newCols )
{
	if ( newRows != newCols )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	//  if ( !fBand ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	if ( fElem  )  delete []fElem;
	
	unsigned long size = (2 * fBand + 1);
    size *= newRows;
	if(size>0) {
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
	unsigned long size = (2 * fBand + 1);
    size *= Dim();
	
	TVar *p = fElem,*plast=fElem+size;
	while(p < plast) *p++ = (TVar)(0.0);
	
	this->fDecomposed = 0;
	
	return( 1 );
}


/***************/
/*** SetBand ***/
// DEPENDS ON THE STORAGE FORMAT
template<class TVar>
int
TPZFBMatrix<TVar>::SetBand( long newBand )
{
	if ( newBand >= Dim() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SetBand <the band must be lower than the matrix dimension " );
	
	unsigned long newSize = (2 * newBand + 1);
    newSize *= Dim();
	TVar *newElem = new TVar[ newSize ] ;
	if ( newElem == NULL )
		return TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	TVar *dst = newElem;
	TVar *src = fElem;
	for ( long r = 0; r < Dim(); r++ )
		for ( long i = -newBand; i <= newBand; i++ )
			*dst++ = ( (i >= -fBand) && (i <= fBand) ) ? *src++ : (TVar)0.0;
	
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
	
	long end, begin;
	//REAL *p = fElem;
	for ( long r = 0; r < Dim(); r++ )
    {
		begin = MAX( r - fBand, 0 );
		end   = MIN( r + fBand + 1, Dim() );
		for ( long c = begin; c < end; c++ )
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
TPZFBMatrix<TVar>::Decompose_LU(std::list<long> &singular)
{
    if (  this->fDecomposed && this->fDecomposed == ELU) {
        return ELU;
    } else if(this->fDecomposed) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
    }
    long rows = this->Rows();
    long  min = rows - 1;
    long imax;
    TVar nn;
    
    TVar *kFirstPtr = fElem+fBand+1;
    long k;
    for ( k = 0; k < min ; k++ )
    {
        TVar *iPtr = fElem;
		iPtr += ((unsigned long)fBand*(2*k+1)+k);
        TVar pivot = *iPtr;
        if ( IsZero(pivot))
        {
            (*iPtr) = (*iPtr)+(TVar)1.;
            pivot = pivot+(TVar)1.;
//            std::cout << __PRETTY_FUNCTION__ << " at row " << k << " is singular\n";
			//            TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
        }
        
        //       if(IsZero(pivot)){
        //         if(pivot < 0.) *iPtr = -1.e-10;
        //         else *iPtr = +1.e-10;
        //         pivot = *iPtr;
        //       }
        
        TVar *jFirstPtr = fElem;
		jFirstPtr += (((unsigned long)fBand)*(2*k+3)+k+1);
        imax = k+fBand+1;
        if(imax > rows) imax = rows;
        TVar *kLastPtr = fElem;
		kLastPtr += (((unsigned long)fBand)*(2*k+1)+imax);
        TVar *kPtr, *jPtr;
        for ( long i = k+1; i < imax; i++ )
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
    TVar *iPtr = fElem;
	iPtr += (((unsigned long)fBand)*(2*k+1)+k);
    TVar pivot = *iPtr;
    if ( IsZero(pivot))
    {
        (*iPtr) = (*iPtr)+(TVar)1.;
        pivot = pivot+(TVar)1.;
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
	long rows = this->Rows();
	long  min = rows - 1;
	long imax;
	TVar nn;
	
	TVar *kFirstPtr = fElem+fBand+1;
	for ( long k = 0; k < min ; k++ )
    {
		TVar *iPtr = fElem+((unsigned long)fBand)*(2*k+1)+k;
		TVar pivot = *iPtr;
		if ( IsZero(pivot) )
		{
			TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
		}
		
		TVar *jFirstPtr = fElem+((unsigned long)fBand)*(2*k+3)+k+1;
		imax = k+fBand+1;
		if(imax > rows) imax = rows;
		TVar *kLastPtr = fElem+((unsigned long)fBand)*(2*k+1)+imax;
		TVar *kPtr, *jPtr;
		for ( long i = k+1; i < imax; i++ )
		{
			iPtr += 2*fBand;
			*iPtr /= pivot;
			nn = *iPtr;
			kPtr = kFirstPtr;
			jPtr = jFirstPtr;
			while(kPtr < kLastPtr) *jPtr++ -= nn * *kPtr++;
			jFirstPtr += 2*fBand;
		}
		kFirstPtr += 2*fBand+1;
    }
	this->fDecomposed = ELU;
	return 1;
}

template<class TVar>
int TPZFBMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const{
	
    long rowb = B->Rows();
    long colb = B->Cols();
	TVar *BfElem = &(B->operator ()(0,0));
    if ( rowb != this->Rows() )
		TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SubstitutionLU <incompatible dimensions>" );
	long i;
    for ( i = 0; i < rowb; i++ ) {
        for ( long col = 0; col < colb; col++ ) {
			long firstj = i-fBand;
			if(firstj < 0) firstj=0;
			TVar *jptr = fElem+((unsigned long)fBand)*(2*i+1)+firstj;
			TVar *iiptr = fElem+((unsigned long)fBand)*(2*i+1)+i;
			TVar *Biptr = BfElem+i+((unsigned long)rowb)*col;
			TVar *Bjptr = BfElem+firstj+((unsigned long)rowb)*col;
			while(jptr < iiptr) {
				*Biptr -= *jptr++ * *Bjptr++;
			}
        }
    }
    for (long col=0; col<colb; col++) {
        for ( i = rowb-1; i >= 0; i-- ) {
			long jlast = i+fBand+1;
			if(jlast > rowb) jlast = rowb;
			TVar *jlastptr = fElem + ((unsigned long)fBand)*(2*i+1)+jlast;
			TVar *jptr = fElem + ((unsigned long)fBand)*(2*i+1)+i+1;
			TVar *Biptr = BfElem+i+((unsigned long)rowb)*col;
			long j=i+1;
			TVar *Bjptr = BfElem+j+((unsigned long)rowb)*col;
			while(jptr < jlastptr) {
				*Biptr -= *jptr++ * *(Bjptr++);
			}
            B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
    }
    return( 1 );
}


/************************** Private **************************/

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
template class TPZFBMatrix<float>;
template class TPZFBMatrix<long>;
template class TPZFBMatrix<int>;
template class TPZFBMatrix<std::complex<long double> >;
template class TPZFBMatrix<std::complex<double> >;
template class TPZFBMatrix<std::complex<float> >;
