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

TPZFBMatrix::TPZFBMatrix()
: TPZMatrix( 0, 0 )
{
	fElem = NULL;
	fBand = 0;
}

/********************/
/*** Constructors ***/

TPZFBMatrix::TPZFBMatrix( int dim, int band_width )
: TPZMatrix( dim, dim )
{
	int size;
	fBand = band_width;
	size  = dim * (2*fBand + 1);
	if(size) {
		fElem = new REAL[ size ] ;
		if ( fElem == NULL ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	} else {
		fElem = NULL;
	}
	
	// Zera a Matriz.
	REAL *p   = fElem;
	REAL *end = fElem + size;
	while ( p < end )
		*p++ = 0.0;
}



/*********************************/
/*** Constructor( TPZFBMatrix& ) ***/

TPZFBMatrix::TPZFBMatrix (const TPZFBMatrix & A)
: TPZMatrix( A.Dim(), A.Dim() )
{
	// Philippe 20/10/97
	int size = A.Dim()*(2 * A.fBand + 1);
	fBand = A.fBand;
	fElem = new REAL[ size ] ;
	
	if ( fElem == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Constructor <memory allocation error>." );
	
	// Copia a matriz
	REAL *src = A.fElem;
	REAL *dst = fElem;
	for ( int i = 0; i < size; i++ )
		*dst++ = *src++;
}



/******************/
/*** Destructor ***/

TPZFBMatrix::~TPZFBMatrix ()
{
	if ( fElem != NULL )
		delete []fElem ;
}



/***********/
/*** Put ***/

int
TPZFBMatrix::Put(const int row,const int col,const REAL& value )
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

const REAL&
TPZFBMatrix::Get(const int row,const int col ) const
{
	if ( (row >= Dim()) || (col >= Dim()) )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}



/******** Operacoes com matrizes FULL BAND  ********/

/******************/
/*** Operator = ***/

TPZFBMatrix&
TPZFBMatrix::operator=(const TPZFBMatrix & A )
{
	if(this == &A) return *this;
	//  int size = A.Dim() * A.fBand;
	//Philippe 23/10/97
	int size = A.Dim() * (1+2*A.fBand);
	int oldsize = Dim()+2*fBand;
	if(oldsize != size && fElem) delete [] fElem;
	
	if(oldsize != size && size) fElem = new REAL[size] ;
	else if(size == 0) fElem = 0;
	if ( size && fElem == NULL )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator= <memory allocation error>." );
	
	fRow  = A.fRow;
	fCol  = A.fCol;
	fBand = A.fBand;
	
	// Copia a matriz
	REAL *src = A.fElem;
	REAL *dst = fElem;
	memcpy(dst,src,size*sizeof(REAL));
	//  for ( int i = 0; i < size; i++ )
	//	 *dst++ = *src++;
	
	return *this;
}



/*******************************/
/*** Operator+( TPZFBMatrix & ) ***/
// DEPENDS ON THE STORAGE FORMAT
TPZFBMatrix
TPZFBMatrix::operator+(const TPZFBMatrix & A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
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
	
	TPZFBMatrix res( Dim(), newBand );
	REAL *pm = fElem + (inc_pm >> 1);
	REAL *pa = A.fElem + (inc_pa >> 1);
	REAL *pr = res.fElem + (inc_pr >> 1);
	
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

TPZFBMatrix
TPZFBMatrix::operator-(const TPZFBMatrix & A ) const
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
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
	
	TPZFBMatrix res( Dim(), newBand );
	REAL *pm = fElem + (inc_pm >> 1);
	REAL *pa = A.fElem + (inc_pa >> 1);
	REAL *pr = res.fElem + (inc_pr >> 1);
	
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

TPZFBMatrix &
TPZFBMatrix::operator+=(const TPZFBMatrix & A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator+ <matrixs with different dimensions>" );
	
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
	
	REAL *pm = fElem + (inc_pm >> 1);
	REAL *pa = A.fElem + (inc_pa >> 1);
	
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

TPZFBMatrix &
TPZFBMatrix::operator-=(const TPZFBMatrix & A )
{
	if ( A.Dim() != Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Operator- <matrixs with different dimensions>" );
	
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
	
	REAL *pm = fElem + (inc_pm >> 1);
	REAL *pa = A.fElem + (inc_pa >> 1);
	
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

void TPZFBMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						  const REAL alpha,const REAL beta ,const int opt,const int stride ) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && Cols()*stride != x.Rows()) || Rows()*stride != x.Rows())
		TPZMatrix::Error(__PRETTY_FUNCTION__, "TPZFBMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		TPZMatrix::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::MultAdd incompatible dimensions\n");
	}
	PrepareZ(y,z,beta,opt,stride);
	int rows = Rows();
	int xcols = x.Cols();
	int ic, r;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			int begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBand, 0 );
				end   = MIN( r + fBand + 1, Dim() );
				REAL val = z.GetVal(r*stride,ic);
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
				REAL val = z.GetVal(r*stride,ic);
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

TPZFBMatrix TPZFBMatrix::operator-() const {
	TPZFBMatrix temp(*this);
	temp *= -1.;
	return temp;
}


/******** Operacoes com valores NUMERICOS ********/



/**************************/
/*** Operator*( value ) ***/

TPZFBMatrix
TPZFBMatrix::operator*(const REAL value ) const
{
	TPZFBMatrix res( Dim(), fBand );
	int size = Dim() * (2*fBand + 1);
	
	REAL *dst = res.fElem;
	REAL *src = fElem;
	for ( int i = 0; i < size; i++ )
		*dst++ = (*src++) * value;
	
	return( res );
}





/***************************/
/*** Operator*=( value ) ***/

TPZFBMatrix &
TPZFBMatrix::operator*=(const REAL value )
{
	if ( value != 1.0 )
    {
		int size = Dim() * (2*fBand + 1);
		REAL *dst = fElem;
		for ( int i = 0; i < size; i++ )
			*dst++ *= value;
    }
	
	return( *this );
}



/**************/
/*** Resize ***/
// DEPENDS ON THE STORAGE FORMAT

int
TPZFBMatrix::Resize(const int newRows,const int newCols)
{
	if ( newRows != newCols )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	if ( !fBand )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	
	if ( newRows == Dim() )
		return( 1 );
	
	int bandSize = 2 * fBand + 1;
	REAL *newElem = new REAL[ newRows * bandSize ] ;
	if ( !newElem )
		return TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	int minDim = ( Dim() < newRows ? Dim() : newRows );
	REAL *src = fElem;
	REAL *dst = newElem;
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
	fRow  = newRows;
	fCol  = newCols;
	return( 1 );
}



/*************/
/*** Redim ***/
int
TPZFBMatrix::Redim(const int newRows,const int newCols )
{
	if ( newRows != newCols )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	//  if ( !fBand ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	
	if ( fElem  )  delete []fElem;
	
	int size = newRows * (2 * fBand + 1);
	if(size) {
		fElem = new REAL[ size ] ;
		if ( fElem == NULL ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	} else {
		fElem = NULL;
	}
	
	/*  REAL *dst = fElem;
	 for ( int i = 0; i < size; i++ )
	 *dst++ = 0.0;   */
	
	fRow = fCol = newRows;
	Zero();
	
	return( 1 );
}


/***************/
/**** Zero ****/
int
TPZFBMatrix::Zero()
{
	int size = Dim() * (2 * fBand + 1);
	
	REAL *p = fElem,*plast=fElem+size;
	while(p < plast) *p++ = 0.0;
	
	fDecomposed = 0;
	
	return( 1 );
}


/***************/
/*** SetBand ***/
// DEPENDS ON THE STORAGE FORMAT
int
TPZFBMatrix::SetBand( int newBand )
{
	if ( newBand >= Dim() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "SetBand <the band must be lower than the matrix dimension " );
	
	int newSize = Dim() * (2 * newBand + 1);
	REAL *newElem = new REAL[ newSize ] ;
	if ( newElem == NULL )
		return TPZMatrix::Error(__PRETTY_FUNCTION__, "Resize <memory allocation error>." );
	
	REAL *dst = newElem;
	REAL *src = fElem;
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
void
TPZFBMatrix::Transpose (TPZMatrix *const T) const
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
int
TPZFBMatrix::Decompose_LU(std::list<int> &singular)
{
    if (  fDecomposed && fDecomposed == ELU) {
        return ELU;
    } else if(fDecomposed) {
        TPZMatrix::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
    }
    int rows = Rows();
    int  min = rows - 1;
    int imax;
    REAL nn;
    
    REAL *kFirstPtr = fElem+fBand+1;
    int k;
    for ( k = 0; k < min ; k++ )
    {
        REAL *iPtr = fElem+fBand*(2*k+1)+k;
        REAL pivot = *iPtr;
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
        
        REAL *jFirstPtr = fElem+fBand*(2*k+3)+k+1;
        imax = k+fBand+1;
        if(imax > rows) imax = rows;
        REAL *kLastPtr = fElem+fBand*(2*k+1)+imax;
        REAL *kPtr, *jPtr;
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
    REAL *iPtr = fElem+fBand*(2*k+1)+k;
    REAL pivot = *iPtr;
    if ( IsZero(pivot))
    {
        (*iPtr)++;
        pivot++;
        std::cout << __PRETTY_FUNCTION__ << " at row " << k << " is singular\n";
        //            TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
    }
	
    fDecomposed = ELU;
    return 1;
	
}

int
TPZFBMatrix::Decompose_LU()
{
	if (  fDecomposed && fDecomposed == ELU) {
		return ELU;
	} else if(fDecomposed) {
		TPZMatrix::Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
	}
	int rows = Rows();
	int  min = rows - 1;
	int imax;
	REAL nn;
	
	REAL *kFirstPtr = fElem+fBand+1;
	for ( int k = 0; k < min ; k++ )
    {
		REAL *iPtr = fElem+fBand*(2*k+1)+k;
		REAL pivot = *iPtr;
		if ( IsZero(pivot) )
		{
			TPZMatrix::Error(__PRETTY_FUNCTION__, "Decompose_LU <matrix is singular>" );
		}
		
		//       if(IsZero(pivot)){
		//         if(pivot < 0.) *iPtr = -1.e-10;
		//         else *iPtr = +1.e-10;
		//         pivot = *iPtr;
		//       }
		
		REAL *jFirstPtr = fElem+fBand*(2*k+3)+k+1;
		imax = k+fBand+1;
		if(imax > rows) imax = rows;
		REAL *kLastPtr = fElem+fBand*(2*k+1)+imax;
		REAL *kPtr, *jPtr;
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
	fDecomposed = ELU;
	return 1;
}


int TPZFBMatrix::Substitution( TPZFMatrix *B ) const{
	
    int rowb = B->Rows();
    int colb = B->Cols();
	REAL *BfElem = &(B->operator ()(0,0));
    if ( rowb != Rows() )
		TPZMatrix::Error(__PRETTY_FUNCTION__, "SubstitutionLU <incompatible dimensions>" );
	int i;
    for ( i = 0; i < rowb; i++ ) {
        for ( int col = 0; col < colb; col++ ) {
			int firstj = i-fBand;
			if(firstj < 0) firstj=0;
			REAL *jptr = fElem+fBand*(2*i+1)+firstj;
			REAL *iiptr = fElem+fBand*(2*i+1)+i;
			REAL *Biptr = BfElem+i+rowb*col;
			REAL *Bjptr = BfElem+firstj+rowb*col;
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
			REAL *jlastptr = fElem + fBand*(2*i+1)+jlast;
			REAL *jptr = fElem + fBand*(2*i+1)+i+1;
			REAL *Biptr = BfElem+i+rowb*col;
			int j=i+1;
			REAL *Bjptr = BfElem+j+rowb*col;
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

int
TPZFBMatrix::Clear()
{
	if ( fElem != NULL )
		delete( fElem );
	
	fElem = NULL;
	fRow  = fCol = 0;
	fBand = 0;
	return( 1 );
}

#ifdef OOPARLIB

int TPZFBMatrix::Unpack( TReceiveStorage *buf ){
	TPZMatrix::Unpack(buf);
	fElem= new REAL[fBand];    //dim *(2*fBand+1)??
	buf->UpkDouble(fElem,fBand);
	return 1;
}



TSaveable *TPZFBMatrix::Restore(TReceiveStorage *buf) {
	TPZFBMatrix *m = new TPZFBMatrix();
	m->Unpack(buf);
	return m;
}

int TPZFBMatrix::Pack( TSendStorage *buf ) const {
	TPZMatrix::Pack(buf);
	buf->PkDouble(fElem,fBand);
	return 1;
}


int TPZFBMatrix::DerivedFrom(const long Classid) const {
	if(Classid == GetClassID()) return 1;
	return TPZMatrix::DerivedFrom(Classid);
}

int TPZFBMatrix::DerivedFrom(const char *classname) const {
	
	if(!strcmp(ClassName(),classname)) return 1;
	return TPZMatrix::DerivedFrom(classname);
}

#endif

