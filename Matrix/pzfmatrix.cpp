/**
 * @file
 * @brief Contains the implementation of the TPZFMatrix<>methods.
 */

#ifdef __BOORLANDC__
#include <alloc.h>
#endif

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzerror.h"


#include <math.h>
#include <stdlib.h>
#include <sstream>
#include <string>
#include "pzaxestools.h"

#include "pzlog.h"

#ifdef DEBUG
#define DEBUG2
#endif

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzfmatrix"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif


#ifdef USING_ATLAS
extern "C"{
#include <cblas.h>
};
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
void cblas_daxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
#endif
#ifdef USING_BLAS
extern "C"{
#include "cblas.h"
};
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
void cblas_daxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
#endif

using namespace std;

/********************/
/*** Constructors ***/

template <class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZMatrix<TVar> &mat) : TPZMatrix<TVar>(mat), fElem(0),fGiven(0),fSize(0) {
	if(this->fRow*this->fCol) {
		fElem = new TVar[this->fRow*this->fCol];
		TVar * p = fElem;
		int i,j;
		for(j=0; j<this->fCol; j++) {
			for(i=0; i<this->fRow; i++) {
				*p++ = mat.GetVal(i,j);
			}
		}
	}
}

/********************************/
/*** Constructor( TPZFMatrix<TVar> & ) ***/

template<class TVar>
TPZFMatrix<TVar>::TPZFMatrix(const TPZFMatrix<TVar> &A)
: TPZMatrix<TVar>( A.fRow, A.fCol ), fElem(0), fGiven(0), fSize(0) {
    int size = this->fRow * this->fCol;
    if(!size) return;
    fElem = new TVar[ size ] ;
#ifdef DEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    // Copia a matriz
    TVar * src = A.fElem;
    TVar * p = fElem;
    memcpy(p,src,(size_t)size*sizeof(TVar));
}


/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/
template<class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZFMatrix<TVar> &A ) {
	if(this == &A) return *this;
	long size = A.fRow * A.fCol;
	
	TVar * newElem = fElem;
	if(fSize < size && size != this->fRow*this->fCol) {
		newElem = new TVar
		[size] ;
	} else if (fSize >= size) {
		newElem = fGiven;
	}
	
	if ( newElem == NULL && size > 0) Error( "Operator= <memory allocation error>." );
	if (fElem && fElem != newElem && fElem != fGiven) delete[]( fElem );
	this->fRow  = A.fRow;
	this->fCol  = A.fCol;
	fElem = newElem;
	
	// Copia a matriz
	memcpy(fElem,A.fElem,(size_t)size*sizeof(TVar));
	
	return *this;
}

template <class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int> &destination) {
	if(rhs.Cols() != this->Cols()) {
		PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
		DebugStop();
		return;
	}
	int ncol = this->Cols();
	int nrow = rhs.Rows();
	int i,j;
	for(j=0; j<ncol; j++) {
		for(i=0; i<nrow; i++) {
			operator()(destination[i],j) += rhs(i,j);
		}
	}
}

template<class TVar>
void TPZFMatrix<TVar>::AddFel(TPZFMatrix<TVar> &rhs,TPZVec<int> &source, TPZVec<int> &destination) {
	if(rhs.Cols() != this->Cols() && source.NElements()) {
		PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
		DebugStop();
		return;
	}
	int ncol = this->Cols();
	int nrow = source.NElements();
	int i,j;
	for(j=0; j<ncol; j++) {
		for(i=0; i<nrow; i++) {
			operator()(destination[i],j) += rhs(source[i],j);
		}
	}
}
/*******************************/
/*** Operator+( TPZFMatrix>& ) ***/

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator+(const TPZFMatrix<TVar> &A ) const {
	if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
		Error( "Operator+ <matrixs with different dimensions>" );
	
	TPZFMatrix<TVar> res;
	res.Redim( this->Rows(), this->Cols() );
	long size = ((long)this->Rows()) * this->Cols();
	TVar * pm = fElem, *plast = fElem+size;
	TVar * pa = A.fElem;
	TVar * pr = res.fElem;
	
	while(pm < plast) *pr++ = (*pm++) + (*pa++);
	
	return( res );
}

/*******************************/
/*** Operator-( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator-(const TPZFMatrix<TVar> &A ) const {
	if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
		Error( "Operator- <matrixs with different dimensions>" );
	
	TPZFMatrix<TVar> res;
    res.Redim( this->Rows(), this->Cols() );
	long size = ((long)this->Rows()) * this->Cols();
	TVar * pm = fElem;
	TVar * pa = A.fElem;
	TVar * pr = res.fElem, *prlast =pr+size;
	
	while(pr < prlast) *pr++ = (*pm++) - (*pa++);
	return( res );
}

template<>
void TPZFMatrix<int>::GramSchmidt(TPZFMatrix<int> &Orthog, TPZFMatrix<int> &TransfToOrthog)
{
	std::cout << "Nothing to do\n";
}

template <class TVar>
void TPZFMatrix<TVar>::GramSchmidt(TPZFMatrix<TVar> &Orthog, TPZFMatrix<TVar> &TransfToOrthog)
{
#ifdef LOG4CXX2
	{
		std::stringstream sout;
		Print("GrSchmidt Entrada",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
    double scale = 1.;
    for(int j = 0; j < this->Cols(); j++){
		double norm = 0.;
		for(int i = 0; i < this->Rows(); i++){
			norm += fabs(this->GetVal(i,j)*this->GetVal(i,j));
		}
		norm = sqrt(norm);
		if(norm > 1e-10){
			if(1./norm > scale) scale = 1./norm;
		}
    }
	
    this->operator *=( scale );
	
    int QTDcomp = this->Rows();
    int QTDvec = this->Cols();
    Orthog.Resize(QTDcomp,QTDvec);
    Orthog.Zero();
    /// Making a copy of *this (Ortog = *this)
    for(int r = 0; r < QTDcomp; r++)
    {
        for(int c = 0; c < QTDvec; c++)
        {
            Orthog(r,c) = GetVal(r,c);
        }
    }
	
#ifdef DEBUG
    int check = 0;
    for(int c = 0; c < QTDvec; c++)
    {
        double summ = 0.;
        for(int r = 0; r < QTDcomp; r++)
        {
            summ += fabs(GetVal(r,c));
        }
        if(fabs(summ) < 0.00001)
        {
			std::stringstream sout;
            sout << "Null Vector on Gram-Schmidt Method! Col = " << c << "\n";
			LOGPZ_ERROR(logger,sout.str())
            check = 1;
        }
    }
#endif
	
    TVar dotUp, dotDown;
    for(int c = 1; c < QTDvec; c++)
    {
        for(int stop = 0; stop < c; stop++)
        {
            dotUp = 0.;
            dotDown = 0.;
            for(int r = 0; r < QTDcomp; r++)
            {
                dotUp += GetVal(r,c)*Orthog(r,stop);
                dotDown += Orthog(r,stop)*Orthog(r,stop);
            }
            if(fabs(dotDown) < 1.E-8) 
            { 
#ifdef DEBUG
                if(check == 0)
                {
					std::stringstream sout;
                    sout << "Parallel Vectors on Gram-Schmidt Method! Col = " << stop << "\n";
					LOGPZ_ERROR(logger,sout.str())
                }
#endif
				
                for(int r = 0; r < QTDcomp; r++) 
                { 
                    Orthog(r,stop) = 0.; 
                }
            }
            else
            {
#ifdef LOG4CXX2
				{
					std::stringstream sout;
					sout << "dotdown = " << dotDown << " dotup = " << dotUp;
					LOGPZ_DEBUG(logger,sout.str())
				}
#endif
                for(int r = 0; r < QTDcomp; r++)
                {
                    Orthog(r,c) -= dotUp*Orthog(r,stop)/dotDown;
                }
            }
        }
    }
    for(int c = 0; c < QTDvec; c++)
    {
        dotUp = 0.; 
        for(int r = 0; r < QTDcomp; r++)
        {
            dotUp += Orthog(r,c)*Orthog(r,c);
        }
        if(fabs(dotUp) > 1.e-8)
        {
            for(int r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = Orthog(r,c)/sqrt(dotUp);
            }
        }
		else {
#ifdef LOG4CXX
			std::stringstream sout;
			sout << "Linearly dependent columns dotUp = " << dotUp;
			LOGPZ_ERROR(logger,sout.str())
#endif
            for(int r = 0; r < QTDcomp; r++)
            {
                Orthog(r,c) = 0.;
            }			
		}
		
    }
    Orthog.Multiply(*this,TransfToOrthog,1);
	
    this->operator *=( 1./scale );
    TransfToOrthog.operator *=( 1./scale );
	
#ifdef LOG4CXX2
	{
		std::stringstream sout;
		sout << "Output this = ";
		Print("Output GS",sout);
		Orthog.Print("Orthog matrix",sout);
		TransfToOrthog.Print("TransfToOrthog matrix",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
#ifdef DEBUG
	TPZFNMatrix<9, TVar> OrthogT;
	Orthog.Transpose(&OrthogT);
	TPZAxesTools<TVar>::VerifyAxes(OrthogT);
#endif
}

template <class TVar>
void TPZFMatrix<TVar>::DeterminantInverse(TVar &determinant, TPZFMatrix<TVar> &inverse)
{
	TPZFNMatrix<100, TVar> copy(*this);
	inverse.Redim(this->Rows(),this->Rows());
	int r;
	for(r=0; r<this->Rows(); r++) inverse(r,r) = 1.;
	copy.Solve_LU(&inverse);
	determinant = 1.;
	for(r=0; r<this->Rows(); r++) determinant *= copy(r,r);
}


template <class TVar>
void TPZFMatrix<TVar>::ConstMultiply(const TPZFMatrix<TVar> & x,TPZFMatrix<TVar> & B,const int opt) const{
	if (!opt){
		if (this->Cols() != x.Rows()){
			Error( "Error in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied" );
			std::cout << "\nError in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied\n";
		}
		const int BRows = this->Rows();
		const int BCols = x.Cols();
		const int KSize = this->Cols();
		TVar sum;
		B.Resize(BRows, BCols);
		int i, j, k;
		for(i = 0; i < BRows; i++){
			for(j = 0; j < BCols; j++){
				sum = 0.;
				for(k = 0; k < KSize; k++){
					sum += this->g(i,k) * x.g(k,j);
				}
				B.s(i,j) = sum;
			}//for j
		}//for i
	}
	else{
		if (this->Rows() != x.Rows()){
			Error( "Error in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied" );
			std::cout << "\nError in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied\n";
		}
		const int BRows = this->Cols();
		const int BCols = x.Cols();
		const int KSize = this->Rows();
		TVar sum;
		B.Resize(BRows, BCols);
		int i, j, k;
		for(i = 0; i < BRows; i++){
			for(j = 0; j < BCols; j++){
				sum = 0.;
				for(k = 0; k < KSize; k++){
					sum += this->g(k,i) * x.g(k,j);
				}
				B.s(i,j) = sum;
			}//for j
		}//for i
	}
}//void

template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TVar *ptr, int rows, int cols, const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
							   const TVar alpha,const TVar beta ,const int opt ,const int stride) 
{
	if ((!opt && cols*stride != x.Rows()) || (opt && rows*stride != x.Rows())) {
		Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
		return;
	}
	if(beta != (TVar)0. && ((!opt && rows*stride != y.Rows()) || (opt && cols*stride != y.Rows()) || y.Cols() != x.Cols())) {
		Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
		return;
	}
	if(!opt) {
		if(z.Cols() != x.Cols() || z.Rows() != rows*stride) {
			z.Redim(rows*stride,x.Cols());
		}
	} else {
		if(z.Cols() != x.Cols() || z.Rows() != cols*stride) {
			z.Redim(cols*stride,x.Cols());
		}
	}
	unsigned numeq = opt ? cols : rows;
	long xcols = x.Cols();
	int ic, c;
	if(!(rows*cols)) return;
	for (ic = 0; ic < xcols; ic++) {
		TVar *zp = &z(0,ic), *zlast = zp+numeq*stride;
		if(beta != (TVar)0.) {
			const TVar *yp = &y.g(0,ic);
			if(beta != (TVar)1. || (&z != &y && stride != 1)) {
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			} else if(&z != &y) {
				memcpy(zp,yp,numeq*sizeof(TVar));
			}
		} else {
			while(zp != zlast) {
				*zp = 0.;
				zp += stride;
			}
		}
	}
	
	
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			for ( c = 0; c<cols; c++) {
				TVar * zp = &z(0,ic), *zlast = zp+rows*stride;
				const TVar * fp = ptr +rows*c;
				const TVar * xp = &x.g(c*stride,ic);
				while(zp < zlast) {
					*zp += alpha* *fp++ * *xp;
					zp += stride;
				}
			}
		} else {
			const TVar * fp = ptr;
			TVar *zp = &z(0,ic);
			for (c = 0; c<cols; c++) {
				TVar val = 0.;
				// bug correction philippe 5/2/97
				//					 REAL * xp = &x(0,ic), xlast = xp + numeq*stride;
				const TVar *xp = &x.g(0,ic);
				const TVar *xlast = xp + rows*stride;
				while(xp < xlast) {
					val += *fp++ * *xp;
					xp += stride;
				}
				*zp += alpha *val;
				zp += stride;
			}
		}
	}
}

template <class TVar>
void TPZFMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
							   const TVar alpha,const TVar beta,const int opt,const int stride) const {
	if ((!opt && this->Cols()*stride != x.Rows()) || (opt && this->Rows()*stride != x.Rows())) {
		Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
		return;
	}
	if(beta != (TVar)0. && ((!opt && this->Rows()*stride != y.Rows()) || (opt && this->Cols()*stride != y.Rows()) || y.Cols() != x.Cols())) {
		Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
		return;
	}
	if(!opt) {
		if(z.Cols() != x.Cols() || z.Rows() != this->Rows()*stride) {
			z.Redim(this->Rows()*stride,x.Cols());
		}
	} else {
		if(z.Cols() != x.Cols() || z.Rows() != this->Cols()*stride) {
			z.Redim(this->Cols()*stride,x.Cols());
		}
	}
	if(this->Cols() == 0)
	{
		z.Zero();
	}
	unsigned numeq = opt ? this->Cols() : this->Rows();
	long rows = this->Rows();
	long cols = this->Cols();
	long xcols = x.Cols();
	int ic, c;
    if (numeq)
    {
        for (ic = 0; ic < xcols; ic++) {
            TVar *zp = &z(0,ic), *zlast = zp+numeq*stride;
            if(beta != (TVar)0.) {
                const TVar *yp = &y.g(0,ic);
                if(beta != (TVar)1. || (&z != &y && stride != 1)) {
                    while(zp < zlast) {
                        *zp = beta * (*yp);
                        zp += stride;
                        yp += stride;
                    }
                } else if(&z != &y) {
                    memcpy(zp,yp,numeq*sizeof(TVar));
                }
            } else {
                while(zp != zlast) {
                    *zp = 0.;
                    zp += stride;
                }
            }
        }
    }
	
	if(!(rows*cols)) return;
	
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			for ( c = 0; c<cols; c++) {
				TVar * zp = &z(0,ic), *zlast = zp+rows*stride;
				TVar * fp = fElem +rows*c;
				const TVar * xp = &x.g(c*stride,ic);
				while(zp < zlast) {
					*zp += alpha* *fp++ * *xp;
					zp += stride;
				}
			}
		} else {
			TVar * fp = fElem,  *zp = &z(0,ic);
			for (c = 0; c<cols; c++) {
				TVar val = 0.;
				// bug correction philippe 5/2/97
				//					 REAL * xp = &x(0,ic), xlast = xp + numeq*stride;
				const TVar *xp = &x.g(0,ic);
				const TVar *xlast = xp + rows*stride;
				while(xp < xlast) {
					val += *fp++ * *xp;
					xp += stride;
				}
				*zp += alpha *val;
				zp += stride;
			}
		}
	}
}

/********************************/
/*** Operator+=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> & TPZFMatrix<TVar>::operator+=(const TPZFMatrix<TVar> &A ) {
	if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
		Error( "Operator+= <matrixs with different dimensions>" );
	
	long size = ((long)this->Rows()) * this->Cols();
	TVar * pm = fElem, *pmlast=pm+size;
	TVar * pa = A.fElem;
	while(pm < pmlast) (*pm++) += (*pa++);
	return( *this );
}

/*******************************/
/*** Operator-=( TPZFMatrix<>& ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator-=(const TPZFMatrix<TVar> &A ) {
	if ( (A.Rows() != this->Rows())  ||  (A.Cols() != this->Cols()) )
		Error( "Operator-= <matrixs with different dimensions>" );
	
	long size = ((long)this->Rows()) * this->Cols();
	TVar * pm = fElem;
	TVar * pa = A.fElem;
	
	for ( long i = 0; i < size; i++ ) *pm++ -= *pa++;
	
	return( *this );
}

template <class TVar>
void TPZFMatrix<TVar>::ZAXPY(const TVar alpha,const TPZFMatrix<TVar> &p) {
	
#ifndef USING_ATLAS
#ifndef USING_BLAS
	TVar * pt = fElem;
	TVar * pp = p.fElem;
	TVar * ptlast = fElem + this->fRow*this->fCol;
	while(pt < ptlast) *pt++ += alpha * *pp++;
#endif
#endif
#ifdef USING_ATLAS
	//Como definir o tamanho dos vetores
	int size  = (fRow*fCol) ;
	cblas_daxpy(size, alpha, &p.fElem[0], 1, &fElem[0], 1);
#endif
#ifdef USING_BLAS
	//Como definir o tamanho dos vetores
	int size  = (fRow*fCol) ;
	cblas_daxpy(size, alpha, &p.fElem[0], 1, &fElem[0], 1);
#endif
}

template<class TVar>
void TPZFMatrix<TVar>::TimesBetaPlusZ(const TVar beta,const TPZFMatrix<TVar> &z) {
#ifdef USING_ATLAS
	int size = fRow*fCol;
	cblas_dscal(size,beta,fElem,1);
	cblas_daxpy(size,1.,z.fElem,1,fElem,1);
#else
	
	TVar * pt = fElem,  *ptlast = fElem + this->fRow*this->fCol;
	TVar * pz = z.fElem;
	//	while(pt < ptlast) *pt++ *= (beta) + *pz++;
	while(pt < ptlast) {
		*pt *= (beta);
		*pt++ += *pz++;
	}
#endif
}

/******** Operacoes com MATRIZES GENERICAS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator=(const TPZMatrix<TVar> &A ) {
	int arows  = A.Rows();
	int acols  = A.Cols();
	int size = arows * acols;
	if(fElem != fGiven) {
		delete []fElem;
		fElem = 0;
	}
	this->fRow  =  arows;
	this->fCol  = acols;
	if(fSize < size) {
		fElem = new TVar[ arows * acols ] ;
	} else {
		fElem = fGiven;
	}
	TVar * dst = fElem;
	for ( int c = 0; c < this->fCol; c++ )
		for ( int r = 0; r < this->fRow; r++ ) 
			*dst++ = A.Get( r, c );
	return( *this );
}


/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/
template <class TVar>
TPZFMatrix<TVar>& TPZFMatrix<TVar>::operator=(const TVar value ) {
	long size = ((long)this->fRow) * this->fCol;
	TVar * dst   = fElem;
	for ( long i = 0; i < size; i++ )
		*dst++ = value;
	this->fDecomposed = 0;
	return *this;
}



/***************************/
/*** Operator+=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator+=(const TVar value ) {
	long size = ((long)this->Rows()) * this->Cols();
	
	TVar * dst = fElem, *dstlast = dst+size;
	while ( dst < dstlast ) *dst++ += value;
	return( *this );
}



/**************************/
/*** Operator+( value ) ***/
template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator+(const TVar value ) const {
	TPZFMatrix<TVar> res( *this );
	long size = ((long)this->Rows()) * this->Cols();
	
	TVar * dst = res.fElem,  *dstlast = dst+size;
	while ( dst < dstlast )
		*dst++ += value;
	
	return( res );
}

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator-  (const TVar val ) const {
	return operator+( -val ); 
}


/**************************/
/*** Operator*( value ) ***/

template <class TVar>
TPZFMatrix<TVar> TPZFMatrix<TVar>::operator*(const TVar value ) const
{
	TPZFMatrix<TVar> res( *this );
	res *= value;
	return( res );
}

/***************************/
/*** Operator*=( value ) ***/
template <class TVar>
TPZFMatrix<TVar> &TPZFMatrix<TVar>::operator*=( const TVar value ) {
	long size = ((long)this->Rows()) * this->Cols();
	TVar * dst = fElem, *dstlast = dst+size;
	while ( dst < dstlast ) *dst++ *= value;
	return( *this );
}

/**************/
/*** Resize ***/
template <class TVar>
int TPZFMatrix<TVar>::Resize(const int newRows,const int newCols) {
	if ( newRows == this->Rows() && newCols == this->Cols() ) return( 1 );
	long newsize = ((long)newRows)*newCols;
	TVar * newElem;
	if(fGiven && fElem != fGiven && newsize <= fSize) 
	{
		newElem = fGiven;
	} else 
	{
		newElem = new TVar[ newRows * newCols ] ;
	}
	if ( newElem == NULL )
		Error( "Resize <memory allocation error>." );
	
	long minRow  = ( this->fRow < newRows ? this->fRow : newRows );
	long minCol  = ( this->fCol < newCols ? this->fCol : newCols );
	TVar * src;
	TVar * dst;
	long r, c;
	
	for ( c = 0; c < minCol; c++ ) {
		// Copia as linhas da matriz antiga para a nova.
		// Copia os elementos de uma linha.
		dst = newElem + c*newRows;
		src = fElem + c*this->fRow;
		for ( r = 0; r < minRow; r++ ) *dst++ = *src++;
		
		// Se a nova linha for maior (mais colunas), preenche o
		//  resto da linha com ZEROS.
		for ( ; r < newRows; r++ ) *dst++ = 0.0;
	}
	
	// Preenche as linha que sobrarem (se sobrarem) com ZEROS.
	for ( ;c < newCols; c++ ) 
	{
		dst = newElem + c*newRows;
		for (r = 0 ; r < newRows; r++ ) *dst++ = 0.0;
	}
	
	if (fElem && fElem != fGiven )delete[]( fElem );
	fElem = newElem;
	this->fRow  = newRows;
	this->fCol  = newCols;
	return( 1 );
}

template <class TVar>
int TPZFMatrix<TVar>::Remodel(const int newRows,const int newCols) {
	if(newRows*newCols != this->fRow*this->fCol) return -1;
	this->fRow = newRows;
	this->fCol = newCols;
	return 1;
}

/********************/
/*** Transpose () ***/
template <class TVar>
void TPZFMatrix<TVar>::Transpose(TPZMatrix<TVar> *const T) const{
	T->Resize( this->Cols(), this->Rows() );
	//Transposta por filas
	TVar * p = fElem;
	for ( int c = 0; c < this->Cols(); c++ ) {
		for ( int r = 0; r < this->Rows(); r++ ) {
			T->PutVal( c, r, *p++ );
			//            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		}
	}
}

template <class TVar>
void TPZFMatrix<TVar>::Transpose() {
	TPZFMatrix<TVar> temp;
	Transpose(&temp);
	*this = temp;
}

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(TPZVec<int> &index) {
	
	if (this->fDecomposed) return 0;
	
	if ( this->Rows() != this->Cols() ) {
		cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
		return 0;
	}
	
	int i,j,k;
	TVar sum = 0.;
	int nRows = this->Rows();
	int nCols = this->Cols();
	
	index.Resize(nRows);
	//inicializo o vetor de índices para o caso de pivotamento  
	for (i=0;i<nRows;i++) index[i] = i;
	
	//LU
	for (j=0;j<nCols;j++){
		// cout << "line..." << j << endl;
		for (i=0;i<=j;i++){
			sum = 0.;
			for (k=0;k<i;k++){
				sum += this->Get(i,k) * this->Get(k,j);
			}
			TVar aux = this->Get(i,j);
			PutVal(i,j,(aux - sum));
			//cout << "0_A[" << i << "," << j << "]= " << Get(i,j) << endl;
		}
		//Print(cout);
		TVar piv = this->Get(j,j);
		//  cout << "Pivo 1 =" << piv << endl;
		int row = j;
		for (i=j+1;i<nRows;i++){
			sum = 0.;
			for (k=0;k<(j);k++){
				sum += this->Get(i,k) * this->Get(k,j);
			}
			TVar aux = this->Get(i,j);
			PutVal(i,j,(aux - sum));
			//cout << "1_A[" << i << "," << j << "]= " << Get(i,j) << endl;
			
			if (fabs(this->Get(i,j)) > fabs(piv)){
				piv = this->Get(i,j);
				//  cout << "Pivo 2 =" << piv << endl;
				row = i;
			}
		}
		//    Print(cout);
		if (row > j){
			for (k=0;k<nCols;k++){
				TVar aux = this->Get(j,k);
				PutVal(j,k,this->Get(row,k));
				//cout << "2_A[" << j << "," << k << "]= " << Get(j,k) << endl;
				PutVal(row,k,aux);
				//cout << "3_A[" << row << "," << k << "]= " << Get(row,k) << endl;
			}
			k = index[j];
			index[j] = index[row];
			index[row] = k;
		}
		//    cout << "Pivo = " << piv << endl;
		for (i=j+1;i<nRows;i++){
			if (fabs(piv) < 1e-12) { 
				cout << "Pivot < 1e-12. Probably matrix is singular." << endl;
				DebugStop();
			}
			TVar aux = this->Get(i,j) / piv;
			PutVal(i,j,aux);
			//cout << "4_A[" << i << "," << j << "]= " << Get(i,j) << endl;
		}
		//Print(cout);
	}    
	this->fDecomposed = ELUPivot;
	return 1;
}


/*****************/
/*** DecomposeLU ***/
template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU(std::list<int> &singular) {
	return Decompose_LU();
}

template <class TVar>
int TPZFMatrix<TVar>::Decompose_LU() {
	
	if (  this->fDecomposed && this->fDecomposed != ELU)  Error( "Decompose_LU <TPZFMatrix<>already Decomposed with other scheme>" );
	if (this->fDecomposed) return 1;
	TVar nn;
	TVar * ptrpivot,*pik, *pij,*pkj;
	
	int i,j,k,rows=this->Rows(),cols=this->Cols();
	int  min = ( cols < (rows) ) ? cols : rows;
	
	ptrpivot=&fElem[0];
	for (  k = 0; k < min ; k++ )
	{
		if ( IsZero( *ptrpivot ) ){
			//Power plus...
			if (fabs(*ptrpivot) > 0){
				for (j=k+1;j<rows;j++){
					if (fabs(*(ptrpivot + j - k) - *(ptrpivot)) > 1e-12){
						Error( "DecomposeLU <matrix is singular> even after Power Plus..." );
						cout << "DecomposeLU <matrix is singular> even after Power Plus...\n" ;
					}
				}
			}
			else{
				Error( "DecomposeLU <matrix is singular>" );
				cout << "DecomposeLU <matrix is singular>\n";
			}
		}
		pik=ptrpivot;
		for ( i = k+1; i < rows; i++ )
		{
			pik+=1;
			nn = (*pik)/(*ptrpivot);
			(*pik)=nn;
			pkj=&fElem[k*cols+k];
			pij=&fElem[k*cols+i];
			for ( j = k+1; j < this->Cols(); j++ )
			{
				pkj+=cols;pij+=cols;
				(*pij)-=nn*(*pkj);
			}
		}
		ptrpivot+=this->Rows()+1;
	}
	
	this->fDecomposed=1;
	return 1; 
}

template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int rows, TPZFMatrix<TVar> *B)
{
	int rowb = B->Rows();
	int colb = B->Cols();
	if ( rowb != rows ) Error( "static::SubstitutionLU <incompatible dimensions>" );
	int i,j;
	for ( i = 0; i < rowb; i++ ) {
		for ( int col = 0; col < colb; col++ )
			for (j = 0; j < i; j++ )
				//B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
				PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
	}
	
	for (int col=0; col<colb; col++){
		for ( i = rowb-1; i >= 0; i-- ) {
			for (j = i+1; j < rowb ; j++ )
				//B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
				PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, j) * GETVAL(B, rowb, j, col));
			if ( IsZero( SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/ ) ) {
				if (fabs(SELECTEL(ptr, rows, i, i)/*GetVal(i, i)*/) > 0.){
					if (fabs(GETVAL(B, rowb, i, col) - SELECTEL(ptr, rows, i, i)/*B->GetVal(i, col) - GetVal(i, i)*/) > 1e-12){
						Error( "static::BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
					}
				}else  Error( "static::BackSub(SubstitutionLU) <Matrix is singular" );
			}
			PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/SELECTEL(ptr, rows, i, i));
			//B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
	}
	return( 1 );
}

/****************/
/*** Substitution ***/

template <class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const {
	
#ifndef DEBUG2
	if(this->fDecomposed != ELU) {
		Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
	}
	int rowb = B->Rows();
	int colb = B->Cols();
	int row = this->Rows();
	if ( rowb != this->Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
	
	
	int i,j;
	for ( i = 0; i < rowb; i++ ) {
		for ( int col = 0; col < colb; col++ )
			for (j = 0; j < i; j++ )
				//B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
				PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
	}
	
	for (int col=0; col<colb; col++){
		for ( i = rowb-1; i >= 0; i-- ) {
			for (j = i+1; j < rowb ; j++ )
				//B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
				PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col) - GETVAL(this, row, i, j) * GETVAL(B, rowb, j, col));
			if ( IsZero( GETVAL(this, row, i, i)/*GetVal(i, i)*/ ) ) {
				if (fabs(GETVAL(this, row, i, i)/*GetVal(i, i)*/) > 0.){
					if (fabs(GETVAL(B, rowb, i, col) - GETVAL(this, row, i, i)/*B->GetVal(i, col) - GetVal(i, i)*/) > 1e-12){
						Error( "BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
					}
				}else  Error( "BackSub(SubstitutionLU) <Matrix is singular" );
			}
			PUTVAL(B, rowb, i, col, GETVAL(B, rowb, i, col)/GETVAL(this, row, i, i));
			//B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
	}
	return( 1 );
	
#else
	
	if(this->fDecomposed != ELU) {
		Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
	}
	int rowb = B->Rows();
	int colb = B->Cols();
	if ( rowb != this->Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
	
	
	int i,j;
	for ( i = 0; i < rowb; i++ ) {
		for ( int col = 0; col < colb; col++ )
			for (j = 0; j < i; j++ )
				B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
	}
	
	for (int col=0; col<colb; col++){
		for ( i = rowb-1; i >= 0; i-- ) {
			for (j = i+1; j < rowb ; j++ )
				B->PutVal( i, col, B->GetVal(i, col) - GetVal(i, j) * B->GetVal(j, col) );
			if ( IsZero( GetVal(i, i) ) ) {
				if (fabs(GetVal(i, i)) > 0.){
					if (fabs(B->GetVal(i, col) - GetVal(i, i)) > 1e-12){
						Error( "BackSub(SubstitutionLU) <Matrix is singular even after Power Plus..." );
					}
				}else  Error( "BackSub(SubstitutionLU) <Matrix is singular" );
			}
			B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
	}
	return( 1 );
	
#endif
}

template<class TVar>
int TPZFMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B, TPZVec<int> &index ) const{
	
	if(!B){
		PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
		return 0;
	}
	
	TPZFMatrix<TVar> &b = *B;
	
	if (!this->fDecomposed){
		PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
		return 0;
	}
	
	if (this->fDecomposed != ELUPivot){
		PZError << __PRETTY_FUNCTION__ << "\nfDecomposed != ELUPivot" << endl;
	}
	
	int nRows = this->Rows();
	
	if (index.NElements() != nRows || b.Rows() != nRows)
	{
		cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
		<< "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
		return 0;
	}
	
	int i,j;
	TVar sum = 0;
	
	TPZVec<TVar> v(nRows);
	
	
	for (i=0;i<nRows;i++)
	{
		v[i] = b(index[i]);
	}
	
	//Ly=b
	for (i=0;i<nRows;i++)
	{
		sum = 0.;
		for (j=0;j<(i);j++) sum +=this->Get(i,j) * v[j];
		v[i] -= sum;
	}
	
	//Ux=y
	for (i=(nRows-1);i>-1;i--)
	{
		sum = 0.;
		for (j=(i+1);j<nRows;j++) sum += this->Get(i,j) * v[j];
		v[i] = (v[i] - sum) / this->Get(i,i);
	}
	
	for (i=0;i<nRows;i++) b(i) = v[i];
	return 1;
}

template <class TVar>
int TPZFMatrix<TVar>::Substitution(const TVar *ptr, int rows, TPZFMatrix<TVar> *B, TPZVec<int> &index )
{
	
	if(!B){
		PZError << __PRETTY_FUNCTION__ << "TPZFMatrix<>*B eh nulo" << endl;
		return 0;
	}
	
	TPZFMatrix<TVar> &b = *B;
	
	
	
	if (index.NElements() != rows || b.Rows() != rows)
	{
		cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
		<< "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
		return 0;
	}
	
	int i,j;
	TVar sum = 0;
	
	TPZVec<TVar> v(rows);
	
	
	for (i=0;i<rows;i++)
	{
		v[i] = b(index[i]);
	}
	
	//Ly=b
	for (i=0;i<rows;i++)
	{
		sum = 0.;
		for (j=0;j<(i);j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
		v[i] -= sum;
	}
	
	//Ux=y
	for (i=(rows-1);i>-1;i--)
	{
		sum = 0.;
		for (j=(i+1);j<rows;j++) sum += SELECTEL(ptr,rows,i,j) * v[j];
		v[i] = (v[i] - sum) / SELECTEL(ptr,rows,i,i);
	}
	
	for (i=0;i<rows;i++) b(i) = v[i];
	return 1;
}

/** @brief Implement dot product for matrices */
template<class TVar>
TVar Dot(const TPZFMatrix<TVar> &A, const TPZFMatrix<TVar> &B) {
	int size = (A.Rows())*A.Cols();
	TVar result = 0.;
	if(!size) return result;
#ifdef USING_ATLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;
	
#elif USING_BLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;
	
#else
	const TVar *fpA = &A.g(0,0), *fpB = &B.g(0,0);
	const TVar *fpLast = fpA+size;
	while(fpA < fpLast) 
	{
		result += (*fpA++ * *fpB++); 
	}
	return result;
#endif
}

template
std::complex<float> Dot(const TPZFMatrix< std::complex<float> > &A, const TPZFMatrix< std::complex<float> > &B);

template
std::complex<double> Dot(const TPZFMatrix< std::complex<double> > &A, const TPZFMatrix< std::complex<double> > &B);

template
std::complex<long double> Dot(const TPZFMatrix< std::complex<long double> > &A, const TPZFMatrix< std::complex<long double> > &B);

template
long double Dot(const TPZFMatrix<long double> &A, const TPZFMatrix<long double> &B);

template
double Dot(const TPZFMatrix<double> &A, const TPZFMatrix<double> &B);

template
float Dot(const TPZFMatrix<float> &A, const TPZFMatrix<float> &B);

template
int Dot(const TPZFMatrix<int> &A, const TPZFMatrix<int> &B);

/** @brief Increments value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator+(const TVar value, const TPZFMatrix<TVar> &A ) {
	return( A + value );
}

/** @brief Decrements value over all entries of the matrix A. */
template <class TVar>
TPZFMatrix<TVar> operator-(const TVar value, const TPZFMatrix<TVar> &A ) {
	return( A - value );
}

/************************** Private **************************/

/*************/
/*** Error ***/
template<class TVar>
int TPZFMatrix<TVar>::Error(const char *msg1,const char *msg2 ) {
	ostringstream out;
	out << "TPZFMatrix::" << msg1;
	if(msg2) out << msg2;
	out << ".\n";
	LOGPZ_ERROR (logger, out.str().c_str());
	DebugStop();
	return 0;
}


/*************/
/*** Clear ***/
template <class TVar>
int TPZFMatrix<TVar>::Clear() {
	if(fElem && fElem != fGiven) delete[]( fElem );
	fElem = NULL;
	this->fRow  = this->fCol = 0;
	return( 1 );
}

template <class TVar>
void TPZFMatrix<TVar>::Read( TPZStream &buf, void *context ){
	TPZMatrix<TVar>::Read(buf,context);
	int row = this->fRow;
	int col = this->fCol;
	this->fRow = this->fCol = 0;
	Resize(row,col);
	buf.Read(fElem,this->fRow*this->fCol);
}

template <class TVar>
void TPZFMatrix<TVar>::Write( TPZStream &buf, int withclassid ) {
	TPZMatrix<TVar>::Write(buf,withclassid);
	buf.Write(fElem,this->fRow*this->fCol);
}

template <>
void TPZFMatrix<float>::Read( TPZStream &buf, void *context ){
    DebugStop();
}

template <>
void TPZFMatrix<float>::Write( TPZStream &buf, int withclassid ) {
    DebugStop();
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSaveable *copy, bool override)
{
	TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
	if(!fmat) return false;
	
	bool matresult = TPZMatrix<TVar>::Compare(copy,false);
	int nel = this->fRow*this->fCol;
	TVar diff=0.;
	int numdif = 0;
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		if(fElem[iel] != fmat->fElem[iel]) 
		{
			matresult = false;
			numdif++;
		}
		diff += fabs(fElem[iel]-fmat->fElem[iel]);
	}
	if(!matresult)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " did not compare ";
		sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
		sout << " difference in norm L1 " << diff;
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(!matresult && override) 
	{
		this->operator=(*fmat);
	}
	return matresult;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template<class TVar>
bool TPZFMatrix<TVar>::Compare(TPZSaveable *copy, bool override) const
{
	TPZFMatrix<TVar> *fmat = dynamic_cast<TPZFMatrix<TVar> *> (copy);
	if(!fmat) return false;
    
	bool matresult = TPZMatrix<TVar>::Compare(copy,false);
	int nel = this->fRow*this->fCol;
	TVar diff=0.;
	int numdif = 0;
	int iel;
	for(iel=0; iel<nel; iel++)
	{
		if(fElem[iel] != fmat->fElem[iel]) 
		{
			matresult = false;
			numdif++;
		}
		diff += fabs(fElem[iel]-fmat->fElem[iel]);
	}
	if(!matresult)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " did not compare ";
		sout << " number different terms " << numdif << " number terms " << this->fRow*this->fCol;
		sout << " difference in norm L1 " << diff;
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(!matresult && override) 
	{
		DebugStop();
	}
	return matresult;
}

template <class TVar>
void TPZFMatrix<TVar>::PrintStatic(const TVar *ptr, int rows, int cols, const char *name, std::ostream& out,const MatrixOutputFormat form){
	
	if(form == EFormatted) {
		out << "Writing matrix '";
		if(name) out << name;
		out << "' (" << rows << " x " << cols << "):\n";
		
		for ( int row = 0; row < rows; row++) {
			out << "\t";
			for ( int col = 0; col < cols; col++ ) {
				out << SELECTEL(ptr,rows, row, col) << "  ";
			}
			out << "\n";
		}
    	out << "\n";
	} else if (form == EInputFormat) {
		out << rows << " " << cols << endl;
		for ( int row = 0; row < rows; row++) {
			for ( int col = 0; col < cols; col++ ) {
				TVar val = SELECTEL(ptr,rows,row, col);
				if(val != (TVar)0.) out << row << ' ' << col << ' ' << val << std::endl;
			}
		}
		out << "-1 -1 0.\n";
	} else if( form == EMathematicaInput)
	{
		char number[32];
		out << name << "\n{ ";
		for ( int row = 0; row < rows; row++) {
			out << "\n{ ";
			for ( int col = 0; col < cols; col++ ) {
				TVar val = SELECTEL(ptr,rows,row, col);
				sprintf(number, "%16.16lf", (REAL)fabs(val));
				out << number;
				if(col < cols-1)
					out << ", ";
				if((col+1) % 6 == 0)out << std::endl;
			}
			out << " }";
			if(row < rows-1)
				out << ",";
		}
		
		out << " }\n";
		
	}
	
}


template<class TVar>
int TPZFMatrix<TVar>::ClassId() const   
{ 
	return TPZFMATRIXID; 
}


template <class TVar>
int TPZFMatrix<TVar>::SetSize(const int newRows,const int newCols) {
	long newsize = ((long)newRows)*newCols;
	long oldsize = this->fRow*this->fCol;
	if(newsize == oldsize) return 1;
	if(fElem && fElem != fGiven)
	{
		delete []fElem;
		fElem = 0;
	}
	if(fGiven && newsize <= fSize) 
	{
		fElem = fGiven;
	} else 
	{
		fElem = new TVar[ newRows * newCols ] ;
	}
	if (newsize && fElem == NULL )
		Error( "Resize <memory allocation error>." );
	
	return( 1 );
}

#include <complex>
template class TPZFMatrix< std::complex<float> >;
template class TPZFMatrix< std::complex<double> >;
template class TPZFMatrix< std::complex<long double> >;

template class TPZFMatrix<long double>;
template class TPZFMatrix<double>;
template class TPZFMatrix<int>;
template class TPZFMatrix<float>;
template class TPZRestoreClass< TPZFMatrix<REAL> , TPZFMATRIXID>;
