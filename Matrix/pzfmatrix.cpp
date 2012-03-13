/**
 * @file
 * @brief Contains the implementation of the TPZFMatrix methods.
 */
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tfullmat.c
//
// Class:  TPZFMatrix
//
// Obs.:   Implementa matrizes cheias (normais).
//
// Versao: 04 / 1996.
//


#ifdef __BOORLANDC__
#include <alloc.h>
#endif

#include "pzfmatrix.h"
//#include "pztempmat.h"
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




TPZFMatrix::TPZFMatrix(const TPZMatrix &mat) : TPZMatrix(mat), fElem(0),fGiven(0),fSize(0) {
	if(fRow*fCol) {
		fElem = new REAL[fRow*fCol];
		REAL * p = fElem;
		int i,j;
		for(j=0; j<fCol; j++) {
			for(i=0; i<fRow; i++) {
				*p++ = mat.GetVal(i,j);
			}
		}
	}
}




/********************************/
/*** Constructor( TPZFMatrix& ) ***/

TPZFMatrix::TPZFMatrix (const TPZFMatrix & A)
: TPZMatrix( A.fRow, A.fCol ), fElem(0), fGiven(0), fSize(0) {
    int size = fRow * fCol;
    if(!size) return;
    fElem = new REAL[ size ] ;
#ifdef DEBUG2
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    // Copia a matriz
    REAL * src = A.fElem;
    REAL * p = fElem;
    memcpy(p,src,(size_t)size*sizeof(REAL));
}







/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/
/*
 TPZFMatrix &TPZFMatrix::operator=( TPZTempFMatrix Atemp ) {
 TPZFMatrix &A = Atemp.Object();
 if(this == &A) return *this;
 
 if (fElem && fElem != fGiven) delete[]( fElem );
 //free(fElem); //
 
 if(A.fElem != A.fGiven) {
 fRow  = A.fRow;
 fCol  = A.fCol;
 fDecomposed = A.fDecomposed;
 fDefPositive = A.fDefPositive;
 fElem = A.fElem;
 fSize = A.fSize;
 
 A.fRow = 0;
 A.fCol = 0;
 A.fDecomposed = 0;
 A.fDefPositive = 0;
 A.fElem = 0;
 A.fSize = 0;
 } else {
 *this = A;
 }
 
 return *this;
 }
 */
/******************/
/*** Operator = ***/

/*
 TPZFMatrix::TPZFMatrix( TPZTempFMatrix Atemp ) : fGiven(Atemp.Object().fGiven) {
 TPZFMatrix &A = Atemp.Object();
 
 fRow  = A.fRow;
 fCol  = A.fCol;
 fDecomposed = A.fDecomposed;
 fDefPositive = A.fDefPositive;
 fElem = A.fElem;
 fSize = A.fSize;
 
 A.fRow = 0;
 A.fCol = 0;
 A.fDecomposed = 0;
 A.fDefPositive = 0;
 A.fElem = 0;
 A.fSize = 0;
 
 }
 */
TPZFMatrix &TPZFMatrix::operator=(const TPZFMatrix &A ) {
	if(this == &A) return *this;
	long size = A.fRow * A.fCol;
	
	REAL * newElem = fElem;
	if(fSize < size && size != fRow*fCol) {
		newElem = new REAL[size] ;
	} else if (fSize >= size) {
		newElem = fGiven;
	}
	
	if ( newElem == NULL && size > 0) Error( "Operator= <memory allocation error>." );
	if (fElem && fElem != newElem && fElem != fGiven) delete[]( fElem );
	fRow  = A.fRow;
	fCol  = A.fCol;
	fElem = newElem;
	
	// Copia a matriz
	memcpy(fElem,A.fElem,(size_t)size*sizeof(REAL));
	
	return *this;
}

void TPZFMatrix::AddFel(TPZFMatrix &rhs,TPZVec<int> &destination) {
	if(rhs.Cols() != Cols()) {
		PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
		DebugStop();
		return;
	}
	int ncol = Cols();
	int nrow = rhs.Rows();
	int i,j;
	for(j=0; j<ncol; j++) {
		for(i=0; i<nrow; i++) {
			operator()(destination[i],j) += rhs(i,j);
		}
	}
}

void TPZFMatrix::AddFel(TPZFMatrix &rhs,TPZVec<int> &source, TPZVec<int> &destination) {
	if(rhs.Cols() != Cols() && source.NElements()) {
		PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
		DebugStop();
		return;
	}
	int ncol = Cols();
	int nrow = source.NElements();
	int i,j;
	for(j=0; j<ncol; j++) {
		for(i=0; i<nrow; i++) {
			operator()(destination[i],j) += rhs(source[i],j);
		}
	}
}
/*******************************/
/*** Operator+( TPZFMatrix & ) ***/

TPZFMatrix TPZFMatrix::operator+(const TPZFMatrix &A ) const {
	if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		Error( "Operator+ <matrixs with different dimensions>" );
	
	TPZFMatrix res;
	res.Redim( Rows(), Cols() );
	long size = ((long)Rows()) * Cols();
	REAL * pm = fElem, *plast = fElem+size;
	REAL * pa = A.fElem;
	REAL * pr = res.fElem;
	
	while(pm < plast) *pr++ = (*pm++) + (*pa++);
	
	return( res );
}



/*******************************/
/*** Operator-( TPZFMatrix & ) ***/

TPZFMatrix TPZFMatrix::operator-(const TPZFMatrix &A ) const {
	if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		Error( "Operator- <matrixs with different dimensions>" );
	
	TPZFMatrix res;
    res.Redim( Rows(), Cols() );
	long size = ((long)Rows()) * Cols();
	REAL * pm = fElem;
	REAL * pa = A.fElem;
	REAL * pr = res.fElem, *prlast =pr+size;
	
	while(pr < prlast) *pr++ = (*pm++) - (*pa++);
	return( res );
}

void TPZFMatrix::GramSchmidt(TPZFMatrix &Orthog, TPZFMatrix &TransfToOrthog)
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
			norm += this->GetVal(i,j)*this->GetVal(i,j);
		}//for i
		norm = sqrt(norm);
		if(norm > 1e-10){
			if(1./norm > scale) scale = 1./norm;
		}
    }//for j
	
    this->operator *=( scale );
	
    int QTDcomp = Rows();
    int QTDvec = Cols();
    Orthog.Resize(QTDcomp,QTDvec);
    Orthog.Zero();
    // Making a copy of *this (Ortog = *this)
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
	
    double dotUp, dotDown;
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
            if(dotDown < 1.E-8) 
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
        if(dotUp > 1.e-8)
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
	TPZFNMatrix<9> OrthogT;
	Orthog.Transpose(&OrthogT);
	TPZAxesTools::VerifyAxes(OrthogT);
#endif
}

void TPZFMatrix::DeterminantInverse(REAL &determinant, TPZFMatrix &inverse)
{
	TPZFNMatrix<100> copy(*this);
	inverse.Redim(Rows(),Rows());
	int r;
	for(r=0; r<Rows(); r++) inverse(r,r) = 1.;
	copy.Solve_LU(&inverse);
	determinant = 1.;
	for(r=0; r<Rows(); r++) determinant *= copy(r,r);
}



void TPZFMatrix::ConstMultiply(const TPZFMatrix & x,TPZFMatrix & B,const int opt) const{
	if (!opt){
		if (this->Cols() != x.Rows()){
			Error( "Error in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied" );
			std::cout << "\nError in TPZFMatrix::ConstMultiply() - matrices have wrong sizes to be multiplied\n";
		}
		const int BRows = this->Rows();
		const int BCols = x.Cols();
		const int KSize = this->Cols();
		REAL sum;
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
		REAL sum;
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

void TPZFMatrix::MultAdd(const REAL *ptr, int rows, int cols, const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
                         const REAL alpha,const REAL beta ,const int opt ,const int stride) 
{
	if ((!opt && cols*stride != x.Rows()) || (opt && rows*stride != x.Rows())) {
		Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
		return;
	}
	if(beta != 0. && ((!opt && rows*stride != y.Rows()) || (opt && cols*stride != y.Rows()) || y.Cols() != x.Cols())) {
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
		REAL *zp = &z(0,ic), *zlast = zp+numeq*stride;
		if(beta != 0.) {
			const REAL *yp = &y.g(0,ic);
			if(beta != 1. || (&z != &y && stride != 1)) {
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			} else if(&z != &y) {
				memcpy(zp,yp,numeq*sizeof(REAL));
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
				REAL * zp = &z(0,ic), *zlast = zp+rows*stride;
				const REAL * fp = ptr +rows*c;
				const REAL * xp = &x.g(c*stride,ic);
				while(zp < zlast) {
					*zp += alpha* *fp++ * *xp;
					zp += stride;
				}
			}
		} else {
			const REAL * fp = ptr;
			REAL *zp = &z(0,ic);
			for (c = 0; c<cols; c++) {
				REAL val = 0.;
				// bug correction philippe 5/2/97
				//					 REAL * xp = &x(0,ic), xlast = xp + numeq*stride;
				const REAL *xp = &x.g(0,ic);
				const REAL *xlast = xp + rows*stride;
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


void TPZFMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
						 const REAL alpha,const REAL beta,const int opt,const int stride) const {
	if ((!opt && Cols()*stride != x.Rows()) || (opt && Rows()*stride != x.Rows())) {
		Error( "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" );
		return;
	}
	if(beta != 0. && ((!opt && Rows()*stride != y.Rows()) || (opt && Cols()*stride != y.Rows()) || y.Cols() != x.Cols())) {
		Error( "TPZFMatrix::MultAdd matrix y with incompatible dimensions>" );
		return;
	}
	if(!opt) {
		if(z.Cols() != x.Cols() || z.Rows() != Rows()*stride) {
			z.Redim(Rows()*stride,x.Cols());
		}
	} else {
		if(z.Cols() != x.Cols() || z.Rows() != Cols()*stride) {
			z.Redim(Cols()*stride,x.Cols());
		}
	}
	if(Cols() == 0)
	{
		z.Zero();
	}
	unsigned numeq = opt ? Cols() : Rows();
	long rows = Rows();
	long cols = Cols();
	long xcols = x.Cols();
	int ic, c;
	if(!(rows*cols)) return;
	for (ic = 0; ic < xcols; ic++) {
		REAL *zp = &z(0,ic), *zlast = zp+numeq*stride;
		if(beta != 0.) {
			const REAL *yp = &y.g(0,ic);
			if(beta != 1. || (&z != &y && stride != 1)) {
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			} else if(&z != &y) {
				memcpy(zp,yp,numeq*sizeof(REAL));
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
				REAL * zp = &z(0,ic), *zlast = zp+rows*stride;
				REAL * fp = fElem +rows*c;
				const REAL * xp = &x.g(c*stride,ic);
				while(zp < zlast) {
					*zp += alpha* *fp++ * *xp;
					zp += stride;
				}
			}
		} else {
			REAL * fp = fElem,  *zp = &z(0,ic);
			for (c = 0; c<cols; c++) {
				REAL val = 0.;
				// bug correction philippe 5/2/97
				//					 REAL * xp = &x(0,ic), xlast = xp + numeq*stride;
				const REAL *xp = &x.g(0,ic);
				const REAL *xlast = xp + rows*stride;
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
/*** Operator+=( TPZFMatrix & ) ***/

TPZFMatrix & TPZFMatrix::operator+=(const TPZFMatrix &A ) {
	if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		Error( "Operator+= <matrixs with different dimensions>" );
	
	long size = ((long)Rows()) * Cols();
	REAL * pm = fElem, *pmlast=pm+size;
	REAL * pa = A.fElem;
	while(pm < pmlast) (*pm++) += (*pa++);
	return( *this );
}



/*******************************/
/*** Operator-=( TPZFMatrix & ) ***/

TPZFMatrix &TPZFMatrix::operator-=(const TPZFMatrix &A ) {
	if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		Error( "Operator-= <matrixs with different dimensions>" );
	
	long size = ((long)Rows()) * Cols();
	REAL * pm = fElem;
	REAL * pa = A.fElem;
	
	for ( long i = 0; i < size; i++ ) *pm++ -= *pa++;
	
	return( *this );
}

void TPZFMatrix::ZAXPY(const REAL alpha,const TPZFMatrix &p) {
	
#ifndef USING_ATLAS
#ifndef USING_BLAS
	REAL * pt = fElem;
	REAL * pp = p.fElem;
	REAL * ptlast = fElem + fRow*fCol;
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

void TPZFMatrix::TimesBetaPlusZ(const REAL beta,const TPZFMatrix &z) {
#ifdef USING_ATLAS
	int size = fRow*fCol;
	cblas_dscal(size,beta,fElem,1);
	cblas_daxpy(size,1.,z.fElem,1,fElem,1);
#else
	
	REAL * pt = fElem,  *ptlast = fElem + fRow*fCol;
	REAL * pz = z.fElem;
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

TPZFMatrix &TPZFMatrix::operator=(const TPZMatrix &A ) {
	int arows  = A.Rows();
	int acols  = A.Cols();
	int size = arows * acols;
	if(fElem != fGiven) {
		delete []fElem;
		fElem = 0;
	}
	fRow  =  arows;
	fCol  = acols;
	if(fSize < size) {
		fElem = new REAL[ arows * acols ] ;
	} else {
		fElem = fGiven;
	}
	REAL * dst = fElem;
	for ( int c = 0; c < fCol; c++ )
		for ( int r = 0; r < fRow; r++ ) 
			*dst++ = A.Get( r, c );
	return( *this );
}



/******************/
/*** Operator + ***/
/*
 TPZFMatrix TPZFMatrix::operator+(const TPZMatrix &A ) const {
 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
 Error( "Operator+ (TPZMatrix &) <different dimensions>" );
 
 TPZFMatrix res;
 res.Redim( Rows(), Cols() );
 REAL * src = fElem;
 REAL * dst = res.fElem;
 long row = Rows();
 long col = Cols();
 for ( int c = 0; c < col; c++ )
 for ( int r = 0; r < row; r++ )
 *dst++ = (*src++) + A.Get( r, c );
 
 return( res );
 }
 */


/******************/
/*** Operator - ***/
/*
 TPZFMatrix TPZFMatrix::operator-(const TPZMatrix &A ) const {
 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
 Error( "Operator+ (TPZMatrix &) <different dimensions>" );
 
 TPZFMatrix res(*this);
 res -= A;
 return( res );
 }
 */


/******************/
/*** Operator * ***/
/*
 TPZFMatrix TPZFMatrix::operator*(const TPZMatrix &A ) const {
 
 if ( Cols() != A.Rows() )
 Error( "Operator* (TPZMatrix &) <incompatible dimensions>" );
 
 
 TPZFMatrix res;
 res.Redim( Rows(), A.Cols() );
 
 int     r,c,i,
 acols=A.Cols(),
 rows=Rows(),
 cols=Cols();
 
 REAL * pr = res.fElem, 
 *pm,
 *pm_aux;
 
 for (  c = 0; c < acols; c++) {
 pm = fElem;     // Percorre os elementos desta matriz.
 for (  r = 0; r < rows; r++) {
 pm_aux = pm;
 
 // Calcula um elemento da resposta.
 for ( i = 0; i < cols; i++, pm_aux += rows ) {
 *pr += (*pm_aux) * A.Get( i, c );
 }
 pm++;
 pr++;
 }
 }
 
 return( res );
 }
 */


/*******************/
/*** Operator += ***/
/*
 TPZFMatrix &TPZFMatrix::operator+=(const TPZMatrix &A ) {
 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
 Error( "Operator+ (TPZMatrix &) <different dimensions>" );
 
 REAL * pm = fElem;
 long cols = Cols();
 long rows = Rows();
 for ( int c = 0; c < cols; c++ )
 for ( int r = 0; r < rows; r++ )
 *pm++ += A.Get( r, c );
 
 return( *this );
 }
 */


/*******************/
/*** Operator -= ***/
/*
 TPZFMatrix &TPZFMatrix::operator-=(const TPZMatrix &A ) {
 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
 Error( "Operator+ (TPZMatrix &) <different dimensions>" );
 
 REAL * pm = fElem;
 int cols = Cols();
 int rows = Rows();
 for ( int c = 0; c < cols; c++ )
 for ( int r = 0; r < rows; r++ )
 *pm++ -= A.Get( r, c );
 
 return( *this );
 }
 */


/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/

TPZFMatrix& TPZFMatrix::operator=(const REAL value ) {
	long size = ((long)fRow) * fCol;
	REAL * dst   = fElem;
	for ( long i = 0; i < size; i++ )
		*dst++ = value;
	fDecomposed = 0;
	return *this;
}



/***************************/
/*** Operator+=( value ) ***/

TPZFMatrix &TPZFMatrix::operator+=(const REAL value ) {
	long size = ((long)Rows()) * Cols();
	
	REAL * dst = fElem, *dstlast = dst+size;
	while ( dst < dstlast ) *dst++ += value;
	return( *this );
}



/**************************/
/*** Operator+( value ) ***/

TPZFMatrix TPZFMatrix::operator+(const REAL value ) const {
	TPZFMatrix res( *this );
	long size = ((long)Rows()) * Cols();
	
	REAL * dst = res.fElem,  *dstlast = dst+size;
	while ( dst < dstlast )
		*dst++ += value;
	
	return( res );
}

TPZFMatrix TPZFMatrix::operator-  (const REAL val ) const {
	return operator+( -val ); 
}


/**************************/
/*** Operator*( value ) ***/

TPZFMatrix
TPZFMatrix::operator*(const REAL value ) const
{
	TPZFMatrix res( *this );
	res *= value;
	return( res );
}



/***************************/
/*** Operator*=( value ) ***/

TPZFMatrix &TPZFMatrix::operator*=( const REAL value ) {
	long size = ((long)Rows()) * Cols();
	REAL * dst = fElem, *dstlast = dst+size;
	while ( dst < dstlast ) *dst++ *= value;
	return( *this );
}



/**************/
/*** Resize ***/


int TPZFMatrix::Resize(const int newRows,const int newCols) {
	if ( newRows == Rows() && newCols == Cols() ) return( 1 );
	long newsize = ((long)newRows)*newCols;
	REAL * newElem;
	if(fGiven && fElem != fGiven && newsize <= fSize) 
	{
		newElem = fGiven;
	} else 
	{
		newElem = new REAL[ newRows * newCols ] ;
	}
	if ( newElem == NULL )
		Error( "Resize <memory allocation error>." );
	
	long minRow  = ( fRow < newRows ? fRow : newRows );
	long minCol  = ( fCol < newCols ? fCol : newCols );
	REAL * src;
	REAL * dst;
	long r, c;
	
	for ( c = 0; c < minCol; c++ ) {
		// Copia as linhas da matriz antiga para a nova.
		// Copia os elementos de uma linha.
		dst = newElem + c*newRows;
		src = fElem + c*fRow;
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
	fRow  = newRows;
	fCol  = newCols;
	return( 1 );
}

int TPZFMatrix::Remodel(const int newRows,const int newCols) {
	if(newRows*newCols != fRow*fCol) return -1;
	fRow = newRows;
	fCol = newCols;
	return 1;
}


/*************/
/*** Redim ***/

/*int TPZFMatrix::Redim(const int newRows,const int newCols) {
 if ( newRows == Rows() && newCols == Cols() ) {
 Zero();
 return( 1 );
 }
 int newsize = newRows*newCols;
 #ifdef __BOORLANDC__
 if(fElem && fElem != fGiven)  delete []fElem;
 // farfree(fElem);//
 
 #else
 if(fElem && fElem != fGiven) delete []fElem;
 //free(fElem);//
 
 #endif
 if(fGiven && newsize <= fSize) {
 fElem = fGiven;
 } else if(newsize == 0) {
 fElem = NULL;
 } else {
 #ifdef __BOORLANDC__
 //fElem = (REAL *) farcalloc(newsize,sizeof(REAL));//
 fElem = new( REAL[ newsize ] );
 #else
 //fElem = (REAL *) calloc(newsize,sizeof(REAL));//
 fElem = new( REAL[ newsize ] );
 #endif
 }
 if (newsize && fElem == NULL )
 Error( "Resize <memory allocation error>." );
 
 fRow  = newRows;
 fCol  = newCols;
 
 Zero();
 
 return( 1 );
 }
 */
/***************/
/****Zero*******/

/*int TPZFMatrix::Zero() {
 int size = fRow * fCol * sizeof(REAL);
 memset(fElem,'\0',size);
 //	 REAL * p = fElem, plast = p+size;
 //	 while(p < plast) *p++ = 0.0;
 fDecomposed = 0;
 return( 1 );
 }*/



/********************/
/*** Transpose () ***/


void TPZFMatrix::Transpose(TPZMatrix *const T) const{
	T->Resize( Cols(), Rows() );
	//Transposta por filas
	REAL * p = fElem;
	for ( int c = 0; c < Cols(); c++ ) {
		for ( int r = 0; r < Rows(); r++ ) {
			T->PutVal( c, r, *p++ );
			//            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		}
	}
}

void TPZFMatrix::Transpose() {
	TPZFMatrix temp;
	Transpose(&temp);
	*this = temp;
}

int TPZFMatrix::Decompose_LU(TPZVec<int> &index) {
	
	if (fDecomposed) return 0;
	
	if ( this->Rows() != this->Cols() ) {
		cout << "TPZFPivotMatrix::DecomposeLU ERRO : A Matriz não é quadrada" << endl;
		return 0;
	}
	
	int i,j,k;
	REAL sum = 0.;
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
				sum += Get(i,k) * Get(k,j);
			}
			REAL aux = Get(i,j);
			PutVal(i,j,(aux - sum));
			//cout << "0_A[" << i << "," << j << "]= " << Get(i,j) << endl;
		}
		//Print(cout);
		REAL piv = Get(j,j);
		//  cout << "Pivo 1 =" << piv << endl;
		int row = j;
		for (i=j+1;i<nRows;i++){
			sum = 0.;
			for (k=0;k<(j);k++){
				sum += Get(i,k) * Get(k,j);
			}
			REAL aux = Get(i,j);
			PutVal(i,j,(aux - sum));
			//cout << "1_A[" << i << "," << j << "]= " << Get(i,j) << endl;
			
			if (fabs(Get(i,j)) > fabs(piv)){
				piv = Get(i,j);
				//  cout << "Pivo 2 =" << piv << endl;
				row = i;
			}
		}
		//    Print(cout);
		if (row > j){
			for (k=0;k<nCols;k++){
				REAL aux = Get(j,k);
				PutVal(j,k,Get(row,k));
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
			REAL aux = Get(i,j) / piv;
			PutVal(i,j,aux);
			//cout << "4_A[" << i << "," << j << "]= " << Get(i,j) << endl;
		}
		//Print(cout);
	}    
	fDecomposed = ELUPivot;
	return 1;
}


/*****************/
/*** DecomposeLU ***/
int TPZFMatrix::Decompose_LU(std::list<int> &singular) {
	return Decompose_LU();
}


int TPZFMatrix::Decompose_LU() {
	
	if (  fDecomposed && fDecomposed != ELU)  Error( "Decompose_LU <TPZFMatrix already Decomposed with other scheme>" );
	if (fDecomposed) return 1;
	REAL nn;
	REAL * ptrpivot,*pik, *pij,*pkj;
	
	int i,j,k,rows=Rows(),cols=Cols();
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
			for ( j = k+1; j < Cols(); j++ )
			{
				pkj+=cols;pij+=cols;
				(*pij)-=nn*(*pkj);
			}
		}
		ptrpivot+=Rows()+1;
	}
	
	fDecomposed=1;
	return 1; 
}

int TPZFMatrix::Substitution(const REAL *ptr, int rows, TPZFMatrix *B)
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

int TPZFMatrix::Substitution( TPZFMatrix *B ) const {
	
#ifndef DEBUG2
	if(fDecomposed != ELU) {
		Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
	}
	int rowb = B->Rows();
	int colb = B->Cols();
	int row = Rows();
	if ( rowb != Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
	
	
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
	
	if(fDecomposed != ELU) {
		Error("TPZFMatrix::Decompose_LU substitution called for a wrongly decomposed matrix");
	}
	int rowb = B->Rows();
	int colb = B->Cols();
	if ( rowb != Rows() ) Error( "SubstitutionLU <incompatible dimensions>" );
	
	
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


int TPZFMatrix::Substitution( TPZFMatrix *B, TPZVec<int> &index ) const{
	
	if(!B){
		PZError << __PRETTY_FUNCTION__ << "TPZFMatrix *B eh nulo" << endl;
		return 0;
	}
	
	TPZFMatrix &b = *B;
	
	if (!fDecomposed){
		PZError <<  __PRETTY_FUNCTION__ << "Matriz não decomposta" << endl;
		return 0;
	}
	
	if (fDecomposed != ELUPivot){
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
	REAL sum = 0;
	
	TPZVec<REAL> v(nRows);
	
	
	for (i=0;i<nRows;i++)
	{
		v[i] = b(index[i]);
	}
	
	//Ly=b
	for (i=0;i<nRows;i++)
	{
		sum = 0.;
		for (j=0;j<(i);j++) sum += Get(i,j) * v[j];
		v[i] -= sum;
	}
	
	//Ux=y
	for (i=(nRows-1);i>-1;i--)
	{
		sum = 0.;
		for (j=(i+1);j<nRows;j++) sum += Get(i,j) * v[j];
		v[i] = (v[i] - sum) / Get(i,i);
	}
	
	for (i=0;i<nRows;i++) b(i) = v[i];
	return 1;
}


int TPZFMatrix::Substitution(const REAL *ptr, int rows, TPZFMatrix *B, TPZVec<int> &index )
{
	
	if(!B){
		PZError << __PRETTY_FUNCTION__ << "TPZFMatrix *B eh nulo" << endl;
		return 0;
	}
	
	TPZFMatrix &b = *B;
	
	
	
	if (index.NElements() != rows || b.Rows() != rows)
	{
		cout << "TMatrix::Substituicao ERRO : vetores com dimensões incompativeis:\n"
		<< "this->fIndex = " << index.NElements() << "  b = " << b.Rows() << endl;
		return 0;
	}
	
	int i,j;
	REAL sum = 0;
	
	TPZVec<REAL> v(rows);
	
	
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
REAL Dot(const TPZFMatrix &A, const TPZFMatrix &B) {
	int size = (A.Rows())*A.Cols();
	REAL result = 0.;
	if(!size) return result;
#ifdef USING_ATLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;
	
#elif USING_BLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;
	
#else
	const REAL *fpA = &A.g(0,0), *fpB = &B.g(0,0);
	const REAL *fpLast = fpA+size;
	while(fpA < fpLast) 
	{
		result += (*fpA++ * *fpB++); 
	}
	return result;
#endif
}

/** @brief Increments value over all entries of the matrix A. */
TPZFMatrix operator+(const REAL value, const TPZFMatrix &A ) {
	return( A + value );
}

/** @brief Decrements value over all entries of the matrix A. */
TPZFMatrix operator-(const REAL value, const TPZFMatrix &A ) {
	return( A - value );
}

/*TPZTempFMatrix operator*(const REAL value, const TPZFMatrix &A ) {
 return( A * value );
 }
 */


/************************** Private **************************/

/*************/
/*** Error ***/

int TPZFMatrix::Error(const char *msg1,const char *msg2 ) {
	ostringstream out;
	out << "TPZFMatrix::" << msg1;
	if(msg2) out << msg2;
	out << ".\n";
	LOGPZ_ERROR (logger, out.str().c_str());
	DebugStop();
	// int temp;//para testes
	// cin >> temp;//para testes
	//DebugStop();//para testes
	return 0;
}



/*************/
/*** Clear ***/

int TPZFMatrix::Clear() {
	if(fElem && fElem != fGiven) delete[]( fElem );
	fElem = NULL;
	fRow  = fCol = 0;
	return( 1 );
}

void TPZFMatrix::Read( TPZStream &buf, void *context ){
	TPZMatrix::Read(buf,context);
	int row = fRow;
	int col = fCol;
	fRow = fCol = 0;
	Resize(row,col);
	buf.Read(fElem,fRow*fCol);
}

void TPZFMatrix::Write( TPZStream &buf, int withclassid ) {
	TPZMatrix::Write(buf,withclassid);
	buf.Write(fElem,fRow*fCol);
}

// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZFMatrix::Compare(TPZSaveable *copy, bool override)
{
	TPZFMatrix *fmat = dynamic_cast<TPZFMatrix *> (copy);
	if(!fmat) return false;
	
	bool matresult = TPZMatrix::Compare(copy,false);
	int nel = fRow*fCol;
	REAL diff=0.;
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
		sout << " number different terms " << numdif << " number terms " << fRow*fCol;
		sout << " difference in norm L1 " << diff;
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(!matresult && override) 
	{
		this->operator=(*fmat);
	}
	return matresult;
}

// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZFMatrix::Compare(TPZSaveable *copy, bool override) const
{
	TPZFMatrix *fmat = dynamic_cast<TPZFMatrix *> (copy);
	if(!fmat) return false;
    
	bool matresult = TPZMatrix::Compare(copy,false);
	int nel = fRow*fCol;
	REAL diff=0.;
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
		sout << " number different terms " << numdif << " number terms " << fRow*fCol;
		sout << " difference in norm L1 " << diff;
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(!matresult && override) 
	{
		DebugStop();
	}
	return matresult;
}

void TPZFMatrix::PrintStatic(const REAL *ptr, int rows, int cols, const char *name, std::ostream& out,const MatrixOutputFormat form){
	
	//  out.width( 8 );
	//  out.precision( 4 );
	
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
				REAL val = SELECTEL(ptr,rows,row, col);
				if(val != 0.) out << row << ' ' << col << ' ' << val << std::endl;
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
				REAL val = SELECTEL(ptr,rows,row, col);
				sprintf(number, "%16.16lf", (double)val);
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

int TPZFMatrix::ClassId() const   
{ 
	return TPZFMATRIXID; 
}

template class TPZRestoreClass< TPZFMatrix, TPZFMATRIXID>;



int TPZFMatrix::SetSize(const int newRows,const int newCols) {
	long newsize = ((long)newRows)*newCols;
	long oldsize = fRow*fCol;
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
		fElem = new REAL[ newRows * newCols ] ;
	}
	if (newsize && fElem == NULL )
		Error( "Resize <memory allocation error>." );
	
	return( 1 );
}

