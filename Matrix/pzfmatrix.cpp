
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


#include <string.h>
#include <math.h>
#include <stdlib.h>

#ifdef __BOORLANDC__
#include <alloc.h>
#endif

#include "pzfmatrix.h"
#include "pztempmat.h"
#include "pzvec.h"
#include "pzerror.h"


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
     #include <cblas.h>
     };
double cblas_ddot(const int N, const double *X, const int incX,
                  const double *Y, const int incY);
void cblas_daxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);
#endif


#define IsZero( a )  ( (a) > -1.e-20 && (a) < 1.e-20 )




/********************/
/*** Constructors ***/




TPZFMatrix::TPZFMatrix(const TPZMatrix &mat) : TPZMatrix(mat), fElem(0),fGiven(0),fSize(0) {
  if(fRow*fCol) {
    fElem = new REAL[fRow*fCol];
    REALPtr p = fElem;
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
    fElem = new( REAL[ size ] );
#ifdef DEBUG
    if ( size && fElem == NULL ) Error( "Constructor <memory allocation error>." );
#endif
    // Copia a matriz
    REALPtr src = A.fElem;
    REALPtr p = fElem;
    memcpy(p,src,(size_t)size*sizeof(REAL));
  }







/******** Operacoes com matrizes FULL  ********/

/******************/
/*** Operator = ***/

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

/******************/
/*** Operator = ***/

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

TPZFMatrix &TPZFMatrix::operator=(const TPZFMatrix &A ) {
	if(this == &A) return *this;
	 long size = A.fRow * A.fCol;

	 REALPtr newElem = fElem;
	 if(fSize < size && size != fRow*fCol) {
		 newElem = new( REAL[size] );
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
	if(rhs.Cols() != Cols()) {
   	PZError << "TPZFMatrix::AddFel number of columns does not correspond\n";
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

TPZTempFMatrix TPZFMatrix::operator+(const TPZFMatrix &A ) const {
	 if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		  Error( "Operator+ <matrixs with different dimensions>" );

  TPZTempFMatrix res;
  res.Object().Redim( Rows(), Cols() );
  long size = ((long)Rows()) * Cols();
  REALPtr pm = fElem, plast = fElem+size;
  REALPtr pa = A.fElem;
  REALPtr pr = res.Object().fElem;

  while(pm < plast) *pr++ = (*pm++) + (*pa++);

  return( res );
}



/*******************************/
/*** Operator-( TPZFMatrix & ) ***/

TPZTempFMatrix TPZFMatrix::operator-(const TPZFMatrix &A ) const {
	 if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		  Error( "Operator- <matrixs with different dimensions>" );

	 TPZTempFMatrix res;
    res.Object().Redim( Rows(), Cols() );
	 long size = ((long)Rows()) * Cols();
	 REALPtr pm = fElem;
	 REALPtr pa = A.fElem;
	 REALPtr pr = res.Object().fElem, prlast =pr+size;

	 while(pr < prlast) *pr++ = (*pm++) - (*pa++);
	 return( res );
}

void TPZFMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
			 const REAL alpha,const REAL beta,const int opt,const int stride) const{
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
	REALPtr zp = &z(0,ic), zlast = zp+rows*stride;
	REALPtr fp = fElem +rows*c;
	const REAL * xp = &x.g(c*stride,ic);
	while(zp < zlast) {
	  *zp += alpha* *fp++ * *xp;
	  zp += stride;
	}
      }
    } else {
      REALPtr fp = fElem,  zp = &z(0,ic);
      for (c = 0; c<cols; c++) {
	REAL val = 0.;
	// bug correction philippe 5/2/97
	//					 REALPtr xp = &x(0,ic), xlast = xp + numeq*stride;
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


/*******************************/
/*** Operator*( TPZFMatrix & ) ***/

TPZTempFMatrix TPZFMatrix::operator*(const TPZFMatrix &A )const {
	 if ( Cols() != A.Rows() )
		  Error( "Operator* <matrixs with incompatible dimensions>" );
	 TPZTempFMatrix res;
    res.Object().Redim( Rows(), A.Cols() );
	 MultAdd(A,A,res.Object(),1.,0.,0);
	 return( res );
}



/********************************/
/*** Operator+=( TPZFMatrix & ) ***/

TPZFMatrix & TPZFMatrix::operator+=(const TPZFMatrix &A ) {
	 if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		  Error( "Operator+= <matrixs with different dimensions>" );

	 long size = ((long)Rows()) * Cols();
	 REALPtr pm = fElem, pmlast=pm+size;
	 REALPtr pa = A.fElem;
	 while(pm < pmlast) (*pm++) += (*pa++);
	 return( *this );
}



/*******************************/
/*** Operator-=( TPZFMatrix & ) ***/

TPZFMatrix &TPZFMatrix::operator-=(const TPZFMatrix &A ) {
	 if ( (A.Rows() != Rows())  ||  (A.Cols() != Cols()) )
		  Error( "Operator-= <matrixs with different dimensions>" );

	 long size = ((long)Rows()) * Cols();
	 REALPtr pm = fElem;
	 REALPtr pa = A.fElem;

	 for ( long i = 0; i < size; i++ ) *pm++ -= *pa++;

	 return( *this );
}

void TPZFMatrix::ZAXPY(const REAL alpha,const TPZFMatrix &p) {

#ifndef USING_ATLAS
#ifndef USING_BLAS
	REALPtr pt = fElem;
	REALPtr pp = p.fElem;
	REALPtr ptlast = fElem + fRow*fCol;
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
	REALPtr pt = fElem,  ptlast = fElem + fRow*fCol;
	REALPtr pz = z.fElem;
//	while(pt < ptlast) *pt++ *= (beta) + *pz++;
	while(pt < ptlast) {
		*pt *= (beta);
		*pt++ += *pz++;
	}

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
    fElem = new( REAL[ arows * acols ] );
  } else {
    fElem = fGiven;
  }
  REALPtr dst = fElem;
  for ( int c = 0; c < fCol; c++ )
    for ( int r = 0; r < fRow; r++ ) 
      *dst++ = A.Get( r, c );
  return( *this );
}



/******************/
/*** Operator + ***/

TPZTempFMatrix TPZFMatrix::operator+(const TPZMatrix &A ) const {
	 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
		  Error( "Operator+ (TPZMatrix &) <different dimensions>" );

	 TPZTempFMatrix res;
    res.Object().Redim( Rows(), Cols() );
	 REALPtr src = fElem;
	 REALPtr dst = res.Object().fElem;
	 long row = Rows();
	 long col = Cols();
	 for ( int c = 0; c < col; c++ )
		  for ( int r = 0; r < row; r++ )
				*dst++ = (*src++) + A.Get( r, c );

	 return( res );
}



/******************/
/*** Operator - ***/

TPZTempFMatrix TPZFMatrix::operator-(const TPZMatrix &A ) const {
	 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
		  Error( "Operator+ (TPZMatrix &) <different dimensions>" );

	 TPZTempFMatrix res(*this);
	 res.Object() -= A;
	 return( res );
}



/******************/
/*** Operator * ***/

TPZTempFMatrix TPZFMatrix::operator*(const TPZMatrix &A ) const {

	 if ( Cols() != A.Rows() )
		  Error( "Operator* (TPZMatrix &) <incompatible dimensions>" );


	 TPZTempFMatrix res;
    res.Object().Redim( Rows(), A.Cols() );

	 int     r,c,i,
				acols=A.Cols(),
				rows=Rows(),
				cols=Cols();

	 REALPtr pr = res.Object().fElem, /*Percorre os elementos da matriz 'res'*/
				pm,
				pm_aux;

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



/*******************/
/*** Operator += ***/

TPZFMatrix &TPZFMatrix::operator+=(const TPZMatrix &A ) {
	 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
		  Error( "Operator+ (TPZMatrix &) <different dimensions>" );

	 REALPtr pm = fElem;
	 long cols = Cols();
	 long rows = Rows();
	 for ( int c = 0; c < cols; c++ )
		  for ( int r = 0; r < rows; r++ )
				*pm++ += A.Get( r, c );

	 return( *this );
}



/*******************/
/*** Operator -= ***/

TPZFMatrix &TPZFMatrix::operator-=(const TPZMatrix &A ) {
	 if ( (Rows() != A.Rows()) || (Cols() != A.Cols()) )
		  Error( "Operator+ (TPZMatrix &) <different dimensions>" );

	 REALPtr pm = fElem;
	 int cols = Cols();
	 int rows = Rows();
	 for ( int c = 0; c < cols; c++ )
		  for ( int r = 0; r < rows; r++ )
				*pm++ -= A.Get( r, c );

	 return( *this );
}



/******** Operacoes com valores NUMERICOS ********/

/******************/
/*** Operator = ***/

TPZFMatrix& TPZFMatrix::operator=(const REAL value ) {
	 long size = ((long)fRow) * fCol;
	 REALPtr dst   = fElem;
	 for ( long i = 0; i < size; i++ )
		  *dst++ = value;
	 fDecomposed = 0;
	 return *this;
}



/***************************/
/*** Operator+=( value ) ***/

TPZFMatrix &TPZFMatrix::operator+=(const REAL value ) {
	 long size = ((long)Rows()) * Cols();

	 REALPtr dst = fElem, dstlast = dst+size;
	 while ( dst < dstlast ) *dst++ += value;
	 return( *this );
}



/**************************/
/*** Operator+( value ) ***/

TPZTempFMatrix TPZFMatrix::operator+(const REAL value ) const {
	 TPZTempFMatrix res( *this );
	 long size = ((long)Rows()) * Cols();

	 REALPtr dst = res.Object().fElem,  dstlast = dst+size;
	 while ( dst < dstlast )
		  *dst++ += value;

	 return( res );
}

TPZTempFMatrix TPZFMatrix::operator-  (const REAL val ) const {
	return operator+( -val ); 
}


/**************************/
/*** Operator*( value ) ***/

TPZTempFMatrix
TPZFMatrix::operator*(const REAL value ) const
{
	 TPZTempFMatrix res( *this );
	 res.Object() *= value;
	 return( res );
}



/***************************/
/*** Operator*=( value ) ***/

TPZFMatrix &TPZFMatrix::operator*=( const REAL value ) {
  long size = ((long)Rows()) * Cols();
  REALPtr dst = fElem, dstlast = dst+size;
  while ( dst < dstlast ) *dst++ *= value;
  return( *this );
}



/**************/
/*** Resize ***/

int TPZFMatrix::Resize(const int newRows,const int newCols) {
	 if ( newRows == Rows() && newCols == Cols() )
		  return( 1 );
	 long newsize = ((long)newRows)*newCols;
	 REALPtr newElem;
	 if(fGiven && fElem != fGiven && newsize <= fSize) {
		  newElem = fGiven;
	 } else {
		 // newElem = (REALPtr) calloc(newRows*newCols,sizeof(REAL));//
       newElem = new( REAL[ newRows * newCols ] );
	 }
	 if ( newElem == NULL )
		  Error( "Resize <memory allocation error>." );

	 long minRow  = ( fRow < newRows ? fRow : newRows );
	 long minCol  = ( fCol < newCols ? fCol : newCols );
	 REALPtr src;
	 REALPtr dst;
	 long r, c;

	 for ( c = 0; c < minCol; c++ ) {
		  // Copia as linhas da matriz antiga para a nova.
		  // Copia os elementos de uma linha.
		  dst = newElem + c*newRows;
		  src = fElem + c*fRow;
		  for ( r = 0; r < minRow; r++ )
				*dst++ = *src++;

		  // Se a nova linha for maior (mais colunas), preenche o
		  //  resto da linha com ZEROS.
		  for ( ; r < newRows; r++ )
				*dst++ = 0.0;
	 }

	 // Preenche as linha que sobrarem (se sobrarem) com ZEROS.
	 for ( ;c < newCols; c++ ) {
		  dst = newElem + c*newRows;
		  for (r = 0 ; r < newRows; r++ ) *dst++ = 0.0;
	 }

	 if (fElem && fElem != fGiven )delete[]( fElem );
	 fElem = newElem;
	 fRow  = newRows;
	 fCol  = newCols;
	 return( 1 );
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
		  //fElem = (REALPtr) farcalloc(newsize,sizeof(REAL));//
        fElem = new( REAL[ newsize ] );
#else
		  //fElem = (REALPtr) calloc(newsize,sizeof(REAL));//
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
//	 REALPtr p = fElem, plast = p+size;
//	 while(p < plast) *p++ = 0.0;
	 fDecomposed = 0;
	 return( 1 );
}*/



/********************/
/*** Transpose () ***/


void TPZFMatrix::Transpose(TPZMatrix *const T) const{
	  T->Resize( Cols(), Rows() );
//Transposta por filas
	 REALPtr p = fElem;
	 for ( int c = 0; c < Cols(); c++ ) {
		 for ( int r = 0; r < Rows(); r++ ) {
				T->PutVal( c, r, *p++ );
//            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		  }
	 }
}

void TPZFMatrix::Transpose() {
	  Resize( Cols(), Rows() );
//Transposta por filas
    REAL val;
	 for ( int c = 0; c < Cols(); c++ ) {
		 for ( int r = 0; r < Rows(); r++ ) {
       	val = GetVal(r,c);
         PutVal(r,c,GetVal(c,r));
         PutVal(c,r,val);
//            cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		  }
	 }
    int row = Rows();
	 fRow = fCol;
    fCol = row;
}


/*****************/
/*** DecomposeLU ***/

int TPZFMatrix::Decompose_LU() {

  if (  fDecomposed && fDecomposed != ELU)  Error( "Decompose_LU <TPZFMatrix already Decomposed with other scheme>" );
  if (fDecomposed) return 1;
  REAL nn;
  REALPtr ptrpivot,pik, pij,pkj;
  
  int i,j,k,rows=Rows(),cols=Cols();
  int  min = ( cols < (rows-1) ) ? cols : rows - 1;
  
  ptrpivot=&fElem[0];
  for (  k = 0; k < min ; k++ )
  {
    if ( IsZero( *ptrpivot ) ){
      //Power plus...
      if (fabs(*ptrpivot) > 0){
        for (j=k+1;j<rows;j++){
          if (fabs(*(ptrpivot + j - k) - *(ptrpivot)) > 1e-12)
            Error( "DecomposeLU <matrix is singular> even after Power Plus..." );
        }
      }
      else Error( "DecomposeLU <matrix is singular>" );
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


/****************/
/*** Substitution ***/

int TPZFMatrix::Substitution( TPZFMatrix *B ) const {

#ifndef DEBUG
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

REAL Dot(const TPZFMatrix &A,const TPZFMatrix &B) {
	int size = (A.Rows())*A.Cols();
	double result = 0.;
#ifdef USING_ATLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;

#elif USING_BLAS
	result = cblas_ddot(size, &A.g(0,0), 1, &B.g(0,0), 1);
	return result;

#else
	const REAL *fpA = &A.g(0,0), *fpB = &B.g(0,0);
	const REAL *fpLast = fpA+size;
	while(fpA < fpLast) result += *fpA++ * *fpB++;
	return result;
#endif
}

TPZTempFMatrix operator+(const REAL value, const TPZFMatrix &A ) {
	 return( A + value );
}



TPZTempFMatrix operator-(const REAL value, const TPZFMatrix &A ) {
	 return( A - value );
}




TPZTempFMatrix operator*(const REAL value, const TPZFMatrix &A ) {
	 return( A * value );
}



/************************** Private **************************/

/*************/
/*** Error ***/

int TPZFMatrix::Error(const char *msg1,const char *msg2 ) const {
  cout << "TPZFMatrix::" << msg1;
  if(msg2) cout << msg2;
  cout << ".\n";
  int temp;//para testes
  cin >> temp;//para testes
  //exit( 1 );//para testes
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

