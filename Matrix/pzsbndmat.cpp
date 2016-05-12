/**
 * @file
 * @brief Contains the implementation of the TPZSBMatrix methods.
 */

#include <math.h>
#include <stdlib.h>

#include "pzfmatrix.h"
#include "pzsbndmat.h"

#ifdef USING_LAPACK
/** CBlas Math Library */
#ifdef MACOSX
#include <Accelerate/Accelerate.h>
#else
#include "cblas.h"
#define BLAS_MULT
#endif
#endif

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
TPZSBMatrix<TVar>::TPZSBMatrix( long dim, long band )
: TPZMatrix<TVar>( dim, dim )
{
    fBand = ( band > (dim - 1) ? (dim - 1) : band );
    fDiag.resize(Size());
    
    Zero();
}

/** Fill the matrix with random values (non singular matrix) */
template <class TVar>
void TPZSBMatrix<TVar>::AutoFill(long nrow, long ncol, int symmetric) {
    
    if (nrow != ncol || symmetric == 0) {
        DebugStop();
    }
    fBand = nrow/10;
    if (fBand == 0) {
        fBand = nrow-1;
    }
    Resize(nrow, ncol);
    
    long i, j;
    TVar val, sum;
    /** Fill data */
    for(i=0;i<this->Rows();i++) {
        sum = 0.0;
        long jmax = i+fBand+1;
        if (jmax >= this->Rows()) {
            jmax = this->Rows();
        }
        for (j=0; j<i; j++) {
            sum += fabs(GetVal(i, j));
        }
        for(j=i;j<jmax;j++) {
            val = ((TVar)rand())/((TVar)RAND_MAX);
            if(!PutVal(i,j,val))
            {
                std::cout << "AutoFill (TPZMatrix) failed.";
                DebugStop();
            }
            if(i!=j) sum += fabs(val);
        }
        if (this->Rows() == this->Cols()) {
            /** Making diagonally dominant and non zero in diagonal */
            if(fabs(sum) > fabs(GetVal(i,i)))            // Deve satisfazer:  |Aii| > SUM( |Aij| )  sobre j != i
                PutVal(i,i,sum+(TVar)1.);
            // To sure diagonal is not zero.
            if(IsZero(sum) && IsZero(GetVal(i,i)))
                PutVal(i,i,1.);
        }
    }
    
}
/**************/
/*** PutVal ***/

template <class TVar>
int
TPZSBMatrix<TVar>::PutVal(const long r,const long c,const TVar& value )
{
    // inicializando row e col para trabalhar com a triangular superior
    long row(r),col(c);
    if ( row > col )
        this->Swap( &row, &col );
    
    long index;
    if ( (index = col-row) > fBand )
    {
#ifdef PZDEBUG
        if (value != this->gZero) {
            DebugStop();
        }
#endif
        return( 0 );        // O elemento esta fora da banda.
    }
    fDiag[ Index(row,col) ] = value;
    this->fDecomposed = 0;
    return( 1 );
}

/**************/
/*** GetVal ***/
template<class TVar>
const TVar
&TPZSBMatrix<TVar>::GetVal(const long r,const long c ) const
{
    
    // inicializando row e col para trabalhar com a triangular superior
    long row(r),col(c);
    if ( row > col )
        this->Swap( &row, &col );
    
    long index;
    if ( (index = col-row) > fBand )
        return( this->gZero );        // O elemento esta fora da banda.
    
    return( fDiag[ Index(row,col) ] );
}

template<class TVar>
TVar &TPZSBMatrix<TVar>::operator()(long row, long col)
{
    // inicializando row e col para trabalhar com a triangular superior
    if ( row > col )
        this->Swap( &row, &col );
    
    long index;
    if ( (index = col-row) > fBand )
        return( this->gZero );        // O elemento esta fora da banda.
    
    return( fDiag[ Index(row,col) ] );

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
    TPZMatrix<TVar>::Print(name,out,form);
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
    
    for ( long row = 0; row < A.Rows(); row++)
    {
        out << "\t";
        for ( long col = 0; col < A.Cols(); col++ )
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
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
    
    // Define os ponteiros e tamanhos para os elementos da maior e da
    //  menor banda.
    TPZSBMatrix Res(*this);
    Res += A;
    return Res;
}

/******************/
/*** Operator - ***/

template<class TVar>
TPZSBMatrix<TVar>
TPZSBMatrix<TVar>::operator-(const TPZSBMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
    {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-( TPZSBMatrix ) <incompatible dimensions>" );
    }
    
    TPZSBMatrix<TVar> res(*this);
    res -= A;
    
    return( res );
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator+=(const TPZSBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
    {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator+=( TPZSBMatrix ) <incompatible dimensions>" );
    }
    long siz= fDiag.size();
    for (long i=0; i<siz; i++) {
        fDiag[i] += A.fDiag[i];
    }
    return( *this );
}

/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator-=(const TPZSBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
    {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "operator-=( TPZSBMatrix ) <incompatible dimensions>" );
    }
    long siz= fDiag.size();
    for (long i=0; i<siz; i++) {
        fDiag[i] -= A.fDiag[i];
    }
    
    return( *this );
}

template<class TVar>
void TPZSBMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                const TVar alpha,const TVar beta ,const int opt) const {
    // Computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot overlap in memory
    if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "TPZSBMatrix::MultAdd <matrixs with incompatible dimensions>" );
    if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"TPZSBMatrix::MultAdd incompatible dimensions\n");
    }
    this->PrepareZ(y,z,beta,opt);
    long rows = this->Rows();
    long xcols = x.Cols();
    long ic, r;
    for (ic = 0; ic < xcols; ic++) {
        long begin, end;
        for ( r = 0; r < rows; r++ ) {
            begin = MAX( r - fBand, 0 );
            end   = MIN( r + fBand + 1, rows );
            TVar val = 0.;
            // Calcula um elemento da resposta.
            for ( long i = begin ; i < end; i++ ) val += GetVal( r, i ) * x.GetVal(i, ic );
            val *= alpha;
            val += z.GetVal(r,ic);
            z.PutVal( r , ic, val );
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
    TPZSBMatrix<TVar> res( *this );
    
    long siz= fDiag.size();
    for (long i=0; i<siz; i++) {
        res.fDiag[i] *= value;
    }
    
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
    long siz= fDiag.size();
    for (long i=0; i<siz; i++) {
        fDiag[i] *= value;
    }
    
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
TPZSBMatrix<TVar>::Resize(const long newDim ,const long)
{
    if ( newDim == this->Dim() )
        return( 1 );
    
    Redim(newDim,newDim);
    return( 1 );
}

/*************/
/*** Redim ***/
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//
template<class TVar>
int
TPZSBMatrix<TVar>::Redim(const long newDim ,const long otherDim)
{
    if (newDim != otherDim) {
        DebugStop();
    }
    
    if ( newDim != this->Dim() )
    {
        TPZMatrix<TVar>::Redim(newDim,newDim);
        if (fBand > newDim-1) {
            fBand = newDim-1;
        }
        fDiag.resize(Size());
    }
    
    Zero();
    this->fDecomposed = 0;
    this->fDefPositive = 0;
    return( 1 );
}

template<class TVar>
int
TPZSBMatrix<TVar>::Zero()
{
    long siz= fDiag.size();
    for (long i=0; i<siz; i++) {
        fDiag[i] = (TVar)0.;
    }
    
    
    this->fDecomposed = 0;
    this->fDefPositive = 0;
    return( 1 );
}


/****************/
/*** Set Band ***/
template<class TVar>
int
TPZSBMatrix<TVar>::SetBand(long newBand )
{
    if ( this->fBand == newBand )
        return( 1 );
    
    if ( newBand > (this->Dim() - 1) )
    {
        newBand = this->Dim()-1;
    }
    fBand = newBand;
    fDiag.resize(Size());
    Zero();
    
    return( 1 );
}

/********************* Resolucao de sistemas *********************/


#ifdef USING_LAPACK
/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky(std::list<long> &singular)
{
    return Decompose_Cholesky();
}

template<>
int TPZSBMatrix<float>::Decompose_Cholesky()
{
    if (fDecomposed == ECholesky) {
        return 1;
    }
    if (fDecomposed != ENoDecompose) {
        DebugStop();
    }
    char uplo[]="Upper";
    int n = Dim();
    int kd = fBand;
    int nrhs = 0;
    float *ab = &fDiag[0];
    int ldab = fBand+1;
    float b = 0;
    int info;
    
    //    spbsv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__kd#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    spbsv_(uplo, &n, &kd, &nrhs, ab, &ldab, &b, &n, &info);
    
    if (info != 0) {
        DebugStop();
    }
    fDecomposed = ECholesky;
    return 1;
}

template<>
int TPZSBMatrix<double>::Decompose_Cholesky()
{
    if (fDecomposed == ECholesky) {
        return 1;
    }
    if (fDecomposed != ENoDecompose) {
        DebugStop();
    }
    char uplo[]="Upper";
    int n = Dim();
    int kd = fBand;
    int nrhs = 0;
    double *ab = &fDiag[0];
    int ldab = fBand+1;
    double b = 0;
    int info;
    
    //    spbsv_(<#char *__uplo#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__kd#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    dpbsv_(uplo, &n, &kd, &nrhs, ab, &ldab, &b, &n, &info);
    
    if (info != 0) {
        DebugStop();
    }
    fDecomposed = ECholesky;
    return 1;
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky()
{
    return TPZMatrix<TVar>::Decompose_Cholesky();
}
#endif

/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_LDLt(std::list<long> &singular)
{
    return Decompose_LDLt();
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_LDLt()
{
    
    if (  this->fDecomposed )  TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
    
    long j,k,l, begin,end;
    TVar sum;
    
    for ( j = 0; j < this->Dim(); j++ )
    {
        //Print("curernt");
        sum=0.;
        
        begin = MAX( long(j - fBand), 0 );
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
            end   = MIN( long(k + fBand )+1, this->Dim() );
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
        end  = MIN( long(j + fBand )+1, this->Dim() );
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
#ifdef USING_LAPACK
template<>
int
TPZSBMatrix<float>::Subst_Forward( TPZFMatrix<float>*B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    long bcols = B->Cols();
    for(long ic=0; ic<bcols; ic++)
    {
        //    cblas_stbsv(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_UPLO __Uplo#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_DIAG __Diag#>, <#const int __N#>, <#const int __K#>, <#const float *__A#>, <#const int __lda#>, <#float *__X#>, <#const int __incX#>)
        float *bptr = &(*B)(0,ic);
        cblas_stbsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<double>::Subst_Forward( TPZFMatrix<double>*B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    long bcols = B->Cols();
    for(long ic=0; ic<bcols; ic++)
    {
        //    cblas_stbsv(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_UPLO __Uplo#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_DIAG __Diag#>, <#const int __N#>, <#const int __K#>, <#const float *__A#>, <#const int __lda#>, <#float *__X#>, <#const int __incX#>)
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<class TVar>
int
TPZSBMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
    {
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    return TPZMatrix<TVar>::Subst_Forward(B);
}

template<>
int
TPZSBMatrix<float>::Subst_Backward( TPZFMatrix<float>*B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    long bcols = B->Cols();
    for(long ic=0; ic<bcols; ic++)
    {
        //    cblas_stbsv(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_UPLO __Uplo#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_DIAG __Diag#>, <#const int __N#>, <#const int __K#>, <#const float *__A#>, <#const int __lda#>, <#float *__X#>, <#const int __incX#>)
        float *bptr = &(*B)(0,ic);
        cblas_stbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<>
int
TPZSBMatrix<double>::Subst_Backward( TPZFMatrix<double>*B ) const
{
    if ( (B->Rows() != this->Dim()) || this->fDecomposed != ECholesky)
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    long bcols = B->Cols();
    for(long ic=0; ic<bcols; ic++)
    {
        //    cblas_stbsv(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_UPLO __Uplo#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_DIAG __Diag#>, <#const int __N#>, <#const int __K#>, <#const float *__A#>, <#const int __lda#>, <#float *__X#>, <#const int __incX#>)
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<class TVar>
int TPZSBMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    
    return TPZMatrix<TVar>::Subst_Backward(B);
    return ( 1 ) ;
    
}

#endif
/***********************/
/*** Subst L Forward ***/
//
//  Faz a "Forward substitution" assumindo que os elementos
//   da diagonal sao todos iguais a 1.
//
template<class TVar>
int TPZSBMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
    
    long size = fBand + 1;
    long i,j,k;
    for ( k = 0; k < this->Dim(); k++ )
    {
        for ( j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            
            long imin = k-fBand;
            if(imin < 0) imin = 0;
            
            long end=(k-fBand>0)? fBand:k;  //misael
            TVar sum = 0.0;
            for ( i = imin; i < k ; i++ )//misael
            {
                sum += GetVal(k,i) * B->GetVal(i,j);
            }
            // Faz b[k] = (b[k] - sum).
            //
            B->PutVal( k, j, (B->GetVal( k, j ) - sum) );
            
        }
    }
    return( 1 );
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template<class TVar>
int TPZSBMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const
{
    
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
    
    
    long size = fBand + 1;
    for ( long k = 0; k < this->Dim(); k++ )
        for ( long j = 0; j < B->Cols(); j++ )
            B->PutVal( k, j, B->GetVal( k, j) / GetVal(k,k) );
    
    return( 1 );
}

template<class TVar>
int TPZSBMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
        TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
    
    long k,j,i,jmax,stepcol=fBand+2;
    
    for(k=0; k<B->Cols() ; k++)
    {
        for(long i=this->Rows()-1; i>=0; i--)
        {
            jmax=( (i+fBand+1)>this->Rows())? this->Rows() : i+fBand+1;
            TVar sum = 0.;
            for(j=i+1;j<jmax;j++)
            {
                TVar el = GetVal(i,j);
                sum += B->GetVal(j,k)*el;
            }
            B->operator()(i,k) -= sum;
        }
        
    }
    
    return 1;
    
}

/**************************** PRIVATE ****************************/

/*************/
/*** CLear ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Clear()
{
    this->fRow = this->fCol = 0;
    fDiag.resize(0);
    this->fDecomposed = 0;
    return( 1 );
}

/************/
/*** Copy ***/

template<class TVar>
void
TPZSBMatrix<TVar>::Copy(const TPZSBMatrix<TVar> &A )
{
    TPZMatrix<TVar>::operator=(A);
    this->fBand = A.fBand;
    this->fDiag = A.fDiag;
}

#ifdef USING_LAPACK
/*** @name Solve eigenvalues ***/
/** @{ */
template< class TVar>
int
TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenVectors)
{
    TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}

template< class TVar>
int
TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < complex<double> > &w)
{
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SolveEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}

template<>
int
TPZSBMatrix<double>::SolveEigenProblem(TPZVec < complex<double> > &eigenvalues)
{
    
#ifdef USING_LAPACK
    char jobz = 'n'; //compute eigenvalues only
    char uplo = 'u';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec<double> w(0,0.);
    w.Resize( this->Dim() );
    TPZVec <double> z( this->Dim() *this->Dim() );
    int ldz = this->Dim();
    TPZVec <double> work( 3 *this->Dim() );
    int info = -666;
    
    dsbev_(&jobz, &uplo, &n, &kd, fDiag.begin(), &ldab, w.begin(), z.begin(), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
#endif
    eigenvalues.Resize( this->Dim() );
    for (int i = 0 ; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    return( 1 );
}

template<>
int
TPZSBMatrix<complex <double> >::SolveEigenProblem(TPZVec <complex<double> > &eigenvalues)
{
#ifdef USING_LAPACK
    
    char jobz = 'n'; //compute eigenvalues only
    char uplo = 'u';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec<double> w(this->Dim() , 0.);
    w.Resize( this->Dim() );
    TPZVec <complex <double> > z( this->Dim() *this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <double> > work( this->Dim() );
    TPZVec < double > rwork( 3 *this->Dim() );
    int info = -666;
    
    zhbev_(&jobz, &uplo, &n, &kd, (__CLPK_doublecomplex *)fDiag.begin(), &ldab,  w.begin(), (__CLPK_doublecomplex *)z.begin(), &ldz, (__CLPK_doublecomplex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    
#endif
    eigenvalues.Resize( this->Dim() );
    for (int i = 0 ; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
    return( 1 );
}

template<>
int
TPZSBMatrix<double>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenVectors)
{
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec < double > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix<double> z( this->Dim() ,this->Dim() );
    int ldz = this->Dim();
    TPZVec <double> work( 3 *this->Dim() );
    int info = -666;
    
    dsbev_(&jobz, &uplo, &n, &kd, fDiag.begin(), &ldab, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for(int i = 0 ; i < this->Dim() ; i++){
        eigenvalues[i] = w[i];
    }
    eigenVectors.Redim(this->Dim(), this->Dim());
    for (int iVec = 0 ; iVec < this->Dim(); iVec++) {
        for (int iCol = 0; iCol < this->Dim(); iCol++) {
            eigenVectors( iVec , iCol) = z(iVec, iCol);
        }
    }
    
#endif
    
    return( 1 );
}

template<>
int
TPZSBMatrix<complex <double> >::SolveEigenProblem(TPZVec <complex<double> > &eigenvalues, TPZFMatrix <complex <double> > &eigenVectors)
{
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec < double > w;
    w.Resize( this->Dim() );
    TPZFMatrix <complex <double> > z( this->Dim(),this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <double> > work( this->Dim() );
    TPZVec < double > rwork( 3 *this->Dim() );
    int info = -666;
    
    zhbev_(&jobz, &uplo, &n, &kd, (__CLPK_doublecomplex *)fDiag.begin(), &ldab, w.begin(), (__CLPK_doublecomplex *)&z(0,0), &ldz, (__CLPK_doublecomplex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim());
    for(int i = 0; i < this->Dim(); i++){
        eigenvalues[i] = w[i];
    }
    eigenVectors.Redim(this->Dim(), this->Dim());
    for (int iVec = 0 ; iVec < this->Dim(); iVec++) {
        for (int iCol = 0; iCol < this->Dim(); iCol++) {
            eigenVectors( iVec , iCol) = z(iVec,iCol);
        }
    }
    
#endif
    
    return( 1 );
}

template<>
int
TPZSBMatrix<float>::SolveEigenProblem(TPZVec < complex<double> > &eigenvalues)
{
    
#ifdef USING_LAPACK
    char jobz = 'n'; //compute eigenvalues only
    char uplo = 'u';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec<float> w(0,0.);
    w.Resize( this->Dim() );
    TPZVec <float> z( this->Dim() *this->Dim() );
    int ldz = this->Dim();
    TPZVec <float> work( 3 *this->Dim() );
    int info = -666;
    
    ssbev_(&jobz, &uplo, &n, &kd, fDiag.begin(), &ldab, w.begin(), z.begin(), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
#endif
    eigenvalues.Resize( this->Dim() );
    for (int i = 0 ; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    return( 1 );
}

template<>
int
TPZSBMatrix<complex <float> >::SolveEigenProblem(TPZVec <complex<double> > &eigenvalues)
{
#ifdef USING_LAPACK
    
    char jobz = 'n'; //compute eigenvalues only
    char uplo = 'u';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec<float> w(this->Dim() , 0.);
    w.Resize( this->Dim() );
    TPZVec <complex <float> > z( this->Dim() *this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <float> > work( this->Dim() );
    TPZVec < float > rwork( 3 *this->Dim() );
    int info = -666;
    
    chbev_(&jobz, &uplo, &n, &kd, (__CLPK_complex *)fDiag.begin(), &ldab,  w.begin(), (__CLPK_complex *)z.begin(), &ldz, (__CLPK_complex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    
#endif
    eigenvalues.Resize( this->Dim() );
    for (int i = 0 ; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
    return( 1 );
}

template<>
int
TPZSBMatrix<float>::SolveEigenProblem(TPZVec < std::complex<double> > &eigenvalues, TPZFMatrix < std::complex<double> > &eigenVectors)
{
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec < float > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix <float> z( this->Dim(),this->Dim() );
    int ldz = this->Dim();
    TPZVec <float> work( 3 *this->Dim() );
    int info = -666;
    
    ssbev_(&jobz, &uplo, &n, &kd, fDiag.begin(), &ldab, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for(int i = 0 ; i < this->Dim() ; i++){
        eigenvalues[i] = w[i];
    }
    eigenVectors.Redim(this->Dim(), this->Dim());
    for (int iVec = 0 ; iVec < this->Dim(); iVec++) {
        for (int iCol = 0; iCol < this->Dim(); iCol++) {
            eigenVectors( iVec , iCol) = z(iVec,iCol);
        }
    }
    
#endif
    
    return( 1 );
}

template<>
int
TPZSBMatrix<complex <float> >::SolveEigenProblem(TPZVec <complex<double> > &eigenvalues, TPZFMatrix <complex <double> > &eigenVectors)
{
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int kd = this->fBand;
    int ldab = this->fBand + 1;
    TPZVec < float > w;
    w.Resize( this->Dim() );
    TPZFMatrix <complex <float> > z( this->Dim(), this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <float> > work( this->Dim() );
    TPZVec < float > rwork( 3 *this->Dim() );
    int info = -666;
    
    chbev_(&jobz, &uplo, &n, &kd, (__CLPK_complex *)fDiag.begin(), &ldab, w.begin(), (__CLPK_complex *)&z(0,0), &ldz, (__CLPK_complex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"SolveEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim());
    for(int i = 0; i < this->Dim(); i++){
        eigenvalues[i] = w[i];
    }
    eigenVectors.Redim(this->Dim(), this->Dim());
    for (int iVec = 0 ; iVec < this->Dim(); iVec++) {
        for (int iCol = 0; iCol < this->Dim(); iCol++) {
            eigenVectors( iVec , iCol) = z(iVec,iCol);
        }
    }
    
#endif
    
    return( 1 );
}

template< class TVar>
int
TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &B , TPZVec < complex<double> > &w, TPZFMatrix < complex<double> > &eigenVectors)
{
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}
template< class TVar>
int
TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &B , TPZVec < complex<double> > &w)
{
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <LAPACK does not support this specific data type>" );
    return( 0 );
}
template<>
int
TPZSBMatrix<float>::SolveGeneralisedEigenProblem(TPZSBMatrix<float> &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenVectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec< float > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix<float> z( this->Dim(), this->Dim() );
    int ldz = this->Dim();
    TPZVec <float> work( 3 *this->Dim() );
    int info = -666;
    
    ssbgv_(&jobz, &uplo, &n, &ka, &kb, fDiag.begin(), &ldab, B.fDiag.begin(), &ldbb, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    eigenVectors.Resize( this->Dim() , this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
        for (int j = 0 ; j < this->Dim() ; j++) {
            eigenVectors(i,j) = z(i,j);
        }
    }
    
#endif
    
    return( 1 );
}


template<>
int
TPZSBMatrix<float>::SolveGeneralisedEigenProblem(TPZSBMatrix<float> &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'N'; //do NOT compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec< float > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix<float> z( this->Dim(), this->Dim() );
    int ldz = this->Dim();
    TPZVec <float> work( 3 *this->Dim() );
    int info = -666;
    
    ssbgv_(&jobz, &uplo, &n, &ka, &kb, fDiag.begin(), &ldab, B.fDiag.begin(), &ldbb, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
#endif
    
    return( 1 );
}

template<>
int
TPZSBMatrix<complex<float> >::SolveGeneralisedEigenProblem(TPZSBMatrix<complex<float> > &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenVectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec<float> w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix <complex <float> > z( this->Dim() ,this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <float> > work( this->Dim() );
    TPZVec < float > rwork( 3 *this->Dim() );
    int info = -666;
    
    chbgv_(&jobz, &uplo, &n, &ka, &kb, (__CLPK_complex *)fDiag.begin(), &ldab,  (__CLPK_complex *)B.fDiag.begin(), &ldbb, w.begin(), (__CLPK_complex *)&z(0,0), &ldz, (__CLPK_complex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    eigenVectors.Resize( this->Dim() , this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
        for (int j = 0 ; j < this->Dim() ; j++) {
            eigenVectors(i,j) = z(i,j);
        }
    }
    
#endif
    
    return( 1 );
}


template<>
int
TPZSBMatrix<complex<float> >::SolveGeneralisedEigenProblem(TPZSBMatrix<complex<float> > &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<float>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'N'; //do NOT compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec<float> w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix <complex <float> > z( this->Dim() ,this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <float> > work( this->Dim() );
    TPZVec < float > rwork( 3 *this->Dim() );
    int info = -666;
    
    chbgv_(&jobz, &uplo, &n, &ka, &kb, (__CLPK_complex *)fDiag.begin(), &ldab,  (__CLPK_complex *)B.fDiag.begin(), &ldbb, w.begin(), (__CLPK_complex *)&z(0,0), &ldz, (__CLPK_complex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<float> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
#endif
    
    return( 1 );
}


template<>
int
TPZSBMatrix<double>::SolveGeneralisedEigenProblem(TPZSBMatrix<double> &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenVectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec< double > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix<double> z( this->Dim(), this->Dim() );
    int ldz = this->Dim();
    TPZVec <double> work( 3 *this->Dim() );
    int info = -666;
    
    dsbgv_(&jobz, &uplo, &n, &ka, &kb, fDiag.begin(), &ldab, B.fDiag.begin(), &ldbb, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    eigenVectors.Resize( this->Dim() , this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
        for (int j = 0 ; j < this->Dim() ; j++) {
            eigenVectors(i,j) = z(i,j);
        }
    }
    
#endif
    
    return( 1 );
}


template<>
int
TPZSBMatrix<double>::SolveGeneralisedEigenProblem(TPZSBMatrix<double> &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'N'; //do NOT compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec< double > w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix<double> z( this->Dim(), this->Dim() );
    int ldz = this->Dim();
    TPZVec <double> work( 3 *this->Dim() );
    int info = -666;
    
    dsbgv_(&jobz, &uplo, &n, &ka, &kb, fDiag.begin(), &ldab, B.fDiag.begin(), &ldbb, w.begin(), &z(0,0), &ldz, work.begin(), &info);
    if( info > 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__,"SolveGeneralisedEigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
#endif
    
    return( 1 );
}
template<>
int
TPZSBMatrix<complex<double> >::SolveGeneralisedEigenProblem(TPZSBMatrix<complex<double> > &B , TPZVec <complex<double> > &eigenvalues, TPZFMatrix < complex<double> > &eigenVectors)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'V'; //compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec<double> w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix <complex <double> > z( this->Dim() ,this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <double> > work( this->Dim() );
    TPZVec < double > rwork( 3 *this->Dim() );
    int info = -666;
    
    zhbgv_(&jobz, &uplo, &n, &ka, &kb, (__CLPK_doublecomplex *)fDiag.begin(), &ldab,  (__CLPK_doublecomplex *)B.fDiag.begin(), &ldbb, w.begin(), (__CLPK_doublecomplex *)&z(0,0), &ldz, (__CLPK_doublecomplex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    eigenVectors.Resize( this->Dim() , this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
        for (int j = 0 ; j < this->Dim() ; j++) {
            eigenVectors(i,j) = z(i,j);
        }
    }
    
#endif
    
    return( 1 );
}


template<>
int
TPZSBMatrix<complex<double> >::SolveGeneralisedEigenProblem(TPZSBMatrix<complex<double> > &B , TPZVec <complex<double> > &eigenvalues)
{
    if (  this->fRow != B.Rows() && this->fCol != B.Cols() )
    {
        TPZMatrix<double>::Error(__PRETTY_FUNCTION__, "SolveGeneralisedEigenProblem <Uncompatible Dimensions>" );
    }
    
#ifdef USING_LAPACK
    char jobz = 'N'; //do NOT compute eigenvectors
    char uplo = 'U';//assume upper triangular
    int n = this->Dim();
    int ka = this->fBand;
    int kb = B.fBand;
    int ldab = this->fBand + 1;
    int ldbb = this->fBand + 1;
    TPZVec<double> w(0,0.);
    w.Resize( this->Dim() );
    TPZFMatrix <complex <double> > z( this->Dim() ,this->Dim() );
    int ldz = this->Dim();
    TPZVec <complex <double> > work( this->Dim() );
    TPZVec < double > rwork( 3 *this->Dim() );
    int info = -666;
    
    zhbgv_(&jobz, &uplo, &n, &ka, &kb, (__CLPK_doublecomplex *)fDiag.begin(), &ldab,  (__CLPK_doublecomplex *)B.fDiag.begin(), &ldbb, w.begin(), (__CLPK_doublecomplex *)&z(0,0), &ldz, (__CLPK_doublecomplex *)work.begin(),rwork.begin(), &info);
    if( info > 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <The algorithm failed to converge>");
    }
    else if( info < 0){
        TPZMatrix<complex<double> >::Error(__PRETTY_FUNCTION__,"Solve_EigenProblem <Invalid argument. Check info value for more information>");
    }
    eigenvalues.Resize( this->Dim() );
    for (int i = 0; i < this->Dim() ; i++) {
        eigenvalues[i] = w[i];
    }
    
#endif
    
    return( 1 );
}

/** @} */
#endif


// Inicializando os templates
template class TPZSBMatrix<float>;
template class TPZSBMatrix<double>;
template class TPZSBMatrix< complex<float> >;
template class TPZSBMatrix< complex<double> >;
template class TPZSBMatrix<long double>;

