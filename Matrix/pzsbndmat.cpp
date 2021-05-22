/**
 * @file
 * @brief Contains the implementation of the TPZSBMatrix methods.
 */

#include <math.h>
#include <stdlib.h>
#include <random>
#include "pzfmatrix.h"
#include "pzsbndmat.h"

#ifdef USING_LAPACK
/** CBlas Math Library */
#include "TPZLapack.h"
#include "TPZLapackEigenSolver.h"
#define BLAS_MULT
#endif

#include <sstream>
#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzsbmatrix");
#endif

using namespace std;

/*******************/
/*** TPZSBMatrix ***/

/**************************** PUBLIC ****************************/

/*****************************/
/*** Construtor (int) ***/

template<class TVar>
TPZSBMatrix<TVar>::TPZSBMatrix( int64_t dim, int64_t band )
: TPZRegisterClassId(&TPZSBMatrix::ClassId),
TPZMatrix<TVar>( dim, dim )
{
    fBand = ( band > (dim - 1) ? (dim - 1) : band );
    fDiag.resize(Size());
    
    Zero();
}



/** Fill the matrix with random values (non singular matrix) */
template <class TVar>
void TPZSBMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric) {
    if (nrow != ncol || symmetric == 0) {
        DebugStop();
    }
    fBand = nrow/10;
    if (fBand == 0) {
        fBand = nrow-1;
    }
    Resize(nrow, ncol);
    
    int64_t i, j;
    TVar val = 0, sum;
    /** Fill data */
    for(i=0;i<this->Rows();i++) {
        sum = 0.0;
        int64_t jmax = i+fBand+1;
        if (jmax >= this->Rows()) {
            jmax = this->Rows();
        }
        for (j=0; j<i; j++) {
            sum += fabs(GetVal(i, j));
        }
        for(j=i;j<jmax;j++) {
            val = this->GetRandomVal();
            if constexpr (is_complex<TVar>::value){
                if(j==i) val = fabs(val);
            }
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
TPZSBMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar& value )
{
    // inicializando row e col para trabalhar com a triangular superior
    auto val = value;
    int64_t row(r),col(c);
    if ( row > col ){
        //hermitian matrix
        if constexpr (is_complex<TVar>::value){
            val = std::conj(value);
        }
        this->Swap( &row, &col );
    }
    
    int64_t index;
    if ( (index = col-row) > fBand )
    {
#ifdef PZDEBUG
        if (value != (TVar) 0) {
            DebugStop();
        }
#endif
        return( 0 );        // O elemento esta fora da banda.
    }
    fDiag[ Index(row,col) ] = val;
    this->fDecomposed = 0;
    return( 1 );
}

template<class TVar>
const TVar
TPZSBMatrix<TVar>::GetVal(const int64_t row,const int64_t col ) const
{
    if ( row > col ){
        if (auto index = row-col; index > fBand ){
            return (TVar) 0;//out of band
        }else if constexpr(is_complex<TVar>::value){
            return( std::conj(fDiag[Index(col,row)]) );
        }else{
            return( fDiag[ Index(col,row) ] );
        }
    }
    
    if (auto index = col-row; index > fBand ){
        return (TVar) 0;//out of band
    }else{
        return( fDiag[ Index(row,col) ] );
    }
}

template<>
std::complex<float>
&TPZSBMatrix< std::complex< float> >::operator()(const int64_t r,const int64_t c )
{
    
    // inicializando row e col para trabalhar com a triangular superior
    int64_t row(r),col(c);
    bool mustConj = false;
    if ( row > col ){
        this->Swap( &row, &col );
        mustConj = true;
    }
    
    int64_t index;
    if ( (index = col-row) > fBand )
        return( this->gZero );        // O elemento esta fora da banda.
    if( mustConj ){
        static std::complex<float> cpVal;
        cpVal = std::conj( fDiag[ Index(row,col) ] );
        return cpVal;
    }
    else
        return( fDiag[ Index(row,col) ] );
}

template<>
std::complex<double>
&TPZSBMatrix< std::complex<double> >::operator()(const int64_t r,const int64_t c )
{
    
    // inicializando row e col para trabalhar com a triangular superior
    int64_t row(r),col(c);
    bool mustConj = false;
    if ( row > col ){
        this->Swap( &row, &col );
        mustConj = true;
    }
    
    int64_t index;
    if ( (index = col-row) > fBand )
        return( this->gZero );        // O elemento esta fora da banda.
    if( mustConj ){
        static std::complex<double> cpVal;
        cpVal = std::conj( fDiag[ Index(row,col) ] );
        return cpVal;
    }
    else
        return( fDiag[ Index(row,col) ] );
}


template<class TVar>
TVar &TPZSBMatrix<TVar>::operator()(int64_t row, int64_t col)
{
    // row and col are set to work with upper triang mat
    if ( row > col )
        this->Swap( &row, &col );
    
    int64_t index;
    if ( (index = col-row) > fBand )
        return this->gZero;        // element out of band
    
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
    
    for ( int64_t row = 0; row < A.Rows(); row++)
    {
        out << "\t";
        for ( int64_t col = 0; col < A.Cols(); col++ )
            out << A.Get( row, col) << "  ";
        out << "\n";
    }
    
    return  out << "\n";
}

/******** Operacoes com matrizes BANDA SIMETRICA  ********/

/******************/
/*** Operator + ***/

template<class TVar>
TPZSBMatrix<TVar>
TPZSBMatrix<TVar>::operator+(const TPZSBMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
    auto res(*this);
    const auto size = res.fDiag.size();
    for(auto i = 0; i < size; i++) res.fDiag[i] += A.fDiag[i];
    return res;
}

/******************/
/*** Operator - ***/

template<class TVar>
TPZSBMatrix<TVar> 
TPZSBMatrix<TVar>::operator-(const TPZSBMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
    auto res(*this);
    const auto size = res.fDiag.size();
    for(auto i = 0; i < size; i++) res.fDiag[i] -= A.fDiag[i];
    return res;
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator+=(const TPZSBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
    *this = *this+A;
    return *this;
}

/*******************/
/*** Operator -= ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator-=(const TPZSBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() || fBand != A.fBand)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZSBMatrix ) <incompatible dimensions>" );
    *this = *this-A;
    return *this;
}

template<class TVar>
void TPZSBMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
                                const TVar alpha,const TVar beta ,const int opt) const {
    // Computes z = beta * y + alpha * opt(this)*x
    //          z and x cannot overlap in memory
    if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
       this->Error(__PRETTY_FUNCTION__, "TPZSBMatrix::MultAdd <matrixs with incompatible dimensions>" );
    if(x.Cols() != y.Cols() ||x.Rows() != y.Rows()) {
       this->Error(__PRETTY_FUNCTION__,"TPZSBMatrix::MultAdd incompatible dimensions\n");
    }
    this->PrepareZ(y,z,beta,opt);
    const int64_t rows = this->Rows();
    const int64_t xcols = x.Cols();
    auto MyGetVal = [this](const int row, const int col, const bool opt){
        if constexpr (is_complex<TVar>::value){
            if(!opt) return this->GetVal(row,col);
            else return this->GetVal(col,row);
        }else{
            return this->GetVal(row,col);
        }
    };
    
    for (auto ic = 0; ic < xcols; ic++) {
        for (auto r = 0; r < rows; r++ ) {
            const auto begin = MAX( r - fBand, 0 );
            const auto end   = MIN( r + fBand + 1, rows );
            TVar val = z.GetVal(r,ic);
            for ( int64_t i = begin ; i < end; i++ ){
                val += MyGetVal( r, i, opt ) * x.GetVal(i, ic );
            }
            val *= alpha;
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
    auto res(*this);
    for(auto &el : res.fDiag) el*=value;
    return res;
}

/******************************/
/*** Operator += ( REAL ) ***/

/******************************/
/*** Operator *= ( REAL ) ***/

template<class TVar>
TPZSBMatrix<TVar> &
TPZSBMatrix<TVar>::operator*=(const TVar value )
{
    int64_t siz= fDiag.size();
    for (int64_t i=0; i<siz; i++) {
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
TPZSBMatrix<TVar>::Resize(const int64_t newDim ,const int64_t)
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
TPZSBMatrix<TVar>::Redim(const int64_t newDim ,const int64_t otherDim)
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
    int64_t siz= fDiag.size();
    for (int64_t i=0; i<siz; i++) {
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
TPZSBMatrix<TVar>::SetBand(int64_t newBand )
{
    if ( this->fBand == newBand )
        return( 1 );
    // AQUI!!!
    int64_t nB = newBand;
    if ( newBand > (this->Dim() - 1) )
    {
        //newBand = this->Dim()-1;
      nB = this->Dim()-1;
    }
    //fBand = newBand;
    fBand = nB;
    fDiag.resize(Size());
    Zero();
    
    return( 1 );
}

/********************* Resolucao de sistemas *********************/



/**************************/
/*** Decompose Cholesky ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    return Decompose_Cholesky();
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky()
{
    return TPZMatrix<TVar>::Decompose_Cholesky();
}

#ifdef USING_LAPACK
template<>
int
TPZSBMatrix<std::complex< float > >::Decompose_Cholesky()
{
    if (fDecomposed == ECholesky) {
        return 1;
    }
    if (fDecomposed != ENoDecompose) {
        DebugStop();
    }
    
    char uplo[] = "Upper";
    int n = this->Dim();
    int lda = this->fBand + 1;
    int kd = this->fBand;
    int info = -666;
    
    cpbtrf_(uplo, &n, &kd , (varfloatcomplex*) fDiag.begin(), &lda, &info);
    if( info > 0){
       this->Error(__PRETTY_FUNCTION__,"Decompose_Cholesky <The matrix is not positive definite>");
    }
    else if ( info < 0){
       this->Error(__PRETTY_FUNCTION__,"Decompose_Cholesky <Invalid argument. Check info value for more information>");
    }
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
}

//final_ok
template<>
int
TPZSBMatrix<std::complex< double > >::Decompose_Cholesky()
{
    if (fDecomposed == ECholesky) {
        return 1;
    }
    if (fDecomposed != ENoDecompose) {
        DebugStop();
    }
    
    char uplo[] = "Upper";
    int n = this->Dim();
    int lda = this->fBand + 1;
    int kd = this->fBand;
    int info = -666;
    zpbtrf_(uplo, &n, &kd, (vardoublecomplex *) fDiag.begin(), &lda, &info);
    if( info > 0){
       this->Error(__PRETTY_FUNCTION__,"Decompose_Cholesky <The matrix is not positive definite>");
    }
    else if ( info < 0){
       this->Error(__PRETTY_FUNCTION__,"Decompose_Cholesky <Invalid argument. Check info value for more information>");
    }
    
    this->fDecomposed  = ECholesky;
    this->fDefPositive = 1;
    return( 1 );
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
#endif

/**********************/
/*** Decompose LDLt ***/
template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    return Decompose_LDLt();
}

template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_LDLt()
{
    
    if (  this->fDecomposed ) this->Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
    
    int64_t j,k,l, begin,end;
    TVar sum;
    
    for ( j = 0; j < this->Dim(); j++ )
    {
        //Print("curernt");
        sum=0.;
        
        begin = MAX( int64_t(j - fBand), 0 );
        //cout<<"begin="<<begin<<"\n";
        for ( k=begin; k<j; k++)
        {
            if constexpr(is_complex<TVar>::value){
                sum=sum-GetVal(k,k)*GetVal(k,j)*std::conj(GetVal(k,j));
            }else{
                sum=sum-GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
            }
            //cout<<"(k,j)"<<k<<" "<<j<<"\n";
        }
        
        
        //	 operator()(j,j)=GetVal(j,j)+sum;
        PutVal(j,j,GetVal(j,j)+sum);
        //cout<<"\n(j,j)"<<j<<" "<<j<<"\n\n";
        for ( k=0; k<j; k++)
        {
            end   = MIN( int64_t(k + fBand )+1, this->Dim() );
            for( l=j+1; l<end;l++)
            {
                if constexpr(is_complex<TVar>::value){
                    PutVal(l,j, GetVal(l,j)-
                           GetVal(k,k)*GetVal(j,k)*std::conj(GetVal(l,k)));
                }else{
                    PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
                }
                /*cout<<"end="<<end<<"\n";
                 cout<<"(l,j)"<<l<<" "<<j<<"\n";
                 cout<<"(j,k)"<<j<<" "<<k<<"\n";
                 cout<<"(l,k)"<<l<<" "<<k<<"\n\n";
                 */
            }
        }
        
        if ( IsZero(GetVal(j,j)) )this->Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Zero on diagonal>" );
        end  = MIN( int64_t(j + fBand )+1, this->Dim() );
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
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
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
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        //    cblas_stbsv(<#const enum CBLAS_ORDER __Order#>, <#const enum CBLAS_UPLO __Uplo#>, <#const enum CBLAS_TRANSPOSE __TransA#>, <#const enum CBLAS_DIAG __Diag#>, <#const int __N#>, <#const int __K#>, <#const float *__A#>, <#const int __lda#>, <#float *__X#>, <#const int __incX#>)
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<float> >::Subst_Forward( TPZFMatrix<std::complex<float> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<float> *bptr = &(*B)(0,ic);
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<double> >::Subst_Forward( TPZFMatrix<std::complex<double> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<double> *bptr = &(*B)(0,ic);
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<float>::Subst_Backward( TPZFMatrix<float>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        float *bptr = &(*B)(0,ic);
        cblas_stbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<>
int
TPZSBMatrix<double>::Subst_Backward( TPZFMatrix<double>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<float> >::Subst_Backward( TPZFMatrix<std::complex<float> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<float> *bptr = &(*B)(0,ic);
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<>
int
TPZSBMatrix<std::complex<double> >::Subst_Backward( TPZFMatrix<std::complex<double> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<double>  *bptr = &(*B)(0,ic);
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<float>::Subst_LForward( TPZFMatrix<float>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        float *bptr = &(*B)(0,ic);
        cblas_stbsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<double>::Subst_LForward( TPZFMatrix<double>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<float> >::Subst_LForward( TPZFMatrix<std::complex<float> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<float> *bptr = &(*B)(0,ic);
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<double> >::Subst_LForward( TPZFMatrix<std::complex<double> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<double> *bptr = &(*B)(0,ic);
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<float>::Subst_LBackward( TPZFMatrix<float>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        float *bptr = &(*B)(0,ic);
        cblas_stbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<>
int
TPZSBMatrix<double>::Subst_LBackward( TPZFMatrix<double>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        double *bptr = &(*B)(0,ic);
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

template<>
int
TPZSBMatrix<std::complex<float> >::Subst_LBackward( TPZFMatrix<std::complex<float> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<float> *bptr = &(*B)(0,ic);
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}


template<>
int
TPZSBMatrix<std::complex<double> >::Subst_LBackward( TPZFMatrix<std::complex<double> >*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> uncompatible matrices") ;
    }
    int n = Rows();
    int kd = fBand;
    int lda = 1+fBand;
    int64_t bcols = B->Cols();
    for(int64_t ic=0; ic<bcols; ic++)
    {
        std::complex<double>  *bptr = &(*B)(0,ic);
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fDiag[0], lda, bptr , 1);
    }
    return 1;
}

#endif

template<class TVar>
int
TPZSBMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>*B ) const
{
    if ( (B->Rows() != this->Dim()) || ! this->fDecomposed)
    {
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    }
    return TPZMatrix<TVar>::Subst_Forward(B);
}
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
       this->Error(__PRETTY_FUNCTION__,"Subst_LForward-> uncompatible matrices") ;
    
    int64_t size = fBand + 1;
    int64_t i,j,k;
    for ( k = 0; k < this->Dim(); k++ )
    {
        for ( j = 0; j < B->Cols(); j++ )
        {
            // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
            
            int64_t imin = k-fBand;
            if(imin < 0) imin = 0;
            
            int64_t end=(k-fBand>0)? fBand:k;  //misael
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
       this->Error(__PRETTY_FUNCTION__,"Subst_Diag-> uncompatible matrices") ;
    
    
    int64_t size = fBand + 1;
    for ( int64_t k = 0; k < this->Dim(); k++ )
        for ( int64_t j = 0; j < B->Cols(); j++ )
            B->PutVal( k, j, B->GetVal( k, j) / GetVal(k,k) );
    
    return( 1 );
}

template<class TVar>
int TPZSBMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
       this->Error(__PRETTY_FUNCTION__,"Subst_Forward-> uncompatible matrices") ;
    
    return TPZMatrix<TVar>::Subst_Backward(B);
    return ( 1 ) ;
    
}

template<class TVar>
int TPZSBMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const
{
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed )
       this->Error(__PRETTY_FUNCTION__,"Subst_LBackward-> uncompatible matrices") ;
    
    int64_t k,j,jmax,stepcol=fBand+2;
    
    for(k=0; k<B->Cols() ; k++)
    {
        for(int64_t i=this->Rows()-1; i>=0; i--)
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
TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &eigenvalues, TPZFMatrix < CTVar > &eigenVectors)
{
   if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                 && is_arithmetic_pz<TVar>::value){
       TPZLapackEigenSolver<TVar> solver;
       return solver.SolveEigenProblem(*this,eigenvalues,eigenVectors);
   }else{
       PZError<<__PRETTY_FUNCTION__;
       PZError<<"\nERROR: Incompatible types.\nAborting...\n";
       DebugStop();
       return -1;
   }
}

template< class TVar>
int
TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &w)
{
   if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                 && is_arithmetic_pz<TVar>::value){
       TPZLapackEigenSolver<TVar> solver;
       return solver.SolveEigenProblem(*this,w);
   }else{
       PZError<<__PRETTY_FUNCTION__;
       PZError<<"\nERROR: Incompatible types.\nAborting...\n";
       DebugStop();
       return -1;
   }
}

template< class TVar>
int
TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &B , TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors)
{
    if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveGeneralisedEigenProblem(*this,B,w,eigenVectors);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}
template< class TVar>
int
TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix<TVar> &B , TPZVec < CTVar > &w)
{
   if constexpr ((std::is_same_v<RTVar,float> || std::is_same_v<RTVar,double>)
                  && is_arithmetic_pz<TVar>::value){
        TPZLapackEigenSolver<TVar> solver;
        return solver.SolveGeneralisedEigenProblem(*this,B,w);
    }else{
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nERROR: Incompatible types.\nAborting...\n";
        DebugStop();
        return -1;
    }
}

/** @} */
#else
#define NON_LAPACK \
  PZError<<__PRETTY_FUNCTION__<<" requires Lapack\n";\
  PZError<<" Set either USING_LAPACK=ON or USING_MKL=ON on CMake ";\
  PZError<<" when configuring NeoPZ library"<<std::endl;\
  DebugStop();\
  return -1;


template<class TVar>
int TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors){NON_LAPACK}

template<class TVar>
int TPZSBMatrix<TVar>::SolveEigenProblem(TPZVec < CTVar > &w){NON_LAPACK}

template<class TVar>
int TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < CTVar > &w, TPZFMatrix < CTVar > &eigenVectors){NON_LAPACK}

template<class TVar>
int TPZSBMatrix<TVar>::SolveGeneralisedEigenProblem(TPZSBMatrix< TVar > &B , TPZVec < CTVar > &w){NON_LAPACK}
#undef NON_LAPACK
#endif

template<class TVar>
int TPZSBMatrix<TVar>::ClassId() const{
    return Hash("TPZSBMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

// Inicializando os templates
template class TPZSBMatrix<float>;
template class TPZSBMatrix<double>;
template class TPZSBMatrix<long double>;
template class TPZSBMatrix<std::complex<float>>;
template class TPZSBMatrix<std::complex<double>>;
template class TPZSBMatrix<std::complex<long double>>;
