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
    this->fSymProp = SymProp::Herm;
    fBand = ( band > (dim - 1) ? (dim - 1) : band );
    fStorage.resize(Size());
    
    Zero();
}


template<class TVar>
void TPZSBMatrix<TVar>::SetSymmetry (SymProp sp){
    if(sp == SymProp::NonSym){
        PZError<<__PRETTY_FUNCTION__
               <<"\nTrying to set matrix with symmetric storage as non symmetric\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
    TPZBaseMatrix::SetSymmetry(sp);
}

/** Fill the matrix with random values (non singular matrix) */
template <class TVar>
void TPZSBMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, SymProp sym) {
    if (nrow != ncol || sym == SymProp::NonSym) {
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
    SetSymmetry(sym);
    
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
            if(this->fSymProp == SymProp::Herm) {val = std::conj(value);}
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
    fStorage[ Index(row,col) ] = val;
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
            if(this->fSymProp == SymProp::Herm) {
                return( std::conj(fStorage[Index(col,row)]) );
            } else {
                return fStorage[Index(col,row)];
            }
        }else{
            return( fStorage[ Index(col,row) ] );
        }
    }
    
    if (auto index = col-row; index > fBand ){
        return (TVar) 0;//out of band
    }else{
        return( fStorage[ Index(row,col) ] );
    }
}


template<class TVar>
TVar &TPZSBMatrix<TVar>::operator()(int64_t row, int64_t col)
{
    // row and col are set to work with upper triang mat
    if ( row > col ){
        this->Swap( &row, &col );
        if constexpr (is_complex<TVar>::value){
            if(this->fSymProp == SymProp::Herm){
                PZError<<__PRETTY_FUNCTION__
                       <<"\nTrying to access lower triang hermitian mat by reference\n"
                       <<"Aborting..."
                       <<std::endl;
                DebugStop();
            }
        }
    }
    
    int64_t index;
    if ( (index = col-row) > fBand )
        return this->gZero;        // element out of band
    
    return( fStorage[ Index(row,col) ] );

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
    const auto size = res.fStorage.size();
    for(auto i = 0; i < size; i++) res.fStorage[i] += A.fStorage[i];
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
    const auto size = res.fStorage.size();
    for(auto i = 0; i < size; i++) res.fStorage[i] -= A.fStorage[i];
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
    this->MultAddChecks(x,y,z,alpha,beta,opt);
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
    for(auto &el : res.fStorage) el*=value;
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
    int64_t siz= fStorage.size();
    for (int64_t i=0; i<siz; i++) {
        fStorage[i] *= value;
    }
    
    return( *this );
}

template<class TVar>
int
TPZSBMatrix<TVar>::Resize(const int64_t newDim ,const int64_t)
{
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Resize should not be called for Band matrices\n";
    DebugStop();
    
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
    // if (newDim != otherDim) {
    //     DebugStop();
    // }
    
    // if ( newDim != this->Dim() )
    // {
    //     TPZMatrix<TVar>::Redim(newDim,newDim);
    //     if (fBand > newDim-1) {
    //         fBand = newDim-1;
    //     }
    //     fStorage.resize(Size());
    // }
    
    // Zero();
    // this->fDecomposed = ENoDecompose;
    // this->fDefPositive = 0;

    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Redim should not be called for Band matrices\n";
    DebugStop();
    return( 1 );
}

template<class TVar>
int
TPZSBMatrix<TVar>::Zero()
{
    int64_t siz= fStorage.size();
    for (int64_t i=0; i<siz; i++) {
        fStorage[i] = (TVar)0.;
    }
    
    
    this->fDecomposed = ENoDecompose;
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
    fStorage.resize(Size());
    Zero();
    
    return( 1 );
}

/********************* Resolucao de sistemas *********************/



template<class TVar>
int
TPZSBMatrix<TVar>::Decompose_Cholesky()
{
    if (  this->fDecomposed && this->fDecomposed != ECholesky) this->Error( "Decompose_Cholesky <Matrix already Decomposed>" );
    if (  this->fDecomposed ) return ECholesky;
    if ( this->Rows()!=this->Cols() ) this->Error( "Decompose_Cholesky <Matrix must be square>" );
    //return 0;

#ifdef PZDEBUG
    const bool cond =
        (is_complex<TVar>::value && this->GetSymmetry() == SymProp::Sym) ||
        this->GetSymmetry() == SymProp::NonSym;
    if (cond){
        PZError<<__PRETTY_FUNCTION__
               <<"\nCalling Cholesky decomposition on non symmetric matrix! Aborting..."
               <<std::endl;
        DebugStop();
    }
#endif
    int64_t dim=this->Rows();
    for (int64_t i=0 ; i<dim; i++) {
        for(int64_t k=0; k<i; k++) {//diagonal elements
            TVar sum = 0;
            if constexpr (is_complex<TVar>::value){
                sum += GetVal(i,k)*std::conj(GetVal(i,k));
            }else{
                sum += GetVal(i,k)*GetVal(i,k);
            }
            PutVal( i,i,GetVal(i,i)-sum );
        }
        TVar tmp = sqrt(GetVal(i,i));
        PutVal( i,i,tmp );
        for (int64_t j=i+1;j<dim; j++) {//off-diagonal elements
            for(int64_t k=0; k<i; k++) {
                TVar sum = 0.;
                if constexpr (is_complex<TVar>::value){
                    sum += GetVal(i,k)*std::conj(GetVal(j,k));
                }else{
                    sum += GetVal(i,k)*GetVal(j,k);
                }
                PutVal( i,j,GetVal(i,j)-sum);
            }
            TVar tmp2 = GetVal(i,i);
            if ( IsZero(tmp2) ) {
                this->Error( "Decompose_Cholesky <Zero on diagonal>" );
            }
            PutVal(i,j,GetVal(i,j)/GetVal(i,i) );
            if constexpr (is_complex<TVar>::value){
                PutVal(j,i,std::conj(GetVal(i,j)));
            }else{
                PutVal(j,i,GetVal(i,j));
            }
        }
    }
    this->fDecomposed = ECholesky;
    return ECholesky;
    
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
    lapack_int n = this->Dim();
    lapack_int lda = this->fBand + 1;
    lapack_int kd = this->fBand;
    lapack_int info = -666;
    
    cpbtrf_(uplo, &n, &kd , (varfloatcomplex*) fStorage.begin(), &lda, &info);
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
    lapack_int n = this->Dim();
    lapack_int lda = this->fBand + 1;
    lapack_int kd = this->fBand;
    lapack_int info = -666;
    zpbtrf_(uplo, &n, &kd, (vardoublecomplex *) fStorage.begin(), &lda, &info);
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
    lapack_int n = Dim();
    lapack_int kd = fBand;
    lapack_int nrhs = 0;
    float *ab = &fStorage[0];
    lapack_int ldab = fBand+1;
    float b = 0;
    lapack_int info;
    
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
    lapack_int n = Dim();
    lapack_int kd = fBand;
    lapack_int nrhs = 0;
    double *ab = &fStorage[0];
    lapack_int ldab = fBand+1;
    double b = 0;
    lapack_int info;
    
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
TPZSBMatrix<TVar>::Decompose_LDLt()
{
    
    if (  this->fDecomposed ) this->Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Matrix already Decomposed>" );
    
    if constexpr(is_complex<TVar>::value) {
        std::cout << "LDLt decomposition for complex numbers not implemented\n";
        DebugStop();
    }
    int64_t j,k,l, begin,end;
    TVar sum;
    
    for ( j = 0; j < this->Rows(); j++ )
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
        
        if ( IsZero(GetVal(j,j)) )
        {
            this->Error(__PRETTY_FUNCTION__, "Decompose_LDLt <Zero on diagonal>" );
        }
        end  = MIN( int64_t(j + fBand )+1, this->Dim() );
        //cout<<"end="<<end<<"\n";
        for( l=j+1; l<end;l++)
        {
            //cout<<"(l,j)"<<l<<" "<<j<<"\n";
            PutVal( l,j,GetVal(l,j)/GetVal(j,j) ) ;
        }
    }
    this->fDecomposed  = ELDLt;
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
        cblas_stbsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_stbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_stbsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasConjTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_stbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_dtbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ctbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
        cblas_ztbsv(CblasColMajor, CblasUpper, CblasNoTrans, CblasUnit, n, kd, &fStorage[0], lda, bptr , 1);
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
    for ( int64_t r = 0; r < this->Dim(); r++ ) {
        TVar pivot = GetVal( r, r );
        for ( int64_t c = 0; c < B->Cols();  c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
            //
            TVar sum = 0.0;
            for ( int64_t i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
            
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, (B->GetVal(r, c) - sum) / pivot );
        }
    }
    return( 1 );
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
    
    if ( (B->Rows() != this->Dim()) || !this->fDecomposed || this->fDecomposed != ECholesky) {
        this->Error(__PRETTY_FUNCTION__,"Subst_Backward-> wrong parameters") ;
        return( -1 );
    }
    for ( int64_t r = this->Dim()-1;  r >= 0;  r-- ) {
        TVar pivot = GetVal( r, r );
        for ( int64_t c = 0; c < B->Cols(); c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
            //
            TVar sum = 0.0;
            for ( int64_t i = this->Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, (B->GetVal(r, c) - sum) / pivot );
        }
    }
    return( 1 );

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
    fStorage.resize(0);
    this->fDecomposed = ENoDecompose;
    return( 1 );
}

/************/
/*** Copy ***/

template<class TVar>
void
TPZSBMatrix<TVar>::Copy(const TPZSBMatrix<TVar> &A )
{
    TPZMatrix<TVar>::operator=(A);
    std::cout << __PRETTY_FUNCTION__ << " Please implement me!\n";
    DebugStop();
    this->fBand = A.fBand;
    this->fStorage = A.fStorage;
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
