/**
 * @file
 * @brief Contains the implementation of the TPZFBMatrix methods.
 */

#include <math.h>
#include "pzfmatrix.h"
#include "pzbndmat.h"
#ifdef USING_LAPACK
/** CBlas Math Library */
#include "TPZLapack.h"

#endif

#include <random>

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzfbmatrix");
#endif

using namespace std;

/*******************/
/*** Constructor ***/
template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix()
: TPZRegisterClassId(&TPZFBMatrix::ClassId),
TPZMatrix<TVar>( 0, 0 ), fElem()
{
	fBandLower = 0;
    fBandUpper = 0;
}

/********************/
/*** Constructors ***/

template<class TVar>
TPZFBMatrix<TVar>::TPZFBMatrix( int64_t dim, int64_t band_width )
: TPZRegisterClassId(&TPZFBMatrix::ClassId),
TPZMatrix<TVar>( dim, dim ), fElem(dim*(3*band_width+1),0.), fBandLower(band_width), fBandUpper(band_width)
{
}



/*********************************/



/******************/
/*** Destructor ***/
template<class TVar>
TPZFBMatrix<TVar>::~TPZFBMatrix ()
{
}



/***********/
/*** Put ***/

template<class TVar>
int
TPZFBMatrix<TVar>::Put(const int64_t row,const int64_t col,const TVar& value )
{
	if ( (row >= Dim()) || (col >= Dim()) || row<0 || col<0 || row-col > fBandLower || col-row > fBandUpper)
    {
		cout << "TPZFBMatrix::Put: " << row << "," << col << "," << Dim();
		cout << "\n";
        DebugStop();
		return( 0 );
    }
	
	return( PutVal( row, col, value ) );
}



/***********/
/*** Get ***/

template<class TVar>
const TVar
TPZFBMatrix<TVar>::Get(const int64_t row,const int64_t col ) const
{
	if ( (row >= Dim()) || (col >= Dim()) )
		this->Error(__PRETTY_FUNCTION__, "Get <indices out of band matrix range>" );
	
	return( GetVal( row, col ) );
}



/******** Operacoes com matrizes FULL BAND  ********/

/******************/

 template<class TVar>
TPZFBMatrix<TVar>
TPZFBMatrix<TVar>::operator+(const TPZFBMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() ||
         fBandUpper != A.fBandUpper ||
         fBandLower != A.fBandLower)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZFBMatrix ) <incompatible dimensions>" );
    auto res(*this);
    const auto size = res.fElem.size();
    for(auto i = 0; i < size; i++) res.fElem[i] += A.fElem[i];
    return res;
}

/******************/
/*** Operator - ***/

template<class TVar>
TPZFBMatrix<TVar> 
TPZFBMatrix<TVar>::operator-(const TPZFBMatrix<TVar> &A ) const
{
    if ( this->Dim() != A.Dim() ||
         fBandUpper != A.fBandUpper ||
         fBandLower != A.fBandLower)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZFBMatrix ) <incompatible dimensions>" );
    auto res(*this);
    const auto size = res.fElem.size();
    for(auto i = 0; i < size; i++) res.fElem[i] -= A.fElem[i];
    return res;
}

/*******************/
/*** Operator += ***/

template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator+=(const TPZFBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() ||
         fBandUpper != A.fBandUpper ||
         fBandLower != A.fBandLower)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZFBMatrix ) <incompatible dimensions>" );
    *this = *this+A;
    return *this;
}

/*******************/
/*** Operator -= ***/

template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator-=(const TPZFBMatrix<TVar> &A )
{
    if ( this->Dim() != A.Dim() ||
         fBandUpper != A.fBandUpper ||
         fBandLower != A.fBandLower)
       this->Error(__PRETTY_FUNCTION__,"operator+( TPZFBMatrix ) <incompatible dimensions>" );
    *this = *this-A;
    return *this;
}


/******** Operacoes com MATRIZES GENERICAS ********/

/*******************/
/*** MultiplyAdd ***/
//
//  perform a multiply add operation to be used by iterative solvers
//

template<class TVar>
void TPZFBMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						  const TVar alpha,const TVar beta ,const int opt) const {
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
		this->Error(__PRETTY_FUNCTION__, "TPZFBMatrix::MultAdd <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Rows() != y.Rows()) {
		this->Error(__PRETTY_FUNCTION__,"TPZFBMatrix::MultAdd incompatible dimensions\n");
	}
	this->PrepareZ(y,z,beta,opt);
	int64_t rows = this->Rows();
	int64_t xcols = x.Cols();
	int64_t ic, r;
	if(opt == 0) {
		for (ic = 0; ic < xcols; ic++) {
			int64_t begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBandLower, 0 );
				end   = MIN( r + fBandUpper + 1, Dim() );
				TVar val = z.GetVal(r,ic);
				// Calcula um elemento da resposta.
				for ( int64_t i = begin ; i < end; i++ ) val += alpha * GetVal( r, i ) * x.GetVal( i, ic );
				z.PutVal( r, ic, val );
			}
		}
	} else {
		for (ic = 0; ic < xcols; ic++) {
			int64_t begin, end;
			for ( r = 0; r < rows; r++ ) {
				begin = MAX( r - fBandLower, 0 );
				end   = MIN( r + fBandUpper + 1, Dim() );
				TVar val = z.GetVal(r,ic);
				// Calcula um elemento da resposta.
				for ( int64_t i = begin ; i < end; i++ ) val += alpha * GetVal( i, r ) * x.GetVal( i, ic );
				z.PutVal( r, ic, val );
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
	auto res(*this);
  for(auto &el : res.fElem) el *=value;
	return res;
}

template<class TVar>
TPZFBMatrix<TVar> &
TPZFBMatrix<TVar>::operator*=(const TVar value )
{
	if ( value != (TVar)1.0 )
    {
        int64_t sz = fElem.size();
        for (int64_t i=0; i<sz; i++) {
            fElem[i] *= value;
        }
    }
	
	return( *this );
}

template <class TVar>
void TPZFBMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, SymProp sp) {
    if (nrow != ncol) {
        DebugStop();
    }
    int Band = nrow/3;
    if (Band == 0) {
        Band = nrow;
    }
    TPZFBMatrix A(nrow, Band);
    *this = A;
    
    for (int64_t i=0; i<nrow; i++) {
        int64_t jmin = i-Band;
        if (jmin< 0) {
            jmin = 0;
        }
        int64_t jmax = i+Band+1;
        if (jmax > nrow) {
            jmax = nrow;
        }
        TVar sum = 0.;
        int64_t j = jmin;
        if (sp != SymProp::NonSym) {
            for (; j<i; j++) {
                if constexpr (is_complex<TVar>::value){
                    if(sp == SymProp::Herm){PutVal(i, j, std::conj(GetVal(j,i)));}
                    else {PutVal(i, j, GetVal(j,i));}
                }else{
                    PutVal(i, j, GetVal(j,i));
                }
                sum += fabs(GetVal(i, j));
            }
        }
        for (; j<jmax; j++) {
            TVar val = this->GetRandomVal();
            if constexpr (is_complex<TVar>::value){
                if(j==i) val = fabs(val);
            }
            PutVal(i, j, val);
            if(i!=j) sum += fabs(val);
        }
        PutVal(i, i, sum+(TVar)1.);
    }
    
//    this->Print("AutoFill = ",std::cout,EMathematicaInput);
}
/**************/
/*** Resize ***/
// DEPENDS ON THE STORAGE FORMAT

template<class TVar>
int
TPZFBMatrix<TVar>::Resize(const int64_t newRows,const int64_t newCols)
{
	if ( newRows != newCols )
		this->Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	
    Redim(newRows,newRows);
	return( 1 );
}



/*************/
/*** Redim ***/
template<class TVar>
int
TPZFBMatrix<TVar>::Redim(const int64_t newRows,const int64_t newCols )
{
	if ( newRows != newCols )
		this->Error(__PRETTY_FUNCTION__, "Resize <Band matrix must be NxN>" );
	
	//  if ( !fBand ) TPZMatrix::Error(__PRETTY_FUNCTION__, "Bandwith = NULL" );
	

    if (fBandLower > newRows-1) {
        fBandLower = newRows-1;
    }
    if (fBandUpper > newRows-1) {
        fBandUpper = newRows-1;
    }
    TPZMatrix<TVar>::Redim(newRows,newRows);
	uint64_t size = newRows*(2*fBandLower+fBandUpper + 1);
    fElem.Resize(size);
	Zero();
	
	return( 1 );
}


/***************/
/**** Zero ****/
template<class TVar>
int
TPZFBMatrix<TVar>::Zero()
{
    uint64_t size = fElem.size();

    for (int64_t i=0; i<size; i++) {
        fElem[i] = (TVar)0.;
    }
	
	this->fDecomposed = ENoDecompose;
	
	return( 1 );
}


/***************/
/*** SetBand ***/
// DEPENDS ON THE STORAGE FORMAT
template<class TVar>
int
TPZFBMatrix<TVar>::SetBand( int64_t newBand )
{
	if ( newBand >= Dim() )
		this->Error(__PRETTY_FUNCTION__, "SetBand <the band must be lower than the matrix dimension " );
	
	uint64_t newSize = Dim()*(3 * newBand + 1);
    fBandLower = newBand;
    fBandUpper = newBand;
    fElem.resize(newSize);
    Zero();
	return( 1 );
}



/********************/
/*** Transpose () ***/
template<class TVar>
void
TPZFBMatrix<TVar>::Transpose (TPZMatrix<TVar> *const T) const
{
	T->Resize( Dim(), Dim() );
	
	int64_t end, begin;
	//REAL *p = fElem;
	for ( int64_t r = 0; r < Dim(); r++ )
    {
		begin = MAX( r - fBandLower, 0 );
		end   = MIN( r + fBandUpper + 1, Dim() );
		for ( int64_t c = begin; c < end; c++ )
		{
			T->PutVal( c, r, GetVal( r, c ) );
			//			cout<<"(r,c)= "<<r<<"  "<<c<<"\n";
		}
    }
}



template<class TVar>
int TPZFBMatrix<TVar>::Decompose_LU()
{
    if (  this->fDecomposed && this->fDecomposed == ELU) {
        return ELU;
    } else if(this->fDecomposed) {
        this->Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
    }
    
    TVar nn, pivot;
    int64_t  min = ( this->Cols() < (this->Rows()) ) ? this->Cols() : this->Rows();
    
    for ( int64_t k = 0; k < min ; k++ ) {
        if (IsZero( pivot = GetVal(k, k))) this->Error( "Decompose_LU <matrix is singular>" );
        for ( int64_t i = k+1; i < this->Rows(); i++ ) {
            nn = GetVal( i, k ) / pivot;
            PutVal( i, k, nn );
            for ( int64_t j = k+1; j < this->Cols(); j++ ) PutVal(i,j,GetVal(i,j)-nn*GetVal(k,j));
        }
    }
    this->fDecomposed=ELU;
    return 1;
}

#ifdef USING_LAPACK
template<>
int
TPZFBMatrix<float>::Decompose_LU()
{
	if (  this->fDecomposed && this->fDecomposed == ELU) {
		return ELU;
	} else if(this->fDecomposed) {
		this->Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
	}
    lapack_int rows = Rows();
    lapack_int bandlower = fBandLower;
    lapack_int bandupper = fBandUpper;
    lapack_int nrhs = 0;
    lapack_int ldab = 1+2*fBandLower+fBandUpper;
    fPivot.Resize(rows, 0);
    float B;
    lapack_int info;
    
//    sgbsv_(<#__CLPK_integer *__n#>, <#__CLPK_integer *__kl#>, <#__CLPK_integer *__ku#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
//    lapack_int LAPACKE_sgbsv( int matrix_layout, lapack_int n, lapack_int kl,
//                             lapack_int ku, lapack_int nrhs, float* ab,
//                             lapack_int ldab, lapack_int* ipiv, float* b,
//                             lapack_int ldb );

    TPZManVector<lapack_int,2000> pivot(rows);
    for(int i = 0; i < rows; i++){pivot[i]=fPivot[i];}
    sgbsv_(&rows, &bandlower, &bandupper, &nrhs, &fElem[0], &ldab,&pivot[0], &B,&rows, &info);
    for(int i = 0; i < rows; i++){
        fPivot[i] = pivot[i];
    }
    //int matrix_layout = 0;
//    LAPACKE_sgbsv(matrix_layout,rows, bandlower, bandupper, nrhs, &fElem[0], ldab,&fPivot[0], &B,rows);
    
    if (info != 0) {
        DebugStop();
    }
	this->fDecomposed = ELU;
	return 1;
}
template<>
int
TPZFBMatrix<double>::Decompose_LU()
{
    if (  this->fDecomposed && this->fDecomposed == ELU) {
        return ELU;
    } else if(this->fDecomposed) {
        this->Error(__PRETTY_FUNCTION__,"TPZFBMatrix::Decompose_LU is already decomposed with other scheme");
    }
    lapack_int rows = Rows();
    lapack_int bandlower = fBandLower;
    lapack_int bandupper = fBandUpper;
    lapack_int nrhs = 0;
    lapack_int ldab = 1+2*fBandLower+fBandUpper;
    fPivot.Resize(rows, 0);
    double B;
    lapack_int info;
    
    //    sgbsv_(<#__CLPK_integer *__n#>, <#__CLPK_integer *__kl#>, <#__CLPK_integer *__ku#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    TPZManVector<lapack_int,2000> pivot(rows);
    for(int i = 0; i < rows; i++){pivot[i]=fPivot[i];}
    dgbsv_(&rows, &bandlower, &bandupper, &nrhs, &fElem[0], &ldab,&pivot[0], &B,&rows, &info);

    for(int i = 0; i < rows; i++){
        fPivot[i] = pivot[i];
    }
    if (info != 0) {
        DebugStop();
    }
    this->fDecomposed = ELU;
    return 1;
}


template<>
int TPZFBMatrix<float>::Substitution( TPZFMatrix<float> *B ) const{
    
    if (this->fDecomposed != ELU) {
        DebugStop();
    }
    lapack_int rows = Rows();
    lapack_int bandlower = fBandLower;
    lapack_int bandupper = fBandUpper;
    lapack_int nrhs = B->Cols();
    lapack_int ldab = 1+2*fBandLower+fBandUpper;
    lapack_int info;
    char notrans[]="No Transpose";

    
//    sgbtrs_(<#char *__trans#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__kl#>, <#__CLPK_integer *__ku#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    TPZManVector<lapack_int,2000> pivot(rows);
    for(int i = 0; i < rows; i++){pivot[i]=fPivot[i];}
    sgbtrs_(notrans, &rows, &bandlower, &bandupper, &nrhs, &fElem[0], &ldab, &pivot[0], &B->s(0,0), &rows, &info);
    for(int i = 0; i < rows; i++){fPivot[i] = pivot[i];}
    return( 1 );
}

template<>
int TPZFBMatrix<double>::Substitution( TPZFMatrix<double> *B ) const{
    
    if (this->fDecomposed != ELU) {
        DebugStop();
    }
    lapack_int rows = Rows();
    lapack_int bandlower = fBandLower;
    lapack_int bandupper = fBandUpper;
    lapack_int nrhs = B->Cols();
    lapack_int ldab = 1+2*fBandLower+fBandUpper;
    lapack_int info;
    char notrans[]="No Transpose";
    
    
    //    sgbtrs_(<#char *__trans#>, <#__CLPK_integer *__n#>, <#__CLPK_integer *__kl#>, <#__CLPK_integer *__ku#>, <#__CLPK_integer *__nrhs#>, <#__CLPK_real *__ab#>, <#__CLPK_integer *__ldab#>, <#__CLPK_integer *__ipiv#>, <#__CLPK_real *__b#>, <#__CLPK_integer *__ldb#>, <#__CLPK_integer *__info#>)
    TPZManVector<lapack_int,2000> pivot(rows);
    for(int i = 0; i < rows; i++){pivot[i]=fPivot[i];}
    dgbtrs_(notrans, &rows, &bandlower, &bandupper, &nrhs, &fElem[0], &ldab, &pivot[0], &B->s(0,0), &rows, &info);
    for(int i = 0; i < rows; i++){
        fPivot[i] = pivot[i];
    }
    return( 1 );
}

#endif

template<class TVar>
int TPZFBMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const{
    
    if (this->fDecomposed != ELU) {
        DebugStop();
    }
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    if ( rowb != this->Rows() )
        this->Error( "SubstitutionLU <incompatible dimensions>" );
    int64_t i;
    for ( i = 0; i < rowb; i++ ) {
        for ( int64_t col = 0; col < colb; col++ ) {
            for ( int64_t j = 0; j < i; j++ ) {
                B->PutVal( i, col, B->GetVal(i, col)-GetVal(i, j) * B->GetVal(j, col) );
            }
        }
    }
    for (int64_t col=0; col<colb; col++) {
        for ( i = rowb-1; i >= 0; i-- ) {
            for ( int64_t j = i+1; j < rowb ; j++ ) {
                B->PutVal( i, col, B->GetVal(i, col) -
                          GetVal(i, j) * B->GetVal(j, col) );
            }
            if ( IsZero( GetVal(i, i) ) ) {
                this->Error( "BackSub( SubstitutionLU ) <Matrix is singular" );
            }
            B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
        }
    }
    return( 1 );
}

template<class TVar>
int TPZFBMatrix<TVar>::ClassId() const{
    return Hash("TPZFBMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}

/************************** Private **************************/

/*************/
/*** Clear ***/

template<class TVar>
int
TPZFBMatrix<TVar>::Clear()
{
	fElem.resize(0);
	this->fRow  = this->fCol = 0;
	fBandLower = 0;
    fBandUpper = 0;
	return( 1 );
}


template class TPZFBMatrix<long double>;
template class TPZFBMatrix<double>;
template class TPZFBMatrix<float>;
template class TPZFBMatrix<int64_t>;
template class TPZFBMatrix<int>;
template class TPZFBMatrix<std::complex<long double> >;
template class TPZFBMatrix<std::complex<double> >;
template class TPZFBMatrix<std::complex<float> >;
