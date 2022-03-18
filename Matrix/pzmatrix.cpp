/**
 * @file
 * @brief Contains the implementation of the TPZMatrix<>methods.
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <sstream>
#include <set>
#include <random>
#include <limits>
//#include "pzfmatrix.h"
#include "pzmatrix.h"
#include "TPZMatrixSolver.h"
#include "TPZFMatrixRef.h"
#include "pzvec.h"
#include "pzextractval.h"

#include "fadType.h"
#include "fad.h"

#include <exception>
#include "pzlog.h"
#include <complex>

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzmatrix");
static TPZLogger loggerCheck("pz.checkconsistency");
#endif

#ifdef PZDEBUG
#define DEBUG2
#endif

using namespace std;
//template <class TVar>
//TVar TPZMatrix<TVar>::gZero = TVar(0);

template<class TVar>
TPZMatrix<TVar> &TPZMatrix<TVar>::operator=(const TPZMatrix<TVar> &cp) {
    
    TPZBaseMatrix::operator=(cp);
    return *this;
}

template<class TVar>
TPZMatrix<TVar> &TPZMatrix<TVar>::operator=(TPZMatrix<TVar> &&rval) {
    TPZBaseMatrix::operator=(rval);
    return *this;
}

template<class TVar>
const TVar TPZMatrix<TVar>::GetVal(const int64_t /*row*/, const int64_t /*col*/ ) const
{
    return (TVar)0;
}

template<class TVar>
TPZFMatrixRef<TVar> TPZMatrix<TVar>::Storage()
{
  const auto size = Size();
  const auto elem = Elem();
	TPZFMatrixRef<TVar> ref(size,Elem());
	return ref;
}
template<class TVar>
const TPZFMatrix<TVar> TPZMatrix<TVar>::Storage() const
{
  const auto size = Size();
  const auto elem = Elem();
	TPZFMatrix<TVar> ref(size,1,const_cast<TVar*>(elem),size);
	return ref;
}

template<class TVar>
void TPZMatrix<TVar>::Simetrize() {
  
  if ( Rows() != Cols() ) {
    Error( "Simetrize only work for square matrices" );
  }
  
  int64_t row,col;
  int64_t fDim1 = Rows();
  for(row=0; row<fDim1; row++) {
    for(col=row+1; col<fDim1; col++) {
      this->s(col,row) = this->s(row,col);
    }
  }
  
}

/** @brief Implements product of matrices: \f$ A*B \f$ */
template<class TVar>
TPZFMatrix<TVar> TPZMatrix<TVar>::operator*(const TPZFMatrix<TVar> &B ) const{
  TPZFMatrix<TVar> res;
  res.Redim( this->Rows(), B.Cols() );
	this->Multiply(B,res);
	return res;
}

template<class TVar>
void TPZMatrix<TVar>::CheckTypeCompatibility(const TPZMatrix<TVar>*A,
                                             const TPZMatrix<TVar>*B) const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: not implemented for your matrix type.\nAborting...\n";
  DebugStop();
}

template<class TVar>
void TPZMatrix<TVar>::PrepareZ(const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,const TVar beta,const int opt) const
{
	int64_t numeq = (opt) ? Cols() : Rows();
	int64_t xcols = y.Cols();
	int64_t ic;
	z.Resize(numeq, xcols);
	for (ic = 0; ic < xcols; ic++)
	{
		TVar *zp = &z(0,ic), *zlast = zp+numeq;
		if(beta != (TVar)0)
		{
			const TVar *yp = &y.g(0,ic);
            if(&z != &y)
			{
				memcpy((void *)(zp),(void *)(yp),numeq*sizeof(TVar));
			}
            for (int64_t i=0; i<numeq; i++) {
                z(i,ic) *= beta;
            }
		} else
		{
			while(zp != zlast)
			{
				*zp = 0.;
				zp ++;
			}
		}
	}
}

template<class TVar>
void TPZMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z, const TVar alpha,const TVar beta,const int opt) const {
	if ((!opt && Cols() != x.Rows()) || Rows() != x.Rows())
		Error( "Operator* <matrices with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Rows() != y.Rows()) {
		Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	int64_t rows = Rows();
	int64_t cols = Cols();
	int64_t xcols = x.Cols();
	int64_t ic, c, r;
	PrepareZ(y,z,beta,opt);
	TVar val = 0.;
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			for ( c = 0; c<cols; c++) {
				for ( r = 0; r < rows; r++ ) {
					val = z(r,ic) + alpha * GetVal(r,c) * x.GetVal(c,ic);
					z.PutVal(r,ic,val);
				}
			}
		} else {
			for (r = 0; r<rows; r++) {
            	val = 0.;
				for(c = 0; c<cols; c++) {
					val += GetVal(c,r)* x.GetVal(c,ic);
				}
				z.PutVal(r,ic,alpha*val);
			}
		}
	}
}

template<class TVar>
void TPZMatrix<TVar>::Identity() {
	
    if ( Cols() != Rows() ) {
        Error( "identity (TPZMatrix<>*) <TPZMatrix<>must be square>" );
    }
    for ( int64_t row = 0; row < Rows(); row++) {
        for ( int64_t col = 0; col < Cols(); col++ ) {
            (row == col)? PutVal(row,col,1.):PutVal(row,col,0.);
        }
    }
}

template<class TVar>
void TPZMatrix<TVar>::Diagonal(TVar val) {
    
    if ( Cols() != Rows() ) {
        Error( "Diagonal (TPZMatrix<>*) <TPZMatrix<>must be square>" );
    }
    for ( int64_t row = 0; row < Rows(); row++) {
        for ( int64_t col = 0; col < Cols(); col++ ) {
            (row == col)? PutVal(row,col,val):PutVal(row,col,0.);
        }
    }
}


template<class TVar>
void TPZMatrix<TVar>::Input(std::istream& in )
{
  if constexpr (std::is_same<TVar, TFad<6, REAL>>::value ||
                std::is_same<TVar, Fad<float>>::value ||
                std::is_same<TVar, Fad<double>>::value ||
                std::is_same<TVar, Fad<long double>>::value ||
                std::is_same<TVar, TPZFlopCounter>::value) {
    PZError << __PRETTY_FUNCTION__;
    PZError << " not implemented for this type\n";
    PZError << "Aborting...";
    DebugStop();
  }
    int64_t newRow, newCol;
	in >> newRow;
	in >> newCol;
	Redim( newRow, newCol );
	int64_t i,j;
	TVar elem;
	for(i=0;i<Rows();i++)
		for(j=0;j<Cols();j++)
		{
			in >> elem;
			Put( i,j, elem );
		}
}

/** @brief Overload >> operator to input data of the matrix ***/
template<class TVar>
std::istream & operator>>(std::istream& in,TPZMatrix<TVar> &A)
{
    A.Input(in);
    return in;
}


/*************/
/*** Print ***/

template<class TVar>
void TPZMatrix<TVar>::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const {
    if constexpr (std::is_same<TVar, TFad<6, REAL>>::value ||
                  std::is_same<TVar, Fad<float>>::value ||
                  std::is_same<TVar, Fad<double>>::value ||
                  std::is_same<TVar, Fad<long double>>::value) {
        out << name << std::endl;
        for (int64_t i = 0; i < fRow; i++) {
            for (int64_t j = 0; j < fCol; j++) {
                out << "i = " << i << " j = " << j << " val " << Get(i, j) << std::endl;
            }
        }
        return;
   
    }
    //getting current settings of std::ostream instance
    std::ios cout_state(nullptr);
    cout_state.copyfmt(out);

    typedef std::numeric_limits< RTVar > typlim;
    out << std::setprecision(typlim::max_digits10);

    const auto nrows = Rows();
    const auto ncols = Cols();
    
    if(form == EFormatted) {
        out << "Writing matrix '";
        if(name) out << name;
        out << "' (" << nrows << " x " << ncols << "):\n";
		
        for ( int64_t row = 0; row < nrows; row++) {
            out << "\t";
            for ( int64_t col = 0; col < ncols; col++ ) {
                out << GetVal( row, col) << "  ";//complex numbers are (a,b)
            }
            out << "\n";
        }
        out << "\n";
    } else if (form == EFixedColumn) {
        out << "Writing matrix '";
        if(name) out << name;
        out << "' (" << nrows << " x " << ncols << "):\n";
        TVar value;
        constexpr int columnlength = typlim::max_digits10;

        auto PrintVal = [columnlength](std::ostream &out, RTVar val){
            if(val<0.) {
                out << std::setw(columnlength+1) << std::left << val << " ";
            }else{
                out << " ";
                out << std::setw(columnlength)   << std::left << val << " ";
            }
        };
        
        for ( int64_t row = 0; row < nrows; row++) {
            for ( int64_t col = 0; col < ncols; col++ ) {
                value = Get( row, col);
                if constexpr (std::is_same_v<RTVar,TVar>){
                    PrintVal(out,value);
                }else{
                    out<<"(";
                    PrintVal(out,value.real());
                    out<<",";
                    PrintVal(out,value.imag());
                    out<<")";
                }
            }
            out << "\n";
        }
        out << "\n";
    } else if (form == EInputFormat) {
        out << nrows << " " << ncols << endl;
        for ( int64_t row = 0; row < nrows; row++) {
            for ( int64_t col = 0; col < ncols; col++ ) {
                const TVar val = GetVal(row, col);
                if(val != (TVar)0.) out << row << ' ' << col << ' ' << val << std::endl;
            }
        }
        out << "-1 -1 0.\n";
    } else if( form == EMathematicaInput)
        {
            
            out << name << "\n{ ";
            for ( int64_t row = 0; row < nrows; row++) {
                out << "\n{ ";
                for ( int64_t col = 0; col < ncols; col++ ) {
                    const TVar val = GetVal(row, col);
                    if constexpr (std::is_same_v<TVar,RTVar>){
                        out << val;
                    }else{
                        out << val.real() << "+I"<<val.imag();
                    }
                    if(col < ncols-1)
                        out << ", ";
                    if((col+1) % 6 == 0) out << '\n';
                }
                out << " }";
                if(row < nrows-1)
                    out << ",";
            }
		
            out << " };"<<std::endl;
		
        }else if( form == ECSV)
        {
         
            for ( int64_t row = 0; row < nrows; row++) {
                for ( int64_t col = 0; col < ncols; col++ ) {
                    const TVar val = GetVal(row, col);
                    if constexpr (std::is_same_v<TVar,RTVar>){
                        out << val;
                    }else{
                        const auto signchar = val.imag() > 0 ? '+' : '-';
                        out << val.real() <<signchar<<std::fabs(val.imag())<<"j";
                    }
                    if(col < ncols-1)  out << " , ";
                }
                if(row < nrows-1)
                    out << "\n";
            }

        } else if( form == EMatlabNonZeros)
        {
            out << name;
            for ( int64_t row = 0; row < nrows; row++) {
                out << "\n|";
                for ( int64_t col = 0; col < ncols; col++ )
                    if(IsZero(GetVal(row, col)) ){
                        out << ".";
                    }else{
                        out << "#";
                    }
                out << "|";
            }
            out << "\n";
        }
    else if( form == EMatrixMarket)
        {
            bool sym = IsSymmetric();
            int64_t numzero = 0;
            int64_t nrow = nrows;
            for ( int64_t row = 0; row < nrows; row++) {
                int64_t colmax = nrow;
                if (sym) colmax = row+1;
                for (int64_t col = 0; col < colmax; col++ )
                    {
                        TVar val = GetVal(row, col);
                        if (val != (TVar)0.) {
                            numzero++;
                        }
                    }
            }
            out << "%%"<< name << std::endl;
            out << nrows << ' ' << ncols << ' ' << numzero << std::endl;
            for ( int64_t row = 0; row < nrows; row++) {
                int64_t colmax = nrow;
                if (sym) colmax = row+1;
                for (int64_t col = 0; col < colmax; col++ )
                    {
                        const TVar val = GetVal(row, col);
                        if (val != (TVar)0.) {
                            out << row+1 << ' ' << col+1 << ' ' << val << '\n';
                        }
                    }
            }
        }
    //restore precision
    out.copyfmt(cout_state);
}

template<class TVar>
void TPZMatrix<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &destinationindex) {
	
	int64_t nelem = elmat.Rows();
  	int64_t icoef,jcoef,ieq,jeq;
	if(IsSymmetric()) {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			for(jcoef=icoef; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				TVar prevval = GetVal(ieq,jeq);
				prevval += elmat(icoef,jcoef);
				PutVal(ieq,jeq,prevval);
			}
		}
	} else {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			for(jcoef=0; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				TVar prevval = GetVal(ieq,jeq);
				prevval += elmat(icoef,jcoef);
				PutVal(ieq,jeq,prevval);
			}
		}
	}
}

template<class TVar>
void TPZMatrix<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &source, TPZVec<int64_t> &destinationindex) {
	
	int64_t nelem = source.NElements();
  	int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    TVar prevval;
	if(IsSymmetric()) {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			ieqs = source[icoef];
			for(jcoef=icoef; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				jeqs = source[jcoef];
				prevval = GetVal(ieq,jeq);
				prevval += elmat(ieqs,jeqs);
				PutVal(ieq,jeq,prevval);
			}
		}
	} else {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			ieqs = source[icoef];
			for(jcoef=0; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				jeqs = source[jcoef];
				prevval = GetVal(ieq,jeq);
				prevval += elmat(ieqs,jeqs);
				PutVal(ieq,jeq,prevval);
			}
		}
	}
}


/***************/
/*** Put Sub ***/
//
// Put submatriz A como sub matriz a partir do ponto
//  (sRow, sCol).
//
template<class TVar>
int TPZMatrix<TVar>::PutSub(const int64_t sRow,const int64_t sCol,const TPZFMatrix<TVar> &A ) {
    if(sRow >= this->Rows()){
        LOGPZ_ERROR(logger,"incompatible dimensions")
        DebugStop();
    }

    int64_t minRow = MIN( A.Rows(), Rows() - sRow );
    int64_t minCol = MIN( A.Cols(), Cols() - sCol );
	
    int64_t row = sRow;
    for ( int64_t r = 0; r < minRow; r++, row++ ) {
        int64_t col = sCol;
        for ( int64_t c = 0; c < minCol; c++, col++ ) {
            PutVal( row, col, A.GetVal( r, c ) );
        }
    }
    return( 1 );
}


template<class TVar>
int TPZMatrix<TVar>::GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
							const int64_t colSize, TPZFMatrix<TVar> & A ) const {
    if ( ((sRow + rowSize) > Rows()) || ((sCol + colSize) > Cols()) ) {
        return( Error( "GetSub <t.he sub-matrix is too big>" ) );
    }
    A.Resize( rowSize, colSize );
    int64_t row = sRow;
    for ( int64_t r = 0; r < rowSize; r++, row++ ) {
        int64_t col = sCol;
        for ( int64_t c = 0; c < colSize; c++, col++ ) {
            A.PutVal( r, c, GetVal( row, col ) );
        }
    }
    return( 1 );
}


/***************/
/*** Add Sub ***/
//
// Adds a matriz A to a sub-matriz which starts at (row, col).
//
template<class TVar>
int TPZMatrix<TVar>::AddSub(const int64_t sRow,const int64_t sCol,const TPZFMatrix<TVar> &A ) {
	
    int64_t minRow = MIN( A.Rows(), Rows() - sRow );
    int64_t minCol = MIN( A.Cols(), Cols() - sCol );
	//  REAL v;
	
	int64_t row = sRow;
    for ( int64_t r = 0; r < minRow; r++, row++ ) {
        int64_t col = sCol;
        for ( int64_t c = 0; c < minCol; c++, col++ ) {
            PutVal( row, col, GetVal( row, col ) + A.GetVal( r, c ) );
        }
    }
    return( 1 );
}

/***************/
/*** InsertSub ***/
//
//    Inserts a submatriz with the current matrix without changing the dimensions
//
template<class TVar>
int TPZMatrix<TVar>::InsertSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
							   const int64_t colSize,const int64_t pRow,const int64_t pCol, TPZMatrix<TVar> *pA ) const {
	
	
    if ( ((pRow + rowSize) > pA->Rows()) || ((pCol + colSize) > pA->Cols())) {
        return( Error( "InsertSub <the sub-matrix is too big that target>" ) );
    }
	
    int64_t NewRowSize = rowSize+pRow;
    int64_t NewColSize = colSize+pCol;
	
	
    int64_t row = sRow;
    for ( int64_t r = pRow; r < NewRowSize; r++, row++ ) {
        int64_t col = sCol;
        for ( int64_t c = pCol ; c < NewColSize; c++, col++ ) {
            pA->PutVal( r, c, GetVal( row, col ) );
        }
    }
    return( 1 );
}


/***************/
/*** AddSub ***/
//
// Adds the submatrix to *pA at the given place
//
template<class TVar>
int TPZMatrix<TVar>::AddSub(const int64_t sRow, const int64_t sCol, const int64_t rowSize,
							const int64_t colSize,const int64_t pRow,const int64_t pCol, TPZMatrix<TVar> *pA ) const {
    if ( ((pRow + rowSize) > pA->Rows()) || ((pCol + colSize) > pA->Cols())) {
        Error( "AddSub <the sub-matrix is too big that target>" );
    }
    int64_t NewRowSize = rowSize+pRow;
    int64_t NewColSize = colSize+pCol;
	
    int64_t row = sRow;
    for ( int64_t r = pRow; r < NewRowSize; r++, row++ ) {
        int64_t col = sCol;
        for ( int64_t c = pCol ; c < NewColSize; c++, col++ ) {
            pA->PutVal( r, c, GetVal( row, col )+pA->GetVal(r,c));
        }
    }
    return( 1 );
}

/*****************/
/*** Transpose ***/
template<class TVar>
void TPZMatrix<TVar>::Transpose(TPZMatrix<TVar> *T) const {
	T->Resize( Cols(), Rows() );
	
	for ( int64_t r = 0; r < Rows(); r++ ) {
        for ( int64_t c = 0; c < Cols(); c++ ) {
            T->PutVal( c, r, GetVal( r, c ) );
        }
    }
}

/*************/
/*** Solve ***/
template<class TVar>
int TPZMatrix<TVar>::SolveDirect( TPZFMatrix<TVar> &B , DecomposeType dt, std::list<int64_t> &singular) {
	
	switch ( dt ) {
		case ELU:
			return( Solve_LU( &B ,singular)  );
		case ECholesky:
			return( Solve_Cholesky( &B , singular)  );
		case ELDLt:
			return( Solve_LDLt( &B, singular )  );
		default:
			Error( "Solve  < Unknown decomposition type >" );
			break;
	}
	return ( 0 );
}
template<class TVar>
int TPZMatrix<TVar>::SolveDirect( TPZFMatrix<TVar> &B , DecomposeType dt) {
	
	switch ( dt ) {
		case ELU:
			return( Solve_LU( &B)  );
		case ECholesky:
			return( Solve_Cholesky( &B )  );
		case ELDLt:
			return( Solve_LDLt( &B )  );
		default:
			Error( "Solve  < Unknown decomposition type >" );
			break;
	}
	return ( 0 );
}

template<class TVar>
void TPZMatrix<TVar>::SolveJacobi(int64_t &numiterations,const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
								  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch, REAL &tol,const int FromCurrent) {
	
	
	if(FromCurrent) {
		Residual(result,F,scratch);
	} else {
		scratch = F;
		result.Zero();
	}
	REAL res;
	res = TPZExtractVal::val(Norm(scratch));
	int64_t r = Dim();
	int64_t c = F.Cols();
	for(int64_t it=0; it<numiterations && (fabs(res)) > tol; it++) {
		for(int64_t ic=0; ic<c; ic++) {
			for(int64_t i=0; i<r; i++) {
				result(i,ic) += (scratch)(i,ic)/GetVal(i,i);
			}
		}
		Residual(result,F,scratch);
		res = TPZExtractVal::val(Norm(scratch));
	}
	if(residual) *residual = scratch;
}

template<>
void TPZMatrix<TFad<6,REAL> >::SolveSOR(int64_t & numiterations, const TPZFMatrix<TFad<6,REAL> > &F,
							   TPZFMatrix<TFad<6,REAL> > &result, TPZFMatrix<TFad<6,REAL> > *residual, TPZFMatrix<TFad<6,REAL> > &/*scratch*/, const REAL overrelax,
							   REAL &tol,const int FromCurrent,const int direction) {
    DebugStop();
}
template<>
void TPZMatrix<Fad<REAL> >::SolveSOR(int64_t & numiterations, const TPZFMatrix<Fad<REAL> > &F,
                                        TPZFMatrix<Fad<REAL> > &result, TPZFMatrix<Fad<REAL> > *residual, TPZFMatrix<Fad<REAL> > &/*scratch*/, const REAL overrelax,
                                        REAL &tol,const int FromCurrent,const int direction) {
  DebugStop();
}

template<class TVar>
void TPZMatrix<TVar>::SolveSOR(int64_t & numiterations, const TPZFMatrix<TVar> &F,
							   TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &/*scratch*/, const REAL overrelax,
							   REAL &tol,const int FromCurrent,const int direction) {
	
	if(residual == &F) {
		cout << "TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
		return;
	}
	TVar res = (TVar)2*(TVar)tol+(TVar)1.;
	if(residual) res = Norm(*residual);
	if(!FromCurrent) {
		result.Zero();
	}
	int64_t r = Dim();
	int64_t c = F.Cols();
	int64_t i,ifirst = 0, ilast = r, iinc = 1;
	int64_t it;
	if(direction == -1) {
		ifirst = r-1;  //misael
		ilast = -1;
		iinc = -1;
	}
	TVar eqres;
	for(it=0; it<numiterations && (REAL)(fabs(res)) > fabs(tol); it++) {
		res = 0.;
		for(int64_t ic=0; ic<c; ic++) {
			for(i=ifirst; i!=ilast; i+= iinc) {
				eqres = F.GetVal(i,ic);
				for(int64_t j=0; j<r; j++) {
					eqres -= GetVal(i,j)*result(j,ic);
				}
				res += eqres*eqres;
				result(i,ic) += (TVar)overrelax*eqres/GetVal(i,i);
			}
		}
		res = sqrt(res);
	}
	if(residual) Residual(result,F,*residual);
	numiterations = it;
	tol = fabs(res);
}

template <class TVar>
void TPZMatrix<TVar>::SolveSSOR(int64_t &numiterations, const TPZFMatrix<TVar> &F,
						  TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &scratch, const REAL overrelax,
						  REAL &tol,const int FromCurrent) {
	REAL res = tol*2.+1.;
	int64_t i, one = 1;
	int fromcurrent = FromCurrent;
	for(i=0; i<numiterations && fabs(res) > fabs(tol); i++) {
		one = 1;
		res = tol;
		SolveSOR(one,F,result,residual,scratch,overrelax,res,fromcurrent,1);
		one = 1;
		res = tol;
		fromcurrent = 1;
		SolveSOR(one,F,result,residual,scratch,overrelax,res,fromcurrent,-1);
		//                cout << "SSOR iter = " << i <<  " res = " << res << endl;
	}
	numiterations = i;
	tol = res;
}




#include "cg.h"
template <class TVar>
void TPZMatrix<TVar>::SolveCG(int64_t &numiterations, TPZSolver &preconditioner,
						const	TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
						TPZFMatrix<TVar> *residual, REAL &tol, const int FromCurrent) {
    if constexpr(std::is_floating_point<TVar>::value){
        TPZMatrixSolver<TVar> &precond = dynamic_cast<TPZMatrixSolver<TVar> &>(preconditioner);
        CG(*this, result, F, precond, residual, numiterations, tol, FromCurrent);
    }else{
        PZError<<__PRETTY_FUNCTION__<<" is currently not implemented ";
        PZError<<"for this type\nAborting...";
        DebugStop();
    }
}

#include "gmres.h"
#include "pzstepsolver.h"
template <class TVar>
void TPZMatrix<TVar>::SolveGMRES(int64_t &numiterations, TPZSolver &preconditioner,
                                 TPZFMatrix<TVar> &H, int &numvectors,
                                 const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
                                 TPZFMatrix<TVar> *residual, REAL &tol,const int FromCurrent)  
{
    if constexpr(std::is_floating_point<TVar>::value){
        TPZMatrixSolver<TVar> &precond = dynamic_cast<TPZMatrixSolver<TVar> &>(preconditioner);
      if (F.Cols() > 1) {
        int64_t locnumiter = numiterations;
        TVar loctol = tol;
        int64_t nrow = F.Rows();
        int64_t ncol = F.Cols();
        int64_t col;
        // preconditioner.Solve(F, result);
        for (col = 0; col < ncol; col++) {
          //            std::cout << "Column " << col << std::endl;
          numiterations = locnumiter;
          tol = TPZExtractVal::val(loctol);
          TPZFMatrix<TVar> FCol(nrow, 1);
          F.GetSub(0, col, F.Rows(), 1, FCol);
          TPZFMatrix<TVar> resultCol(nrow, 1, &result(0, col), nrow);
          
          GMRES(*this, resultCol, FCol, precond, H, numvectors,
                numiterations, tol, residual, FromCurrent);
        }
      } else {
        GMRES(*this, result, F, precond, H, numvectors, numiterations,
              tol, residual, FromCurrent);
      }
    }else{
        PZError<<__PRETTY_FUNCTION__<<" is currently not implemented ";
        PZError<<"for this type\nAborting...";
        DebugStop();
    }
}
#include "bicg.h"
template<class TVar>
void TPZMatrix<TVar>::SolveBICG(int64_t &numiterations, TPZSolver &preconditioner,
						  const TPZFMatrix<TVar> &F,
						  TPZFMatrix<TVar> &result,
						  REAL &tol)  {
    if constexpr(std::is_floating_point<TVar>::value){
        TPZMatrixSolver<TVar> &precond = dynamic_cast<TPZMatrixSolver<TVar> &>(preconditioner);
        BiCG (*this, result,F,precond,numiterations,tol);
    }else{
        PZError<<__PRETTY_FUNCTION__<<" is currently not implemented ";
        PZError<<"for this type\nAborting...";
        DebugStop();
    }
}

#include "bicgstab.h"
template <class TVar>
void TPZMatrix<TVar>::SolveBICGStab(int64_t &numiterations, TPZSolver &preconditioner,
							  const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
							  TPZFMatrix<TVar> *residual, REAL &tol,const int FromCurrent)  {
    if constexpr(std::is_floating_point<TVar>::value){
        TPZMatrixSolver<TVar> &precond = dynamic_cast<TPZMatrixSolver<TVar> &>(preconditioner);
        BiCGSTAB(*this,result,F,precond,numiterations,tol,residual,FromCurrent);
    }else{
        PZError<<__PRETTY_FUNCTION__<<" is currently not implemented ";
        PZError<<"for this type\nAborting...";
        DebugStop();
    }    
}

#include "ir.h"
template <class TVar>
void TPZMatrix<TVar>::SolveIR(int64_t &numiterations, TPZSolver &preconditioner,
						const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
						TPZFMatrix<TVar> *residual, REAL &tol,
						const int FromCurrent)  {
    if constexpr(std::is_floating_point<TVar>::value){
        TPZMatrixSolver<TVar> &precond = dynamic_cast<TPZMatrixSolver<TVar> &>(preconditioner);
        IR(*this,result,F,precond,residual,numiterations,tol,FromCurrent);
    }else{
        PZError<<__PRETTY_FUNCTION__<<" is currently not implemented ";
        PZError<<"for this type\nAborting...";
        DebugStop();
    }
}

/*******************************************************************************/

/*****************/
/*** Decompose_LU ***/
template <class TVar>
int TPZMatrix<TVar>::Decompose_LU(std::list<int64_t> &singular) {
	return Decompose_LU();
}
template <class TVar>
int TPZMatrix<TVar>::Decompose_LU() {
	
	if (  fDecomposed && fDecomposed != ELU)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
	if (fDecomposed) return 1;
	
	TVar nn, pivot;
	int64_t  min = ( Cols() < (Rows()) ) ? Cols() : Rows();
	
	for ( int64_t k = 0; k < min ; k++ ) {
		if (IsZero( pivot = GetVal(k, k))) Error( "Decompose_LU <matrix is singular>" );
		for ( int64_t i = k+1; i < Rows(); i++ ) {
			nn = GetVal( i, k ) / pivot;
			PutVal( i, k, nn );
			for ( int64_t j = k+1; j < Cols(); j++ ) PutVal(i,j,GetVal(i,j)-nn*GetVal(k,j));
		}
	}
	fDecomposed=ELU;
	return 1;
}



/****************/
/*** Substitution ***/
template <class TVar>
int TPZMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const{
	
    int64_t rowb = B->Rows();
    int64_t colb = B->Cols();
    if ( rowb != Rows() )
		Error( "SubstitutionLU <incompatible dimensions>" );
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
				Error( "BackSub( SubstitutionLU ) <Matrix is singular" );
            }
            B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
    }
    return( 1 );
}
template <class TVar>
int TPZMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular) {
	return Decompose_LDLt();
}
template <class TVar>
int TPZMatrix<TVar>::Decompose_LDLt() {
	
	if (  fDecomposed && fDecomposed != ELDLt) {
		Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
	} else if(fDecomposed ) {
		return ELDLt;
	}
	if ( Rows()!=Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
	
	const int dim=Rows();
	
	for (auto j = 0; j < dim; j++ ) {
        TVar sum = 0;
		for (auto k=0; k<j; k++) {
            if constexpr(is_complex<TVar>::value){
                sum +=GetVal(k,k)*std::conj(GetVal(k,j))*GetVal(k,j);
            }else{
                sum +=GetVal(k,k)*GetVal(k,j)*GetVal(k,j);
            }   
		}
        PutVal( j,j,GetVal(j,j) - sum );
		TVar tmp = GetVal(j,j);
		if ( IsZero(tmp) ) Error( "Decompose_LDLt <Zero on diagonal>" );
        for(auto l=j+1; l<dim;l++) {
            TVar sum = 0;
            for (auto k=0; k<j; k++) {
                if constexpr(is_complex<TVar>::value){
                    sum += GetVal(k,k)*std::conj(GetVal(j,k))*GetVal(l,k);
                }else{
                    sum += GetVal(k,k)*GetVal(j,k)*GetVal(l,k);
                }
			}
            TVar val = (GetVal(l,j) - sum)/tmp;
            PutVal(l,j,val);
            if constexpr(is_complex<TVar>::value){val=std::conj(val);}
            PutVal(j,l,val);
		}
	}
	fDecomposed  = ELDLt;
	fDefPositive = 0;
	return( 1 );
}
template <class TVar>
int TPZMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular) {
	if (  fDecomposed && fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
	if (  fDecomposed ) return ECholesky;
	if ( Rows()!=Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
	//return 0;
	
	int64_t dim=Dim();
	for (int64_t i=0 ; i<dim; i++) {
		for(int64_t k=0; k<i; k++) {             //elementos da diagonal
			PutVal( i,i,GetVal(i,i)-GetVal(i,k)*GetVal(i,k) );
		}
		if((fabs(GetVal(i,i))) <= fabs((TVar)1.e-12))
		{
			singular.push_back(i);
			PutVal(i,i,1.);
		}
		TVar tmp = sqrt(GetVal(i,i));
        PutVal( i,i,tmp );
        for (int64_t j=i+1;j<dim; j++) {           //elementos fora da diagonal
            for(int64_t k=0; k<i; k++) {
                PutVal( i,j,GetVal(i,j)-GetVal(i,k)*GetVal(k,j) );
            }
            TVar tmp2 = GetVal(i,i);
            if ( IsZero(tmp2) ) {
				Error( "Decompose_Cholesky <Zero on diagonal>" );
            }
            PutVal(i,j,GetVal(i,j)/GetVal(i,i) );
            PutVal(j,i,GetVal(i,j));
			
        }
    }
	fDecomposed = ECholesky;
	return ECholesky;
}
template <class TVar>
int TPZMatrix<TVar>::Decompose_Cholesky() {
	if (  fDecomposed && fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
	if (  fDecomposed ) return ECholesky;
	if ( Rows()!=Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
	//return 0;
	
	int64_t dim=Dim();
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
				Error( "Decompose_Cholesky <Zero on diagonal>" );
            }
            PutVal(i,j,GetVal(i,j)/GetVal(i,i) );
            if constexpr (is_complex<TVar>::value){
                PutVal(j,i,std::conj(GetVal(i,j)));
            }else{
                PutVal(j,i,GetVal(i,j));
            }
        }
    }
	fDecomposed = ECholesky;
	return ECholesky;
	
}

template <class TVar>
int TPZMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar> *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ECholesky)
		return( 0 );
	for ( int64_t r = 0; r < Dim(); r++ ) {
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

/**********************/
/*** Subst Backward ***/
//
//  Faz Ax = b, onde A(NxN) e' triangular superior.
//
template <class TVar>
int TPZMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar> *B ) const {
	
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ECholesky) return( 0 );
	for ( int64_t r = Dim()-1;  r >= 0;  r-- ) {
		TVar pivot = GetVal( r, r );
		for ( int64_t c = 0; c < B->Cols(); c++ ) {
			// Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
			//
			TVar sum = 0.0;
			for ( int64_t i = Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
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
//  Faz Ax = b, onde A e' triangular inferior e A(i,i) = 1.
//
template <class TVar>
int TPZMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar> *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ELDLt) {
		Error("TPZMatrix::Subst_LForward incompatible dimensions\n");
	}
    for ( int64_t r = 0; r < Dim(); r++ ) {
        for ( int64_t c = 0; c < B->Cols();  c++ )    {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
            //
            TVar sum = 0.0;
            for ( int64_t i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
			
            // Faz B[r,c] = (B[r,c] - sum) / A[r,r].
            //
            B->PutVal( r, c, B->GetVal(r, c) - sum );
        }
    }
    return( 1 );
}

/************************/
/*** Subst L Backward ***/
//
//  Faz Ax = b, onde A e' triangular superior e A(i,i) = 1.
//
template <class TVar>
int TPZMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar> *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ELDLt){
		Error("TPZMatrix::Subst_LBackward incompatible dimensions \n");
	}
	
    for ( int64_t r = Dim()-1;  r >= 0;  r-- ) {
        for ( int64_t c = 0; c < B->Cols(); c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
            //
            TVar sum = 0.0;
            for ( int64_t i = Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
			
            // Faz B[r,c] = B[r,c] - sum.
            //
            B->PutVal( r, c, B->GetVal(r, c) - sum );
        }
    }
    return( 1 );
}

/******************/
/*** Subst Diag ***/
//
//  Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//
template <class TVar>
int TPZMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar> *B ) const {
    if ( (B->Rows() != Dim())) {
        Error("TPZMatrix::Subst_Diag incompatible dimensions\n");
	}
	for ( int64_t r = 0; r < Dim(); r++ ) {
		TVar pivot = GetVal( r, r );
		for ( int64_t c = 0; c < B->Cols(); c++ ) {
			B->PutVal( r, c, B->GetVal( r, c ) / pivot );
		}
	}
	return( 1 );
}

/************************** Private **************************/

/*************/
/*** Error ***/
template <class TVar>
int TPZMatrix<TVar>::Error(const char *msg ,const char *msg2) {
  ostringstream out;
  out << "TPZMatrix::" << msg;
  if(msg2) out << msg2;
  out << ".\n";
  std::cout<<out.str()<<std::endl;
  LOGPZ_ERROR (logger, out.str().c_str());
  DebugStop();
	std::bad_exception myex;
	throw myex;
}

template<class TVar>
int TPZMatrix<TVar>::ClassId() const{
    return Hash("TPZMatrix") ^ ClassIdOrHash<TVar>()<<1 ^ TPZBaseMatrix::ClassId();
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template <class TVar>
bool TPZMatrix<TVar>::Compare(TPZSavable *copy, bool override)
{
	TPZMatrix<TVar> *copmat = dynamic_cast<TPZMatrix<TVar> *> (copy);
	if(!copmat) return false;
	bool result = true;
	if(fRow != copmat->fRow || fCol != copmat->fCol || fDecomposed != copmat->fDecomposed) result = false;
	if(!result)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " did not compare ";
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(override && !result)
	{
		this->operator=(*copmat);
	}
	return result;
}

/// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
template <class TVar>
bool TPZMatrix<TVar>::Compare(TPZSavable *copy, bool override) const
{
	TPZMatrix<TVar> *copmat = dynamic_cast<TPZMatrix<TVar> *> (copy);
	if(!copmat) return false;
	bool result = true;
	if(fRow != copmat->fRow || fCol != copmat->fCol || fDecomposed != copmat->fDecomposed) result = false;
	if(!result)
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " did not compare ";
		LOGPZ_ERROR(loggerCheck,sout.str())
	}
	if(override && !result)
	{
		DebugStop();
	}
	return result;
}

template <class TVar>
void TPZMatrix<TVar>::GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar> &block) const
{
    int64_t nel = indices.NElements();
    int64_t iel,jel;
    for(iel=0; iel<nel; iel++)
    {
		for(jel=0; jel<nel; jel++)
		{
			block(iel,jel) = GetVal(indices[iel],indices[jel]);
		}
    }
	
}

template<>
int TPZMatrix<TFad<6,REAL> >::VerifySymmetry(REAL tol) const{
    DebugStop();
    return -1;
}

template <class TVar>
int TPZMatrix<TVar>::VerifySymmetry(REAL tol) const{
	int64_t nrows = this->Rows();
	int64_t ncols = this->Cols();
	if (nrows != ncols) return 0;
	
	for( int64_t i = 0; i < nrows; i++){
		for(int64_t j = 0; j <= i; j++){
            TVar exp = this->Get(i,j) - this->Get(j,i);
			if ( (REAL)(fabs( exp )) > tol ) {
			  	#ifdef STATE_COMPLEX
				cout << "Elemento: " << i << ", " << j << "  -> " << fabs( exp ) << "/" <<
				this->Get(i,j) << endl;
				#else
				cout << "Elemento: " << i << ", " << j << "  -> " << exp << "/" <<
				this->Get(i,j) << endl;
				#endif
				return 0;
			}
		}
	}
	return 1;
}

template<>
bool TPZMatrix<TFad<6,REAL> >::CompareValues(TPZMatrix<TFad<6,REAL> > &M, TFad<6,REAL> tol){
    DebugStop();
    return false;
}
template <class TVar>
bool TPZMatrix<TVar>::CompareValues(TPZMatrix<TVar> &M, TVar tol){
	
	int64_t nrows = this->Rows();
	int64_t ncols = this->Cols();
	if ( (M.Rows() != nrows) || (M.Cols() != ncols) ) return false;
	
	REAL diff;
	for( int64_t i = 0; i < nrows; i++)
		for( int64_t j = 0; j < ncols; j++){
            TVar exps = this->Get(i,j) - M.Get(i,j);
			diff = fabs( exps );
			if (diff > fabs(tol)) return false;
		}
	
	return true;
}

template<>
TFad<6,REAL> TPZMatrix<TFad<6,REAL> >::ReturnNearestValue(TFad<6,REAL> val, TPZVec<TFad<6,REAL> >& Vec, TFad<6,REAL> tol)
{
    DebugStop();
    TFad<6,REAL> res;
    return res;
}

template <class TVar>
TVar TPZMatrix<TVar>::ReturnNearestValue(TVar val, TPZVec<TVar>& Vec, TVar tol)
{
    TVar exps = val - (TVar)Vec[0];
    TVar diff0 = fabs(exps) >= fabs(tol) ?  (val - (TVar)Vec[0]) : (TVar)1.E10;
    TVar diff1, res = (TVar)Vec[0];
    for(int64_t i = 1; i < Vec.NElements(); i++)
    {
        diff1 = val - (TVar)Vec[i];
        diff0 = ( fabs(diff1) >= fabs(tol) && fabs(diff1) < fabs(diff0) ) ? res = Vec[i],diff1 : diff0;
    }
    return res;
}

template<>
bool TPZMatrix<TFad<6,REAL> >::SolveEigensystemJacobi(int64_t &numiterations, REAL & tol, TPZVec<TFad<6,REAL> > & Eigenvalues, TPZFMatrix<TFad<6,REAL> > & Eigenvectors) const{
    DebugStop();
    return false;
}


template <class TVar>
bool TPZMatrix<TVar>::SolveEigensystemJacobi(int64_t &numiterations, REAL & tol, TPZVec<TVar> & Eigenvalues, TPZFMatrix<TVar> & Eigenvectors) const{
	
	int64_t NumIt = numiterations;
	REAL tolerance = tol;
	
#ifdef PZDEBUG2
	if (this->Rows() != this->Cols())
	{
		PZError << __PRETTY_FUNCTION__ <<
		" - Jacobi method of computing eigensystem requires a symmetric square matrix. this->Rows = " << this->Rows() << " - this->Cols() = " << this->Cols() << endl;
		return false;
	}
	
	if (this->VerifySymmetry(1.e-8) == false)
	{
		PZError << __PRETTY_FUNCTION__ <<
		" - Jacobi method of computing eigensystem requires a symmetric square matrix. This matrix is not symmetric." << endl;
		return false;
	}
#endif
	
	const int64_t size = this->Rows();
	
	/** Making a copy of this */
	TPZFNMatrix<9,TVar> Matrix(size,size); //fast constructor in case of this is a stress or strain tensor.
	for(int64_t i = 0; i < size; i++) for(int64_t j = 0; j < size; j++) Matrix(i,j) = this->Get(i,j);
	
	/** Compute Eigenvalues *//////////////////////////////////////
	bool result = Matrix.SolveEigenvaluesJacobi(numiterations, tol, &Eigenvalues);
	if (result == false) return false;
	
	/** Compute Eigenvectors *//////////////////////////////////////
	TPZFNMatrix<3, TVar> VecIni(size,1,0.), VecIni_cp(size,1,0.);
	
	Eigenvectors.Resize(size, size);
	Eigenvectors.Zero();
	for(int64_t eigen = 0; eigen < size; eigen++)
	{
        for(int64_t i = 0; i < size; i++) VecIni.PutVal(i,0,rand());
		
        TPZFNMatrix<9,TVar> Matrix(*this);
		
        TVar answ = ReturnNearestValue(Eigenvalues[eigen], Eigenvalues,((TVar)1.E-5));
        TVar exp = answ - Eigenvalues[eigen];
        if((REAL)(fabs(exp)) > 1.E-5)
        {
            for(int64_t i = 0; i < size; i++) Matrix.PutVal(i,i, this->GetVal(i,i) - (TVar)(Eigenvalues[eigen] - (TVar)(0.01 * fabs(exp))) );
        }
        else
        {
            for(int64_t i = 0; i < size; i++) Matrix.PutVal(i,i, this->GetVal(i,i) - (Eigenvalues[eigen] - TVar(0.01)) );
        }
		
        /** Normalizing Initial Eigenvec */
        REAL norm1 = 0.;
        for(int64_t i = 0; i < size; i++) norm1 += fabs(VecIni(i,0)) * fabs(VecIni(i,0)); norm1 = sqrt(norm1);
        for(int64_t i = 0; i < size; i++) VecIni(i,0) = VecIni(i,0)/(TVar)norm1;
		
        int64_t count = 0;
        double dif = 10., difTemp = 0.;
		
        while(dif > tolerance && count <= NumIt)
        {
			for(int64_t i = 0; i < size; i++) VecIni_cp(i,0) = VecIni(i,0);
			
			/** Estimating Eigenvec */
			Matrix.Solve_LU(&VecIni);
			
			/** Normalizing Final Eigenvec */
			REAL norm2 = 0.;
			for(int64_t i = 0; i < size; i++) norm2 += fabs(VecIni(i,0)) * fabs(VecIni(i,0));
			norm2 = sqrt(norm2);
			
			difTemp = 0.;
			for(int64_t i = 0; i < size; i++)
			{
                VecIni(i,0) = VecIni(i,0)/(TVar)norm2;
                TVar exp = (VecIni_cp(i,0) - VecIni(i,0)) * (VecIni_cp(i,0) - VecIni(i,0));
                difTemp += fabs(exp);
			}
			dif = sqrt(difTemp);
			count++;
        }
		
        /** Copy values from AuxVector to Eigenvectors */
        for(int64_t i = 0; i < size; i++)
        {
			TVar val = VecIni(i,0);
			if((REAL)(fabs(val)) < 1.E-5) val = 0.;
			Eigenvectors(eigen,i) = val;
        }
		
#ifdef PZDEBUG2
        double norm = 0.;
        for(int64_t i = 0; i < size; i++) norm += fabs(VecIni(i,0)) * fabs(VecIni(i,0));
        if (fabs(norm - 1.) > 1.e-10)
        {
            PZError << __PRETTY_FUNCTION__ << endl;
        }
        if(count > NumIt-1)// O metodo nao convergiu !!!
        {
            PZError << __PRETTY_FUNCTION__ << endl;
#ifdef PZ_LOG
            {
                std::stringstream sout;
                Print("Matrix for SolveEigensystemJacobi did not converge",sout);
                LOGPZ_DEBUG(logger,sout.str());
            }
#endif
        }
#endif
	}
	
	return true;
	
}//method

template <>
bool TPZMatrix< std::complex< float > >::SolveEigenvaluesJacobi(int64_t &numiterations, REAL & tol, TPZVec< std::complex< float > > * Sort){
  DebugStop(); // Does not work with complex numbers. To be implemented in the future.
    return false;
}

template <>
bool TPZMatrix< std::complex< double > >::SolveEigenvaluesJacobi(int64_t &numiterations, REAL & tol, TPZVec< std::complex< double > > * Sort){
  DebugStop(); // Does not work with complex numbers. To be implemented in the future.
    return false;
}

template <>
bool TPZMatrix< std::complex< long double > >::SolveEigenvaluesJacobi(int64_t &numiterations, REAL & tol, TPZVec< std::complex< long double > > * Sort){
  DebugStop(); // Does not work with complex numbers. To be implemented in the future.
    return false;
}

template<>
bool TPZMatrix<TFad<6,REAL> >::SolveEigenvaluesJacobi(int64_t &numiterations, REAL & tol, TPZVec<TFad<6,REAL> > * Sort){
    DebugStop();
    return false;
}

template <class TVar>
bool TPZMatrix<TVar>::SolveEigenvaluesJacobi(int64_t &numiterations, REAL & tol, TPZVec<TVar> * Sort){
	
#ifdef PZDEBUG2
	if (this->Rows() != this->Cols()){
		PZError << __PRETTY_FUNCTION__ << " - Jacobi method of computing eigenvalues requires a symmetric square matrix. this->Rows = " << this->Rows() << " - this->Cols() = " << this->Cols() << endl;
		return false;
	}
	
	if (this->VerifySymmetry(1.e-8) == false){
		PZError << __PRETTY_FUNCTION__ << " - Jacobi method of computing eigenvalues requires a symmetric square matrix. This matrix is not symmetric." << endl;
		return false;
	}
#endif
	
	int64_t iter = 0;
	TVar res = (TVar)2. * (TVar)tol;
	const int64_t size = this->Rows();
	int64_t i, j;
	int64_t p = -1, q = -1;
	TVar maxval, theta, cost, sint, aux, aux2, Spp, Sqq, Spq;
	while (iter < numiterations){
		/** First of all find the max value off diagonal */
		maxval = 0.;
		for(i = 0; i < size; i++){
			for(j = 0; j < i; j++) {
				if( fabs(this->operator ( )(i,j) ) > fabs(maxval) ) {
					p = i;
					q = j;
					maxval = fabs( this->operator ( )(i,j) );
				}//if
			}//for j
		}//for i
		
		/** Check if max value off diagonal is lesser than required tolerance */
		res = maxval;
		if ((REAL)(fabs(res)) < tol) break;
		
		/** Compute angle of rotation */
		theta = (TVar)0.5 * atan((TVar)2. * this->operator ( )(p,q) / (this->operator ( )(q,q) - this->operator ( )(p,p) ) );
		cost = cos(theta);
		sint = sin(theta);
		
		/** Apply rotation */
		for(i = 0; i < size; i++){
			if (i != p && i != q){
				
				aux = this->operator ( )(i,p) * (TVar)cost - this->operator ( )(i,q) * (TVar)sint;
				aux2 = this->operator ( )(i,p) * (TVar)sint + this->operator ( )(i,q) * (TVar)cost;
				
				this->operator ( )(i,p) = aux;
				this->operator ( )(p,i) = aux;
				
				this->operator ( )(i,q) = aux2;
				this->operator ( )(q,i) = aux2;
				
			}//if
		}//for i
		
		Spp = this->operator ( )(p,p) * cost * cost -(TVar)2. * this->operator ( )(p,q) * sint * cost + this->operator ( )(q,q) * sint * sint;
		Sqq = this->operator ( )(p,p) * sint * sint +(TVar)2. * this->operator ( )(p,q) * sint * cost + this->operator ( )(q,q) * cost * cost;
		Spq = ( this->operator ( )(p,p) - this->operator ( )(q,q) ) * cost * sint + this->operator ( )(p,q)*( cost*cost - sint*sint );
		
		this->operator ( )(p,p) = Spp;
		this->operator ( )(q,q) = Sqq;
		this->operator ( )(p,q) = Spq;
		this->operator ( )(q,p) = Spq;
		
		iter++;
	}//while
	
	/** Sorting */
	if (Sort){
		multiset< REAL > myset;
		for(i = 0; i < size; i++)
        {   TVar exps = this->operator ( )(i,i);
            myset.insert( (TPZExtractVal::val(exps)) );
        }
		
#ifdef PZDEBUG2
		if ((int64_t)myset.size() != size) PZError << __PRETTY_FUNCTION__ << " - ERROR!" << endl;
#endif
		
		Sort->Resize(size);
		multiset< REAL >::iterator w, e = myset.end();
		for(i = size - 1, w = myset.begin(); w != e; w++, i--){
			Sort->operator [ ](i) = *w;
		}//for
	}//if (Sort)
	
	
	if ((REAL)(fabs(res)) < tol){
		tol = fabs(res);
		numiterations = iter;
		return true;
	}
	
	tol = fabs(res);
	numiterations = iter;
	return false;
	
}//method

template<>
TFad<6,REAL> TPZMatrix<TFad<6,REAL> >::MatrixNorm(int p, int64_t numiter, REAL tol) const{
    DebugStop();
    TFad<6,REAL> res;
    return res;
}

template <class TVar>
TVar TPZMatrix<TVar>::MatrixNorm(int p, int64_t numiter, REAL tol) const{
	const int64_t n = this->Rows();
	if (!n) return 0.;
	if (n != this->Cols()){
		PZError << __PRETTY_FUNCTION__
		<< " - matrix must be square - Rows() = "
		<< this->Rows() << " - Cols() = "
		<< this->Cols() << std::endl;
	}
	switch(p){
		case 0:{
			TVar max = 0.;
			int64_t i, j;
			for(i = 0; i < n; i++){
				REAL sum = 0.;
				for(j = 0; j < n; j++) sum += fabs( this->Get(i,j) );
				if (sum > fabs(max)) max = sum;
			}
			return max;
		}
		case 1:{
			TVar max = 0.;
			int64_t i, j;
			for(i = 0; i < n; i++){
				REAL sum = 0.;
				for(j = 0; j < n; j++) sum += fabs( this->Get(j,i) );
				if (sum > fabs(max)) max = sum;
			}
			return max;
		}
		case 2:{
			TPZFMatrix<TVar> transp(n,n);
			int64_t i, j, k;
			//TRANSPOSE
			for(i = 0; i < n; i++) for(j = 0; j < n; j++) transp(i,j) = this->Get(j,i);
			//MULTIPLY transp = transp.this
			TPZVec<TVar> ROW(n);
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++) ROW[j] = transp(i,j);
				for(j = 0; j < n; j++){
					TVar sum = 0.;
					for(k = 0; k < n; k++){
						sum += ROW[k] * this->Get(k,j);
					}//for k
					transp(i,j) = sum;
				}//for j
			}//for i
			
			TPZVec<TVar> EigenVal;
			bool result = transp.SolveEigenvaluesJacobi(numiter, tol, &EigenVal);
			if (result == false) PZError << __PRETTY_FUNCTION__
				<< " - it was not possible to find Eigenvalues. Iterations = "
				<< numiter << " - error found = " << tol << std::endl;
			return sqrt(EigenVal[0]);
		}//case 2
		default:{
			PZError << __PRETTY_FUNCTION__ << " p = " << p << " is not a correct option" << std::endl;
		}
	}//switch
	return 0.;
}//method

template <class TVar>
TVar TPZMatrix<TVar>::ConditionNumber(int p, int64_t numiter, REAL tol){
	int64_t localnumiter = numiter;
	REAL localtol = tol;
	TPZFMatrix<TVar> Inv;
	TVar thisnorm = this->MatrixNorm(p, localnumiter, localtol);
	if (!this->Inverse(Inv,ENoDecompose)){
		PZError << __PRETTY_FUNCTION__ << " - it was not possible to compute the inverse matrix." << std:: endl;
		return 0.;
	}
	TVar invnorm = Inv.MatrixNorm(p, numiter, tol);
	return thisnorm * invnorm;
}

template<class TVar>
void TPZMatrix<TVar>::Add(const TPZMatrix<TVar>&A,TPZMatrix<TVar>&res) const {
	if ((Rows() != A.Rows()) || (Cols() != A.Cols()) ) {
		Error( "Add(TPZMatrix<>&, TPZMatrix) <different dimensions>" );
	}
  if(this->Size()!=A.Size()){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"ERROR\nDifferent number of entries.\nAborting...\n";
    DebugStop();
  }
  res.CopyFrom(this);
  auto thisFmat =
    dynamic_cast<const TPZFMatrix<TVar>*>(this);
  auto aFmat =
    dynamic_cast<const TPZFMatrix<TVar>*>(&A);
  auto resFmat =
    dynamic_cast<TPZFMatrix<TVar>*>(&res);
  if(resFmat && (!thisFmat || !aFmat)){
    for ( int64_t r = 0; r < Rows(); r++ )
      for ( int64_t c = 0; c < Cols(); c++ )
        res.PutVal( r, c, GetVal(r,c)+A.GetVal(r,c) );
  }else{
    CheckTypeCompatibility(&res,&A);
    res.Storage() += A.Storage();
  }
  
}
template<class TVar>
void TPZMatrix<TVar>::Subtract(const TPZMatrix<TVar> &A,TPZMatrix<TVar> &res) const {
	
	if ((Rows() != A.Rows()) || (Cols() != A.Cols()) ) {
		Error( "Add(TPZMatrix<>&, TPZMatrix<>&) <different dimensions>" );
	}
	if(this->Size()!=A.Size()){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"ERROR\nDifferent number of entries.\nAborting...\n";
    DebugStop();
  }
  res.CopyFrom(this);
  auto thisFmat =
    dynamic_cast<const TPZFMatrix<TVar>*>(this);
  auto aFmat =
    dynamic_cast<const TPZFMatrix<TVar>*>(&A);
  auto resFmat =
    dynamic_cast<TPZFMatrix<TVar>*>(&res);
  if(resFmat && (!thisFmat || !aFmat)){
    for ( int64_t r = 0; r < Rows(); r++ )
      for ( int64_t c = 0; c < Cols(); c++ )
        res.PutVal( r, c, GetVal(r,c)-A.GetVal(r,c) );
  }else{
    CheckTypeCompatibility(&res,&A);
    res.Storage() -= A.Storage();
  }
}

template<class TVar>
void TPZMatrix<TVar>::MultiplyByScalar(const TVar alpha, TPZMatrix<TVar>&res) const
{
  res.CopyFrom(this);
  res.Storage() *= alpha;
}

template<class TVar>
TPZMatrix<TVar> &TPZMatrix<TVar>::operator*=(const TVar val)
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"\nERROR: not implemented.\nAborting...\n";
  DebugStop();
  return *this;
}

template<class TVar>
void TPZMatrix<TVar>::Multiply(const TPZFMatrix<TVar> &A, TPZFMatrix<TVar>&B, int opt) const {
  
	
	if(!opt && (B.Rows() != Rows() || B.Cols() != A.Cols())) {
		B.Redim(Rows(),A.Cols());
	}
	else if (opt && (B.Rows() != Cols() || B.Cols() != A.Cols())) {
		B.Redim(Cols(),A.Cols());
	}else if(opt && this->fDecomposed){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Cannot multiply with opt when matrix is decomposed"
           <<"\nAborting...\n";
    DebugStop();
  }
  
  switch(this->fDecomposed){
  case ENoDecompose:
    MultAdd( A, B, B, 1.0, 0.0, opt);
    break;
  case ELU:
    B=A;
    Substitution(&B);
    break;
  case ECholesky:
    B=A;
    Subst_Forward(&B);
    Subst_Backward(&B);
    break;
  case ELDLt:
    B=A;
    Subst_LForward(&B);
    Subst_Diag(&B);
    Subst_LBackward(&B);
    break;
  default:
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: Cannot multiply with fDecomposed:"<<fDecomposed
           <<"\nAborting...\n";
    DebugStop();
  }
	
}

template<class TVar> 
int TPZMatrix<TVar>::Inverse(TPZFMatrix<TVar>&Inv, DecomposeType dec){
	const int64_t n = this->Rows();
	if (!n) return 0;
	if (n != this->Cols()){
		PZError << __PRETTY_FUNCTION__
		<< " - matrix must be square - Rows() = "
		<< this->Rows() << " - Cols() = "
		<< this->Cols() << std::endl;
		return 0;
	}
	
	Inv.Redim(n,n);
	Inv.Zero();
	int64_t i;
	for(i = 0; i < n; i++) Inv(i,i) = 1.;
	
    if(dec != ENoDecompose)
    {
        this->SolveDirect(Inv,dec);
    }
    else
    {
        const int issimetric = this->IsSymmetric();
        if (issimetric)  return this->SolveDirect(Inv, ELDLt);
        if (!issimetric) return this->SolveDirect(Inv, ELU);
    }
	return 0;
}//method

/** Fill the matrix with random values (non singular matrix) */
template <class TVar>
void TPZMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric) {
    Resize(nrow,ncol);
	int64_t i, j;
	TVar val, sum;
	/** Fill data */
	for(i=0;i<Rows();i++) {
		sum = 0.0;
        j=0;
        if (symmetric) {
            for (; j<i; j++) {
                if constexpr (is_complex<TVar>::value){
                    //hermitian matrices
                    PutVal(i, j, std::conj(GetVal(j,i)));
                }else{
                    PutVal(i, j, GetVal(j,i));
                }
                sum += fabs(GetVal(i, j));
            }
        }
		for(;j<Cols();j++) {
			val = GetRandomVal();
            if constexpr(is_complex<TVar>::value){
                if(j==i) val = fabs(val);
            }
            if(i!=j) sum += fabs(val);
			if(!PutVal(i,j,val))
				Error("AutoFill (TPZMatrix) failed.");
			
		}
		if (Rows() == Cols()) {
		  /** Making diagonally dominant and non zero in diagonal */
		  if(fabs(sum) > fabs(GetVal(i,i)))            // Deve satisfazer:  |Aii| > SUM( |Aij| )  sobre j != i
			  PutVal(i,i,sum+(TVar)1.);
		  // To sure diagonal is not zero.
		  if(IsZero(sum) && IsZero(GetVal(i,i)))
			  PutVal(i,i,1.);
		}
	}
}

template<class TVar>
int TPZMatrix<TVar>::Solve_LDLt( TPZFMatrix<TVar>* B, std::list<int64_t> &singular ) {
    
    int result = Decompose_LDLt(singular);
    if (result == 0) {
        return result;
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B->Print("On input " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Subst_LForward( B );
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B->Print("Only forward " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Subst_Diag( B );
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B->Print("After forward and diagonal " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    result = Subst_LBackward( B );
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B->Print("Final result " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return result;
}

/** @brief Overload << operator to output entries of the matrix ***/
template<class TVar>
std::ostream &operator<<(std::ostream& out,const TPZMatrix<TVar> &A) {
    A.Print("operator << ",out);
    return  out;
}


template<class TVar>
TVar TPZMatrix<TVar>::GetRandomVal() const{

    if constexpr (std::is_integral_v<TVar> ||
                  std::is_same_v<TVar,TPZFlopCounter>){
        return  ((TVar) rand())/((TVar)RAND_MAX);
    }
    
    else if constexpr (is_fad<TVar>::value){
      constexpr typename TVar::value_type lower_bound = 0;
      constexpr typename TVar::value_type upper_bound = 1;
      static std::uniform_real_distribution<typename TVar::value_type>
        unif(lower_bound,upper_bound);
      static std::default_random_engine re;
      TVar a_random = (TVar)unif(re);
      if constexpr (is_complex<TVar>::value){
        a_random+= (TVar)1i * (TVar)unif(re);
      }
      return a_random;
    }else{
        constexpr RTVar lower_bound = 0;
        constexpr RTVar upper_bound = 1;
        static std::uniform_real_distribution<RTVar> unif(lower_bound,upper_bound);
        static std::default_random_engine re;
        TVar a_random = (TVar)unif(re);
        if constexpr (is_complex<TVar>::value){
            a_random+= (TVar)1i * (TVar)unif(re);
        }
        return a_random;
    }
}

template class TPZMatrix< std::complex<float> >;
template class TPZMatrix< std::complex<double> >;
template class TPZMatrix< std::complex<long double> >;

template class TPZMatrix<long double>;
template class TPZMatrix<double>;
template class TPZMatrix<int64_t>;
template class TPZMatrix<int>;
template class TPZMatrix<float>;
template class TPZMatrix<TPZFlopCounter>;

template class TPZMatrix<TFad<6,REAL> >;
template class TPZMatrix<Fad<float> >;
template class TPZMatrix<Fad<double> >;
template class TPZMatrix<Fad<long double> >;



template std::ostream &operator<<(std::ostream &out, const TPZMatrix<int> &A);
template std::ostream &operator<<(std::ostream &out, const TPZMatrix<float> &A);
template std::ostream &operator<<(std::ostream &out, const TPZMatrix<double> &A);
template std::ostream &operator<<(std::ostream &out, const TPZMatrix<long double> &A);
template std::ostream &operator<<(std::ostream &out, const TPZMatrix<std::complex<float> > &A);
template std::ostream &operator<<(std::ostream &out, const TPZMatrix<std::complex<double> > &A);
