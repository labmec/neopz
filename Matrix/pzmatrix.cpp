/**
 * @file
 * @brief Contains the implementation of the TPZMatrix methods.
 */
//
// Aut.hor: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tmatrix.c
//
// Class:  TPZMatrix
//
// Obs.:   Implementa uma classe base para as matrizes:
//
//         Sky Line         (TSkylMatrix),
//         Sparse Simetric  (TSSMatrix),
//         Band Simetric    (TBSMatrix),
//         Full             (TPZFMatrix),
//         Band             (TBMatrix),
//         Sparse           (TSpMatrix).
//
// Versao: 04 / 1996.
//

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <memory.h>
#include <sstream>
#include <set>

#include "pzmatrix.h"
#include "pzfmatrix.h"
//#include "pztempmat.h"
#include "pzsolve.h"
#include "pzvec.h"

#include <sstream>
#include <exception>
#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzmatrix"));
static LoggerPtr loggerCheck(Logger::getLogger("pz.checkconsistency"));
#endif

#ifdef DEBUG
#define DEBUG2
#endif

using namespace std;

REAL TPZMatrix::gZero = 0.;


TPZMatrix::~TPZMatrix()
{
	fDecomposed = 0;
	fDefPositive = 0;
	fRow = 0;
	fCol = 0;
}

void
TPZMatrix::Add(const TPZMatrix &A,TPZMatrix&B) const {
	if ((Rows() != A.Rows()) || (Cols() != A.Cols()) ) {
		Error( "Add(TPZMatrix &, TPZMatrix) <different dimensions>" );
	}
	
	B.Redim( A.Rows(), A.Cols() );
	
	for ( int r = 0; r < Rows(); r++ )
		for ( int c = 0; c < Cols(); c++ ) B.PutVal( r, c, GetVal(r,c)+A.GetVal(r,c) );
}

void TPZMatrix::Substract(const TPZMatrix &A,TPZMatrix &result) const {
	
	if ((Rows() != A.Rows()) || (Cols() != A.Cols()) ) {
		Error( "Add(TPZMatrix &, TPZMatrix &) <different dimensions>" );
	}
	
    result.Resize( Rows(), Cols() );
    for ( int r = 0; r < Rows(); r++ ) {
        for ( int c = 0; c < Cols(); c++ ) {
            result.PutVal( r, c, GetVal(r,c)-A.GetVal(r,c) );
        }
    }
}

/** @brief Implements sum of matrices: \f$ A+B \f$ */
TPZFMatrix operator+(const TPZMatrix &A, const TPZMatrix &B ) {
	TPZFMatrix temp;
    temp.Redim( A.Rows(), A.Cols() );
    A.Add(B,temp);
    return temp;
}


/** @brief Implements difference of matrices: \f$ A-B \f$ */
TPZFMatrix operator-(const TPZMatrix &A, const TPZMatrix &B ) {
	TPZFMatrix temp;
    TPZFMatrix res;
    res.Redim( A.Rows(), A.Cols() );
    A.Substract(B,res);
    return temp;
}

/** @brief Implements product of matrices: \f$ A*B \f$ */
TPZFMatrix operator*( TPZMatrix &A, const TPZFMatrix &B ) {
    TPZFMatrix res;
    res.Redim( A.Rows(), B.Cols() );
	A.Multiply(B,res);
	return res;
}

void TPZMatrix::PrepareZ(const TPZFMatrix &y, TPZFMatrix &z,const REAL beta,const int opt,const int stride) const
{
	int numeq = (opt) ? Cols()*stride : Rows()*stride;
	int xcols = y.Cols();
	int ic;
	if(!z.Rows()) return;
	for (ic = 0; ic < xcols; ic++)
	{
		REAL *zp = &z(0,ic), *zlast = zp+numeq;
		if(beta != 0)
		{
			const REAL *yp = &y.g(0,ic);
			if(beta != 1. || (beta == 1. && stride != 1 && &z != &y))
			{
				while(zp < zlast)
				{
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			} else if(&z != &y && stride == 1)
			{
				memcpy(zp,yp,numeq*sizeof(REAL));
			}
		} else
		{
			while(zp != zlast)
			{
				*zp = 0.;
				zp += stride;
			}
		}
	}
}

void TPZMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z, const REAL alpha,const REAL beta,const int opt,const int stride) const {
	if ((!opt && Cols() != x.Rows()*stride) || Rows() != x.Rows()*stride)
		Error( "Operator* <matrixs with incompatible dimensions>" );
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows() || x.Rows() != z.Rows()) {
		Error ("TPZFMatrix::MultiplyAdd incompatible dimensions\n");
	}
	int rows = Rows();
	int cols = Cols();
	int xcols = x.Cols();
	int ic, c, r;
	PrepareZ(y,z,beta,opt,stride);
	REAL val = 0.;
	for (ic = 0; ic < xcols; ic++) {
		if(!opt) {
			for ( c = 0; c<cols; c++) {
				for ( r = 0; r < rows; r++ ) {
					val = z(r*stride,ic) + alpha * GetVal(r,c) * x.GetVal(c*stride,ic);
					z.PutVal(r*stride,ic,val);
				}
			}
		} else {
			for (r = 0; r<rows; r++) {
            	val = 0.;
				for(c = 0; c<cols; c++) {
					val += GetVal(c,r)* x.GetVal(c*stride,ic);
				}
				z.PutVal(r*stride,ic,alpha*val);
			}
		}
	}
}


/*
 void TPZMatrix::InnerProd(TPZFMatrix & D) {
 if ( (Cols() != D.Rows()) && ( D.Rows()!=D.Cols() ) ) {
 Error( "InnerProd (TPZMatrix &) <incompatible dimensions>" );
 }
 TPZFMatrix temp( Rows(), D.Cols() );
 int r,c,i;
 for ( r = 0; r < Rows(); r++ ) {
 for ( c = 0; c < D.Cols(); c++ ) {
 REAL val = 0.0;
 for ( i = 0; i < Cols(); i++ ) {
 val += GetVal( r, i ) * D.GetVal( i, c );
 }
 temp.PutVal( r, c, val );
 }
 }
 D.Resize( Rows(),Rows() );
 for ( r = 0; r < temp.Rows(); r++ ) {
 for (  c = 0; c < Rows(); c++ ) {
 REAL val = 0.0;
 for ( i = 0; i < temp.Cols(); i++ ) {
 val += temp.GetVal( r, i ) * GetVal( c, i );
 }
 D.PutVal( r, c, val );
 }
 }
 }
 */

void TPZMatrix::Identity() {
	
    if ( Cols() != Rows() ) {
        Error( "identity (TPZMatrix *) <TPZMatrix must be square>" );
    }
    for ( int row = 0; row < Rows(); row++) {
        for ( int col = 0; col < Cols(); col++ ) {
            (row == col)? PutVal(row,col,1.):PutVal(row,col,0.);
        }
    }
}



/*************/
/*** Input ***/
void
TPZMatrix::Input(std::istream& in )
{
	
	
	// Read a Matriz (RxC) with format:
	//
	//  rows, cols
	//  a00 a01 a02 ... a0C
	//  a10 a11 a12 ... a1C
	//  a20 a21 a22 ... a2C
	//      ...
	//  aR0 aR1 aR2 ... aRC
	//
	int newRow, newCol;
	in >> newRow;
	in >> newCol;
	Redim( newRow, newCol );
	/*   int  row, col;
	 in >> row >> col;
	 while(row >0) {
	 REAL elem;
	 in >> elem;
	 Put( row, col, elem );
	 in >> row >> col;
	 }
	 */
	int i,j;
	REAL elem;
	for(i=0;i<Rows();i++)
		for(j=0;j<Cols();j++)
		{
			in >> elem;
			Put( i,j, elem );
		}
	
}


/** @brief Overload >> operator to input data of the matrix ***/
std::istream & operator>>(std::istream& in,TPZMatrix &A)
{
	
	// Read a Matriz (RxC) with format:
	//
	//  rows cols
	//  a00 a01 a02 ... a0C
	//  a10 a11 a12 ... a1C
	//  a20 a21 a22 ... a2C
	//      ...
	//  aR0 aR1 aR2 ... aRC
	//
    A.Input(in);
    return in;
}


/*************/
/*** Print ***/
void TPZMatrix::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const{
	
	//  out.width( 8 );
	//  out.precision( 4 );
	
	if(form == EFormatted) {
		out << "Writing matrix '";
		if(name) out << name;
		out << "' (" << Rows() << " x " << Cols() << "):\n";
		
		for ( int row = 0; row < Rows(); row++) {
			out << "\t";
			for ( int col = 0; col < Cols(); col++ ) {
				out << Get( row, col) << "  ";
			}
			out << "\n";
		}
    	out << "\n";
	} else if (form == EInputFormat) {
		out << Rows() << " " << Cols() << endl;
		for ( int row = 0; row < Rows(); row++) {
			for ( int col = 0; col < Cols(); col++ ) {
				REAL val = Get (row, col);
				if(val != 0.) out << row << ' ' << col << ' ' << val << std::endl;
			}
		}
		out << "-1 -1 0.\n";
	} else if( form == EMathematicaInput)
	{
		char number[32];
		out << name << "\n{ ";
		for ( int row = 0; row < Rows(); row++) {
			out << "\n{ ";
			for ( int col = 0; col < Cols(); col++ ) {
				REAL val = Get (row, col);
				sprintf(number, "%16.16lf",(double)val);
				out << number;
				if(col < Cols()-1)
					out << ", ";
				if((col+1) % 6 == 0)out << std::endl;
			}
			out << " }";
			if(row < Rows()-1)
				out << ",";
		}
		
		out << " };\n";
		
	}else if( form == EMatlabNonZeros)
	{
		out << name;
		for ( int row = 0; row < Rows(); row++) {
			out << "\n|";
			for ( int col = 0; col < Cols(); col++ )
				if(IsZero(Get (row, col)) ){
					out << "."; 
				}else{
					out << "#";
				}
		    //out << IsZero(Get (row, col)) ? "." : "#";
			out << "|";
		}
		out << "\n";
	}
	
}

/** @brief Overload << operator to output entries of the matrix ***/
std::ostream &operator<<(std::ostream& out,const TPZMatrix &A) {
    A.Print("operator << ",out);
    return  out;
}

void TPZMatrix::AddKel(TPZFMatrix &elmat, TPZVec<int> &destinationindex) {
	
	int nelem = elmat.Rows();
  	int icoef,jcoef,ieq,jeq;
	if(IsSimetric()) {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			for(jcoef=icoef; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				REAL prevval = GetVal(ieq,jeq);
				prevval += elmat(icoef,jcoef);
				PutVal(ieq,jeq,prevval);
			}
		}
	} else {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			for(jcoef=0; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				REAL prevval = GetVal(ieq,jeq);
				prevval += elmat(icoef,jcoef);
				PutVal(ieq,jeq,prevval);
			}
		}
	}
}


void TPZMatrix::AddKel(TPZFMatrix &elmat, TPZVec<int> &source, TPZVec<int> &destinationindex) {
	
	int nelem = source.NElements();
  	int icoef,jcoef,ieq,jeq,ieqs,jeqs;
	if(IsSimetric()) {
		for(icoef=0; icoef<nelem; icoef++) {
			ieq = destinationindex[icoef];
			ieqs = source[icoef];
			for(jcoef=icoef; jcoef<nelem; jcoef++) {
				jeq = destinationindex[jcoef];
				jeqs = source[jcoef];
				REAL prevval = GetVal(ieq,jeq);
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
				REAL prevval = GetVal(ieq,jeq);
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
int TPZMatrix::PutSub(const int sRow,const int sCol,const TPZFMatrix &A ) {
    int minRow = MIN( A.Rows(), Rows() - sRow );
    int minCol = MIN( A.Cols(), Cols() - sCol );
	
    int row = sRow;
    for ( int r = 0; r < minRow; r++, row++ ) {
        int col = sCol;
        for ( int c = 0; c < minCol; c++, col++ ) {
            PutVal( row, col, A.GetVal( r, c ) );
        }
    }
    return( 1 );
}



/***************/
/*** Get Sub ***/
//
//  Le a sub-matriz de dimensoes (colSize, rowSize) para a matriz A.
//  O inicio da sub-matriz e' dado por (sRow, sCol).
//
int TPZMatrix::GetSub(const int sRow,const int sCol,const int rowSize,
					  const int colSize, TPZFMatrix & A ) const {
    if ( ((sRow + rowSize) > Rows()) || ((sCol + colSize) > Cols()) ) {
        return( Error( "GetSub <t.he sub-matrix is too big>" ) );
    }
    A.Resize( rowSize, colSize );
    int row = sRow;
    for ( int r = 0; r < rowSize; r++, row++ ) {
        int col = sCol;
        for ( int c = 0; c < colSize; c++, col++ ) {
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
int TPZMatrix::AddSub(const int sRow,const int sCol,const TPZFMatrix &A ) {
	
    int minRow = MIN( A.Rows(), Rows() - sRow );
    int minCol = MIN( A.Cols(), Cols() - sCol );
	//  REAL v;
	
	int row = sRow;
    for ( int r = 0; r < minRow; r++, row++ ) {
        int col = sCol;
        for ( int c = 0; c < minCol; c++, col++ ) {
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
int TPZMatrix::InsertSub(const int sRow,const int sCol,const int rowSize,
						 const int colSize,const int pRow,const int pCol, TPZMatrix *pA ) const {
	
	
    if ( ((pRow + rowSize) > pA->Rows()) || ((pCol + colSize) > pA->Cols())) {
        return( Error( "InsertSub <the sub-matrix is too big that target>" ) );
    }
	
    int NewRowSize = rowSize+pRow;
    int NewColSize = colSize+pCol;
	
	
    int row = sRow;
    for ( int r = pRow; r < NewRowSize; r++, row++ ) {
        int col = sCol;
        for ( int c = pCol ; c < NewColSize; c++, col++ ) {
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
int TPZMatrix::AddSub(const int sRow, const int sCol, const int rowSize,
					  const int colSize,const int pRow,const int pCol, TPZMatrix *pA ) const {
    if ( ((pRow + rowSize) > pA->Rows()) || ((pCol + colSize) > pA->Cols())) {
        Error( "AddSub <the sub-matrix is too big that target>" );
    }
    int NewRowSize = rowSize+pRow;
    int NewColSize = colSize+pCol;
	
    int row = sRow;
    for ( int r = pRow; r < NewRowSize; r++, row++ ) {
        int col = sCol;
        for ( int c = pCol ; c < NewColSize; c++, col++ ) {
            pA->PutVal( r, c, GetVal( row, col )+pA->GetVal(r,c));
        }
    }
    return( 1 );
}




/*****************/
/*** Transpose ***/
void TPZMatrix::Transpose(TPZMatrix *T) const {
	T->Resize( Cols(), Rows() );
	
	for ( int r = 0; r < Rows(); r++ ) {
        for ( int c = 0; c < Cols(); c++ ) {
            T->PutVal( c, r, GetVal( r, c ) );
        }
    }
}



/*************/
/*** Solve ***/
int TPZMatrix::SolveDirect( TPZFMatrix &B , DecomposeType dt, std::list<int> &singular) {
	
	switch ( dt ) {
		case ELU:
			return( Solve_LU( &B ,singular)  );
		case ECholesky:
			return( Solve_Cholesky( &B , singular)  );
		case ELDLt:
			return( Solve_LDLt( &B, singular )  );
		default:
			Error( "Solve  < Unknow decomposition type >" );
			break;
	}
	return ( 0 );
}

int TPZMatrix::SolveDirect( TPZFMatrix &B , DecomposeType dt) {
	
	switch ( dt ) {
		case ELU:
			return( Solve_LU( &B)  );
		case ECholesky:
			return( Solve_Cholesky( &B )  );
		case ELDLt:
			return( Solve_LDLt( &B )  );
		default:
			Error( "Solve  < Unknow decomposition type >" );
			break;
	}
	return ( 0 );
}


void TPZMatrix::SolveJacobi(int &numiterations,const TPZFMatrix &F, TPZFMatrix &result,
							TPZFMatrix *residual, TPZFMatrix &scratch, REAL &tol,const int FromCurrent) {
	
	
	if(FromCurrent) {
		Residual(result,F,scratch);
	} else {
		scratch = F;
		result.Zero();
	}
	REAL res;
	res = Norm(scratch);
	int r = Dim();
	int c = F.Cols();
	for(int it=0; it<numiterations && res > tol; it++) {
		for(int ic=0; ic<c; ic++) {
			for(int i=0; i<r; i++) {
				result(i,ic) += (scratch)(i,ic)/GetVal(i,i);
			}
		}
		Residual(result,F,scratch);
		res = Norm(scratch);
	}
	if(residual) *residual = scratch;
}


void TPZMatrix::SolveSOR(int & numiterations, const TPZFMatrix &F,
						 TPZFMatrix &result, TPZFMatrix *residual, TPZFMatrix &/*scratch*/, const REAL overrelax,
						 REAL &tol,const int FromCurrent,const int direction) {
	
	if(residual == &F) {
		cout << "TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
		return;
	}
	REAL res = 2*tol+1.;
	if(residual) res = Norm(*residual);
	if(!FromCurrent) {
		result.Zero();
	}
	int r = Dim();
	int c = F.Cols();
	int i,ifirst = 0, ilast = r, iinc = 1;
	int it;
	if(direction == -1) {
		ifirst = r-1;  //misael
		ilast = -1;
		iinc = -1;
	}
	REAL eqres;
	for(it=0; it<numiterations && res > tol; it++) {
		res = 0.;
		for(int ic=0; ic<c; ic++) {
			for(i=ifirst; i!=ilast; i+= iinc) {
				eqres = F.GetVal(i,ic);
				for(int j=0; j<r; j++) {
					eqres -= GetVal(i,j)*result(j,ic);
				}
				res += eqres*eqres;
				result(i,ic) += overrelax*eqres/GetVal(i,i);
			}
		}
		res = sqrt(res);
	}
	if(residual) Residual(result,F,*residual);
	numiterations = it;
	tol = res;
}

void TPZMatrix::SolveSSOR(int &numiterations, const TPZFMatrix &F,
						  TPZFMatrix &result, TPZFMatrix *residual, TPZFMatrix &scratch, const REAL overrelax,
						  REAL &tol,const int FromCurrent) {
	REAL res = (tol*REAL(2.))+REAL(1.);
	int i, one = 1;
	int fromcurrent = FromCurrent;
	for(i=0; i<numiterations && res > tol; i++) {
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

void TPZMatrix::SolveCG(int &numiterations, TPZSolver &preconditioner,
						const	TPZFMatrix &F, TPZFMatrix &result,
						TPZFMatrix *residual, REAL &tol, const int FromCurrent) {
	CG(*this, result, F, preconditioner, residual, numiterations, tol, FromCurrent);
}

#include "gmres.h"
#include "pzstepsolver.h"

void TPZMatrix::SolveGMRES(int &numiterations, TPZSolver &preconditioner,
						   TPZFMatrix &H, int &numvectors,
						   const TPZFMatrix &F, TPZFMatrix &result,
						   TPZFMatrix *residual, REAL &tol,const int FromCurrent)  
{
    if(F.Cols() > 1)
    {
        int locnumiter = numiterations;
        REAL loctol = tol;
        int nrow = F.Rows();
        int ncol = F.Cols();
        int col;
        //preconditioner.Solve(F, result);
        for (col=0; col<ncol; col++) {
            std::cout << "Column " << col << std::endl;
            numiterations = locnumiter;
            tol = loctol;
            TPZFMatrix FCol(nrow,1);
            memcpy(&FCol(0,0), &F.GetVal(0,col), nrow*sizeof(REAL));
            TPZFMatrix resultCol(nrow,1,&result(0,col),nrow);
            GMRES(*this,resultCol,FCol,preconditioner,H,numvectors,numiterations,tol,residual,FromCurrent);
        }
    }
    else
    {
        GMRES(*this,result,F,preconditioner,H,numvectors,numiterations,tol,residual,FromCurrent);
    }
}

#include "bicg.h"

void TPZMatrix::SolveBICG(int &numiterations, TPZSolver &preconditioner,
						  const TPZFMatrix &F,
						  TPZFMatrix &result,
						  REAL &tol)  {
	BiCG<TPZMatrix,TPZFMatrix,TPZSolver,REAL>(*this,result,F,preconditioner,numiterations,tol);
}

#include "bicgstab.h"

void TPZMatrix::SolveBICGStab(int &numiterations, TPZSolver &preconditioner,
							  const TPZFMatrix &F, TPZFMatrix &result,
							  TPZFMatrix *residual, REAL &tol,const int FromCurrent)  {
	BiCGSTAB<TPZMatrix,TPZFMatrix,TPZSolver,REAL>(*this,result,F,preconditioner,numiterations,tol,residual,FromCurrent);
}



#include "ir.h"

void TPZMatrix::SolveIR(int &numiterations, TPZSolver &preconditioner,
						const TPZFMatrix &F, TPZFMatrix &result,
						TPZFMatrix *residual, REAL &tol,
						const int FromCurrent)  {
	IR(*this,result,F,preconditioner,residual,numiterations,tol,FromCurrent);
}


/*****************/
/*** Decompose_LU ***/
int TPZMatrix::Decompose_LU(std::list<int> &singular) {
	return Decompose_LU();
}

int TPZMatrix::Decompose_LU() {
	
	if (  fDecomposed && fDecomposed != ELU)  Error( "Decompose_LU <Matrix already Decomposed with other scheme>" );
	if (fDecomposed) return 1;
	
	REAL nn, pivot;
	int  min = ( Cols() < (Rows()) ) ? Cols() : Rows();
	
	for ( int k = 0; k < min ; k++ ) {
		if (IsZero( pivot = GetVal(k, k))) Error( "Decompose_LU <matrix is singular>" );
		for ( int i = k+1; i < Rows(); i++ ) {
			nn = GetVal( i, k ) / pivot;
			PutVal( i, k, nn );
			for ( int j = k+1; j < Cols(); j++ ) PutVal(i,j,GetVal(i,j)-nn*GetVal(k,j));
		}
	}
	fDecomposed=ELU;
	return 1;
}



/****************/
/*** Substitution ***/

int TPZMatrix::Substitution( TPZFMatrix *B ) const{
	
    int rowb = B->Rows();
    int colb = B->Cols();
    if ( rowb != Rows() )
		Error( "SubstitutionLU <incompatible dimensions>" );
	int i;
    for ( i = 0; i < rowb; i++ ) {
        for ( int col = 0; col < colb; col++ ) {
            for ( int j = 0; j < i; j++ ) {
                B->PutVal( i, col, B->GetVal(i, col)-GetVal(i, j) * B->GetVal(j, col) );
            }
        }
    }
    for (int col=0; col<colb; col++) {
        for ( i = rowb-1; i >= 0; i-- ) {
            for ( int j = i+1; j < rowb ; j++ ) {
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

int TPZMatrix::Decompose_LDLt(std::list<int> &singular) {
	return Decompose_LDLt();
}

int TPZMatrix::Decompose_LDLt() {
	
	if (  fDecomposed && fDecomposed != ELDLt) {
		Error( "Decompose_LDLt <Matrix already Decomposed with other scheme> " );
	} else if(fDecomposed ) {
		return ELDLt;
	}
	if ( Rows()!=Cols() ) Error( "Decompose_LDLt <Matrix must be square>" );
	
	int j,k,l,dim=Rows();
	
	for ( j = 0; j < dim; j++ ) {
		for ( k=0; k<j; k++) {
			PutVal( j,j,GetVal(j,j) - GetVal(k,k)*GetVal(k,j)*GetVal(k,j) );
		}
		for ( k=0; k<j; k++) {
			for( l=j+1; l<dim;l++) {
				PutVal(l,j, GetVal(l,j)-GetVal(k,k)*GetVal(j,k)*GetVal(l,k) );
				PutVal(j,l,GetVal(l,j) );
			}
		}
		REAL tmp = GetVal(j,j);
		if ( IsZero(tmp) ) Error( "Decompose_LDLt <Zero on diagonal>" );
		for( l=j+1; l<dim;l++) {
			PutVal(l,j, GetVal(l,j)/GetVal(j,j) ) ;
			PutVal(j,l, GetVal(l,j) );
		}
	}
	fDecomposed  = ELDLt;
	fDefPositive = 0;
	return( 1 );
}
int TPZMatrix::Decompose_Cholesky(std::list<int> &singular) {
	if (  fDecomposed && fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
	if (  fDecomposed ) return ECholesky;
	if ( Rows()!=Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
	//return 0;
	
	int dim=Dim();
	for (int i=0 ; i<dim; i++) {
		for(int k=0; k<i; k++) {             //elementos da diagonal
			PutVal( i,i,GetVal(i,i)-GetVal(i,k)*GetVal(i,k) );
		}
		if(GetVal(i,i) <= 1.e-12)
		{
			singular.push_back(i);
			PutVal(i,i,1.);
		}
		REAL tmp = sqrt(GetVal(i,i));
        PutVal( i,i,tmp );
        for (int j=i+1;j<dim; j++) {           //elementos fora da diagonal
            for(int k=0; k<i; k++) {
                PutVal( i,j,GetVal(i,j)-GetVal(i,k)*GetVal(k,j) );
            }
            REAL tmp2 = GetVal(i,i);
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

int TPZMatrix::Decompose_Cholesky() {
	if (  fDecomposed && fDecomposed != ECholesky) Error( "Decompose_Cholesky <Matrix already Decomposed>" );
	if (  fDecomposed ) return ECholesky;
	if ( Rows()!=Cols() ) Error( "Decompose_Cholesky <Matrix must be square>" );
	//return 0;
	
	int dim=Dim();
	for (int i=0 ; i<dim; i++) {
		for(int k=0; k<i; k++) {             //elementos da diagonal
			PutVal( i,i,GetVal(i,i)-GetVal(i,k)*GetVal(i,k) );
		}
		REAL tmp = sqrt(GetVal(i,i));
        PutVal( i,i,tmp );
        for (int j=i+1;j<dim; j++) {           //elementos fora da diagonal
            for(int k=0; k<i; k++) {
                PutVal( i,j,GetVal(i,j)-GetVal(i,k)*GetVal(k,j) );
            }
            REAL tmp2 = GetVal(i,i);
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


int TPZMatrix::Subst_Forward( TPZFMatrix *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ECholesky)
		return( 0 );
	for ( int r = 0; r < Dim(); r++ ) {
		REAL pivot = GetVal( r, r );
		for ( int c = 0; c < B->Cols();  c++ ) {
			// Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
			//
			REAL sum = 0.0;
			for ( int i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
			
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
int TPZMatrix::Subst_Backward( TPZFMatrix *B ) const {
	
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ECholesky) return( 0 );
	for ( int r = Dim()-1;  r >= 0;  r-- ) {
		REAL pivot = GetVal( r, r );
		for ( int c = 0; c < B->Cols(); c++ ) {
			// Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
			//
			REAL sum = 0.0;
			for ( int i = Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
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
int TPZMatrix::Subst_LForward( TPZFMatrix *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ELDLt) {
		Error("TPZMatrix::Subst_LForward incompatible dimensions\n");
	}
    for ( int r = 0; r < Dim(); r++ ) {
        for ( int c = 0; c < B->Cols();  c++ )    {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = 0, ..., r-1.
            //
            REAL sum = 0.0;
            for ( int i = 0; i < r; i++ ) sum += GetVal(r, i) * B->GetVal(i, c);
			
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
int TPZMatrix::Subst_LBackward( TPZFMatrix *B ) const {
	if ( (B->Rows() != Dim()) || !fDecomposed || fDecomposed != ELDLt){
		Error("TPZMatrix::Subst_LBackward incompatible dimensions \n");
	}
	
    for ( int r = Dim()-1;  r >= 0;  r-- ) {
        for ( int c = 0; c < B->Cols(); c++ ) {
            // Faz sum = SOMA( A[r,i] * B[i,c] ); i = N, ..., r+1.
            //
            REAL sum = 0.0;
            for ( int i = Dim()-1; i > r; i-- ) sum += GetVal(r, i) * B->GetVal(i, c);
			
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
int TPZMatrix::Subst_Diag( TPZFMatrix *B ) const {
    if ( (B->Rows() != Dim())) {
        Error("TPZMatrix::Subst_Diag incompatible dimensions\n");
	}
	for ( int r = 0; r < Dim(); r++ ) {
		REAL pivot = GetVal( r, r );
		for ( int c = 0; c < B->Cols(); c++ ) {
			B->PutVal( r, c, B->GetVal( r, c ) / pivot );
		}
	}
	return( 1 );
}



/************************** Private **************************/


/*************/
/*** Error ***/
int TPZMatrix::Error(const char *msg ,const char *msg2) {
    ostringstream out;
    out << "TPZMatrix::" << msg;
    if(msg2) out << msg2;
    out << ".\n";
    LOGPZ_ERROR (logger, out.str().c_str());
    DebugStop();
	std::bad_exception myex;
	throw myex;
	//    DebugStop();
	//    return 0;
}

void TPZMatrix::Read( TPZStream &buf, void *context ){
	TPZSaveable::Read(buf,context);
	buf.Read(&fRow,1);
	buf.Read(&fCol,1);
	int tmp;
	buf.Read(&tmp,1);
	fDecomposed = (char) tmp;
}

void TPZMatrix::Write( TPZStream &buf, int withclassid ) {
	TPZSaveable::Write(buf,withclassid);
	buf.Write(&fRow,1);
	buf.Write(&fCol,1);
	int tmp = fDecomposed;
	buf.Write(&tmp,1);
}

// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZMatrix::Compare(TPZSaveable *copy, bool override)
{
	TPZMatrix *copmat = dynamic_cast<TPZMatrix *> (copy);
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

// Compare the object for identity with the object pointed to, eventually copy the object
/**
 * compare both objects bitwise for identity. Put an entry in the log file if different
 * overwrite the calling object if the override flag is true
 */
bool TPZMatrix::Compare(TPZSaveable *copy, bool override) const
{
	TPZMatrix *copmat = dynamic_cast<TPZMatrix *> (copy);
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


void TPZMatrix::GetSub(const TPZVec<int> &indices,TPZFMatrix &block) const
{
    int nel = indices.NElements();
    int iel,jel;
    for(iel=0; iel<nel; iel++)
    {
		for(jel=0; jel<nel; jel++)
		{
			block(iel,jel) = GetVal(indices[iel],indices[jel]);
		}
    }
	
}

int TPZMatrix::VerifySymmetry(REAL tol) const{
	int nrows = this->Rows();
	int ncols = this->Cols();
	if (nrows != ncols) return 0;
	
	for( int i = 0; i < nrows; i++){
		for(int j = 0; j <= i; j++){
			if ( fabs( this->Get(i,j) - this->Get(j,i) ) > tol ) {
				cout << "Elemento: " << i << ", " << j << "  -> " << fabs(this->Get(i,j) - this->Get(j,i) ) << "/" <<
				this->Get(i,j) << endl;
				return 0;
			}
		}
	}
	return 1;
}

bool TPZMatrix::CompareValues(TPZMatrix &M, REAL tol){
	
	int nrows = this->Rows();
	int ncols = this->Cols();
	if ( (M.Rows() != nrows) || (M.Cols() != ncols) ) return false;
	
	REAL diff;
	for( int i = 0; i < nrows; i++)
		for( int j = 0; j < ncols; j++){
			diff = fabs( this->Get(i,j) - M.Get(i,j) );
			if (diff > tol) return false;
		}
	
	return true;
}

REAL TPZMatrix::ReturnNearestValue(REAL val, TPZVec<REAL>& Vec, REAL tol)
{
    REAL diff0 = fabs(val - Vec[0]) >= tol ?  (val - Vec[0]) : 1.E10;
    REAL diff1, res = Vec[0];
    for(int i = 1; i < Vec.NElements(); i++)
    {
        diff1 = val - Vec[i];
        diff0 = ( fabs(diff1) >= tol && fabs(diff1) < fabs(diff0) ) ? res = Vec[i],diff1 : diff0;
    }
    return res;
}

bool TPZMatrix::SolveEigensystemJacobi(int &numiterations, REAL & tol, TPZVec<REAL> & Eigenvalues, TPZFMatrix & Eigenvectors) const{
	
	int NumIt = numiterations;
	REAL tolerance = tol;
	
#ifdef DEBUG2
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
	
	const int size = this->Rows();
	
	/** Making a copy of this */
	TPZFNMatrix<9> Matrix(size,size); //fast constructor in case of this is a stress or strain tensor.
	for(int i = 0; i < size; i++) for(int j = 0; j < size; j++) Matrix(i,j) = this->Get(i,j);
	
	/** Compute Eigenvalues */
	bool result = Matrix.SolveEigenvaluesJacobi(numiterations, tol, &Eigenvalues);
	if (result == false) return false;
	
	/** Compute Eigenvectors */
	TPZFNMatrix<3> VecIni(size,1,0.), VecIni_cp(size,1,0.);
	
	Eigenvectors.Resize(size, size);
	Eigenvectors.Zero();
	for(int eigen = 0; eigen < size; eigen++)
	{
        for(int i = 0; i < size; i++) VecIni.PutVal(i,0,rand());
		
        TPZFNMatrix<9> Matrix(*this);
		
        REAL answ = ReturnNearestValue(Eigenvalues[eigen], Eigenvalues,1.E-5);
        if(fabs(answ - Eigenvalues[eigen]) > 1.E-5)
        {
            for(int i = 0; i < size; i++) Matrix.PutVal(i,i, this->GetVal(i,i) - (Eigenvalues[eigen] - 0.01 * fabs(answ-Eigenvalues[eigen])) );
        }
        else
        {
            for(int i = 0; i < size; i++) Matrix.PutVal(i,i, this->GetVal(i,i) - (Eigenvalues[eigen] - 0.01) );
        }
		
        /** Normalizing Initial Eigenvec */
        REAL norm1 = 0.;
        for(int i = 0; i < size; i++) norm1 += VecIni(i,0) * VecIni(i,0); norm1 = sqrt(norm1);
        for(int i = 0; i < size; i++) VecIni(i,0) = VecIni(i,0)/norm1;
		
        int count = 0;
        double dif = 10., difTemp = 0.;
		
        while(dif > tolerance && count <= NumIt)
        {
			for(int i = 0; i < size; i++) VecIni_cp(i,0) = VecIni(i,0);
			
			/** Estimating Eigenvec */
			Matrix.Solve_LU(&VecIni);
			
			/** Normalizing Final Eigenvec */
			REAL norm2 = 0.;
			for(int i = 0; i < size; i++) norm2 += VecIni(i,0) * VecIni(i,0);
			norm2 = sqrt(norm2);
			
			difTemp = 0.;
			for(int i = 0; i < size; i++)
			{
                VecIni(i,0) = VecIni(i,0)/norm2;
                difTemp += (VecIni_cp(i,0) - VecIni(i,0)) * (VecIni_cp(i,0) - VecIni(i,0));
			}
			dif = sqrt(difTemp);
			count++;
        }
		
        /** Copy values from AuxVector to Eigenvectors */
        for(int i = 0; i < size; i++)
        {
			double val = VecIni(i,0);
			if(fabs(val) < 1.E-5) val = 0.;
			Eigenvectors(eigen,i) = val;
        }
		
#ifdef DEBUG2
        double norm = 0.;
        for(int i = 0; i < size; i++) norm += VecIni(i,0) * VecIni(i,0);
        if (fabs(norm - 1.) > 1.e-10)
        {
            PZError << __PRETTY_FUNCTION__ << endl;
        }
        if(count > NumIt-1)// O metodo nao convergiu !!!
        {
            PZError << __PRETTY_FUNCTION__ << endl;
#ifdef LOG4CXX
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


bool TPZMatrix::SolveEigenvaluesJacobi(int &numiterations, REAL & tol, TPZVec<REAL> * Sort){
	
#ifdef DEBUG2
	if (this->Rows() != this->Cols()){
		PZError << __PRETTY_FUNCTION__ << " - Jacobi method of computing eigenvalues requires a symmetric square matrix. this->Rows = " << this->Rows() << " - this->Cols() = " << this->Cols() << endl;
		return false;
	}
	
	if (this->VerifySymmetry(1.e-8) == false){
		PZError << __PRETTY_FUNCTION__ << " - Jacobi method of computing eigenvalues requires a symmetric square matrix. This matrix is not symmetric." << endl;
		return false;
	}
#endif
	
	int iter = 0;
	REAL res = 2. * tol;
	const int size = this->Rows();
	int i, j;
	int p, q;
	REAL maxval, theta, cost, sint, aux, aux2, Spp, Sqq, Spq;
	while (iter < numiterations){
		/** First of all find the max value off diagonal */
		maxval = 0.;
		for(i = 0; i < size; i++){
			for(j = 0; j < i; j++) {
				if( fabs(this->operator ( )(i,j) ) > maxval ) {
					p = i;
					q = j;
					maxval = fabs( this->operator ( )(i,j) );
				}//if
			}//for j
		}//for i
		
		//    cout << "iter: " << iter << " - " << maxval << endl;
		
		/** Check if max value off diagonal is lesser than required tolerance */
		res = maxval;
		if (res < tol) break;
		//     {
		//       tol = res;
		//       numiterations = iter;
		//       return true;
		//     }//if
		
		/** Compute angle of rotation */
		theta = 0.5 * atan(2. * this->operator ( )(p,q) / (this->operator ( )(q,q) - this->operator ( )(p,p) ) );
		cost = cos(theta);
		sint = sin(theta);
		
		/** Apply rotation */
		for(i = 0; i < size; i++){
			if (i != p && i != q){
				
				aux = this->operator ( )(i,p) * cost - this->operator ( )(i,q) * sint;
				aux2 = this->operator ( )(i,p) * sint + this->operator ( )(i,q) * cost;
				
				this->operator ( )(i,p) = aux;
				this->operator ( )(p,i) = aux;
				
				this->operator ( )(i,q) = aux2;
				this->operator ( )(q,i) = aux2;
				
			}//if
		}//for i
		
		Spp = this->operator ( )(p,p) * cost * cost -2. * this->operator ( )(p,q) * sint * cost + this->operator ( )(q,q) * sint * sint;
		Sqq = this->operator ( )(p,p) * sint * sint +2. * this->operator ( )(p,q) * sint * cost + this->operator ( )(q,q) * cost * cost;
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
		for(i = 0; i < size; i++) myset.insert( this->operator ( )(i,i) );
		
#ifdef DEBUG2
		if ((int)myset.size() != size) PZError << __PRETTY_FUNCTION__ << " - ERROR!" << endl;
#endif
		
		Sort->Resize(size);
		multiset< REAL >::iterator w, e = myset.end();
		for(i = size - 1, w = myset.begin(); w != e; w++, i--){
			Sort->operator [ ](i) = *w;
		}//for
		
	}//if (Sort)
	
	
	if (res < tol){
		tol = res;
		numiterations = iter;
		return true;
	}
	
	tol = res;
	numiterations = iter;
	return false;
	
}//method

REAL TPZMatrix::MatrixNorm(int p, int numiter, REAL tol) const{
	const int n = this->Rows();
	if (!n) return 0.;
	if (n != this->Cols()){
		PZError << __PRETTY_FUNCTION__
		<< " - matrix must be square - Rows() = "
		<< this->Rows() << " - Cols() = "
		<< this->Cols() << std::endl;
	}
	switch(p){
		case 0:{
			REAL max = 0.;
			int i, j;
			for(i = 0; i < n; i++){
				REAL sum = 0.;
				for(j = 0; j < n; j++) sum += fabs( this->Get(i,j) );
				if (sum > max) max = sum;
			}
			return max;
		}
		case 1:{
			REAL max = 0.;
			int i, j;
			for(i = 0; i < n; i++){
				REAL sum = 0.;
				for(j = 0; j < n; j++) sum += fabs( this->Get(j,i) );
				if (sum > max) max = sum;
			}
			return max;
		}
		case 2:{
			TPZFMatrix transp(n,n);
			int i, j, k;
			//TRANSPOSE
			for(i = 0; i < n; i++) for(j = 0; j < n; j++) transp(i,j) = this->Get(j,i);
			//MULTIPLY transp = transp.this
			TPZVec<REAL> ROW(n);
			for(i = 0; i < n; i++){
				for(j = 0; j < n; j++) ROW[j] = transp(i,j);
				for(j = 0; j < n; j++){
					REAL sum = 0.;
					for(k = 0; k < n; k++){
						sum += ROW[k] * this->Get(k,j);
					}//for k
					transp(i,j) = sum;
				}//for j
			}//for i
			
			TPZVec<REAL> EigenVal;
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

int TPZMatrix::Inverse(TPZFMatrix &Inv){
	const int n = this->Rows();
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
	int i;
	for(i = 0; i < n; i++) Inv(i,i) = 1.;
	
	const int issimetric = this->IsSimetric();
	if (issimetric)  return this->SolveDirect(Inv, ELDLt);
	if (!issimetric) return this->SolveDirect(Inv, ELU);
	return 0;
}//method

REAL TPZMatrix::ConditionNumber(int p, int numiter, REAL tol){
	int localnumiter = numiter;
	REAL localtol = tol;
	TPZFMatrix Inv;
	REAL thisnorm = this->MatrixNorm(p, localnumiter, localtol);
	if (!this->Inverse(Inv)){
		PZError << __PRETTY_FUNCTION__ << " - it was not possible to compute the inverse matrix." << std:: endl;
		return 0.;
	}
	REAL invnorm = Inv.MatrixNorm(p, numiter, tol);
	return thisnorm * invnorm;
}

void TPZMatrix::Multiply(const TPZFMatrix &A, TPZFMatrix&B, int opt, int stride) const {
	if ((opt==0 && Cols()*stride != A.Rows()) || (opt ==1 && Rows()*stride != A.Rows()))
		Error( "Multiply (TPZMatrix &,TPZMatrix&) <incompatible dimensions>" );
	if(!opt && (B.Rows() != Rows()*stride || B.Cols() != A.Cols())) {
		B.Redim(Rows()*stride,A.Cols());
	}
	else if (opt && (B.Rows() != Cols()*stride || B.Cols() != A.Cols())) {
		B.Redim(Cols()*stride,A.Cols());
	}
	MultAdd( A, B, B, 1.0, 0.0, opt,stride);
}

/** Fill the matrix with random values (non singular matrix) */
void TPZMatrix::AutoFill() {
	long i, j;
	REAL val, sum;
	/** Fill data */
	for(i=0;i<Rows();i++) {
		sum = 0.0;
		for(j=0;j<Cols();j++) {
			val = ((REAL)rand())/RAND_MAX;
			if(!PutVal(i,j,val))
				Error("AutoFill (TPZMatrix) failed.");
			if(i!=j) sum += fabs(val);
		}
		/** Making diagonally dominant and non zero in diagonal */
		if(sum > GetVal(i,i))            // Deve satisfazer:  |Aii| > SUM( |Aij| )  sobre j != i
			PutVal(i,i,sum);
		// To sure diagonal is not zero.
		if(IsZero(sum) && IsZero(GetVal(i,i)))
			PutVal(i,i,1.);
	}
}
