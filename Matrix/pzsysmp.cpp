/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */

#include <memory.h>

#include "pzsysmp.h"
#include "pzfmatrix.h"

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix() : TPZMatrix<TVar>() {
#ifdef CONSTRUCTOR
    cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix(const long rows,const long cols ) : TPZMatrix<TVar>(rows,cols) {
#ifdef CONSTRUCTOR
	cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSYsmpMatrix<TVar>::~TPZSYsmpMatrix() {
	// Deletes everything associated with a TPZSYsmpMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZSYsmpMatrix()\n";
#endif
}



// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

template<class TVar>
const TVar &TPZSYsmpMatrix<TVar>::GetVal(const long r,const long c ) const {
	// Get the matrix entry at (row,col) without bound checking
    long row(r),col(c);
    if (r > c) {
        long temp = r;
        row = col;
        col = temp;
    }
	for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
		if ( fJA[ic] == col ) return fA[ic];
	}
	return this->gZero;
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
							 TPZFMatrix<TVar> &z,
							 const TVar alpha,const TVar beta,const int opt) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	long  ir, ic;
	long  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	if(beta != 0) {
        z = y*beta;
	} else {
        z.Zero();
	}
	
	// Compute alpha * A * x
    long ncols = x.Cols();
    for (long col=0; col<ncols; col++)
    {
        for(long ir=0; ir<this->Rows(); ir++) {
            for(long ic=fIA[ir]; ic<fIA[ir+1]; ic++) {
                long jc = fJA[ic];
                z(ir,col) += alpha * fA[ic] * x.g(jc,col);
                if(jc != ir)
                {
                    z(jc,col) += alpha * fA[ic] * x.g(ir,col);
                }
            }
        }
    }
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
	if(form != EInputFormat) {
		out << "\nTSYsmpMatrix Print: " << title << '\n'
		<< "\tRows    = " << this->Rows()  << '\n'
		<< "\tColumns = " << this->Cols() << '\n';
		int i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<=this->Rows(); i++) {
			std::cout << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(i=this->Rows()+1; i<fIA[this->Rows()]-1; i++) {
			std::cout << i      << "\t\t"
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
	} else {
		TPZMatrix<TVar>::Print(title,out,form);
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::ComputeDiagonal() {
	if(!fDiag.size()) fDiag.resize(this->Rows());
	for(int ir=0; ir<this->Rows(); ir++) {
		fDiag[ir] = GetVal(ir,ir);
	}
}

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrix<TVar>::AutoFill(long nrow, long ncol, int symmetric)
{
    if (!symmetric || nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    orig.AutoFill(nrow,ncol,symmetric);
    
    TPZVec<long> IA(nrow+1);
    TPZStack<long> JA;
    TPZStack<TVar> A;
    IA[0] = 0;
    TPZVec<std::set<long> > eqs(nrow);
    for (long row=0; row<nrow; row++) {
        eqs[row].insert(row);
        for (long col = 0; col<ncol; col++) {
            REAL test = rand()*1./RAND_MAX;
            if (test > 0.5) {
                eqs[row].insert(col);
                if (symmetric) {
                    eqs[col].insert(row);
                }
            }
        }
    }
    long pos=0;
    for (long row=0; row< nrow; row++) {
        for (std::set<long>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
            if(*col >= row)
            {
                JA.Push(*col);
                A.Push(orig(row,*col));
            }
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}

template class TPZSYsmpMatrix<double>;