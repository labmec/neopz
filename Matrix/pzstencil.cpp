/**
 * @file
 * @brief Contains the implementation of the TPZStencilMatrix methods.
 */
/******************************************************************************
 *
 * Class realization:   TPZStencilMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on sparse matrices stored by
 *                      stencils
 *
 * Operations:          Mult
 *                      MultAdd
 *                      Print
 *
 * Solvers:             SOR
 *                      SSOR
 *
 *****************************************************************************/


#include "pzstencil.h"
#include "pzfmatrix.h"

#include <memory.h>
#include <string>
using namespace std;

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

TPZStencilMatrix::TPZStencilMatrix( int rows, int cols ) {
	// Constructs an empty TPZStencilMatrix
	fRows = rows;
	fCols = cols;
	fMystencils = 0;
	fNumberOfStencilPointers = 0;
	fStencilNumbers = 0;
	fSolver = -1;
	fMaxIterations = 4;
	fSORRelaxation = 1.;
	fDiag = 0;
#ifdef CONSTRUCTOR
	cerr << "TPZStencilMatrix(int rows,int cols)\n";
#endif
}

TPZStencilMatrix::~TPZStencilMatrix() {
	// Deletes everything associated with a TPZStencilMatrix
	for(int i=0; i<fNumberOfStencilPointers; i++) delete fMystencils[i];
	delete fMystencils;
	fMystencils = 0;
	fNumberOfStencilPointers = 0;
	delete fStencilNumbers;
	fStencilNumbers = 0;
	delete fDiag;
	fDiag = 0;
	fSolver = -1;
#ifdef DESTRUCTOR
	cerr << "~TPZStencilMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

const REAL & TPZStencilMatrix::GetVal(const int row,const int col ) const {
	// Get the matrix entry at (row,col) without bound checking
	int stenindex;
	MPStencil *st;
	
	// Walk past each row of the matrix trying to determine where the
	// first stencil in row begins
	static TPZStencilMatrix const *whatsthis = 0;
	static int ix=0;
	static int ir = 0;
	if( whatsthis != this || row < ir ) {
		whatsthis = this;
		for( ir=ix=0; ir<row; ir++) {
			stenindex = fStencilNumbers[ir];
			st = fMystencils[stenindex];
			ix += st->fInc;
		}
	}
	else {
		for( ; ir<row; ir++) {
			stenindex = fStencilNumbers[ir];
			st = fMystencils[stenindex];
			ix += st->fInc;
		}
	}
	
	// Now look through the requested row and see if there is anything
	// in column col
	st = fMystencils[fStencilNumbers[row]];
	int nitems = st->fNumberOfItems;
	int *ia = st->fIA;
	REAL *a = st->fA;
	int it=0;
	
	while(++it<nitems) {
		int numintegers = *ia++;
		int *ialast = ia+numintegers;
		while(ia < ialast) {
			if ( ix+(*ia++) == col ) return *a;
			it++;
		}
		a+=numintegers+1;
	}
	return gZero;
}

// ****************************************************************************
// 
// Multiply and Multiply-Add
// 
// ****************************************************************************

void TPZStencilMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y,
							   TPZFMatrix &z,
							   const REAL alpha,const REAL beta,const int opt,const int stride ) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	int ix=0;
	int r = (opt) ? Cols() : Rows();
	if(beta != 0) {
		REAL *zp = &(z(0,0)), *zlast = zp+r*stride;
		const REAL *yp = &(y.GetVal(0,0));
		if(beta != 1. || (&z != &y && stride != 1)) {
			while(zp < zlast) {
				*zp = beta * (*yp);
				zp += stride;
				yp += stride;
			}
		}
		else if(&z != &y) {
			memcpy(zp,yp,r*sizeof(REAL));
		}
	} else {
		REAL *zp = &(z(0,0)), *zlast = zp+r*stride;
		while(zp != zlast) {
			*zp = 0.;
			zp += stride;
		}
	}
	if(opt == 0) {
		
		for(int ir=0; ir<fRows; ir++) {
			int stenindex = fStencilNumbers[ir];
			MPStencil *st = fMystencils[stenindex];
			int nitems = st->fNumberOfItems;
			int *ia = st->fIA;
			REAL *a = st->fA;
			int it=0;
			while(++it<nitems) {
				int numintegers = *ia++;
				int *ialast = ia+numintegers;
				REAL val;
				val = 0.;
				const REAL *xp = &(x.GetVal(ix*stride,0));
				while(ia < ialast) {
					val += *(xp+((*ia++)*stride));
					it++;
				}
				z(ir*stride,0) += alpha*val*(*a);
				a+=numintegers+1;
			}
			ix += st->fInc;
		}
	} else {
		
		int ir;
		for(ir=0; ir<fRows; ir++) {
			MPStencil *st = fMystencils[fStencilNumbers[ir]];
			int nitems = st->fNumberOfItems;
			int *ia = st->fIA;
			REAL *a = st->fA;
			int it=0;
			while(++it<nitems) {
				int numintegers = *ia++;
				int *ialast = ia + numintegers;
				REAL xval = alpha*x.GetVal(ir*stride,0)*(*a);
				while(ia<ialast) {
					z((ix+(*ia++))*stride,0) += xval;
					it++;
				}
				a += numintegers+1;
			}
			ix += st->fInc;
		}
	}
}

// ****************************************************************************
// 
// Print the matrix
// 
// ****************************************************************************

void TPZStencilMatrix::Print(const char *title, ostream &out,const MatrixOutputFormat form ) const {
	// Print the matrix along with a identification title
	if(form != EInputFormat) {
		out << "\nTStencilMatrix Print: " << title << '\n'
        << "\tRows    = " << fRows  << '\n'
        << "\tColumns = " << fCols << '\n';
		out << '\t';
		for(int ir=0; ir<fNumberOfStencilPointers; ir++) {
			MPStencil *st = fMystencils[ir];
			if (st != 0) {
				int nitems = st->fNumberOfItems;
				int *ia = st->fIA;
				REAL *a = st->fA;
				out << "\nStencil " << ir
				<< ", increment " << st->fInc  << ", "
				<< nitems << " items\n";
				out << "\tIA\tA\n";
				out << "\t--\t-\n";
				int it=0;
				while(++it<nitems) {
					int numintegers = *ia;
					out << '\t' << *ia++ << '\t' << *a << '\n';
					int *ialast = ia+numintegers;
					while(ia < ialast) {
						out << '\t' << *ia++ << '\n';
						it++;
					}
					a+=numintegers+1;
				}
			}
		}
		
		if(fStencilNumbers) {
			out << "\nIndex\tStencil number\n";
			out <<   "-----\t--------------\n";
			for(int i=0; i<fRows; i++) {
				out << i << '\t' << fStencilNumbers[i] << '\n';
			}
		}
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************


void TPZStencilMatrix::ComputeDiagonal() {
	
	if(fDiag) return;
	fDiag = new REAL[Rows()];
	for(int ir=0; ir<fRows; ir++) {
		MPStencil *st = fMystencils[fStencilNumbers[ir]];
		int nitems = st->fNumberOfItems;
		int *ia = st->fIA;
		REAL *a = st->fA;
		if(st->fInc != 1) {
			cout << "TPZStencilMatrix::ComputeDiagonal "
			"Computing the diagonal of a nonsquare matrix?";
		}
		int it=0;
		while(++it<nitems) {
			int numintegers = *ia++;
			int *ialast = ia+numintegers;
			while(ia < ialast) {
				if(*ia++ == 0) {
					fDiag[ir] = *a;
					it=nitems;
					break;
				}
				it++;
			}
			a+=numintegers+1;
		}
	}
}

void TPZStencilMatrix::SolveSOR( int &numiterations,const TPZFMatrix &rhs, TPZFMatrix &x,
								TPZFMatrix *residual, TPZFMatrix &/*scratch*/,
								const REAL overrelax, REAL &tol,
								const int FromCurrent, const int direction ) {
	
	if(!fDiag) ComputeDiagonal();
	int irStart = 0,irLast = fRows,irInc = 1;
	if(direction < 0) {
		irStart = fRows-1;
		irLast = -1;
		irInc = -1;
	}
	if(!FromCurrent) x.Zero();
	REAL eqres = 2.*tol;
	int iteration;
	for(iteration=0; iteration<numiterations && eqres >= tol; iteration++) {
		eqres = 0.;
		int ir=irStart;
		while(ir != irLast) {
			int stenindex = fStencilNumbers[ir];
			MPStencil *st = fMystencils[stenindex];
			int nitems = st->fNumberOfItems;
			int *ia = st->fIA;
			REAL *a = st->fA;
			int it=0;
			REAL xnewval=rhs.GetVal(ir,0);
			while(++it<nitems) {
				int numintegers = *ia++;
				int *ialast = ia+numintegers;
				REAL val =0.;
				REAL *xp = &x(ir,0);
				while(ia < ialast) {
					val -= *(xp+(*ia++));
					it++;
				}
				xnewval += *a * val;
				a+=numintegers+1;
			}
			eqres += xnewval*xnewval;
			x(ir,0) += overrelax*(xnewval/fDiag[ir]);
			ir += irInc;
		}
		eqres = sqrt(eqres);
	}
	tol = eqres;
	numiterations = iteration;
	if(residual) Residual(x,rhs,*residual);
}


// ****************************************************************************
//
// Construct pieces of the stencil matrices
//
// ****************************************************************************

void TPZStencilMatrix::SetStencil( int stencilnumber, int inc,
								  int *IA, REAL *A ) {
	// initiates Stencil number "stencilnumber" with the data
	
	if(stencilnumber < 0) return;
	if(stencilnumber >= fNumberOfStencilPointers)
		IncreaseStencilPointers(stencilnumber);
	if(fMystencils[stencilnumber]) delete fMystencils[stencilnumber];
	fMystencils[stencilnumber] = new MPStencil(inc,IA,A);
}

void TPZStencilMatrix::IncreaseStencilPointers( int numsten ) {
	int newnum = fNumberOfStencilPointers+10;
	if(newnum < numsten) newnum = numsten+1;
	MPStencil **newptr = new MPStencil*[newnum];
	memcpy(newptr,fMystencils,sizeof(MPStencil *)*fNumberOfStencilPointers);
	int i=fNumberOfStencilPointers;
	while(i< newnum) newptr[i++] = 0;
	delete fMystencils;
	fMystencils = newptr;
	fNumberOfStencilPointers = newnum;
}

void TPZStencilMatrix::SetNodeStencils( int *stencilnumber ) {
	// Associates the given stencil number with each row
	if(fStencilNumbers) delete fStencilNumbers;
	fStencilNumbers = new int[fRows];
	memcpy(fStencilNumbers,stencilnumber,sizeof(int)*fRows);
}


// ****************************************************************************
// 
// Subclass MPStencil: Constructors and the destructor
// 
// ****************************************************************************

TPZStencilMatrix::MPStencil::MPStencil( int inc, int *IA, REAL *A ) {
	fInc = inc;
	fNumberOfItems=0;
	int *IAPtr = IA;
	while(*IAPtr) {
		fNumberOfItems += *IAPtr +1;
		IAPtr += *IAPtr +1;
	}
	fNumberOfItems++;
	fIA = new int[fNumberOfItems];
	fA = new REAL[fNumberOfItems];
	memcpy(fIA,IA,sizeof(int)*fNumberOfItems);
	memcpy(fA,A,sizeof(REAL)*fNumberOfItems);
}

TPZStencilMatrix::MPStencil::~MPStencil() {
	fNumberOfItems = 0;
	delete fIA; fIA = 0;
	delete fA; fA = 0;
}
