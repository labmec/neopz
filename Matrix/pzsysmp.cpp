/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */
/******************************************************************************
 *
 * Class realization:   TPZSYsmpMatrix
 *
 * Class type:          Derived from TPZMatrix
 *
 * Purpose:             Define operations on symmetric sparse matrices stored
 *                      in the (old) Yale Sparse Matrix Package format.
 *
 * Operations:          Mult
 *                      MultAdd
 *                      Print
 *
 * Solvers:             SOR
 *
 *****************************************************************************/

#include <memory.h>

#include "pzsysmp.h"
#include "pzfmatrix.h"

// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

TPZSYsmpMatrix::TPZSYsmpMatrix(const int rows,const int cols ) : TPZMatrix(rows,cols) {
	// Constructs an empty TPZSYsmpMatrix
	//    fRows = rows;
	//    fCols = cols;
	//    fSolver = -1;
	//    fSymmetric = 0;
	//    fMaxIterations = 4;
	//    fSORRelaxation = 1.;
	fDiag = 0;
	fIA = 0;
	fJA = 0;
	fA = 0;
#ifdef CONSTRUCTOR
	cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}

TPZSYsmpMatrix::~TPZSYsmpMatrix() {
	// Deletes everything associated with a TPZSYsmpMatrix
	delete fDiag;       fDiag = 0;
	delete fIA;			fIA = 0;
	delete fJA;			fJA = 0;
	delete fA;				fA = 0;
#ifdef DESTRUCTOR
	cerr << "~TPZSYsmpMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

const REAL & TPZSYsmpMatrix::GetVal(const int r,const int c ) const {
	// Get the matrix entry at (row,col) without bound checking
	
	// Now look through the requested row and see if there is anything
	// in column col
	int row(r),col(c);
	if ( col < row ) {
		int t = col;
		col = row;
		row = t;
	}
	col++;
	for(int ic=fIA[row]-1 ; ic < fIA[row+1]-1; ic++ ) {
		if ( fJA[ic] == col ) return fA[ic];
	}
	gZero = 0.;
	return gZero;
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

void TPZSYsmpMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y,
							 TPZFMatrix &z,
							 const REAL alpha,const REAL beta,const int opt,const int stride ) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	int  ir, ic;
	int  r = (opt) ? Cols() : Rows();
	
	// Determine how to initialize z
	if(beta != 0) {
		REAL *zp = &(z(0,0));
		const REAL *yp = &(y.g(0,0));
		REAL *zlast = zp+r*stride;
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
	
	// Compute alpha * A * x
	int jc;
	for(ir=0; ir<Rows(); ir++) {
		REAL xi = x.g(ir,0);
		for(ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
			jc = fJA[ic] - 1;
			z(ir*stride,0) += alpha * fA[ic] * x.g(jc*stride,0);
			if ( jc != ir ) {
				z(jc*stride,0) += alpha * fA[ic] * xi;
			}
		}
	}
	
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

void TPZSYsmpMatrix::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
	if(form != EInputFormat) {
		out << "\nTSYsmpMatrix Print: " << title << '\n'
		<< "\tRows    = " << Rows()  << '\n'
		<< "\tColumns = " << Cols() << '\n';
		int i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<=Rows(); i++) {
			std::cout << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(i=Rows()+1; i<fIA[Rows()]-1; i++) {
			std::cout << i      << "\t\t"
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
	} else {
		TPZMatrix::Print(title,out,form);
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************


void TPZSYsmpMatrix::ComputeDiagonal() {
	if(!fDiag) fDiag = new REAL[Rows()];
	for(int ir=0; ir<Rows(); ir++) {
		fDiag[ir] = GetVal(ir,ir);
	}
}


void TPZSYsmpMatrix::SolveSOR( int &numiterations,const TPZFMatrix &rhs, TPZFMatrix &x,
							  TPZFMatrix *residual, TPZFMatrix &scratch,
							  const REAL overrelax, REAL &tol,
							  const int FromCurrent,const int direction )  {
	
	if(!fDiag) {
		std::cout << "TPZSYsmpMatrix::SolveSOR cannot be called without diagonal\n";
		numiterations = 0;
		if(residual) {
			Residual(x,rhs,*residual);
			tol = sqrt(Norm(*residual));
		}
		return;
	}
	int irStart = 0,irLast = Rows(),irInc = 1;
	if(direction < 0) {
		irStart = Rows()-1;
		irLast = -1;
		irInc = -1;
	}
	if(!FromCurrent) x.Zero();
	REAL eqres = 2.*tol;
	REAL xnewval;
	
	int  iteration, ij, ic;
	REAL xi;
	
	for(iteration=0; iteration<numiterations && eqres >= tol; iteration++) {
		eqres = 0.;
		int ir=irStart;
		if( irInc < 0 ) {
			scratch = rhs;
			for(ir=0; ir<Rows(); ir++) {
				xi = x(ir,0);
				for(ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
					ij = fJA[ic]-1;
					if( ij != ir ) {
						scratch(ij,0) -= fA[ic] * xi;
					}
				}
			}
			for(ir=irStart; ir!=irLast; ir+=irInc) {
				xnewval = 0.0;
				for(ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
					ij = fJA[ic]-1;
					xnewval -= fA[ic] * x(ij,0);
				}
				scratch(ir,0) += xnewval;
				eqres += scratch(ir,0)*scratch(ir,0);
				x(ir,0) += overrelax*(scratch(ir,0)/fDiag[ir]);
			}
		}
		else {
			scratch = rhs;
			while(ir != irLast) {
				xnewval=0.0;
				for(ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
					xnewval -= fA[ic] * x(fJA[ic]-1,0);
				}
				scratch(ir,0) += xnewval;
				eqres += scratch(ir,0)*scratch(ir,0);
				x(ir,0) += overrelax*(scratch(ir,0)/fDiag[ir]);
				xi = x(ir,0);
				for(ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
					ij = fJA[ic]-1;
					if( ij != ir ) {
						scratch(ij,0) -= fA[ic] * xi;
					}
				}
				ir += irInc;
			}
		}
		eqres = sqrt(eqres);
	}
	tol = eqres;
	numiterations = iteration;
	if(residual) Residual(x,rhs,*residual);
}

