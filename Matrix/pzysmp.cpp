#include <memory.h>
#include <string.h>

#include "pzysmp.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#ifdef USING_BLAS
	double cblas_ddoti(const int N, const double *X, const int *indx,
                   const double *Y);
#endif


// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************

void TPZFYsmpMatrix::InitializeData(){}
void TPZFYsmpMatrix::Multiply(TPZFYsmpMatrix & B, TPZFYsmpMatrix & Res){
    int i,j,k;
    if (B.Rows()!=Rows()) return;
    int rows = Rows();
    double aux=0.;
    for(i=0;i<rows;i++){
        for(j=0;j<rows;j++){
            for(k=0;k<rows;k++){
	      // C[i][j] += A[i][k]*B[k][j];
	      aux+=GetVal(i,k)*B.GetVal(i,k);
	    }
	    Res.PutVal(i,j,aux);
	    aux=0.;
	}
    }
}
void TPZFYsmpMatrix::PutVal(const int row, const int col, double Value){
    int k;
    int flag=0;
    for(k=fIA[row];k<fIA[row+1];k++){
      if(fJA[k]==col){
        flag=1;
	    fA[k]=Value;
    	break;
      }
    }
    if(!flag) cout << "TPZFYsmpMatrix::PutVal: Non existing position on sparse matrix: line = " << row << " column " << col << endl;
}
void TPZFYsmpMatrix::AddKel(TPZFMatrix & elmat, TPZVec<int> & destinationindex){
    int i,j,k;
    REAL value=0.;
    int ipos,jpos;
    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(i,j);
            //cout << "j= " << j << endl;
            if(value){
                //cout << "fIA[ipos] " << fIA[ipos] << "     fIA[ipos+1] " << fIA[ipos+1] << endl;
                int flag = 0;
                for(k=fIA[ipos];k<fIA[ipos+1];k++){
                    if(fJA[k]==jpos || fJA[k]==-1){
                        //cout << "fJA[k] " << fJA[k] << " jpos "<< jpos << "   " << value << endl;
                        //cout << "k " << k << "   "<< jpos << "   " << value << endl;
                        flag=1;
                        if(fJA[k]==-1){
                            fJA[k]=jpos;
                            fA[k]=value;
                           // cout << jpos << "   " << value << endl;
                            break;
                        }else{
                            fA[k]+=value;
                            break;
                        }
                    }
                }
                if(!flag) cout << "TPZFYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << endl;         }
        }
    }
}
void TPZFYsmpMatrix::AddKelOld(TPZFMatrix & elmat, TPZVec < int > & destinationindex){
  int i=0;
  int j=0;
  int ilocal=0;
  int jlocal=0;
  int nel = destinationindex.NElements();
  for(i=0;i<nel;i++){
    ilocal = destinationindex[i];
    for(j=0;j<nel;j++){
      jlocal = destinationindex[j];
      //Element(ilocal,jlocal)+=elmat(i,j);
    }
  }
  
}
TPZFYsmpMatrix::TPZFYsmpMatrix(const int rows,const int cols ) : TPZMatrix(rows,cols) {
  // Constructs an empty TPZFYsmpMatrix
  //    fSolver = -1;
  fSymmetric = 0;
  //    fMaxIterations = 4;
  //    fSORRelaxation = 1.;
  fDiag = 0;
#ifdef CONSTRUCTOR
  cerr << "TPZFYsmpMatrix(int rows,int cols)\n";
#endif
}

TPZFYsmpMatrix::~TPZFYsmpMatrix() {
  // Deletes everything associated with a TPZFYsmpMatrix
  delete []fDiag;
  fDiag = 0;
  delete []fIA;
  fIA=0;
  delete []fJA;
  fJA=0;
  delete []fA;
  fA=0;
  //    fSolver = -1;
#ifdef DESTRUCTOR
  cerr << "~TPZFYsmpMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

const REAL & TPZFYsmpMatrix::GetVal(const int row,const int col ) const {
  // Get the matrix entry at (row,col) without bound checking

  // Now look through the requested row and see if there is anything
  // in column col
/*  int loccol = col+1;
  for(int ic=fIA[row]-1 ; ic < fIA[row+1]-1; ic++ ) {
    if ( fJA[ic] == loccol ) return fA[ic];
  }*/
  int loccol = col;
  for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
    if ( fJA[ic] == loccol && fJA[ic] != -1 ) return fA[ic];
  }
  return gZero;
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

void TPZFYsmpMatrix::MultAdd(const TPZFMatrix &x,const TPZFMatrix &y,
			     TPZFMatrix &z,
			     const REAL alpha,const REAL beta,const int opt,const int stride ) const {
  // computes z = beta * y + alpha * opt(this)*x
  //          z and x cannot share storage
  int  ir, ic, icol, xcols;
  xcols = x.Cols();
  REAL sum;
  int  r = (opt) ? Cols() : Rows();

  // Determine how to initialize z
  for(ic=0; ic<xcols; ic++) {
    REAL *zp = &(z(0,ic));
    if(beta != 0) {
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
  }
  // Compute alpha * A * x
  for(ic=0; ic<xcols; ic++) {
    if(opt == 0) {

      for(ir=0; ir<Rows(); ir++) {
#ifndef USING_BLAS
	for(sum = 0.0, icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
	  if(fJA[icol]==-1) break;//Checa a existência de dado ou não
	  sum += fA[icol] * x.g((fJA[icol])*stride,ic);
	}
#endif
#ifdef USING_BLAS
	int size = fIA[ir+1] - fIA[ir];
	sum = cblas_ddoti(size, &fA[fIA[ir]], &fJA[fIA[ir]],
                   &x.g(0,ic));
#endif
	z(ir*stride,ic) += alpha * sum;
      }
    }

  // Compute alpha * A^T * x
    else {

      int jc;
      int icol;
      for(ir=0; ir<Rows(); ir++) {
		for(icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
		  if(fJA[icol]==-1) break; //Checa a existência de dado ou não
		  jc = fJA[icol];
		  z(jc*stride,ic) += alpha * fA[icol] * x.g(jc*stride,ic);
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

void TPZFYsmpMatrix::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
  // Print the matrix along with a identification title
  if(form != EInputFormat) {
    out << "\nTFYsmpMatrix Print: " << title << '\n'
	<< "\tRows    = " << Rows()  << '\n'
	<< "\tColumns = " << Cols() << '\n';
    int i;
    out << "\tIA\tJA\tA\n"
	<< "\t--\t--\t-\n";
    for(i=0; i<=Rows(); i++) {
      out << i      << '\t'
	   << fIA[i] << '\t'
	   << fJA[i] << '\t'
	   << fA[i]  << '\n';
    }
    for(i=Rows(); i<fIA[Rows()]; i++) {
      out << i      << "\t\t"
	   << fJA[i] << '\t'
	   << fA[i]  << '\n';
    }
  }
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************


void TPZFYsmpMatrix::ComputeDiagonal() {
  if(fDiag) return;
  int rows = Rows();
  fDiag = new double [rows];
  for(int ir=0; ir<rows; ir++) {
    fDiag[ir] = GetVal(ir,ir);
  }
}


void TPZFYsmpMatrix::SolveSOR( int &numiterations, const TPZFMatrix &rhs, TPZFMatrix &x,
			       TPZFMatrix *residual, TPZFMatrix &/*scratch*/,
			       const REAL overrelax, REAL &tol,
			       const int FromCurrent,const int direction ) const {

  //    if(!fDiag) ComputeDiagonal();
  int irStart = 0,irLast = Rows(),irInc = 1;
  if(direction < 0) {
    irStart = Rows()-1;
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
      REAL xnewval=rhs.g(ir,0);
      for(int ic=fIA[ir]-1; ic<fIA[ir+1]-1; ic++) {
	xnewval -= fA[ic] * x(fJA[ic]-1,0);
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

int TPZFYsmpMatrix::Zero()
{
   int size = fIA[fRow] * sizeof(REAL);
   int diagSize = min(fRow, fCol) * sizeof(REAL);
   memset(fA,'\0',size);
   memset(fDiag,'\0', diagSize);
   return 1;
}
