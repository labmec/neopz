#include "pztrnsform.h"
#include "pzvec.h"

TPZTransform::TPZTransform(int dim) :
  fMult(dim,dim,fStore,9), fSum(dim,1,fStore+9,3) {
  
  fRow = dim;
  fCol = dim;
  fMult.Zero();
  fSum.Zero();
  int d;
  
  for(d=0; d<dim; d++) {
    fMult(d,d) = 1.;
  }
  
}

TPZTransform::TPZTransform() :
  fMult(), fSum() {
}

 

TPZTransform::TPZTransform(int row,int col) : fMult(row,col,fStore,9)
     ,fSum(row,1,fStore+9,3) {
  fRow = row;
  fCol = col;
  fMult.Zero();
  fSum.Zero();
  int d;
  if (fRow == fCol) {
    for(d=0; d<fRow; d++) {
      fMult(d,d) = 1.;
    }
  }
}

TPZTransform::TPZTransform(const TPZTransform &t) : fMult(t.fRow,t.fCol,fStore,9),
       fSum(t.fRow,1,fStore+9,3) {
  fRow = t.fRow;
  fCol = t.fCol;
  fMult = t.fMult;
  fSum = t.fSum;
}

TPZTransform::~TPZTransform() {
	fRow = 0;
   fCol = 0;
   fMult.Resize(0,0);
   fSum.Resize(0,0);
}

TPZTransform &TPZTransform::operator=(const TPZTransform &t) {
  fMult = t.fMult;
  fSum = t.fSum;
  fRow = t.fRow;
  fCol = t.fCol;
  return *this;
}

void TPZTransform::SetMatrix(TPZFMatrix &mult, TPZFMatrix &sum) {
  fRow = mult.Rows();
  fCol = mult.Cols();
  fMult = mult;
  fSum = sum;
}

TPZTransform TPZTransform::Multiply(TPZTransform &right) {
  TPZTransform res(fRow,right.fCol);
  fMult.Multiply(right.fMult,res.fMult);
  fMult.Multiply(right.fSum,res.fSum);
  res.fSum += fSum;
  return res;
}

void TPZTransform::Apply(TPZVec<REAL> &in, TPZVec<REAL> &out){
  
  int i,j;
  for(i=0; i<fRow; i++) {
    out[i] = fSum(i,0);
    for(j=0; j<fCol; j++) {
      out[i] += fMult(i,j)*in[j];
    }
  }
}

void TPZTransform::PrintInputForm(ostream &out) {
	int i,j;
	out << "{";
	for(j=0; j<3; j++) {
		if(j) out << ',';
		out << "{";
		for(i=0; i<3; i++) {
			if(i) out << ',';
			if(i<fRow && j < fCol) out << fMult(i,j);
			else out << -99;
		}
		out << '}';
	}
	out << ",{";
	for(i=0; i<3; i++) {
		if(i) out << ',';
		if(i<fRow) out << fSum(i,0);
		else out << -99;
	}
	out << "}}";
}
#include <math.h>
int TPZTransform::Compare(TPZTransform &t,REAL tol){

	if(fCol != t.fCol || fRow != t.fRow)
		return 1;
	int i,j;
	for(i=0;i<fRow;i++){
		if(fabs(fSum(i,0) - t.fSum(i,0)) > tol) return 1;
		for(j=0;j<fCol;j++){
			if(fabs(fMult(i,j) - t.fMult(i,j)) > tol) return 1;
		}
	}
	return 0;
}
