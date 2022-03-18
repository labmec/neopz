/**
 * @file
 * @brief Contains the implementation of the TPZTransform<> methods. 
 */

#include "pztrnsform.h"
#include "pzvec.h"

#include "fad.h"

using namespace std;

template<class T>
TPZTransform<T>::TPZTransform(int dim) :
fMult(dim,dim,0.), fSum(dim,1,0.) {
	fRow = dim;
	fCol = dim;
	fMult.Zero();
	fSum.Zero();
	int d;
	
	for(d=0; d<dim; d++) {
		fMult(d,d) = 1.;
	}
	
}

template<class T>
TPZTransform<T>::TPZTransform() :
fMult(0,0,0.), fSum(0,1,0.) {
    fRow = 0;
    fCol = 0;
}



template<class T>
TPZTransform<T>::TPZTransform(int row,int col) : fMult(row,col,0.), fSum(row,1,0.) {
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

template<class T>
TPZTransform<T>::TPZTransform(const TPZTransform<T> &t) : fMult(t.fMult),
fSum(t.fSum)
{
	fRow = t.fRow;
	fCol = t.fCol;
}

template<class T>
void TPZTransform<T>::CopyFrom(const TPZTransform<REAL> &cp)
{
    fRow = cp.fRow;
    fCol = cp.fCol;
    fMult.Resize(fRow,fCol);
    fSum.Resize(fRow,1);
    for (int i=0; i<fRow; i++) {
        fSum(i,0) = cp.fSum.GetVal(i,0);
        for (int j=0; j<fCol; j++) {
            fMult(i,j) = cp.fMult.GetVal(i,j);
        }
    }
}



template<class T>
TPZTransform<T>::~TPZTransform() {
	fRow = 0;
	fCol = 0;
	fMult.Resize(0,0);
	fSum.Resize(0,0);
}

template<class T>
TPZTransform<T> &TPZTransform<T>::operator=(const TPZTransform<T> &t) {
	fMult = t.fMult;
	fSum = t.fSum;
	fRow = t.fRow;
	fCol = t.fCol;
	return *this;
}

template<class T>
void TPZTransform<T>::SetMatrix(TPZFMatrix<T> &mult, TPZFMatrix<T> &sum) {
	fRow = mult.Rows();
	fCol = mult.Cols();

#ifdef PZDEBUG

    if ((sum.Cols()!=1)||(sum.Rows()!=mult.Rows())) {
        DebugStop();
    }
    
#endif
	fMult = mult;
	fSum = sum;
}

template<class T>
TPZTransform<T> TPZTransform<T>::Multiply(TPZTransform<T> &right) {
	TPZTransform<T> res(fRow,right.fCol);
	fMult.Multiply(right.fMult,res.fMult);
	fMult.Multiply(right.fSum,res.fSum);
	res.fSum += fSum;
	return res;
}

template<class T>
void TPZTransform<T>::Apply(const TPZVec<T> &in, TPZVec<T> &out){
#ifdef PZDEBUG
    
    
    if(fCol != in.size() || fRow != out.size())
    {
        DebugStop();
    }
#endif
	int i,j;
	for(i=0; i<fRow; i++) {
		out[i] = fSum(i,0);
		for(j=0; j<fCol; j++) {
			out[i] += fMult(i,j)*in[j];
		}
	}
}

template<class T>
void TPZTransform<T>::PrintInputForm(ostream &out) {
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
template<class T>
int TPZTransform<T>::CompareTransform(TPZTransform<T> &t,REAL tol){
	
	if(fCol != t.fCol || fRow != t.fRow)
		return 1;
	int i,j;
	for(i=0;i<fRow;i++){
        T check = fSum(i,0) - t.fSum(i,0);
		if(fabs(check) > tol) return 1;
		for(j=0;j<fCol;j++){
            T check = fMult(i,j) - t.fMult(i,j);
			if(fabs(check) > tol) return 1;
		}
	}
	return 0;
}

template class TPZTransform<REAL>;
template class TPZTransform<Fad<REAL> >;
