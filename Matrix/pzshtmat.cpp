/**
 * @file
 * @brief Contains the implementation of the TPZGenMatrix methods.
 */

#include "pzshtmat.h"
//#include "pzerror.h"
using namespace std;

template <class TObj>
TPZGenMatrix<TObj>::TPZGenMatrix(){
	fRows = 0;
	this->fCols = 0;
	fMem = NULL;
}

template <class TObj>
TPZGenMatrix<TObj>::TPZGenMatrix(int Rows, int columns) {
	fRows = Rows;
	this->fCols = columns;
	fMem = new int[Rows*columns];
	
	if(fMem){
		int *pbeg,*pfinal;
		pfinal = &(fMem[Rows*columns]);
		for (pbeg = fMem; pbeg<pfinal; *pbeg++=0) ;
	} else {
		cout << "TPZGenMatrix<TObj>.ct-->Cannot create matrix structure\n";
	}
}


template <class TObj>
TPZGenMatrix<TObj>::TPZGenMatrix(const TPZGenMatrix<TObj> & A) {
	
	//***** WARNING *****
	// matrices created with copy initializer are always temporary, eg, they
	// share the same storage with another matrix
	int naloc = A.fRows*A.fCols;
	fMem = new TObj[naloc];
	if(fMem) {
		fRows = A.fRows;
		this->fCols = A.fCols;
		TObj *f = fMem,*l = f+naloc,*fa = A.fMem;
		while(f<l) *f++ = *fa++;
	} else {
		fRows = 0;
		this->fCols = 0;
	}
}

template <class TObj>
void TPZGenMatrix<TObj>::Resize(const int newrow, const int newcol) {
	TObj *sht = new TObj[newrow*newcol];
	if(!sht) return;
	int minrow = (newrow < Rows()) ? newrow : Rows();
	int mincol = (newcol < Cols()) ? newcol : Cols();
	for(int i=0; i<minrow; i++) {
		for(int j=0; j<mincol; j++) {
			sht[i*newcol+j] = (*this)(i,j);
		}
	}
	delete []fMem;
	fMem = sht;
	fRows = newrow;
	this->fCols = newcol;
}

template <class TObj>
TPZGenMatrix<TObj>::~TPZGenMatrix() {
	if(fMem) delete []fMem;
	fMem = 0;
	fRows = 0;
	this->fCols = 0;
}

template <class TObj>
TPZGenMatrix<TObj>& TPZGenMatrix<TObj>::operator= (const TPZGenMatrix<TObj> & rval) {
	
	if(this == &rval) return *this;
	if(fMem) delete fMem;
	fRows = rval.fRows;
	this->fCols = rval.fCols;
	fMem = new int[fRows*this->fCols];
	
	if(fMem) {
		TObj *pbeg,*pfinal;
		pfinal = &(fMem[fRows*this->fCols]);
		TObj *rvalm = rval.fMem;
		for (pbeg = fMem; pbeg<pfinal; pbeg++,rvalm++) *pbeg = *rvalm;
	} else {
		cout << "TPZGenMatrix<TObj>.= -->Cannot allocate matrix storage\n";
	}
	return (*this);
}


template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator+ (const TPZGenAMatrix<TObj> & rval) const {
	if ( (this->fRows != rval.fRows) || (this->fCols != rval.fCols) )
		cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";
	
	
	TPZGenAMatrix<TObj> sum(this->Rows(), this->Cols());
	TObj *pbeg1 = this->fMem, *pfin;
	TObj *pbeg2 = sum.fMem, *pbeg3 = rval.fMem;
	pfin = pbeg1 + (this->fRows*this->fCols);
	for ( ; pbeg1<pfin; pbeg1++, pbeg2++, pbeg3++)
		*pbeg2 = *pbeg1 + *pbeg3;
	return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator+ (const TObj x) const {
	
	TPZGenAMatrix<TObj> sum(this->Rows(), this->Cols());
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		sum.fMem[i] = this->fMem[i] + x;
	return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- () const {
	
	TPZGenAMatrix<TObj> unaryminus(this->Rows(), this->Cols());
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		unaryminus.fMem[i] = - this->fMem[i];
	return unaryminus;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator* (const TPZGenAMatrix<TObj> & rval) const {
	
	if (this->fCols != rval.fRows)
		cout << "ERROR-> unsuitable TPZGenMatrix<TObj> size for multiplication\n";
	cout.flush();
	
	TPZGenAMatrix<TObj> result(this->Rows(), rval.Cols());
	
	TObj *ptr = result.fMem;
	for (int i=0; i<(this->fRows*this->fCols); i+=this->fCols)
		for (int j=0; j<rval.fCols; j++, ptr++)
			for (int k=0; k<this->fCols; k++)
				*ptr += this->fMem[i+k] * rval.fMem[j+k*rval.fCols];
	
	return result;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator* (const TObj x) const {
	
	TPZGenAMatrix<TObj> result(this->Rows(), this->Cols());
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		result.fMem[i] = this->fMem[i] * x;
	return result;
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator+=(const TPZGenAMatrix<TObj> & rval) {
	
	if ( (this->fRows != rval.fRows) || (this->fCols != rval.fCols) )
		cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";
	cout.flush();
	
	int *pbeg1 = this->fMem, *pfin;
	int *pbeg3 = rval.fMem;
	pfin = pbeg1 + (this->fRows*this->fCols);
	for ( ; pbeg1<pfin; pbeg1++, pbeg3++)
		*pbeg1 += *pbeg3;
	return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator+=(const TObj x) {
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		this->fMem[i] += x;
	return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator-=(const TPZGenAMatrix<TObj> & rval) {
	
	if ( (this->fRows != rval.fRows) || (this->fCols != rval.fCols) )
		cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";
	cout.flush();
	
	int *pbeg1 = this->fMem, *pfin;
	int *pbeg3 = rval.fMem;
	pfin = pbeg1 + (this->fRows*this->fCols);
	for ( ; pbeg1<pfin; pbeg1++, pbeg3++)
		*pbeg1 -= *pbeg3;
	return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator-=(const TObj x) {
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		this->fMem[i] -= x;
	return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator*=(const TObj x) {
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		this->fMem[i] *= x;
	return (*this);
}

template <class TObj>
TObj & TPZGenMatrix<TObj>::operator()(const int i,const int j) const{
	
	if (i>=0 && i<this->fRows && j>=0 && j<this->fCols)
		return fMem[i*this->fCols+j];
	else {
		cout << "ERROR-> TPZGenMatrix<TObj> index out of range\n"
		" i = " << i << " j = " << j << " Rows = " <<
		this->fRows << " columns = " << this->fCols << "\n";
		cout.flush();
		return fMem[0];
	}
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::Transpose () const {
	
	TPZGenAMatrix<TObj> transp(this->Cols(),this->Rows());
	
	TObj *ptr  = this->fMem;
	
	for (int i=0; i < this->fRows; i++) {		//loop over columns of Transpose
		TObj *trp = transp.fMem + i;
		for (int j=0; j<this->fCols; j++, ptr++, trp += this->fRows) {
			*trp = *ptr;
		}
	}
	return transp;
}
/*
 void TPZGenMatrix<TObj>::Print (ostream & out) {
 
 if(fMem == NULL) {
 cout << "NULL TPZGenMatrix<TObj>\n";
 return;
 }
 out << "TPZGenMatrix<TObj>  Rows = " << this->fRows << " columns = " << this->fCols << "\n";
 for (int i=0; i<this->fRows; i++) {
 out << "\n row " << i;
 for (int j=0; j<this->fCols; j++) {
 if ( !(j%6) ) out << "\n";
 out << "  " << fMem[(i*this->fCols)+j];
 }
 }
 out << "\n";
 return;
 }
 
 */
template <class TObj>
void TPZGenMatrix<TObj>::Print (const char *c, ostream & out) const {
	
	if(fMem == NULL) {
		cout << "NULL TPZGenMatrix<TObj>\n";
		return;
	}
	out << c << endl;
	out << "TPZGenMatrix<TObj>  Rows = " << this->fRows << " columns = " << this->fCols << endl;
	for (int i=0; i<this->fRows; i++) {
		out << "\n row " << i;
		for (int j=0; j<this->fCols; j++) {
			if ( !(j%6) ) out << "\n";
			out << "  " << fMem[(i*this->fCols)+j];
		}
	}
	out << "\n";
	return;
}


template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- (const TPZGenAMatrix<TObj> & rval) const {
	if ( (this->fRows != rval.fRows) || (this->fCols != rval.fCols) )
		cout << "ERROR-> different TPZGenMatrix<TObj> size for subtraction";
	cout.flush();
	
	
	TPZGenAMatrix<TObj> sum(this->Rows(), this->Cols());
	
	int *pbeg1 = this->fMem, *pfin;
	int *pbeg2 = sum.fMem, *pbeg3 = rval.fMem;
	pfin = pbeg1 + (this->fRows*this->fCols);
	for ( ; pbeg1<pfin; pbeg1++, pbeg2++, pbeg3++)
		*pbeg2 = *pbeg1 - *pbeg3;
	return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- (const TObj x) const {
	
	
	TPZGenAMatrix<TObj> sum(this->Rows(), this->Cols());
	
	for (int i=0; i<this->fRows*this->fCols; i++)
		sum.fMem[i] = this->fMem[i] - x;
	return sum;
}

template class TPZGenMatrix<int>;
template class TPZGenAMatrix<int>;
