
//METHODS DEFINITIONS FOR CLASS MATRIX


#include "pzshtmat.h"
//#include "pzerror.h"

template <class TObj>
TPZGenMatrix<TObj>::TPZGenMatrix(){
  fRows = 0;
  fCols = 0;
  fMem = NULL;
}

template <class TObj>
TPZGenMatrix<TObj>::TPZGenMatrix(int Rows, int columns) {
  fRows = Rows;
  fCols = columns;
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
    fCols = A.fCols;
    TObj *f = fMem,*l = f+naloc,*fa = A.fMem;
    while(f<l) *f++ = *fa++;
  } else {
    fRows = 0;
    fCols = 0;
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
  fCols = newcol;
}

template <class TObj>
TPZGenMatrix<TObj>::~TPZGenMatrix() {
  if(fMem) delete []fMem;
  fMem = 0;
  fRows = 0;
  fCols = 0;
}

template <class TObj>
TPZGenMatrix<TObj>& TPZGenMatrix<TObj>::operator= (const TPZGenMatrix<TObj> & rval) {

  if(this == &rval) return *this;
  if(fMem) delete fMem;
  fRows = rval.fRows;
  fCols = rval.fCols;
  fMem = new int[fRows*fCols];

  if(fMem) {
    TObj *pbeg,*pfinal;
    pfinal = &(fMem[fRows*fCols]);
    TObj *rvalm = rval.fMem;
    for (pbeg = fMem; pbeg<pfinal; pbeg++,rvalm++) *pbeg = *rvalm;
  } else {
    cout << "TPZGenMatrix<TObj>.= -->Cannot allocate matrix storage\n";
  }
  return (*this);
}


template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator+ (const TPZGenAMatrix<TObj> & rval) const {
  if ( (fRows != rval.fRows) || (fCols != rval.fCols) )
    cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";


  TPZGenAMatrix<TObj> sum(Rows(), Cols());
  TObj *pbeg1 = fMem, *pfin;
  TObj *pbeg2 = sum.fMem, *pbeg3 = rval.fMem;
  pfin = pbeg1 + (fRows*fCols);
  for ( ; pbeg1<pfin; pbeg1++, pbeg2++, pbeg3++)
    *pbeg2 = *pbeg1 + *pbeg3;
  return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator+ (const TObj x) const {

  TPZGenAMatrix<TObj> sum(Rows(), Cols());

  for (int i=0; i<fRows*fCols; i++)
    sum.fMem[i] = fMem[i] + x;
  return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- () const {

  TPZGenAMatrix<TObj> unaryminus(Rows(), Cols());

  for (int i=0; i<fRows*fCols; i++)
    unaryminus.fMem[i] = - fMem[i];
  return unaryminus;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator* (const TPZGenAMatrix<TObj> & rval) const {

  if (fCols != rval.fRows)
    cout << "ERROR-> unsuitable TPZGenMatrix<TObj> size for multiplication\n";
  cout.flush();

  TPZGenAMatrix<TObj> result(Rows(), rval.Cols());

  TObj *ptr = result.fMem;
  for (int i=0; i<(fRows*fCols); i+=fCols)
    for (int j=0; j<rval.fCols; j++, ptr++)
      for (int k=0; k<fCols; k++)
	*ptr += fMem[i+k] * rval.fMem[j+k*rval.fCols];

  return result;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator* (const TObj x) const {

  TPZGenAMatrix<TObj> result(Rows(), Cols());

  for (int i=0; i<fRows*fCols; i++)
    result.fMem[i] = fMem[i] * x;
  return result;
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator+=(const TPZGenAMatrix<TObj> & rval) {

  if ( (fRows != rval.fRows) || (fCols != rval.fCols) )
    cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";
  cout.flush();

  int *pbeg1 = fMem, *pfin;
  int *pbeg3 = rval.fMem;
  pfin = pbeg1 + (fRows*fCols);
  for ( ; pbeg1<pfin; pbeg1++, pbeg3++)
    *pbeg1 += *pbeg3;
  return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator+=(const TObj x) {

  for (int i=0; i<fRows*fCols; i++)
    fMem[i] += x;
  return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator-=(const TPZGenAMatrix<TObj> & rval) {

  if ( (fRows != rval.fRows) || (fCols != rval.fCols) )
    cout << "ERROR-> different TPZGenMatrix<TObj> size for addition";
  cout.flush();

  int *pbeg1 = fMem, *pfin;
  int *pbeg3 = rval.fMem;
  pfin = pbeg1 + (fRows*fCols);
  for ( ; pbeg1<pfin; pbeg1++, pbeg3++)
    *pbeg1 -= *pbeg3;
  return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator-=(const TObj x) {

  for (int i=0; i<fRows*fCols; i++)
    fMem[i] -= x;
  return (*this);
}

template <class TObj>
TPZGenAMatrix<TObj>& TPZGenAMatrix<TObj>::operator*=(const TObj x) {

  for (int i=0; i<fRows*fCols; i++)
    fMem[i] *= x;
  return (*this);
}

template <class TObj>
TObj & TPZGenMatrix<TObj>::operator()(const int i,const int j) const{

  if (i>=0 && i<fRows && j>=0 && j<fCols)
    return fMem[i*fCols+j];
  else {
    cout << "ERROR-> TPZGenMatrix<TObj> index out of range\n"
      " i = " << i << " j = " << j << " Rows = " <<
      fRows << " columns = " << fCols << "\n";
    cout.flush();
    return fMem[0];
  }
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::Transpose () const {

  TPZGenAMatrix<TObj> transp(Cols(),Rows());

  TObj *ptr  = fMem;

  for (int i=0; i < fRows; i++) {		//loop over columns of Transpose
    TObj *trp = transp.fMem + i;
    for (int j=0; j<fCols; j++, ptr++, trp += fRows) {
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
	out << "TPZGenMatrix<TObj>  Rows = " << fRows << " columns = " << fCols << "\n";
	for (int i=0; i<fRows; i++) {
		out << "\n row " << i;
		for (int j=0; j<fCols; j++) {
			if ( !(j%6) ) out << "\n";
			out << "  " << fMem[(i*fCols)+j];
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
  out << "TPZGenMatrix<TObj>  Rows = " << fRows << " columns = " << fCols << endl;
  for (int i=0; i<fRows; i++) {
    out << "\n row " << i;
    for (int j=0; j<fCols; j++) {
      if ( !(j%6) ) out << "\n";
      out << "  " << fMem[(i*fCols)+j];
    }
  }
  out << "\n";
  return;
}


template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- (const TPZGenAMatrix<TObj> & rval) const {
  if ( (fRows != rval.fRows) || (fCols != rval.fCols) )
    cout << "ERROR-> different TPZGenMatrix<TObj> size for subtraction";
  cout.flush();


  TPZGenAMatrix<TObj> sum(Rows(), Cols());

  int *pbeg1 = fMem, *pfin;
  int *pbeg2 = sum.fMem, *pbeg3 = rval.fMem;
  pfin = pbeg1 + (fRows*fCols);
  for ( ; pbeg1<pfin; pbeg1++, pbeg2++, pbeg3++)
    *pbeg2 = *pbeg1 - *pbeg3;
  return sum;
}

template <class TObj>
TPZGenAMatrix<TObj> TPZGenAMatrix<TObj>::operator- (const TObj x) const {


  TPZGenAMatrix<TObj> sum(Rows(), Cols());

  for (int i=0; i<fRows*fCols; i++)
    sum.fMem[i] = fMem[i] - x;
  return sum;
}

template class TPZGenMatrix<int>;
template class TPZGenAMatrix<int>;
