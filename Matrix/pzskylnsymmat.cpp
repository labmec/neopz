// ---------------------------------------------------------------------------

#include "pzskylnsymmat.h"

// ---------------------------------------------------------------------------

//
// Author: Nathan Shauer.
//
// File:   pzskylnsymmat
//
// Class:  TPZSkylNSymMatrix
//
// Obs.: This class manages non symmetric skylyne type matrix
//
//
// Versao: 12 / 2011.
//

#include <math.h>
#include <stdlib.h>
#include <random>
//const int templatedepth = 10;

/**
 * Commented out by int64_thin
 * Compilation problem under MaCOSX OS.
 * Used by method TestSpeed which was also commented out.
 */
// #include <sys/timeb.h>

#include "pzfmatrix.h"
#include "pzskylmat.h"
// #include "pzerror.h"

#include <sstream>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.matrix.tpzskylmatrix");
#endif

using namespace std;

/** **************** */
/** * TPZSkylMatrix ** */

/** ************************** PUBLIC *************************** */

/** ************************** */
/** * Construtor (int) ** */

template <class TVar>
TPZSkylNSymMatrix<TVar>::TPZSkylNSymMatrix(const int64_t row, const int64_t col) : TPZRegisterClassId(&TPZSkylNSymMatrix::ClassId),
TPZMatrix<TVar>(row,col),
fElem(row + 1), fElemb(row + 1), fStorage(0), fStorageb(0)
{

    if (row != col) {
        DebugStop();
    }
  // Inicializa a diagonal (vazia).
  fElem.Fill(0);
  fElemb.Fill(0);
}

template <class TVar>
TPZSkylNSymMatrix<TVar>::TPZSkylNSymMatrix(const int64_t dim, const TPZVec<int64_t> &skyline)
		: TPZRegisterClassId(&TPZSkylNSymMatrix::ClassId),
TPZMatrix<TVar>(dim, dim), fElem(dim + 1), fElemb(dim + 1), fStorage(0), fStorageb(0)
{

  // Inicializa a diagonal (vazia).
  fElem.Fill(0);
  fElemb.Fill(0);
  InitializeElem(skyline, fStorage, fElem);
  InitializeElem(skyline, fStorageb, fElemb);
}

template<class TVar>
TPZSkylNSymMatrix<TVar> &
TPZSkylNSymMatrix<TVar>::operator=(const TPZSkylNSymMatrix<TVar> &A )
{
	Clear();
	Copy( A );
	return( *this );
}

/* IMPLEMENTAR
void TPZSkylNSymMatrix::AddSameStruct(TPZSkylMatrix &B, double k)
{
#ifdef PZDEBUG
{
int size = this->fElem.NElements();
if (size != B.fElem.NElements())
{
PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
PZError.flush();
DebugStop();
}
for (int i = 0; i < size; i++)
{
if ((this->fElem[i]-this->fElem[0]) != (B.fElem[i] - B.fElem[0]))
{
PZError << "\nFATAL ERROR at " << __PRETTY_FUNCTION__ << "\n";
PZError.flush();
DebugStop();
}
}
}
#endif

const int n = this->fStorage.NElements();
for (int i = 0; i < n; i++)
this->fStorage[i] += k * B.fStorage[i];

}
 */

template <class TVar>
void TPZSkylNSymMatrix<TVar>::SetSkyline(const TPZVec<int64_t> &skyline)
{
  fElem.Fill(0);
  fElemb.Fill(0);
  InitializeElem(skyline, fStorage, fElem);
  InitializeElem(skyline, fStorageb, fElemb);
}

template <class TVar>
int64_t TPZSkylNSymMatrix<TVar>::NumElements(const TPZVec<int64_t> &skyline)
{
  int64_t dim = skyline.NElements();
  int64_t i, nelem = 0;
  for (i = 0; i < dim; i++)
    nelem += i - skyline[i] + 1;
  return nelem;
}

template <class TVar>
void TPZSkylNSymMatrix<TVar>::InitializeElem(const TPZVec<int64_t> &skyline,
	TPZManVector<TVar> &storage, TPZVec<TVar *> &point)
{
	int64_t dim = skyline.NElements();
	int64_t nel = NumElements(skyline);
	storage.Resize(nel);
	storage.Fill(0.);
	int64_t i;
	point.Resize(dim + 1);
	if (dim)
	{
		point[0] = &storage[0];
		point[dim] = &storage[0] + nel;
	}
	else
	{
		point[0] = 0;
	}
	for (i = 1; i < dim + 1; i++)
		point[i] = point[i - 1] + (i - 1) - skyline[i - 1] + 1;
}

/**
Computes the highest skyline of both objects
 */
template <class TVar>
void TPZSkylNSymMatrix<TVar>::ComputeMaxSkyline(const TPZSkylNSymMatrix &first,
  const TPZSkylNSymMatrix &second, TPZVec<int64_t> &res)
{

  if (first.Rows() != second.Rows())
  {
    cout << "ComputeMaxSkyline : incompatible dimension";
    return;
  }
  int64_t i, dim = first.Rows();
  res.Resize(dim + 1);

  for (i = 1; i < dim + 1; i++)
  {

    int64_t aux = (first.Size(i) > second.Size(i)) ? first.Size(i) : second.Size(i);
    res[i] = i - aux - 1;
  }
}

template <class TVar>
TVar & TPZSkylNSymMatrix<TVar>::operator()(const int64_t r, const int64_t c)
{
  int64_t row(r), col(c);
  if (col >= row)
  {
    // Indice do vetor coluna.
    int64_t index = col - row;
    if (index >= Size(col))
    {
      // Error("TPZSkylMatrix::operator()","Index out of range");
      TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Index out of range");
      DebugStop();
    }
    return fElem[col][index];
  }
  else
  {
    // Indice do vetor coluna.
    int64_t index = row - col;
    if (index >= Size(row))
    {
      // Error("TPZSkylMatrix::operator()","Index out of range");
      TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Index out of range");
      DebugStop();
    }
    return fElemb[row][index];
  }
}

template <class TVar>
TVar &TPZSkylNSymMatrix<TVar>::s(const int64_t row, const int64_t col)
{
  return operator()(row, col);
}

template <class TVar>
TVar & TPZSkylNSymMatrix<TVar>::operator()(const int64_t r)
{
  return operator()(r, r);
}

/** *********** */
/** * PutVal ** */
template <class TVar>
int TPZSkylNSymMatrix<TVar>::PutVal(const int64_t r, const int64_t c, const TVar & value)
{
  // inicializando row e col para trabalhar com a triangular superior
  int64_t row(r), col(c);
  if (col >= row)
  {
    // Indice do vetor coluna.
    int64_t index = col - row;
    // Se precisar redimensionar o vetor.
    if (index >= Size(col) && !IsZero(value))
    {
      cout << "TPZSkylMatrix::PutVal Size" << Size(col);
      cout.flush();
      TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Index out of range");
    }
    else if (index >= Size(col))
      return 1;
    fElem[col][index] = value;
    // delete[]newVet;
  }

  else if (col < row)
  {
    // Indice do vetor coluna.
    int64_t index = row - col;
    // Se precisar redimensionar o vetor.
    if (index >= Size(row) && !IsZero(value))
    {
      cout << "TPZSkylMatrix::PutVal Size" << Size(col);
      cout.flush();
      TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "Index out of range");
    }
    else if (index >= Size(row))
      return 1;
    fElemb[row][index] = value;
    // delete[]newVet;
  }
	this->fDecomposed = 0;
	
  return(1);
}

template <class TVar>
void TPZSkylNSymMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x, const TPZFMatrix<TVar> &y,
  TPZFMatrix<TVar> &z, const TVar alpha, const TVar beta, const int opt)const
{
  // Computes z = beta * y + alpha * opt(this)*x
  // z and x cannot overlap in memory
  if ((!opt && this->Cols() != x.Rows()) || this->Rows() != x.Rows())
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,
    " <matrixs with incompatible dimensions>");
  if (z.Rows() != x.Rows() || z.Cols() != x.Cols())
    z.Redim(x.Rows(), x.Cols());
  if (x.Cols() != y.Cols() || x.Cols() != z.Cols() || x.Rows() != y.Rows()
    || x.Rows() != z.Rows())
  {
    cout << "x.Cols = " << x.Cols() << " y.Cols()" << y.Cols()
        << " z.Cols() " << z.Cols() << " x.Rows() " << x.Rows()
        << " y.Rows() " << y.Rows() << " z.Rows() " << z.Rows() << endl;
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, " incompatible dimensions\n");
  }
  this->PrepareZ(y, z, beta, opt);
  int64_t rows = this->Rows();
  int64_t xcols = x.Cols();
  int64_t ic, r;
  for (ic = 0; ic < xcols; ic++)
  {
    for (r = 0; r < rows; r++)
    {
      int64_t offset = Size(r);
      TVar val = 0.;
      const TVar *p = &x.g((r - offset + 1), ic);
        TVar *diag, *diaglast;
        if (opt == 0)
        {
            diag = fElemb[r] + offset - 1;
            diaglast = fElemb[r];
        }
        else
        {
            diag = fElem[r] + offset - 1;
            diaglast = fElem[r];
        }
      while (diag > diaglast)
      {
        val += *diag--**p;
        p ++;
      }
      if (diag == diaglast)
      {
        diag = fElem[r];
        val += *diag * *p;
      }
      z(r, ic) += val * alpha;
      TVar *zp = &z((r - offset + 1), ic);
      val = x.g(r, ic);
        if(opt == 0)
        {
          diag = fElem[r] + offset - 1;
          diaglast = fElem[r];
        }
        else
        {
            diag = fElemb[r] + offset - 1;
            diaglast = fElemb[r];
        }
      while (diag > diaglast)
      {
        *zp += alpha * *diag--*val;
        zp ++;
      }
      //z.Print("z");
    }
  }
}

/** @brief Updates the values of the matrix based on the values of the matrix */
template <class TVar>
void TPZSkylNSymMatrix<TVar>::UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat)

{
    TPZMatrix<TVar> *matrix = mat.operator->();
    TPZSkylNSymMatrix<TVar> *skylmat = dynamic_cast<TPZSkylNSymMatrix<TVar> *>(matrix);
    if (!skylmat) {
        DebugStop();
    }
    if (fStorage.NElements() != skylmat->fStorage.NElements() || fStorageb.NElements() != skylmat->fStorageb.NElements()) {
        DebugStop();
    }
    memcpy(&fStorage[0], &(skylmat->fStorage[0]) , fStorage.NElements()*sizeof(TVar));
    memcpy(&fStorageb[0], &(skylmat->fStorageb[0]) , fStorageb.NElements()*sizeof(TVar));
    this->fDecomposed = skylmat->fDecomposed;
    this->fDefPositive = skylmat->fDefPositive;
}

/*
void TPZSkylMatrix::SolveSOR(int & numiterations, const TPZFMatrix &F,
TPZFMatrix &result, TPZFMatrix *residual, TPZFMatrix &scratch,
const REAL overrelax, REAL &tol, const int FromCurrent,
const int direction)
{

if (residual == &F)
{
cout <<
"TPZMatrix::SolveSOR called with residual and F equal, no solution\n";
return;
}
REAL res = 2 * tol + 1.; ;
if (residual)
res = Norm(*residual);
if (!FromCurrent)
{
result.Zero();
}
int r = Dim();
int c = F.Cols();
int i, ifirst = 0, ilast = r, iinc = 1;
if (direction == -1)
{
ifirst = r - 1;
ilast = 0;
iinc = -1;
}
int it;
for (it = 0; it < numiterations && res > tol; it++)
{
res = 0.;
scratch = F;
for (int ic = 0; ic < c; ic++)
{
if (direction == 1)
{
//
// compute the upper triangular part first and put into the scractch vector
//
for (i = ifirst; i != ilast; i += iinc)
{
// TPZColuna *mydiag = &fDiag[i];
int offset = Size(i);
REAL val;
REAL *diag;
REAL *diaglast = fElem[i];
REAL *scratchp = &scratch(i - offset + 1, ic);
val = result(i, ic);
diag = fElem[i] + offset - 1;
int lastid = diag - diaglast;
int id;
for (id = 0; id <= lastid; id++)
 * (scratchp + id) -= *(diag - id) * val;
// codeguard fix
while( diag >= diaglast ) *scratchp++ -= *diag-- * val;
//
}
//
// perform the SOR operation
//
for (i = ifirst; i != ilast; i += iinc)
{
// TPZColuna *mydiag = &fDiag[i];
int offset = Size(i);
REAL val = scratch(i, ic);
REAL *p = &result(i - offset + 1, ic);
REAL *diag = fElem[i] + offset - 1;
REAL *diaglast = fElem[i];
while (diag > diaglast)
val -= *diag--**p++;
res += val * val;
result(i, ic) += val * overrelax / *diag;
}
}
else
{
//
// the direction is upward
//
// put the lower triangular part of the multiplication into the scratch vector
//
for (i = ifirst; i != ilast; i += iinc)
{
// TPZColuna *mydiag = &fDiag[i];
int offset = Size(i);
REAL val = scratch(i, ic);
REAL *p = &result(i - offset + 1, ic);
REAL *diag = fElem[i] + offset - 1;
REAL *diaglast = fElem[i];
while (diag > diaglast)
val -= *diag--**p++;
// res += val*val;
scratch(i, ic) = val;
}
//
// perform the SOR operation
//
for (i = ifirst; i != ilast; i += iinc)
{
// TPZColuna *mydiag = &fDiag[i];
int offset = Size(i);
// REAL val = scratch(i,ic);
REAL *diag;
REAL *diaglast = fElem[i];
REAL *scratchp = &scratch(i - offset + 1, ic);
// val= result(i,ic);
REAL val = scratch(i, ic);
val -= *diaglast * result(i, ic);
res += val * val;
val = overrelax * val / *diaglast;
result(i, ic) += val;
val = result(i, ic);
diag = fElem[i] + offset - 1;
while (diag > diaglast)
 * scratchp++ -= *diag--*val;
}
}
}
res = sqrt(res);
}
if (residual)
{
Residual(result, F, *residual);
}
numiterations = it;
tol = res;
}

 */

/** *********** */
/** * GetVal ** */

template <class TVar>
const TVar TPZSkylNSymMatrix<TVar>::GetVal(const int64_t r, const int64_t c)const
{
  // inicializando row e col para trabalhar com a triangular superior
  int64_t row(r), col(c);
  if (col >= row)
  {
    // Indice do vetor coluna.
    int64_t index = col - row;
    // TPZColuna *pCol = &fDiag[col];

    if (index < Size(col)){
      return(fElem[col][index]);
    }else{
      return (TVar)0;
    }
  }
  else
  {
    // Indice do vetor coluna.
    int64_t index = row - col;
    // TPZColuna *pCol = &fDiag[col];

    if (index < Size(row)){
      return(fElemb[row][index]);
    }else{
      return (TVar)0;
    }
  }
}

///

template <class TVar>
const TVar & TPZSkylNSymMatrix<TVar>::GetValSup(const int64_t r, const int64_t c)const
{

  int64_t row(r), col(c);
  int64_t index = col - row;
#ifdef PZDEBUG

  if (row >= this->Dim() || col >= this->Dim() || row < 0 || col < 0)
  {
    cout << "TPZSkylMatrix::GetVal index out of range row = " << row <<
        " col = " << col << endl;
    return this->gZero;
  }
  if (index < Size(col))
    DebugStop();

#endif
  return(fElem[col][index]);
}

template <class TVar>
const TVar & TPZSkylNSymMatrix<TVar>::GetValB(const int64_t r, const int64_t c)const
{

  int64_t row(r), col(c);
  int64_t index = row - col;
#ifdef PZDEBUG

  if (row >= this->Dim() || col >= this->Dim() || row < 0 || col < 0)
  {
    cout << "TPZSkylMatrix::GetVal index out of range row = " << row <<
        " col = " << col << endl;
    return this->gZero;
  }
  if (index < Size(col))
    DebugStop();

#endif
  return(fElem[row][index]);
}

/** ****** Operacoes com matrizes SKY LINE  ******* */

/** *************** */




/** *************** */
/** * Operator + ** */

template<class TVar>
TPZSkylNSymMatrix<TVar>
TPZSkylNSymMatrix<TVar>::operator+(const TPZSkylNSymMatrix<TVar> & A)const
{
  const auto dim = this->Dim();
  if (dim != A.Dim())
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "<incompatible dimensions>");

  TPZVec<int64_t> skylinesize(dim);
  ComputeMaxSkyline(*this, A, skylinesize);
  TPZSkylNSymMatrix<TVar>res(this->Rows(),this->Cols());
  res.SetSkyline(skylinesize);

  TVar *elemMaxA;
  TVar *elemMinA;
  TVar *elemMaxB;
  TVar *elemMinB;
  int sizeMax;
  int sizeMin;

  for (int col = 0; col < dim; col++)
    {
      // Define o tamanho e os elementos da maior e da menor
      // coluna.
      if (Size(col) > A.Size(col))
        {
          sizeMax = Size(col);
          sizeMin = A.Size(col);
          elemMaxA = fElem[col];
          elemMinA = A.fElem[col];
          elemMaxB = fElemb[col];
          elemMinB = A.fElemb[col];
        }
      else
        {
          sizeMax = A.Size(col);
          sizeMin = Size(col);
          elemMaxA = A.fElem[col];
          elemMinA = fElem[col];
          elemMaxB = A.fElem[col];
          elemMinB = fElem[col];
        }

      // Inicializa coluna da matriz resultado.

      // Efetua a SOMA.
      TVar *destA = res.fElem[col];
      TVar *destB = res.fElemb[col];
      int i;
      for (i = 0; i < sizeMin; i++)
        * destA++ = (*elemMaxA++) + (*elemMinA++);
      for (; i < sizeMax; i++)
        * destA++ = *elemMaxA++;
      for (i = 0; i < sizeMin; i++)
        * destB++ = (*elemMaxB++) + (*elemMinB++);
      for (; i < sizeMax; i++)
        * destB++ = *elemMaxB++;
    }

  return res;
}


/** *************** */
/** * Operator - ** */
template<class TVar>
TPZSkylNSymMatrix<TVar>
TPZSkylNSymMatrix<TVar>::operator-(const TPZSkylNSymMatrix<TVar> & A)const
{
  const auto dim = this->Dim();
  if (dim != A.Dim())
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__, "<incompatible dimensions>");

  TPZVec<int64_t> skylinesize(dim);
  ComputeMaxSkyline(*this, A, skylinesize);
  TPZSkylNSymMatrix<TVar> res(this->Rows(),this->Cols());
  res.SetSkyline(skylinesize);

  TVar *elemMaxA;
  TVar *elemMinA;
  TVar *elemMaxB;
  TVar *elemMinB;
  int sizeMax;
  int sizeMin;

  for (int col = 0; col < dim; col++)
    {
      // Define o tamanho e os elementos da maior e da menor
      // coluna.
      if (Size(col) > A.Size(col))
        {
          sizeMax = Size(col);
          sizeMin = A.Size(col);
          elemMaxA = fElem[col];
          elemMinA = A.fElem[col];
          elemMaxB = fElemb[col];
          elemMinB = A.fElemb[col];
        }
      else
        {
          sizeMax = A.Size(col);
          sizeMin = Size(col);
          elemMaxA = A.fElem[col];
          elemMinA = fElem[col];
          elemMaxB = A.fElem[col];
          elemMinB = fElem[col];
        }

      // Inicializa coluna da matriz resultado.

      
      TVar *destA = res.fElem[col];
      TVar *destB = res.fElemb[col];
      int i;
      for (i = 0; i < sizeMin; i++)
        * destA++ = (*elemMaxA++) - (*elemMinA++);
      for (; i < sizeMax; i++)
        * destA++ = -(*elemMaxA++);
      for (i = 0; i < sizeMin; i++)
        * destB++ = (*elemMaxB++) - (*elemMinB++);
      for (; i < sizeMax; i++)
        * destB++ = -(*elemMaxB++);
    }
  return res;
}

/** **************** */
/** * Operator += ** */

template<class TVar>
TPZSkylNSymMatrix<TVar> &
TPZSkylNSymMatrix<TVar>:: operator += (const TPZSkylNSymMatrix<TVar> & A)
{
  TPZSkylNSymMatrix<TVar> res((*this) + A);
  *this = res;
  return *this;
}

 

/** **************** */
/** * Operator -= ** */


template<class TVar>
TPZSkylNSymMatrix<TVar> &
TPZSkylNSymMatrix<TVar>:: operator -= (const TPZSkylNSymMatrix<TVar> & A)
{
  TPZSkylNSymMatrix<TVar> res((*this) - A);
  *this = res;
  return *this;
}


/** ****** Operacoes com valores NUMERICOS ******* */
//
// As operacoes com valores numericos sao efetuadas apenas nos
// elementos alocados. Em especial, as operacoes A = 0.0 e A *= 0.0
// desalocam todos os elementos da matriz.
//

/** ************************** */
/** * Operator * ( REAL ) ** */

template<class TVar>
TPZSkylNSymMatrix<TVar> TPZSkylNSymMatrix<TVar>:: operator*(const TVar value)const
{
  auto res(*this);

  for (int col = 0; col < this->Dim(); col++)
    {
      // Aloca nova coluna para o resultado.
      int colsize = Size(col);
      // Efetua a SOMA.
      TVar *elemRes = res.fElem[col];
      for (int i = 0; i < colsize; i++)
        * elemRes++ *= value;
      elemRes = res.fElemb[col];
      for (int i = 0; i < colsize; i++)
        * elemRes++ *= value;
    }

  return res;
}


template<class TVar>
TVar* &TPZSkylNSymMatrix<TVar>::Elem()
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Not implemented\n.Aborting...\n";
  DebugStop();
  static TVar* t{nullptr};
  return t;
}

template<class TVar>
const TVar* TPZSkylNSymMatrix<TVar>::Elem()const
{
  PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Not implemented\n.Aborting...\n";
  DebugStop();
  return nullptr;
}

/** *************************** */
/** * Operator *= ( REAL ) ** */

template<class TVar>
TPZSkylNSymMatrix<TVar> & TPZSkylNSymMatrix<TVar>:: operator *= (const TVar value)
{
  if (IsZero(value))
    {
      Zero();
      return(*this);
    }

  int col, colmax = this->Dim();
  for (col = 0; col < colmax; col++)
    {
      // Efetua a MULTIPLICACAO.
      TVar *elem = fElem[col];
      TVar *end = fElem[col + 1];
      while (elem < end)
        * elem++ *= value;
    }

  this->fDecomposed = 0;
  return(*this);
}


/** *********** */
/** * Resize ** */
//
// Muda as dimensoes da matriz, mas matem seus valores antigos. Novas
// posicoes sao criadas com ZEROS.
//

/*
int TPZSkylMatrix::Resize(int newDim, int)
{
if (newDim == Dim())
return(1);

fElem.Resize(newDim + 1);
// Cria nova matrix.

// Copia os elementos para a nova matriz.
int min = Min(newDim, Dim());
int i;
for (i = min + 1; i <= newDim; i++)
fElem[i] = fElem[i - 1];

// Zera as posicoes que sobrarem (se sobrarem)
fStorage.Resize(fElem[newDim] - fElem[0]);
fRow = fCol = newDim;
fDecomposed = 0;
return(1);
}

 */

/** ********** */
/** * Redim ** */
//
// Muda as dimensoes da matriz e ZERA seus elementos.
//

/*
int TPZSkylMatrix::Redim(int newDim, int)
{
if (newDim == Dim())
{
Zero();
return(1);
}

Clear();
fElem.Resize(newDim);
fElem.Fill(0);
fRow = fCol = newDim;
fDecomposed = 0;
return(1);
}


 */

/** ****************** */
/** * LU Decomposition ** */
//
// Faz A= LU onde L eh triangular inferior com 1 na diagonal e U eh triangular superior
//
template <class TVar>
int TPZSkylNSymMatrix<TVar>::Decompose_LU()
{
  if (this->fDecomposed == ELU)
    return 1;
  if (this->fDecomposed)
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,
    "Decompose_LU <Matrix already Decomposed>");

  TVar pivot;
  int64_t dimension = this->Dim();
  /* if(Dim() > 100) {
  cout << "\nTPZSkylMatrix Cholesky decomposition Dim = " << Dim() << endl;
  cout.flush();
  } */
  for (int64_t k = 0; k < dimension; k++)
  {
    /* if(!(k%100) && Dim() > 100) {
    cout <<  k << ' ';
    cout.flush();
    }
    if(!(k%1000)) cout << endl; */
    if (Size(k) == 0)
    {
      PZError << "CUIDADO, A MATRIZ NAO TEM POSICOES NA LINHA " << k << "\n";
      return(0);
    }

    // Faz sum = SOMA( A(k,p) * A(k,p) ), p = 1, ..., k-1.
    //
    TVar sum = 0.0;
    TVar *elem_k = fElem[k] + 1;
    TVar *elem_kb = fElemb[k] + 1;
    TVar *end_k = fElem[k] + Size(k);
    for (; elem_k < end_k; elem_k++, elem_kb++)
    {
      sum += (*elem_k) * (*elem_kb);
    }

    // Faz A(k,k) = sqrt( A(k,k) - sum ).
    //
		
    pivot = fElem[k][0] - sum;
    if (fabs(pivot) < 1.e-25)
    {
      cout <<
          "TPZSkylNSymMatrix<TVar>::Decompose_LU a matrix nao e positiva definida"
          << pivot << endl;
      return(0);
    }
    // A matriz nao e' definida positiva.

    fElem[k][0] = pivot;

    // Loop para i = k+1 ... Dim().
    //
    int64_t i = k + 1;
    for (int64_t j = 2; i < dimension; j++, i++)
    {
      // Se tiverem elementos na linha 'i' cuja coluna e'
      // menor do que 'K'...
      if (Size(i) > j)
      {
        // linha
        // Faz sum = SOMA( A(i,p) * A(k,p) ), p = 1,..,k-1.
        sum = 0.0;
        TVar *elem_i = &fElem[i][j];
        TVar *end_i = fElem[i + 1];
        elem_k = &(fElemb[k][1]);
        end_k = fElemb[k] + Size(k);
        while ((elem_i < end_i) && (elem_k < end_k))
          sum += (*elem_i++) * (*elem_k++);

        // Faz A(i,k) = (A(i,k) - sum) / A(k,k)
        fElem[i][j - 1] = (fElem[i][j - 1] - sum);

        // coluna
        sum = 0.0;
        elem_i = &fElemb[i][j];
        end_i = fElemb[i + 1];
        elem_k = &(fElem[k][1]);
        end_k = fElem[k] + Size(k);
        while ((elem_i < end_i) && (elem_k < end_k))
          sum += (*elem_i++) * (*elem_k++);

        // Faz A(i,k) = (A(i,k) - sum) / A(k,k)
        fElemb[i][j - 1] = (fElemb[i][j - 1] - sum) / pivot;
      }
      else if (Size(i) == j)
      {
        fElemb[i][j - 1] /= pivot;
      }

      // Se nao tiverem estes elementos, sum = 0.0.

      // Se nao existir nem o elemento A(i,k), nao faz nada.
    }
    // Print("decomposing");
  }

  this->fDecomposed = ELU;
  this->fDefPositive = 1;
  return(1);
}

/** ****************** */
/** * Subst Forward ** */
//
// Faz Ax = b, onde A e' triangular inferior.
//

/*
int TPZSkylMatrix::Subst_Forward(TPZFMatrix *B)const
{
if ((B->Rows() != Dim()) || fDecomposed != ECholesky)
TPZMatrix::Error(__PRETTY_FUNCTION__,
"TPZSkylMatrix::Subst_Forward not decomposed with cholesky");

// std::cout << "SubstForward this " << (void *) this << " neq " << Dim() << " normb " << Norm(*B) << std::endl;
int dimension = Dim();
for (int j = 0; j < B->Cols(); j++)
{
int k = 0;
while (k < dimension && (*B)(k, j) == 0)
{
k++;
}
// std::cout << "kstart " << k << std::endl;
for (; k < dimension; k++)
{
// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
//
REAL sum = 0.0;
REAL *elem_ki = fElem[k] + 1;
REAL *end_ki = fElem[k + 1];
REAL* BPtr = &(*B)(k, j); // (k-1,j)
// for ( int i = k-1; elem_ki < end_ki; i-- )
// sum += (*elem_ki++) * B->GetVal( i, j );

while (elem_ki < end_ki)
sum += (*elem_ki++) * (*--BPtr); // (*BPtr--)
// Faz B[k,j] = (B[k,j] - sum) / A[k,k].
//
// B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
BPtr = &(*B)(k, j);
 *BPtr -= sum;
 *BPtr /= fElem[k][0];
}
}

return(1);
}

 */
/** ****************** */
/** * Subst Backward ** */
//
// Faz Ax = b, onde A e' triangular superior.
//
template <class TVar>
int TPZSkylNSymMatrix<TVar>::Subst_Backward(TPZFMatrix<TVar> *B)const
{
  // return TSimMatrix::Subst_Backward(B);

  if ((B->Rows() != this->Dim()) || this->fDecomposed != ELU)
    TPZMatrix<TVar>::Error(__PRETTY_FUNCTION__,
    "TPZSkylMatrix::Subst_Backward not decomposed with LU");

  // std::cout << "SubstBackward this " << (void *) this << " neq " << Dim() << " ncols " << B->Cols() << std::endl;

  int64_t Dimension = this->Dim();
  if (!Dimension)
    return 1; // nothing to do
  int64_t j;
  for (j = 0; j < B->Cols(); j++)
  {
    int64_t k = Dimension - 1;
    while (k > 0 && (*B)(k, j) == TVar(0.))
    {
      k--;
    }
    // std::cout << "kstart " << k << std::endl;

    for (; k > 0; k--)
    {
      // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
      //
      TVar val;
      TVar *elem_ki = fElem[k];
      TVar *end_ki = fElem[k + 1];
      TVar *BPtr = &(*B)(k, j);
      *BPtr /= *elem_ki++;
      val = *BPtr;
      // BPtr;
      // substract the column of the skyline matrix from the vector.
      while (elem_ki < end_ki)
        * --BPtr -= (*elem_ki++) * val;
    }
  }
  for (j = 0; j < B->Cols(); j++)
    (*B)(0, j) /= fElem[0][0];
  return(1);
}

/** ******************** */
/** * Subst L Forward ** */
//
// Faz a "Forward substitution" assumindo que os elementos
// da diagonal sao todos iguais a 1.
//
template <class TVar>
int TPZSkylNSymMatrix<TVar>::Subst_LForward(TPZFMatrix<TVar> *B)const
{
  if ((B->Rows() != this->Dim()) || (this->fDecomposed != ELDLt && this->fDecomposed != ELU))
    return(0);

  int64_t dimension = this->Dim();
  for (int64_t k = 0; k < dimension; k++)
  {
    for (int64_t j = 0; j < B->Cols(); j++)
    {
      // Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
      //
      TVar sum = 0.0;
      TVar *elem_ki = fElemb[k] + 1;
      TVar *end_ki = fElemb[k + 1];
      TVar *BPtr = &(*B)(k, j);
      while (elem_ki < end_ki)
        sum += (*elem_ki++) * (*--BPtr);

      // Faz B[k,j] = (B[k,j] - sum) / A[k,k].
      //
      // B->PutVal( k, j, (B->GetVal(k, j) - sum) / row_k->pElem[0] );
      BPtr = &(*B)(k, j);
      *BPtr -= sum;
    }
  }
  return(1);
}

/** *************** */
/** * Subst Diag ** */
//
// Faz Ax = b, sendo que A e' assumida ser uma matriz diagonal.
//

/*
int TPZSkylMatrix::Subst_Diag(TPZFMatrix *B)const
{
if ((B->Rows() != Dim()) || fDecomposed != ELDLt)
return(0);
int dimension = Dim();
for (int j = 0; j < B->Cols(); j++)
{
REAL *BPtr = &(*B)(0, j);
int k = 0;
while (k < dimension)
 * BPtr++ /= *(fElem[k++]);
}
return(1);
}

 */

/** ****************** */
/** * Subst Backward ** */
//
// Faz Ax = b, onde A e' triangular superior.
//

/*
int TPZSkylMatrix::Subst_LBackward(TPZFMatrix *B)const
{
// return TSimMatrix::Subst_Backward(B);

if ((B->Rows() != Dim()) || !fDecomposed || fDecomposed == ECholesky)
TPZMatrix::Error(__PRETTY_FUNCTION__,
"TPZSkylMatrix::Subst_LBackward not decomposed properly");

int Dimension = Dim();
for (int k = Dimension - 1; k > 0; k--)
{
for (int j = 0; j < B->Cols(); j++)
{
// Faz sum = SOMA( A[k,i] * B[i,j] ), para i = 1, ..., k-1.
//
// REAL val = 0.0;
REAL *elem_ki = fElem[k] + 1;
REAL *end_ki = fElem[k + 1];
REAL *BPtr = &(*B)(k, j);
// substract the column of the skyline matrix from the vector.

REAL val = *BPtr;
while (elem_ki < end_ki)
 * --BPtr -= (*elem_ki++) * val;
}
}
return(1);
}

 */

/** ************************** PRIVATE *************************** */

/** ************ */
/** * PutZero ** */

/*
int TPZSkylMatrix::Zero()
{

fStorage.Fill(0.);
fDecomposed = 0;
fDefPositive = 0;
return(1);
}

 */

/** ********** */
/** * Error ** */
/* int
// ESTAVA COMENTADO NO SKYLINE SIMETRICO
TPZSkylMatrix::Error(const char *msg1,const char* msg2 )
{
ostringstream out;
out << "TPZSkylMatrix::" << msg1 << msg2 << ".\n";
//pzerror.Show();
LOGPZ_ERROR (logger, out.str().c_str());
DebugStop();
return 0;
}
 */

/** ********** */
/** * CLear ** */
template <class TVar>
int TPZSkylNSymMatrix<TVar>::Clear()
{
  fStorage.Resize(0);
  fStorageb.Resize(0);
  fStorage.Shrink();
  fStorageb.Shrink();
  fElem.Resize(0);
  fElemb.Resize(0);
  this->fRow = this->fCol = 0;
  this->fDecomposed = 0;
  return(1);
}

/** ********* */
/** * Copy ** */

template <class TVar>
void TPZSkylNSymMatrix<TVar>::Copy(const TPZSkylNSymMatrix &A)
{
  int64_t dimension = A.Dim();
  this->fRow = this->fCol = dimension;
  fElem.Resize(A.fElem.NElements());
  fElemb.Resize(A.fElemb.NElements());
  fStorage = A.fStorage;
  fStorageb = A.fStorageb;
  int64_t i;
  TVar *firstp = 0;

  if (fStorage.NElements())
    firstp = &fStorage[0];
  for (i = 0; i <= dimension; i++)
    fElem[i] = firstp + (A.fElem[i] - A.fElem[0]);

  if (fStorageb.NElements())
    firstp = &fStorageb[0];
  for (i = 0; i <= dimension; i++)
    fElemb[i] = firstp + (A.fElemb[i] - A.fElemb[0]);

  this->fDecomposed = A.fDecomposed;
  this->fDefPositive = A.fDefPositive;
}

/*
TPZSkylMatrix TPZSkylMatrix:: operator -()const
{
return operator*(-1.0);
}

 */

/* OOPAR!

#ifdef OOPARLIB

int TPZSkylMatrix::Unpack(TReceiveStorage *buf)
{
TSimMatrix::Unpack(buf);
int rows;
buf->UpkInt(&rows);
Redim(rows, rows);
int sz;
for (int i = 0; i < rows; i++)
{
sz = fDiag[i].size;
buf->UpkInt(&sz);
for (int j = 0; j < sz; j++)
{
buf->UpkDouble(fDiag[i].pElem);
}
}
return 1;
}

int TPZSkylMatrix::Pack(TSendStorage *buf)
{
TSimMatrix::Pack(buf);
int rows = Rows();
int sz;
for (int i = 0; i < rows; i++)
{
buf->PkInt(&sz);
fDiag[i].size = sz;
fDiag[i].vetSize = sz;
fDiag[i].pElem = new REAL[sz];
for (int j = 0; j < sz; j++)
{
buf->PkDouble(fDiag[i].pElem);
}
}
return 1;
}

TSaveable *TPZSkylMatrix::Restore(TReceiveStorage *buf)
{
TPZSkylMatrix *m = new TPZSkylMatrix(0);
m->Unpack(buf);
return m;
}

int TPZSkylMatrix::DerivedFrom(int64_t Classid)
{
return 1;
return TSimMatrix::DerivedFrom(Classid);
}

int TPZSkylMatrix::DerivedFrom(char *classname)
{
if (!strcmp(ClassName(), classname))
return 1;
return TSimMatrix::DerivedFrom(classname);
}

#endif
 */

/*

void TPZSkylMatrix::DecomposeColumn(int col, int prevcol)
{
REAL *ptrprev; // Pointer to prev column
REAL *ptrcol; // Pointer to col column
int skprev, skcol; // prev and col Skyline height respectively
int minline;

skprev = SkyHeight(prevcol);
skcol = SkyHeight(col);

ptrprev = Diag(prevcol);
ptrcol = Diag(col);

if ((prevcol - skprev) > (col - skcol))
{
minline = prevcol - skprev;
}
else
{
minline = col - skcol;
}
if (minline > prevcol)
{
cout << "error condition\n";
cout.flush();
return;
}
REAL *run1 = ptrprev + (prevcol - minline);
REAL *run2 = ptrcol + (col - minline);
REAL sum = 0;

//  while(run1-ptrprev > templatedepth) {
//  run1-=templatedepth-1;
//  run2-=templatedepth-1;
//  sum += TemplateSum<templatedepth>(run1--,run2--);
//  }


while (run1 != ptrprev)
sum += (*run1--) * (*run2--);
 *run2 -= sum;
if (run1 != run2)
{
 *run2 /= *run1;
}
else
{
 *run2 = sqrt(*run2);
}

}

 */

/*
void TPZSkylMatrix::DecomposeColumn(int col, int prevcol,
std::list<int64_t> &singular)
{
REAL *ptrprev; // Pointer to prev column
REAL *ptrcol; // Pointer to col column
int skprev, skcol; // prev and col Skyline height respectively
int minline;

skprev = SkyHeight(prevcol);
skcol = SkyHeight(col);

ptrprev = Diag(prevcol);
ptrcol = Diag(col);

if ((prevcol - skprev) > (col - skcol))
{
minline = prevcol - skprev;
}
else
{
minline = col - skcol;
}
if (minline > prevcol)
{
cout << "error condition\n";
cout.flush();
return;
}
REAL *run1 = ptrprev + (prevcol - minline);
REAL *run2 = ptrcol + (col - minline);
REAL sum = 0;

//  while(run1-ptrprev > templatedepth) {
//  run1-=templatedepth-1;
//  run2-=templatedepth-1;
//  sum += TemplateSum<templatedepth>(run1--,run2--);
//  }


while (run1 != ptrprev)
sum += (*run1--) * (*run2--);
 *run2 -= sum;
if (run1 != run2)
{
 *run2 /= *run1;
}
else
{
REAL pivot = *run2;
if (pivot < 1.e-10)
{
#ifdef PZ_LOG
std::stringstream sout;
sout << "equation " << col << " is singular pivot " << pivot;
LOGPZ_WARN(logger, sout.str())
#endif
singular.push_back(col);
pivot = 1.;
}

 *run2 = sqrt(pivot);
}

}

 */

/*

void TPZSkylMatrix::DecomposeColumn2(int col, int prevcol)
{

// cout << "lcol " << lcol << " with " << lprevcol << endl;
// cout.flush();
// Cholesky Decomposition
REAL *ptrprev; // Pointer to prev column
REAL *ptrcol; // Pointer to col column
int skprev, skcol; // prev and col Skyline height respectively
int minline;

skprev = SkyHeight(prevcol);
skcol = SkyHeight(col);

ptrprev = Diag(prevcol);
ptrcol = Diag(col);

if ((prevcol - skprev) > (col - skcol))
{
minline = prevcol - skprev;
}
else
{
minline = col - skcol;
}
if (minline > prevcol)
{
cout << "error condition\n";
cout.flush();
return;
}
REAL *run1 = ptrprev + 1;
REAL *run2 = ptrcol + (col - prevcol) + 1;
REAL *lastptr = ptrprev + prevcol - minline + 1;
REAL sum = 0;
REAL *modify = ptrcol + (col - prevcol);

//  while(lastptr-run1 > templatedepth) {
//  sum += TemplateSum<templatedepth>(run1,run2);
//  run1+=templatedepth;
//  run2+=templatedepth;
//  }


while (run1 != lastptr)
sum += (*run1++) * (*run2++);

// cout << "col " << col << " prevcol " << prevcol << " sum " << sum << " modify " << *modify << " ptrprev " << *ptrprev;
 *modify -= sum;
if (col != prevcol)
{
 *modify /= *ptrprev;
}
else
{
if (*modify < 1.e-25)
{
cout <<
"TPZSkylMatrix::DecomposeCholesky a matrix nao e positiva definida"
<< *modify << endl;
 *modify = 1.e-10;
// return;
}

 *modify = sqrt(*modify);
}
// cout << " modified " << *modify << endl;
// cout.flush();

}
*/


template <class TVar>
void TPZSkylNSymMatrix<TVar>::Read(TPZStream &buf, void *context )
{
	TPZMatrix<TVar>::Read(buf, context);
	buf.Read( fStorage);
	buf.Read( fStorage);
	TPZVec<int> skyl(this->Rows()+1,0), skyl2(this->Rows()+1,0);
	buf.Read( skyl);
	buf.Read( skyl2);
	TVar *ptr = 0, *ptr2 = 0;
	if (this->Rows()) {
		ptr = &fStorage[0];
		ptr2 = &fStorageb[0];
	}
	fElem.Resize(this->Rows()+1);
	fElemb.Resize(this->Rows()+1);
	for (int64_t i=0; i<this->Rows()+1; i++) {
		fElem[i] = skyl[i] + ptr;
		fElemb[i] = skyl2[i] + ptr2;
	}
}

template <class TVar>
void TPZSkylNSymMatrix<TVar>::Write( TPZStream &buf, int withclassid ) const
{
	TPZMatrix<TVar>::Write(buf,withclassid);
	buf.Write( fStorage);
	buf.Write( fStorageb);
	TPZVec<int> skyl(this->Rows()+1,0), skyl2(this->Rows()+1,0);
	TVar *ptr = 0, *ptr2 = 0;
	if (this->Rows()) {
		ptr = &fStorage[0];
		ptr2 = &fStorageb[0];
	}
	for (int64_t i=0; i<this->Rows()+1; i++) {
		skyl[i] = fElem[i] - ptr;
		skyl2[i] = fElemb[i] - ptr2;
	}
	buf.Write( skyl);
	buf.Write( skyl2);
}

/** Fill the matrix with random values (non singular matrix) */
template <class TVar>
void TPZSkylNSymMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric) {
    if (nrow != ncol) {
        DebugStop();
    }
    TPZMatrix<TVar>::Resize(nrow,nrow);
    // initialize the skyline
    TPZManVector<int64_t> skyline(this->Rows());
    for (int64_t i=0; i<this->Rows(); i++) {
        int randcol = rand()%(i+1);
        skyline[i] = randcol;
    }
    this->SetSkyline(skyline);
	int64_t i, j;
	TVar val;
	/** Fill data */
	for(i=0;i<this->Rows();i++) {
		for(j=skyline[i];j<=i ;j++) {
          val = this->GetRandomVal();
          if constexpr(is_complex<TVar>::value){
            if(j==i) val = fabs(val);
          }
			if(!PutVal(i,j,val))
            {
				this->Error("AutoFill (TPZMatrix) failed.");
            }
            if (symmetric == 0) {
              val = this->GetRandomVal();
            }else if constexpr(is_complex<TVar>::value){
              val = std::conj(val);
            }
			if(!PutVal(j,i,val))
            {
				this->Error("AutoFill (TPZMatrix) failed.");
            }
		}
    }
    for (i=0; i<this->Rows(); i++) 
    {
        TVar sum = 0.;
        for (j=0; j<this->Rows(); j++) 
        {
          sum += fabs(this->GetVal(i,j));
        }
        sum = fabs(sum);
        /** Making diagonally dominant and non zero in diagonal */
        if(fabs(sum) > fabs(GetVal(i,i))) {           // Deve satisfazer:  |Aii| > SUM( |Aij| )  sobre j != i
            PutVal(i,i,sum+(TVar)1.);
        }
        // To sure diagonal is not zero.
        if(IsZero(sum) && IsZero(GetVal(i,i)))
        {
            PutVal(i,i,1.);
        }
	}
}


template class TPZSkylNSymMatrix<float>;
template class TPZSkylNSymMatrix<double>;
template class TPZSkylNSymMatrix<long double>;

template class TPZSkylNSymMatrix<std::complex<float> >;
template class TPZSkylNSymMatrix<std::complex<double> >;
template class TPZSkylNSymMatrix<std::complex<long double> >;

template class TPZRestoreClass<TPZSkylNSymMatrix<float>>;
template class TPZRestoreClass<TPZSkylNSymMatrix<double>>;
template class TPZRestoreClass<TPZSkylNSymMatrix<long double>>;

template class TPZRestoreClass<TPZSkylNSymMatrix<std::complex<float>>>;
template class TPZRestoreClass<TPZSkylNSymMatrix<std::complex<double>>>;
template class TPZRestoreClass<TPZSkylNSymMatrix<std::complex<long double>>>;

