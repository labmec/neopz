
//
// Author: MISAEL LUIS SANTANA MANDUJANO.
//
// File:   tskylmat.cc
//
// Class:  TPZSkylMatrix
//
// Obs.:   Esta classe gerencia matrizes do tipo SkyLine. Todas
//         matrizes SkyLine sao simetricas (e portanto quadradas).
//
// Versao: 04 / 1996.
//


#ifndef TSKYLMATH
#define TSKYLMATH



#include "pzmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

class TPZFMatrix;

/**@note  Esta classe gerencia matrizes do tipo SkyLine. Todas
          matrizes SkyLine sao simetricas (e portanto quadradas).
*/

/************************************/
class TPZSkylMatrix : public TPZMatrix
{
 public:
  TPZSkylMatrix() : TPZMatrix(0,0),fElem(0),fStorage(0) { }
  TPZSkylMatrix(const int dim);
  /**
     Construct a skyline matrix of dimension dim
     skyline indicates the minimum row number which will be accessed by each equation
  */
  TPZSkylMatrix(const int dim ,const TPZVec<int> &skyline);
  TPZSkylMatrix(const TPZSkylMatrix &A ) : TPZMatrix(A), fElem(0), fStorage(0)  { Copy(A); }
  
  CLONEDEF(TPZSkylMatrix)
  /**
     modify the skyline of the matrix, throwing away its values
     skyline indicates the minimum row number which will be accessed by each equation
  */
  void SetSkyline(const TPZVec<int> &skyline);

  /**
     return the height of the skyline for a given column
  */
  int SkyHeight(int col) { return fElem[col+1]-fElem[col] - 1; }

  /**declare the object as simetric matrix*/
  virtual int IsSimetric() const {return 1;}

  virtual ~TPZSkylMatrix() { Clear(); }

  int    PutVal(const int row,const int col,const REAL &element );
  const REAL &GetVal(const int row,const int col ) const;


  REAL &operator()(const int row, const int col);
  virtual REAL &s(const int row, const int col);


  REAL &operator()(const int row);

  virtual void MultAdd(const TPZFMatrix &x,const TPZFMatrix &y, TPZFMatrix &z,
		       const REAL alpha,const REAL beta ,const int opt = 0,const int stride = 1 ) const ;
  // Operadores com matrizes SKY LINE.
  TPZSkylMatrix &operator= (const TPZSkylMatrix &A );
  //TPZSkylMatrix &operator= (TTempMat<TPZSkylMatrix> A);

  TPZSkylMatrix operator+  (const TPZSkylMatrix &A ) const;
  TPZSkylMatrix operator-  (const TPZSkylMatrix &A ) const;

  TPZSkylMatrix &operator+=(const TPZSkylMatrix &A );
  TPZSkylMatrix &operator-=(const TPZSkylMatrix &A );

  // Operadores com valores NUMERICOS.
  TPZSkylMatrix operator*  (const REAL v ) const;
  TPZSkylMatrix &operator*=( REAL v );

  TPZSkylMatrix operator-() const;// { return operator*(-1.0); }

  // Redimensiona a matriz, mas mantem seus elementos.
  // o segundo parametro � o tamanho das colunas
  int Resize(const int newDim ,const int );

  // Redimensiona a matriz e ZERA seus elementos.
  // o segundo parametro � o tamanho das colunas
  int Redim(const int newDim ,const int );
  int Redim(const int newDim) {return Redim(newDim,newDim);}

  // Zera os Elementos da matriz
  int Zero();


  /*** Resolucao de sistemas ***/

  virtual void SolveSOR(int &numiterations,const TPZFMatrix &F, TPZFMatrix &result,
			TPZFMatrix *residual,TPZFMatrix &scratch,const REAL overrelax, REAL &tol,
			const int FromCurrent = 0,const int direction = 1) ;


  int Decompose_Cholesky();  // Faz A = GGt.
  int Decompose_LDLt    ();  // Faz A = LDLt.
  int Decompose_Cholesky(std::list<int> &singular);  // Faz A = GGt.
  int Decompose_LDLt    (std::list<int> &singular);  // Faz A = LDLt.

  int Subst_Forward  ( TPZFMatrix *b ) const;
  int Subst_Backward ( TPZFMatrix *b ) const;
  int Subst_LForward ( TPZFMatrix *b ) const;
  int Subst_LBackward( TPZFMatrix *b ) const;
  int Subst_Diag     ( TPZFMatrix *b ) const;

  //void TestSpeed(int col, int prevcol);

#ifdef OOPARLIB

  virtual long GetClassID() const    { return TSKYMATRIX_ID; }
  virtual int Unpack( TReceiveStorage *buf );
  static TSaveable *Restore(TReceiveStorage *buf);
  virtual int Pack( TSendStorage *buf ) const;
  virtual char *ClassName() const   { return( "TPZSkylMatrix"); }
  virtual int DerivedFrom(const long Classid) const;
  virtual int DerivedFrom(const char *classname) const; // a class with name classname

#endif

 protected:

  /**
     This method returns a pointer to the diagonal element of the matrix of the col column
  */
  REAL *Diag(int col) { return fElem[col];}

  void DecomposeColumn(int col, int prevcol);
  void DecomposeColumn2(int col, int prevcol);
 private:

  // Aloca uma nova coluna. 'fDiag[col].pElem' deve ser NULL.

//static int  Error(const char *msg1,const char* msg2="" );
  int  Clear();
  void Copy (const TPZSkylMatrix & );
  int Size(const int column) const {return fElem[column+1]-fElem[column];}
  static int NumElements(const TPZVec<int> &skyline);
  static void InitializeElem(const TPZVec<int> &skyline, TPZManVector<REAL> &storage, TPZVec<REAL *> &elem);
  /**
     Computes the highest skyline of both objects
  */
  static void ComputeMaxSkyline(const TPZSkylMatrix &first, const TPZSkylMatrix &second, TPZVec<int> &res);

  /**
     fElem is of size number of equation+1
     fElem[i] is the first element of the skyline of equation i
     fElem[Rows()] is one element beyond the last equation
  */
  TPZVec<REAL *> fElem;
  /**
     fStorage is a unique vector which contains all the data of the skyline matrix
  */
  TPZManVector<REAL> fStorage;
};

template<int N>
inline REAL TemplateSum(const REAL *p1, const REAL *p2){
  return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);

}


template<>
inline REAL TemplateSum<1>(const REAL *p1, const REAL *p2){
  return *p1 * *p2;
}




#endif
