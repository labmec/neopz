/**
 * @file
 * @brief Contains TPZSkyline class which implements a skyline storage format.
 */
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
template<class TVar>
class TPZFMatrix;

/**@note  Esta classe gerencia matrizes do tipo SkyLine. Todas
 matrizes SkyLine sao simetricas (e portanto quadradas).
 */

/**
 * @brief Implements a skyline storage format. A Skyline matrix is symmetric so square. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSkylMatrix : public TPZMatrix<TVar>
{
public:
	TPZSkylMatrix() : TPZMatrix<TVar>(0,0),fElem(0),fStorage(0) { }
	TPZSkylMatrix(const int dim);
	/**
     @brief Construct a skyline matrix of dimension dim
     skyline indicates the minimum row number which will be accessed by each equation
	 */
	TPZSkylMatrix(const int dim ,const TPZVec<int> &skyline);
	TPZSkylMatrix(const TPZSkylMatrix<TVar> &A ) : TPZMatrix<TVar>(A), fElem(0), fStorage(0)  { Copy(A); }
	
	CLONEDEF(TPZSkylMatrix)
	/**
     @brief modify the skyline of the matrix, throwing away its values
     skyline indicates the minimum row number which will be accessed by each equation
	 */
	void SetSkyline(const TPZVec<int> &skyline);
	
	/**
     @brief return the height of the skyline for a given column
	 */
	int SkyHeight(int col) { return fElem[col+1]-fElem[col] - 1; }
	
	/** @brief Add a skyline matrix B with same structure of this
	 *  It makes this += k * B
	 */
	void AddSameStruct(TPZSkylMatrix<TVar> &B, double k = 1.);
	
	/** @brief declare the object as simetric matrix*/
	virtual int IsSimetric() const {return 1;}
	
    /** @brief destructor of the skyline matrix */
	virtual ~TPZSkylMatrix() { Clear(); }
    
    /**
	 * @brief Updates the values of the matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix<TVar> > mat);

	
	int    PutVal(const int row,const int col,const TVar &element );
	const TVar &GetVal(const int row,const int col ) const;
	
	
	TVar &operator()(const int row, const int col);
	virtual TVar &s(const int row, const int col);
	
	
	TVar &operator()(const int row);
	
	virtual void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
						 const TVar alpha,const TVar beta ,const int opt = 0,const int stride = 1 ) const ;
	// Operadores com matrizes SKY LINE.
	TPZSkylMatrix &operator= (const TPZSkylMatrix<TVar> &A );
	//TPZSkylMatrix &operator= (TTempMat<TPZSkylMatrix> A);
	
	TPZSkylMatrix operator+  (const TPZSkylMatrix<TVar> &A ) const;
	TPZSkylMatrix operator-  (const TPZSkylMatrix<TVar> &A ) const;
	
	TPZSkylMatrix &operator+=(const TPZSkylMatrix<TVar> &A );
	TPZSkylMatrix &operator-=(const TPZSkylMatrix<TVar> &A );
	
	// Operadores com valores NUMERICOS.
	TPZSkylMatrix operator*  (const TVar v ) const;
	TPZSkylMatrix &operator*=( TVar value );
	
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
	
	
	/*** @brief To Solve Linear Equations ***/
	// @{
	virtual void SolveSOR(int &numiterations,const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result,
						  TPZFMatrix<TVar> *residual,TPZFMatrix<TVar> &scratch,const TVar overrelax, TVar &tol,
						  const int FromCurrent = 0,const int direction = 1) ;
	
	
	int Decompose_Cholesky();  // Faz A = GGt.
	int Decompose_LDLt    ();  // Faz A = LDLt.
	int Decompose_Cholesky(std::list<int> &singular);  // Faz A = GGt.
	int Decompose_LDLt    (std::list<int> &singular);  // Faz A = LDLt.
	
	int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
	int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LBackward( TPZFMatrix<TVar> *b ) const;
	int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
	// @}
	
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
     @brief This method returns a pointer to the diagonal element of the matrix of the col column
	 */
	TVar *Diag(int col) { return fElem[col];}
	
	void DecomposeColumn(int col, int prevcol);
	void DecomposeColumn(int col, int prevcol, std::list<int> &singular);
	
	void DecomposeColumn2(int col, int prevcol);
private:
	
	// Aloca uma nova coluna. 'fDiag[col].pElem' deve ser NULL.
	
	//static int  Error(const char *msg1,const char* msg2="" );
	int  Clear();
	void Copy (const TPZSkylMatrix<TVar> & );
	int Size(const int column) const {return fElem[column+1]-fElem[column];}
	static int NumElements(const TPZVec<int> &skyline);
	static void InitializeElem(const TPZVec<int> &skyline, TPZManVector<REAL> &storage, TPZVec<REAL *> &elem);
	/**
     @brief Computes the highest skyline of both objects
	 */
	static void ComputeMaxSkyline(const TPZSkylMatrix<TVar> &first, const TPZSkylMatrix<TVar> &second, TPZVec<int> &res);
	
protected:
	/** @brief Storage to keep the first elements to each equation
	 *
	 * fElem is of size number of equation+1
	 * fElem[i] is the first element of the skyline of equation i
	 * fElem[Rows()] is one element beyond the last equation
	 */
	TPZVec<REAL *> fElem;
private:
	/**
     @brief fStorage is a unique vector which contains all the data of the skyline matrix
	 */
	TPZManVector<REAL> fStorage;
};

/** @brief Implements iterative sum over N steps */
template<int N,class TVar>
inline TVar TemplateSum(const TVar *p1, const TVar *p2){
	return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);
	
}
/** @brief Implements product of the values into p1 and p2 */
template<>
inline REAL TemplateSum<1>(const REAL *p1, const REAL *p2){
	return *p1 * *p2;
}

#endif
