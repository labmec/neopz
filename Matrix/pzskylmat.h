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

class TPZFMatrix;

/**@note  Esta classe gerencia matrizes do tipo SkyLine. Todas
 matrizes SkyLine sao simetricas (e portanto quadradas).
 */

/**
 * @brief Implements a skyline storage format. A Skyline matrix is symmetric so square. \ref matrix "Matrix"
 * @ingroup matrix
 */
class TPZSkylMatrix : public TPZMatrix
{
public:
	TPZSkylMatrix() : TPZMatrix(0,0),fElem(0),fStorage(0) { }
	TPZSkylMatrix(const int dim);
	/**
     @brief Construct a skyline matrix of dimension dim
     skyline indicates the minimum row number which will be accessed by each equation
	 */
	TPZSkylMatrix(const int dim ,const TPZVec<int> &skyline);
	TPZSkylMatrix(const TPZSkylMatrix &A ) : TPZMatrix(A), fElem(0), fStorage(0)  { Copy(A); }
	
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
	void AddSameStruct(TPZSkylMatrix &B, double k = 1.);
	
	/** @brief declare the object as simetric matrix*/
	virtual int IsSimetric() const {return 1;}
	
    /** @brief destructor of the skyline matrix */
	virtual ~TPZSkylMatrix() { Clear(); }
    
    /**
	 * @brief Updates the values of the matrix based on the values of the matrix
	 */
	virtual void UpdateFrom(TPZAutoPointer<TPZMatrix> mat);

	
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
	
	
	/*** @brief To Solve Linear Equations ***/
	// @{
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
	REAL *Diag(int col) { return fElem[col];}
	
	void DecomposeColumn(int col, int prevcol);
	void DecomposeColumn(int col, int prevcol, std::list<int> &singular);
	
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
     @brief Computes the highest skyline of both objects
	 */
	static void ComputeMaxSkyline(const TPZSkylMatrix &first, const TPZSkylMatrix &second, TPZVec<int> &res);
	
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
template<int N>
inline REAL TemplateSum(const REAL *p1, const REAL *p2){
	return *p1* *p2 + TemplateSum<N-1>(p1+1,p2+1);
	
}
/** @brief Implements product of the values into p1 and p2 */
template<>
inline REAL TemplateSum<1>(const REAL *p1, const REAL *p2){
	return *p1 * *p2;
}

#endif
