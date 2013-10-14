/**
 * @file
 * @brief Contains TPZSBMatrix class which implements symmetric band matrices.
 */

#ifndef TSBNDMATH
#define TSBNDMATH

#include "pzmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

template<class TVar>
class TPZFMatrix;

/**
 * @brief Implements symmetric band matrices. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSBMatrix : public TPZMatrix<TVar>
{
public:
	TPZSBMatrix() : TPZMatrix<TVar>(0,0) { fDiag = NULL; fBand = 0; }
	TPZSBMatrix(const long dim,const long band );
	TPZSBMatrix(const TPZSBMatrix<TVar> &A ) : TPZMatrix<TVar>(A)  { Copy(A); }
	
	CLONEDEF(TPZSBMatrix)
	
	~TPZSBMatrix() { Clear(); }
	
	int    PutVal(const long row,const long col,const TVar& element );
	const TVar &GetVal(const long row,const long col ) const;
	
	/** @brief Computes z = beta * y + alpha * opt(this)*x */
	/** z and x cannot overlap in memory */
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1.,const TVar beta = 0.,const int opt = 0,const int stride = 1 ) const;
	
	void Print(const char *name = NULL, std::ostream &out = std::cout ,const MatrixOutputFormat form = EFormatted) const;
	//friend std::ostream & operator<< <>(std::ostream& out,const TPZSBMatrix<TVar>  &A); Leonardo removendo o '<>' antes do (std...
	template<class TT>friend std::ostream & operator<< (std::ostream& out,const TPZSBMatrix<TT>  &A); 
	
	/// Operadores com matrizes SKY LINE.
	// @{
	TPZSBMatrix &operator= (const TPZSBMatrix<TVar> &A );
	TPZSBMatrix operator+  (const TPZSBMatrix<TVar> &A ) const;
	TPZSBMatrix operator-  (const TPZSBMatrix<TVar> &A ) const;
	TPZSBMatrix &operator+=(const TPZSBMatrix<TVar> &A );
	TPZSBMatrix &operator-=(const TPZSBMatrix<TVar> &A );
	// @}
	TPZSBMatrix<TVar> operator*  (const TVar v ) const;
	TPZSBMatrix<TVar> &operator*=(const TVar v );
	
	TPZSBMatrix<TVar> operator-() const { return operator*(-1.0); }
	
	/// Redimension the matrix keeping original elements.
	int Resize(const long newDim ,const long);
	
	/// Redimension the matrix and zeroes its elements
	int Redim(const long newDim) {return Redim(newDim,newDim);}
	int Redim(const long newRows ,const long newCols);
	
	/// Zeroes the elements of the matrix
	int Zero();
	
	long GetBand() const { return fBand; }
	int   SetBand(const long newBand );
	
	
	/// To solve linear systems
	// @{
	int Decompose_Cholesky();  // Faz A = GGt.
	int Decompose_LDLt    ();  // Faz A = LDLt.
	int Decompose_Cholesky(std::list<long> &singular);  // Faz A = GGt.
	int Decompose_LDLt    (std::list<long> &singular);  // Faz A = LDLt.
	
	int Subst_Forward  ( TPZFMatrix<TVar> *b ) const;
	int Subst_Backward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LForward ( TPZFMatrix<TVar> *b ) const;
	int Subst_LBackward( TPZFMatrix<TVar> *b ) const;
	int Subst_Diag     ( TPZFMatrix<TVar> *b ) const;
	// @}
	
#ifdef OOPARLIB
	
	virtual long GetClassID() const        { return TSBMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSBMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const; // a class with name classname
	
#endif
	
private:
	
	long  Size() const    { return( this->Dim() * (fBand + 1) ); }
	int  PutZero();
	//static int  Error(const char *msg1,const char* msg2="" ) ;
	int  Clear();
	void Copy (const TPZSBMatrix<TVar> & );
	
	
	TVar *fDiag;
	long  fBand;
};

#endif
