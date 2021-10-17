/**
 * @file
 * @brief Contains TPZSFMatrix class which implements a symmetric full matrix.
 */

#ifndef _TSFULMATHH_
#define _TSFULMATHH_

#include "pzmatrix.h"
#include "pzfmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

/**
 * @brief Implements a symmetric full matrix. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSFMatrix : public TPZMatrix<TVar> {
	
public:
	TPZSFMatrix () : TPZRegisterClassId(&TPZSFMatrix::ClassId),
    TPZMatrix<TVar>( 0,0 )  { fElem = NULL; }
	TPZSFMatrix (const int64_t dim );
	TPZSFMatrix (const TPZSFMatrix<TVar> & );
  TPZSFMatrix (TPZSFMatrix<TVar> && );
	// Usa o maior bloco quadrado possivel, comecado em (0,0).
	// E inicializa com a parte triangular inferior do bloco.
	TPZSFMatrix (const TPZMatrix<TVar> & );
	TPZSFMatrix &operator= (const TPZSFMatrix<TVar> &A );
  TPZSFMatrix &operator= (TPZSFMatrix<TVar> &&A );
  inline TPZSFMatrix<TVar>*NewMatrix() const override {return new TPZSFMatrix<TVar>{};}
	CLONEDEF(TPZSFMatrix)
	
	~TPZSFMatrix();
	
    /** @brief Checks if the current matrix is symmetric */
    virtual int IsSymmetric() const  override {
        return 1;
    }

    friend class TPZSFMatrix<float>;
    friend class TPZSFMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFromDiffPrecision(TPZSFMatrix<TVar2> &orig)
    {
        Resize(orig.Rows(), orig.Cols());
        TPZMatrix<TVar>::CopyFromDiffPrecision(orig);
        int64_t nel = (this->Rows()*(this->Rows()+1))/2;
        for (int64_t el=0; el<nel; el++) {
            fElem[el] = orig.fElem[el];
        }
    }
    
  /** @brief Creates a copy from another TPZSFMatrix*/
  void CopyFrom(const TPZMatrix<TVar> *  mat) override
  {                                                           
    auto *from = dynamic_cast<const TPZSFMatrix<TVar> *>(mat);                
    if (from) {                                               
      *this = *from;                                          
    }                                                         
    else                                                      
      {                                                       
        PZError<<__PRETTY_FUNCTION__;                         
        PZError<<"\nERROR: Called with incompatible type\n."; 
        PZError<<"Aborting...\n";                             
        DebugStop();                                          
      }                                                       
  }
  
	int PutVal(const int64_t row,const int64_t col,const TVar &value ) override;
	const TVar GetVal(const int64_t row,const int64_t col ) const override;
	
	/**
	 * @name Operators with Full simmetric matrices.
	 * @{
	 */
	TPZSFMatrix operator+  (const TPZSFMatrix<TVar> &A ) const;
	TPZSFMatrix operator-  (const TPZSFMatrix<TVar> &A ) const;
	TPZSFMatrix &operator+=(const TPZSFMatrix<TVar> &A );
	TPZSFMatrix &operator-=(const TPZSFMatrix<TVar> &A );
	/** @} */
	
	/**
	 * @name Operators with generic matrices.
	 * @{
	 */
	TPZSFMatrix &operator= (const TPZMatrix<TVar> &A );
	TPZSFMatrix operator+  (const TPZMatrix<TVar> &A ) const;
	TPZSFMatrix operator-  (const TPZMatrix<TVar> &A ) const;
	TPZSFMatrix &operator+=(const TPZMatrix<TVar> &A );
	TPZSFMatrix &operator-=(const TPZMatrix<TVar> &A );
	/** @} */

	/**
	 * @name Operators with numeric values.
	 * @{
	 */
	TPZSFMatrix &operator= (const TVar val );
	TPZSFMatrix operator+  (const TVar val ) const;
	TPZSFMatrix operator-  (const TVar val ) const { return operator+( -val ); }
	TPZSFMatrix operator*  (const TVar val ) const;
	TPZSFMatrix &operator+=(const TVar val );
	TPZSFMatrix &operator-=(const TVar val )  { return operator+=( -val ); }
	TPZSFMatrix &operator*=(const TVar val ) override;
	/** @} */
	
	TPZSFMatrix operator-() const  { return operator*( -1.0 ); }
	
	/** @brief Resize the array but keeps its entirety. */
	int Resize(const int64_t newDim, const int64_t ) override;
	
	/** @brief Resize the array and resets ist entirety. */
	int Redim(const int64_t newRows ,const int64_t) override;
	
	int Redim(const int64_t newDim) {return Redim(newDim,newDim);}
	
	/** @brief Resets all elements. */
	int Zero() override;
	
	
	/**
	 * @name To solve linear systems
	 * @{
	 */
	virtual int Decompose_Cholesky() override;
	virtual int Decompose_LDLt() override;
	virtual int Decompose_Cholesky(std::list<int64_t> &singular) override;
	virtual int Decompose_LDLt(std::list<int64_t> &singular) override;
	
	virtual int Subst_Forward  ( TPZFMatrix<TVar>  *B ) const override;
	virtual int Subst_Backward ( TPZFMatrix<TVar>  *B ) const override;
	virtual int Subst_LForward ( TPZFMatrix<TVar>  *B ) const override;
	virtual int Subst_LBackward( TPZFMatrix<TVar>  *B ) const override;
	virtual int Subst_Diag     ( TPZFMatrix<TVar>  *B ) const override;
	/** @} */
	
#ifdef OOPARLIB
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSFMatrix"); }
	virtual int DerivedFrom(const int64_t Classid) const;
	virtual int DerivedFrom(const char *classname) const;
	
#endif

  int ClassId() const override;
protected:
  /** @brief Checks compatibility of matrices before Add/Subtract operations*/
  inline void CheckTypeCompatibility(const TPZMatrix<TVar>*A,
                                     const TPZMatrix<TVar>*B)const override
  {
    auto aPtr = dynamic_cast<const TPZSFMatrix<TVar>*>(A);
    auto bPtr = dynamic_cast<const TPZSFMatrix<TVar>*>(B);
    if(!aPtr || !bPtr){
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nERROR: incompatible matrices.Aborting...\n";
      DebugStop();
    }
  }
  inline TVar *&Elem() override
  {
    return fElem;
  }
  inline const TVar *Elem() const override
  {
    return fElem;
  }
	
	int64_t Size() const override { return (this->Dim() * (this->Dim()+1)) >> 1; }
private:
	int Clear() override;
	
	TVar   *fElem;
};

/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZSFMatrix<TVar> ::PutVal(const int64_t row,const int64_t col,const TVar &value )
{
	int64_t locrow = row, loccol = col;
	if ( locrow < loccol )
		this->Swap( &locrow, &loccol );
	
	this->fDecomposed = 0;
	this->fElem[ ((locrow*(locrow+1)) >> 1) + loccol ] = value;
	return( 1 );
}

/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar 
TPZSFMatrix<TVar> ::GetVal(const int64_t row,const int64_t col ) const
{
	if ( row < col ){
        if constexpr (is_complex<TVar>::value){
            return( std::conj(fElem[ ((col*(col+1)) >> 1) + row ]) );
        }else{
            return( fElem[ ((col*(col+1)) >> 1) + row ] );
        }
    }
	
	return( fElem[ ((row*(row+1)) >> 1) + col ] );
}

/**************************/
/*** @name Operadores Globais ***/
/** @{ */

/** @brief Increments value for all entries of the A matrix */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator+( const TVar value, const TPZSFMatrix<TVar>  &A )
{
	return( A + value );
}

/** @brief Decrements value for all entries of the A matrix */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator-( const TVar value,const  TPZSFMatrix<TVar>  &A )
{
	return( A - value );
}

/** @brief Implements the scalar product value*A */
template<class TVar>
inline TPZSFMatrix<TVar> 
operator*( const TVar value, const TPZSFMatrix<TVar>  &A )
{
	return( A * value );
}

/** @} */

#endif
