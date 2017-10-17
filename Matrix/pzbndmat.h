/**
 * @file
 * @brief Contains TPZFBMatrix class which defines a non symmetric banded matrix.
 */

#ifndef _TBNDMATHH_
#define _TBNDMATHH_

#include "pzmatrix.h"

#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif

/**
 * @brief Defines a non symmetric banded matrix. \ref matrix "Matrix"
 * @ingroup matrix
 * @author MISAEL LUIS SANTANA MANDUJANO.
 * @since 12/1994
 */
template<class TVar>
class TPZFBMatrix : public TPZMatrix<TVar>
{
	
public:
	virtual int Substitution(TPZFMatrix<TVar> *B) const;
	/** @brief Simple constructor */
	TPZFBMatrix ();
	/**
     * @brief Constructor
     * @param dim Initial dimension of banded matrix
     * @param band_width Initial band width of banded matrix
	 */
	TPZFBMatrix (const long dim,const long band_width = 0 );
	/** @brief Copy constructor */
	TPZFBMatrix (const TPZFBMatrix<TVar> & );
	
	CLONEDEF(TPZFBMatrix)
	/** @brief Simple destructor */
	~TPZFBMatrix();

    friend class TPZFBMatrix<float>;
    friend class TPZFBMatrix<double>;
    
    /// copy the values from a matrix with a different precision
    template<class TVar2>
    void CopyFrom(TPZFBMatrix<TVar2> &orig)
    {
        TPZMatrix<TVar>::CopyFrom(orig);
        fBandLower = orig.fBandLower;
        fBandUpper = orig.fBandUpper;
        fElem.resize(orig.fElem.size());
        long nel = fElem.size();
        for (long el=0; el<nel; el++) {
            fElem[el] = orig.fElem[el];
        }
#ifdef USING_LAPACK
        fPivot = orig.fPivot;
        long nwork = orig.fWork.size();
        for (long el=0; el<nwork; el++) {
            fWork[el] = orig.fWork[el];
        }
#endif
        
    }
    

    
    void AutoFill(long nrow, long ncol, int symmetric);

    
	int    Put(const long row,const long col,const TVar& value );
	const TVar &Get(const long row,const long col ) const;
	
	TVar &operator()(const long row, const long col);
	virtual TVar &s(const long row, const long col);

	inline int    PutVal(const long row,const long col,const TVar& value );
	inline const TVar &GetVal(const long row,const long col ) const;
	
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1,const TVar beta = 0,const int opt = 0) const;
	// Computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot overlap in memory
	
	// Peforms the product (*this)T x D x (*this).
	//  TPZFBMatrix  InnerProd(TPZFBMatrix &D );
	
	TPZFBMatrix &operator= (const TPZFBMatrix<TVar> & A );
	TPZFBMatrix operator+  (const TPZFBMatrix<TVar> & A ) const;
	TPZFBMatrix operator-  (const TPZFBMatrix<TVar> & A ) const;
	TPZFBMatrix &operator+=(const TPZFBMatrix<TVar> & A );
	TPZFBMatrix &operator-=(const TPZFBMatrix<TVar> & A );
	
	TPZFBMatrix operator*  (const TVar val ) const;
	TPZFBMatrix &operator*=(const TVar val );
	
	TPZFBMatrix operator-() const;
	
	long Dim() const     { return this->Rows(); }
	/** @brief Returns band size */
	long GetBandLower() const
    {
        return fBandLower;
    }
    long GetBandUpper() const
    {
        return fBandUpper;
    }
    long GetBand() const
    {
        if (fBandLower != fBandUpper) {
            DebugStop();
        }
        return fBandUpper;
    }
	/**
     * @brief Sets band size
     * @param newBand New band size
	 */
	int SetBand(const long newBand );
	
	/// Redimension the matrix preserving its elements
	int Resize(const long newRows,const long newCols );
	
	/// Redimension the matrix and make zero its elements
	int Redim(const long newRows,const long newCols );
	
	// Zeroes the elements of the matrix
	int Zero();
	
	void Transpose(TPZMatrix<TVar> *const T) const;
    
#ifdef USING_LAPACK
	int       Decompose_LU(std::list<long> &singular);
	int       Decompose_LU();
#endif
	
    public:
virtual int ClassId() const;

    
private:
	
    long Index(long i, long j) const
    {
        return fBandLower+fBandUpper+i-j+(fBandUpper+2*fBandLower+1)*j;
    }
	int Clear();
	
	TPZVec<TVar> fElem;
	long  fBandLower, fBandUpper;
#ifdef USING_LAPACK
    TPZManVector<int,5> fPivot;
    
    TPZVec<TVar> fWork;
#endif

};



/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZFBMatrix<TVar>::PutVal(const long row,const long col,const TVar& value )
{
	if ( (col-row <=fBandUpper) && (row-col <= fBandLower) )
    {
        long index = Index(row,col);
		fElem[index] = value;
    }
    else if(!IsZero(value))
    {
        DebugStop();
    }
	return( 1 );
}



/**************/
/*** GetVal ***/
template<class TVar>
inline const TVar &
TPZFBMatrix<TVar>::GetVal(const long row,const long col ) const {
#ifdef PZDEBUG
    if (row <0 || row > this->fRow || col < 0 || col >= this->fCol) {
        DebugStop();
    }
#endif
    if ( (col-row <=fBandUpper) && (row-col <= fBandLower) )
    {
        return fElem[Index(row,col)];
    }
    static TVar Zero;
    Zero = TVar(0);
	return( Zero );
}

template<class TVar>
inline TVar &TPZFBMatrix<TVar>::operator()(const long row, const long col){  
    if ( (col-row <=fBandUpper) && (row-col <= fBandLower) )
	{
		return( fElem[Index(row,col)] );
	}
    DebugStop();
    static TVar Zero = (TVar)(0);
	return( Zero );
}
template<class TVar>
inline TVar &TPZFBMatrix<TVar>::s(const long row, const long col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

#endif


