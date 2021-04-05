/**
 * @file
 * @brief Contains TPZFBMatrix class which defines a non symmetric banded matrix. Some functionalities depend on LAPACK.

 The functionalities that depend on LAPACK will result in runtime error if LAPACK is not linked to NeoPZ. Search for LAPACK in this header to 
 know which functions are affected by this dependency.
 LAPACK can be linked by setting USING_LAPACK=ON or USING_MKL=ON on CMake
when configuring the library.
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
	virtual int Substitution(TPZFMatrix<TVar> *B) const override;
	/** @brief Simple constructor */
	TPZFBMatrix ();
	/**
     * @brief Constructor
     * @param dim Initial dimension of banded matrix
     * @param band_width Initial band width of banded matrix
	 */
	TPZFBMatrix (const int64_t dim,const int64_t band_width = 0 );
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
        int64_t nel = fElem.size();
        for (int64_t el=0; el<nel; el++) {
            fElem[el] = orig.fElem[el];
        }
        fPivot = orig.fPivot;
        
    }
    

    
    void AutoFill(int64_t nrow, int64_t ncol, int symmetric) override;

    
	int    Put(const int64_t row,const int64_t col,const TVar& value ) override;
	const TVar &Get(const int64_t row,const int64_t col ) const override;
	
	TVar &operator()(const int64_t row, const int64_t col);
	virtual TVar &s(const int64_t row, const int64_t col) override;

	inline int    PutVal(const int64_t row,const int64_t col,const TVar& value ) override;
	inline const TVar &GetVal(const int64_t row,const int64_t col ) const override;
	
	void MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y, TPZFMatrix<TVar> &z,
				 const TVar alpha=1,const TVar beta = 0,const int opt = 0) const override;
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
	
	int64_t Dim() const   override { return this->Rows(); }
	/** @brief Returns band size */
	int64_t GetBandLower() const
    {
        return fBandLower;
    }
    int64_t GetBandUpper() const
    {
        return fBandUpper;
    }
    int64_t GetBand() const
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
	int SetBand(const int64_t newBand );
	
	/// Redimension the matrix preserving its elements
	int Resize(const int64_t newRows,const int64_t newCols ) override;
	
	/// Redimension the matrix and make zero its elements
	int Redim(const int64_t newRows,const int64_t newCols ) override;
	
	// Zeroes the elements of the matrix
	int Zero() override;
	
	void Transpose(TPZMatrix<TVar> *const T) const override;
    //If LAPACK is available, it will use its implementation.
	int       Decompose_LU(std::list<int64_t> &singular) override;
    //If LAPACK is available, it will use its implementation.
	int       Decompose_LU() override;
	
    public:
int ClassId() const override;

    
private:
	
    int64_t Index(int64_t i, int64_t j) const
    {
        return fBandLower+fBandUpper+i-j+(fBandUpper+2*fBandLower+1)*j;
    }
	int Clear() override;
	
	TPZVec<TVar> fElem;
	int64_t  fBandLower, fBandUpper;

    TPZManVector<int,5> fPivot;

};



/**************/
/*** PutVal ***/
template<class TVar>
inline int
TPZFBMatrix<TVar>::PutVal(const int64_t row,const int64_t col,const TVar& value )
{
	if ( (col-row <=fBandUpper) && (row-col <= fBandLower) )
    {
        int64_t index = Index(row,col);
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
TPZFBMatrix<TVar>::GetVal(const int64_t row,const int64_t col ) const {
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
inline TVar &TPZFBMatrix<TVar>::operator()(const int64_t row, const int64_t col){  
    if ( (col-row <=fBandUpper) && (row-col <= fBandLower) )
	{
		return( fElem[Index(row,col)] );
	}
    DebugStop();
    static TVar Zero = (TVar)(0);
	return( Zero );
}
template<class TVar>
inline TVar &TPZFBMatrix<TVar>::s(const int64_t row, const int64_t col) {
	// verificando se o elemento a inserir esta dentro da matriz
	return operator()(row,col);
}

#endif


