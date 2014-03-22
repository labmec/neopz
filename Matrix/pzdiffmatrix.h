/**
 * @file
 * @brief Contains TPZDiffMatrix class which to hold the flux derivatives A B C and diffusive matrix coefficients.
 */

#ifndef PZDIFFMATRIX_H
#define PZDIFFMATRIX_H

#include <ostream>
#include "pzmatrix.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

/** @ingroup matrix */
 enum EStatus {EOk = 0, EIncompDim, EZeroPivot};

/**
 * @ingroup matrix
 * @brief Matrix class to hold the flux derivatives A B C and diffusive matrix coefficients. \ref matrix "Matrix"
 * @author Erick Slis
 * @author Cedric Ayala
 * @since June 1, 2003.
 */
template <class T>
class TPZDiffMatrix
{
public:
    TPZDiffMatrix();
    TPZDiffMatrix(const long rows, const long cols);
    ~TPZDiffMatrix();
    
    TPZDiffMatrix(const TPZDiffMatrix &copy) : fRows(copy.fRows), fCols(copy.fCols),fStore(0), fDecomposed(copy.fDecomposed)
    {
        if(fRows *fCols)
        {
            fStore = new T[fRows*fCols];
            for (long i=0; i<fRows*fCols; i++) {
                fStore[i] = copy.fStore[i];
            }
        }
    }
    
//    TPZDiffMatrix &operator=(const TPZDiffMatrix &copy)
//    {
//        fRows = copy.fRows;
//        fCols = copy.fCols;
//        fDecomposed = copy.fDecomposed;
//        if (fStore) {
//            delete []fStore;
//            fStore = 0;
//        }
//        if (fRows*fCols) {
//            fStore = new T[fRows*fCols];
//            for (int i=0; i<fRows*fCols; i++) {
//                fStore[i] = copy.fStore[i];
//            }            
//        }
//    }
	
    /** @brief Resizes and zeroes the matrix. */
    void Redim(const long rows, const long cols);
	
    /** @brief Multiplies the matrix by a correspondent TPZVec vector. Dimensions are checked. */
    void Multiply(TPZVec<T> & In, TPZVec<T> & Out, const T & scale = T(1.));
	
	
    /** @brief Matrix multiplication. Dimensions are checked. */
    void Multiply(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale = T(1.));
	
    /**
     * @brief Matrix multiplication
	 */
	/** 
	 * Dimensions are checked. \n
     * Results are additively contributed to Out
     */
    void MultiplyAdd(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale = T(1.));
	
	
    /** @brief Copies the matrix, reallocating all coefficients. */
    TPZDiffMatrix<T> & operator=(const TPZDiffMatrix<T> & source);
	
    /** @brief Adds element by element. */
    void Add(TPZDiffMatrix<T> & matrix, const T & scale = T(1.));
	
    /** @brief Matrix data access */
    T & operator()(const long i, const long j = 0);
	
    void PutVal(const long row,const long col,const T & value );
	
	const T &GetVal(const long row,const long col ) const;
	
    /** @brief Transposes the matrix onto the parameter object. */
	/** Resizes it if necessary. */
    TPZDiffMatrix<T> & Transpose(TPZDiffMatrix<T> & matrix);
	
    /**
     * @brief Performs a specific diffusive divergence operation.
     * @param gradx: example: {dU0dx0 dU0dx1 dU0dx2 dU1dx0... dU5dx2}
     * @param varOffset: shall lie between 0 and dim-1. \n
     * Represents the index of spatial derivative to multiply
     * the matrix by.
     * @param dim
     * @param Divergent: vetor to which the operation shall contribute to. \n
     * Must be explicitly zeroed before calling this function.
     */
	/**
	 * The object gradx contain the tangent matrix
     * of the solutions with respect to the dim spatial dimensions.
	 */
    void AddAlignDiv(TPZVec<T> & gradx,
                     const long varOffset,const long dim,
					 TPZVec<T> & Divergent);
	
    /**
     * @brief Computes the divergent for diffusion purposes
     * @param dPhi: the dim number of derivatives of test function
     * @param U the dim+2 solutions at the given point
     * @param dim
     * @param Divergent Vector result containing the computed divergent
     * The operation is additive, so zero it first
     * @param dDivergent Computes an approximate derivative of the divergent
     * with respect to the Ui coefficients.
     */
    void AddDiv(T & dPhi, TPZVec<T> & U,
                const long dim, TPZVec<T> & Divergent,
				TPZDiffMatrix<T> &  dDivergent);
	
    long Cols()const;
	
    long Rows()const;
	
    EStatus Decompose_LU();
	
    EStatus Substitution( TPZDiffMatrix<T> *B ) const;
	
	void Reset();
	
private:
	
    long index(const long i, const long j)const;
	
    long fRows, fCols;
	
    T * fStore;
	
	int fDecomposed;
	
};

/** @brief Re-implements << operator to output of matrices */
template <class T>
std::ostream & operator<<(std::ostream & out, TPZDiffMatrix<T> & A)
{
	long i, j;
	out << "\nTPZDiffMatrix<> " << A.Rows() << " * " << A.Cols();
	for(i=0;i<A.Rows();i++)
	{
		out << "\n\t";
		for(j=0;j<A.Cols();j++)
		{
			out << " " << A(i,j);
		}
	}
    out << std::endl;
	
	return out;
}

template <class T>
inline TPZDiffMatrix<T>::TPZDiffMatrix():fRows(0), fCols(0), fStore(NULL), fDecomposed(ENoDecompose)
{
}

template <class T>
inline TPZDiffMatrix<T>::TPZDiffMatrix(const long rows, const long cols):fRows(0), fCols(0), fStore(NULL), fDecomposed(ENoDecompose)
{
	Redim(rows, cols);
}

template <class T>
inline TPZDiffMatrix<T>::~TPZDiffMatrix()
{
	if(fStore)delete[] fStore;
}


template <class T>
inline void TPZDiffMatrix<T>::Redim(const long rows, const long cols)
{
	if(rows!=fRows || cols!=fCols){
		if(fStore)delete[] fStore;
		fStore = NULL;
		
		fRows = rows;
		fCols = cols;
		fStore = new T[fRows*fCols];
	}
	
	long i=fRows*fCols -1;
	for(;i>=0;i--)fStore[i]=T(0.);
}

template <class T>
inline  TPZDiffMatrix<T>& TPZDiffMatrix<T>::Transpose( TPZDiffMatrix<T> & matrix)
{
	matrix.Redim(fCols, fRows);
	long i, j;
	for(i=0;i<fRows;i++)
		for(j=0;j<fCols;j++)
			matrix(j,i)=operator()(i,j);
	return matrix;
}

template <class T>
inline  long TPZDiffMatrix<T>::index(const long i, const long j)const
{
#ifdef DEBUG
	if(i<0 || i>=fRows)
    {
        PZError << "\nTPZDiffMatrix<T>::index error: row out of bounds\n";
        DebugStop();
    }
	if(j<0 || j>=fCols)
    {
        PZError << "\nTPZDiffMatrix<T>::index error: col out of bounds\n";
        DebugStop();
    }
#endif
	return i*fCols+j;
}

template <class T>
inline  T & TPZDiffMatrix<T>::operator()(const long i, const long j)
{
	return fStore[index(i,j)];
}

template <class T>
inline void TPZDiffMatrix<T>::PutVal(const long row,const long col,const T & value )
{
	fStore[index(row,col)] = value;	
}


template <class T>
inline const T & TPZDiffMatrix<T>::GetVal(const long row,const long col ) const
{
	return fStore[index(row,col)];
}


template <class T>
inline void TPZDiffMatrix<T>::Multiply(TPZVec<T> & In, TPZVec<T> & Out, const T & scale)
{
	if(In.NElements()!=fCols)PZError << "\nTPZDiffMatrix<T>::Multiply error: 'In' vector out of bounds\n";
	if(Out.NElements()!=fRows)Out.Resize(fRows);
	long i, j;
	for(i=0;i<fRows;i++)
	{
		T buff=0.;
		for(j=0;j<fCols;j++)
		{
			buff+=In[j]*operator()(i,j);
		}
		buff*=scale;
		Out[i]=buff;
	}
}

template <class T>
inline void TPZDiffMatrix<T>::Multiply(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale)
{
	if(fCols!=In.fRows)PZError << "\nTPZDiffMatrix<T>::Multiply error: 'In' matrix out of bounds\n";
	Out.Redim(fRows,In.fCols);
	long i, j, k;
	for(i=0;i<fRows;i++)
	{
		for(j=0;j<In.fCols;j++)
		{
			T buff=0.;
			for(k=0;k<fCols;k++)
			{
				buff+=operator()(i,k) * In(k,j);
			}
			buff*=scale;
			Out(i,j)=buff;
		}
	}
}


template <class T>
inline void TPZDiffMatrix<T>::MultiplyAdd(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale)
{
	if(fCols!=In.fRows)PZError << "\nTPZDiffMatrix<T>::MultiplyAdd error: 'In' matrix out of bounds\n";
	if(Out.fCols!=fCols)PZError << "\nTPZDiffMatrix<T>::MultiplyAdd error: 'Out' matrix out of bounds\n";Out.Redim(fRows,In.fCols);
	if(In.fRows!=Out.fRows)PZError << "\nTPZDiffMatrix<T>::MultiplyAdd error: 'Out' matrix out of bounds\n";Out.Redim(fRows,In.fCols);
	long i, j, k;
	for(i=0;i<fRows;i++)
	{
		for(j=0;j<In.fCols;j++)
		{
			T buff=0.;
			for(k=0;k<fCols;k++)
			{
				buff+=operator()(i,k) * In(k,j);
			}
			buff*=scale;
			Out(i,j)+=buff;
		}
	}
}

template <class T>
inline void TPZDiffMatrix<T>::Add(TPZDiffMatrix<T> & matrix, const T & scale)
{
	if(fRows!=matrix.fRows)PZError << "\nTPZDiffMatrix<T>::Add error: row out of bounds\n";
	if(fCols!=matrix.fCols)PZError << "\nTPZDiffMatrix<T>::add error: col out of bounds\n";
	long i = fRows * fCols - 1;
	for(;i>=0;i--)fStore[i]+=matrix.fStore[i]*scale;
}

template <class T>
inline TPZDiffMatrix<T> & TPZDiffMatrix<T>::operator=(const TPZDiffMatrix<T> & source)
{
	Redim(source.fRows, source.fCols);
	long i = fRows * fCols - 1;
	for(;i>=0;i--)fStore[i]=source.fStore[i];
	fDecomposed = source.fDecomposed;
	return *this;
}

template <class T>
inline void TPZDiffMatrix<T>::AddAlignDiv(TPZVec<T> & gradx, const long varOffset,const long dim, TPZVec<T> & Divergent)
{
	long nstate = dim+2;
	if(gradx.NElements()!=(nstate)*dim)PZError << "\nTPZDiffMatrix<T>::AddAlignDiv error: gradx vector of wrong size\n";
	if(Divergent.NElements()!=nstate  )PZError << "\nTPZDiffMatrix<T>::AddAlignDiv error: Divergent vector of wrong size\n";
	//Divergent.Resize(nstate);
	long i, j;
	for(i=0;i<nstate;i++)
	{
		for(j=0;j<nstate;j++)Divergent[i]+=operator()(i,j)*gradx[j*dim+varOffset];
	}
}

template <class T>
inline void TPZDiffMatrix<T>::AddDiv(T & dPhi, TPZVec<T> & U, const long dim, TPZVec<T> & Divergent, TPZDiffMatrix<T> &  dDivergent)
{
	long nstate = dim+2;
	if( U.NElements()!=nstate)PZError << "\nTPZDiffMatrix<T>::AddDiv error: U vector of wrong size\n";
	if(Divergent.NElements()!=nstate)PZError << "\nTPZDiffMatrix<T>::AddDiv error: Divergent vector of wrong size\n";
	
	long i, j;
	for(i=0;i<nstate;i++)
	{
		for(j=0;j<nstate;j++)
		{
			Divergent[i]+=operator()(i,j)*dPhi*U[j];
			dDivergent(i,j)+=operator()(i,j)*dPhi;//?
		}
	}
}

template <class T>
inline long TPZDiffMatrix<T>::Cols()const
{
	return fCols;
}


template <class T>
inline long TPZDiffMatrix<T>::Rows()const
{
	return fRows;
}

template <class T>
inline EStatus TPZDiffMatrix<T>::Decompose_LU() {
	
	if (fDecomposed == ELU)
	{
		PZError << "\npzdiffmatrix.h Attempting to decompose an already decomposed matrix\n";
		return EOk;
	}
	
	T nn, pivot;
	long  min = ( Cols() < (Rows()) ) ? Cols() : Rows();
	
	for ( long k = 0; k < min ; k++ ) {
		if (IsZero( pivot = GetVal(k, k))) return EZeroPivot;
		for ( long i = k+1; i < Rows(); i++ ) {
			nn = GetVal( i, k ) / pivot;
			PutVal( i, k, nn );
			for ( long j = k+1; j < Cols(); j++ ) PutVal(i,j,GetVal(i,j)-nn*GetVal(k,j));
		}
	}
	fDecomposed=ELU;
	return EOk;
}

template <class T>
inline EStatus TPZDiffMatrix<T>::Substitution( TPZDiffMatrix<T> *B ) const{
	
    long rowb = B->Rows();
    long colb = B->Cols();
    if ( rowb != Rows() )
		return EIncompDim;
	long i;
    for ( i = 0; i < rowb; i++ ) {
        for ( long col = 0; col < colb; col++ ) {
            for ( long j = 0; j < i; j++ ) {
                B->PutVal( i, col, B->GetVal(i, col)-GetVal(i, j) * B->GetVal(j, col) );
            }
        }
    }
    for (long col=0; col<colb; col++) {
        for ( i = rowb-1; i >= 0; i-- ) {
            for ( long j = i+1; j < rowb ; j++ ) {
                B->PutVal( i, col, B->GetVal(i, col) -
						  GetVal(i, j) * B->GetVal(j, col) );
            }
            if ( IsZero( GetVal(i, i) ) ) {
				return EZeroPivot;
            }
            B->PutVal( i, col, B->GetVal( i, col) / GetVal(i, i) );
		}
    }
    return EOk;
}

template <class T>
inline void TPZDiffMatrix<T>::Reset()
{
	fDecomposed = ENoDecompose; 
}


#endif
