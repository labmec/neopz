#ifndef PZDIFFMATRIX_H
#define PZDIFFMATRIX_H

#include <ostream>

/**
 * Matrix class to hold the flux derivatives A B C and
 * diffusive matrix coefficients.
 * @author Erick Slis
 * @author Cedric Ayala
 * @since June 1, 2003.
 */

template <class T>
class TPZDiffMatrix
{
  public:
    TPZDiffMatrix();
    TPZDiffMatrix(const int rows, const int cols);
    ~TPZDiffMatrix();

    /**
     * Resizes and zeroes the matrix.
     *
     */
    void Redim(const int rows, const int cols);

    /**
     * Multiplies the matrix by a correspondent TPZVec vector;
     * Dimensions are checked.
     */
    void Multiply(TPZVec<T> & In, TPZVec<T> & Out, const T & scale = T(1.));


    /**
     * Matrix multiplication;
     * Dimensions are checked.
     */
    void Multiply(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale = T(1.));

    /**
     * Matrix multiplication;
     * Dimensions are checked.
     * Results are additively contributed to Out
     */
    void MultiplyAdd(TPZDiffMatrix<T> & In, TPZDiffMatrix<T> & Out, const T & scale = T(1.));


    /**
     * Copies the matrix, reallocating all coefficients.
     *
     */
    TPZDiffMatrix<T> & operator=(const TPZDiffMatrix<T> & source);

    /**
     * Adds element by element.
     *
     */
    void Add(TPZDiffMatrix<T> & matrix, const T & scale = T(1.));

    /**
     * Matrix data access
     */
    T & operator()(const int i, const int j);


    /**
     * Transposes the matrix onto the parameter object.
     * Resizes it if necessary.
     */
    TPZDiffMatrix<T> & Transpose(TPZDiffMatrix<T> & matrix);

    /**
     * performs a specific diffusive divergence operation.
     * The object gradx contain the tangent matrix
     * of the solutions with respect to the dim spatial
     * dimensions.
     *
     * @param gradx: example: {dU0dx0 dU0dx1 dU0dx2 dU1dx0... dU5dx2}
     * @param varOffset: shall lie between 0 and dim-1.
     * Represents the index of spatial derivative to multiply
     * the matrix by.
     * @param dim
     * @param Divergent: vetor to which the operation shall contribute to.
     * Must be explicitly zeroed before calling this function.
     */
    void AddAlignDiv(TPZVec<T> & gradx,
                     const int varOffset,const int dim,
		     TPZVec<T> & Divergent);

    /**
     * Computes the divergent for diffusion purposes
     * @param dPhi: the dim number of derivatives of test function
     * @param U the dim+2 solutions at the given point
     * @param dim
     * @param Divergent Vector result containing the computed divergent
     * The operation is additive, so zero it first
     * @param dDivergent Computes an approximate derivative of the divergent
     * with respect to the Ui coefficients.
     */
    void AddDiv(T & dPhi, TPZVec<T> & U,
                const int dim, TPZVec<T> & Divergent,
		TPZDiffMatrix<T> &  dDivergent);

    int Cols()const;

    int Rows()const;


  private:

    int index(const int i, const int j);

    int fRows, fCols;

    T * fStore;

};

template <class T>
ostream & operator<<(ostream & out, TPZDiffMatrix<T> & A)
{
   int i, j;
   out << "\nTPZDiffMatrix<> " << A.Rows() << " * " << A.Cols();
   for(i=0;i<A.Rows();i++)
   {
      out << "\n\t";
      for(j=0;j<A.Cols();j++)
         {
	    out << " " << A(i,j);
	 }
   }
    out << endl;

   return out;
}

template <class T>
inline TPZDiffMatrix<T>::TPZDiffMatrix():fRows(0), fCols(0), fStore(NULL)
{
}

template <class T>
inline TPZDiffMatrix<T>::TPZDiffMatrix(const int rows, const int cols):fRows(0), fCols(0), fStore(NULL)
{
   Redim(rows, cols);
}

template <class T>
inline TPZDiffMatrix<T>::~TPZDiffMatrix()
{
   if(fStore)delete[] fStore;
}


template <class T>
inline void TPZDiffMatrix<T>::Redim(const int rows, const int cols)
{
   if(rows!=fRows || cols!=fCols){
      if(fStore)delete[] fStore;
      fStore = NULL;

      fRows = rows;
      fCols = cols;
      fStore = new T[fRows*fCols];
      }

   int i=fRows*fCols -1;
   for(;i>=0;i--)fStore[i]=T(0.);
}

template <class T>
inline  TPZDiffMatrix<T>& TPZDiffMatrix<T>::Transpose( TPZDiffMatrix<T> & matrix)
{
   matrix.Redim(fCols, fRows);
   int i, j;
   for(i=0;i<fRows;i++)
      for(j=0;j<fCols;j++)
         matrix(j,i)=operator()(i,j);
   return matrix;
}

template <class T>
inline  int TPZDiffMatrix<T>::index(const int i, const int j)
{
   if(i<0 || i>=fRows)PZError << "\nTPZDiffMatrix<T>::index error: row out of bounds\n";
   if(j<0 || j>=fCols)PZError << "\nTPZDiffMatrix<T>::index error: col out of bounds\n";
   return i*fCols+j;
}

template <class T>
inline  T & TPZDiffMatrix<T>::operator()(const int i, const int j)
{
   return fStore[index(i,j)];
}

template <class T>
inline void TPZDiffMatrix<T>::Multiply(TPZVec<T> & In, TPZVec<T> & Out, const T & scale)
{
   if(In.NElements()!=fCols)PZError << "\nTPZDiffMatrix<T>::Multiply error: 'In' vector out of bounds\n";
   if(Out.NElements()!=fRows)Out.Resize(fRows);
   int i, j;
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
   int i, j, k;
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
   int i, j, k;
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
   int i = fRows * fCols - 1;
   for(;i>=0;i--)fStore[i]+=matrix.fStore[i]*scale;
}

template <class T>
inline TPZDiffMatrix<T> & TPZDiffMatrix<T>::operator=(const TPZDiffMatrix<T> & source)
{
   Redim(source.fRows, source.fCols);
   int i = fRows * fCols - 1;
   for(;i>=0;i--)fStore[i]=source.fStore[i];
   return *this;
}

template <class T>
inline void TPZDiffMatrix<T>::AddAlignDiv(TPZVec<T> & gradx, const int varOffset,const int dim, TPZVec<T> & Divergent)
{
   int nstate = dim+2;
   if(gradx.NElements()!=(nstate)*dim)PZError << "\nTPZDiffMatrix<T>::AddAlignDiv error: gradx vector of wrong size\n";
   if(Divergent.NElements()!=nstate  )PZError << "\nTPZDiffMatrix<T>::AddAlignDiv error: Divergent vector of wrong size\n";
   //Divergent.Resize(nstate);
   int i, j;
   for(i=0;i<nstate;i++)
   {
      for(j=0;j<nstate;j++)Divergent[i]+=operator()(i,j)*gradx[j*dim+varOffset];
   }
}

template <class T>
inline void TPZDiffMatrix<T>::AddDiv(T & dPhi, TPZVec<T> & U, const int dim, TPZVec<T> & Divergent, TPZDiffMatrix<T> &  dDivergent)
{
   int nstate = dim+2;
   if( U.NElements()!=nstate)PZError << "\nTPZDiffMatrix<T>::AddDiv error: U vector of wrong size\n";
   if(Divergent.NElements()!=nstate)PZError << "\nTPZDiffMatrix<T>::AddDiv error: Divergent vector of wrong size\n";

   int i, j;
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
inline int TPZDiffMatrix<T>::Cols()const
{
   return fCols;
}


template <class T>
inline int TPZDiffMatrix<T>::Rows()const
{
   return fRows;
}
/*
template <class T>
inline ostream & TPZDiffMatrix<T>::operator<<(ostream & out)
{
   int i, j;
   out << "\nTPZDiffMatrix<> " << fRows << " * " << fCols;
   for(i=0;i<fRows;i++)
   {
      out << "\n\t";
      for(j=0;j<fCols;j++)
         {
	    out << " " << operator()(i,j);
	 }
   }
    out << endl;

   return out;
}
*/
#endif
