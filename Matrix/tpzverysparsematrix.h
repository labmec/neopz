
#ifndef TPZVERYSPARSEMATRIX_H
#define TPZVERYSPARSEMATRIX_H

#include <iostream>
#include <map>

#include "pzmatrix.h"
#include "pzfmatrix.h"
class TPZFYsmpMatrix;

/**
	@author Agnaldo Monteiro Farias <agnaldo@tuborg>
*/
class TPZVerySparseMatrix: public TPZMatrix {
public:

	friend class TPZFYsmpMatrix; 
	
    TPZVerySparseMatrix();

    virtual ~TPZVerySparseMatrix();
  
    TPZVerySparseMatrix(const TPZVerySparseMatrix &copy) : TPZMatrix(copy), fExtraSparseData(copy.fExtraSparseData)
    {
    }
	
    /** 
	 * Put values checking bounds\n
    */
	
    virtual void PutVal(const int row, const int col, REAL val);
	
	/** 
	 * Get values checking bounds\n
    */
    virtual const REAL &GetVal(const int row, const int col) const;
 
    /** 
     The operators check on the bounds if the DEBUG variable is defined
        @param row Row number.
        @param col Column number.
   */
    virtual REAL &s(const int row, const int col)
    {
        #ifdef DEBUG
        if(row >= Rows() || row<0 || col >= Cols() || col<0)
        {
            Error("TPZFMatrix::operator() "," Index out of bounds");
            DebugStop();
            gZero = 0.;
            return gZero;
        }
        #endif
        return fExtraSparseData[std::pair<int,int>(row,col)];
    }

    
    CLONEDEF(TPZVerySparseMatrix)
    
    /**
    * It computes z = beta * y + alpha * opt(this)*x but z and x can not overlap in memory.
      * @param x Is x on the above operation
      *@param y Is y on the above operation
      * @param z Is z on the above operation
      * @param alpha Is alpha on the above operation
      * @param beta Is beta on the above operation
      * @param opt Indicates if is Transpose or not
      * @param stride Indicates n/N where n is dimension of the right hand side vector and N is matrix dimension
    */
    virtual void MultAdd(const TPZFMatrix & x,const TPZFMatrix & y, TPZFMatrix & z,
                        const REAL alpha=1., const REAL beta = 0., const int opt = 0, const int stride = 1 );
	
	    
protected:

    /**
    Save elements different from zero, of Sparse matrix
    */
    std::map <std::pair<int, int>, REAL> fExtraSparseData;


};

#endif
