/**
 * @file
 * @brief Contains TPZGenMatrix class which implements generic class which holds a matrix of objects.
 */

#ifndef SHTMATRIXHPP
#define SHTMATRIXHPP

#include <iostream>

/** \addtogroup matrix
 * @{
 */
/**
 * @brief Implements generic class which holds a matrix of objects. \ref matrix "Matrix"
 */
template <class TObj>
class TPZGenMatrix {
	
protected:
	/** @brief Pointer to matrix */
	TObj *fMem;
	/** @brief Number of rows and columns */
	long   fRows, fCols;
	
public:
	
	/** @brief Constructor creating Null matrix */
	TPZGenMatrix ();
	/** @brief Constructor creating a rows x columns matrix */
	TPZGenMatrix (const long rows ,const long columns);
	
	/**
	 * @brief Constructor copy initializer
	 * @warning ***** WARNING *****
	 * matrices created with copy initializer are always temporary, eg, they
	 * share the same storage with another matrix
	 */
	TPZGenMatrix (const TPZGenMatrix & A);
	
	~TPZGenMatrix();
	
	void Print (const char *mess,std::ostream & out = std::cout) const;
	
	long Rows() const {return fRows;}
	
	long Cols() const {return fCols;}
	
	void Resize(const long newrow,const long newcol);
	
	//DEFINE OPERATORS
	
	TPZGenMatrix<TObj>& operator=(const TPZGenMatrix<TObj> & rval);	//makes this.p = rval.p
	
	
	TObj & operator()(const long row,const long column = 0) const; //get value, set value
	
};

/** 
 * @brief Implements a generic matrix of objects which implement arithmetic operations. \ref matrix "Matrix"
 */
template <class TObj>
class TPZGenAMatrix : public TPZGenMatrix<TObj> {
	
	TPZGenAMatrix () : TPZGenMatrix<TObj>() {}				//creates NULL matrix
	
	TPZGenAMatrix (const long rows ,const long columns) : TPZGenMatrix<TObj>(rows,columns) {}	//creates a rowsXcolumns matrix
	
	TPZGenAMatrix operator+(const TPZGenAMatrix<TObj> & rval) const;
	
	TPZGenAMatrix operator+(const TObj x) const;
	
	//  friend TPZGenAMatrix operator+(const TObj x, const TPZGenAMatrix<TObj>& rval); // {return (rval+x);}
	
	TPZGenAMatrix operator-(const TPZGenAMatrix<TObj> & rval) const;
	
	TPZGenAMatrix operator-(const TObj rval) const;
	
	//  friend TPZGenAMatrix operator-(const TObj x, const TPZGenAMatrix<TObj>& rval); // {return ((-rval)+x);}
	
	TPZGenAMatrix operator-() const;
	
	TPZGenAMatrix operator*(const TPZGenAMatrix<TObj> & rval) const;
	
	TPZGenAMatrix operator*(const TObj rval) const;
	
	//  friend TPZGenAMatrix operator*(const TObj x, const TPZGenAMatrix<TObj>& rval); // {return (rval*x);}
	
	TPZGenAMatrix<TObj>& operator+=(const TPZGenAMatrix<TObj> & rval);
	
	TPZGenAMatrix<TObj>& operator+=(const TObj x);
	
	TPZGenAMatrix<TObj>& operator-=(const TPZGenAMatrix<TObj> & rval);
	
	TPZGenAMatrix<TObj>& operator-=(const TObj x);
	
	TPZGenAMatrix<TObj>& operator*=(const TObj x);
	
	
	TPZGenAMatrix<TObj> Transpose() const;
	
	
};

/** @} */

#endif

