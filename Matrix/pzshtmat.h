//HEADER FILE FOR CLASS MATRIX


#ifndef SHTMATRIXHPP
#define SHTMATRIXHPP

#include <iostream>

/// generic class which holds a matrix of objects
template <class TObj>
class TPZGenMatrix {

 protected:
  TObj *fMem;  	//pointer to matrix
  int   fRows, fCols;	//number of rows and columns

 public:

  TPZGenMatrix ();				//creates NULL matrix

  TPZGenMatrix (const int rows ,const int columns);	//creates a rowsXcolumns matrix



  //***** WARNING *****
  // matrices created with copy initializer are always temporary, eg, they
  // share the same storage with another matrix

  TPZGenMatrix (const TPZGenMatrix & A);			//copy initializer

  ~TPZGenMatrix();

  void Print (const char *mess,std::ostream & out = std::cout) const;

  int Rows() const {return fRows;}

  int Cols() const {return fCols;}

  void Resize(const int newrow,const int newcol);

  //DEFINE OPERATORS

  TPZGenMatrix<TObj>& operator=(const TPZGenMatrix<TObj> & rval);	//makes this.p = rval.p


  TObj & operator()(const int row,const int column = 0) const; //get value, set value

};

/// A generic matrix of objects which implement arithmetic operations
template <class TObj>
class TPZGenAMatrix : public TPZGenMatrix<TObj> {

  TPZGenAMatrix () : TPZGenMatrix<TObj>() {}				//creates NULL matrix

  TPZGenAMatrix (const int rows ,const int columns) : TPZGenMatrix<TObj>(rows,columns) {}	//creates a rowsXcolumns matrix

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


#endif

