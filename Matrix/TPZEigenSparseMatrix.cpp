/**
 * @file
 * @brief Contains the implementation of the TPZEigenSparseMatrix methods.
 */

#include "TPZEigenSparseMatrix.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <memory.h>
#include <string>
#include <map>
#include <thread>
#include <vector>
#include "tpzverysparsematrix.h"
#include "pzstack.h"

using namespace std;

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************
// ****************************************************************************
//
// Constructor
//
// ****************************************************************************


template<class TVar>
TPZEigenSparseMatrix<TVar>::TPZEigenSparseMatrix() : TPZRegisterClassId(&TPZEigenSparseMatrix::ClassId),
TPZFYsmpMatrix<TVar>()
{
}





template<class TVar>
TPZEigenSparseMatrix<TVar>::TPZEigenSparseMatrix(const int64_t rows,const int64_t cols ) :
TPZRegisterClassId(&TPZEigenSparseMatrix::ClassId),TPZFYsmpMatrix<TVar>(rows,cols) {
#ifdef CONSTRUCTOR
	cerr << "TPZEigenSparseMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZEigenSparseMatrix<TVar>::~TPZEigenSparseMatrix() {
	// Deletes everything associated with a TPZEigenSparseMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZEigenSparseMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************


template<class TVar>
void TPZEigenSparseMatrix<TVar>::CheckTypeCompatibility(const TPZMatrix<TVar>*A, const TPZMatrix<TVar>*B)const
{
  auto incompatSparse = [](){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: incompatible matrices\n.Aborting...\n";
    DebugStop();
  };
	auto aPtr = dynamic_cast<const TPZEigenSparseMatrix<TVar>*>(A);
  auto bPtr = dynamic_cast<const TPZEigenSparseMatrix<TVar>*>(B);
  if(!aPtr || !bPtr){
    incompatSparse();
  }
	bool check{false};
	const auto nIA = aPtr->fIA.size();
	for(auto i = 0; i < nIA; i++){
		check = check || aPtr->fIA[i] != bPtr->fIA[i];
	}

	const auto nJA = aPtr->fJA.size();
	for(auto i = 0; i < nJA; i++){
		check = check || aPtr->fJA[i] != bPtr->fJA[i];
	}
	if(check) incompatSparse();
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator+(const TPZEigenSparseMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
    res += mat;
	return res;
}
template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator-(const TPZEigenSparseMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
    res -= mat;
	return res;
}

template<class TVar>
TPZEigenSparseMatrix<TVar> TPZEigenSparseMatrix<TVar>::operator*(const TVar alpha) const
{
	auto res(*this);
	for(auto &el : res.fA) el *= alpha;
	return res;
}

template<class TVar>
TPZEigenSparseMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator+=(const TPZEigenSparseMatrix<TVar> &A )
{
	TPZEigenSparseMatrix<TVar> res((*this)+A);
	*this = res;
	return *this;
}
template<class TVar>
TPZEigenSparseMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator-=(const TPZEigenSparseMatrix<TVar> &A )
{
	TPZEigenSparseMatrix<TVar> res((*this)-A);
	*this = res;
	return *this;
}
template<class TVar>
TPZMatrix<TVar> &TPZEigenSparseMatrix<TVar>::operator*=(const TVar val)
{
	TPZEigenSparseMatrix<TVar> res((*this)*val);
	*this = res;
	return *this;
}


// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZEigenSparseMatrix<TVar>::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
    TPZFYsmpMatrix<TVar>::Print(title,out,form);
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************




static int64_t  FindCol(int64_t *colf, int64_t *coll, int64_t col)
{
	if(col == *colf) return 0;
	int64_t *begin = colf;
	int64_t *end = coll;
	while (begin != end)
	{
		int64_t dist = (end-begin)/2;
		int64_t *mid = begin+dist;
		if(*mid == col) return (mid-colf);
		else if(*mid > col) end=mid;
		else begin = mid;
	}
	return -1;
}

/// this is a class that doesn't implement direct decompostion
/** @brief decompose the system of equations acording to the decomposition
 * scheme */
template<class TVar>
int TPZEigenSparseMatrix<TVar>::Decompose(const DecomposeType dt)
{
    TPZVec<int> IA;
    for(auto it : fIA.begin(), auto it2:IA.begin(); it != fIA.end(); 
    TPZVec<int> JA;
    int fJA[] = {0};
    if(!fEigenMatrix){
        fEigenMatrix = new EigenSparse(this->fRow,this->fCol,this->fJA.size(),&fIA[0],&fJA[0],&this->fA[0]);
    }
}
/**
 * @brief Solves the linear system using Direct methods
 * @param F The right hand side of the system and where the solution is stored.
 * @param dt Indicates type of decomposition
 */
template<class TVar>
int TPZEigenSparseMatrix<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt)
{
    
}
template<class TVar>
int TPZEigenSparseMatrix<TVar>::SolveDirect ( TPZFMatrix<TVar>& F , const DecomposeType dt) const
{
    
}





template<class TVar>
int TPZEigenSparseMatrix<TVar>::ClassId() const{
    return Hash("TPZEigenSparseMatrix") ^ TPZFYsmpMatrix<TVar>::ClassId() << 1;
}
template class TPZEigenSparseMatrix<long double>;
template class TPZEigenSparseMatrix<double>;
template class TPZEigenSparseMatrix<float>;
template class TPZEigenSparseMatrix<std::complex<long double>>;
template class TPZEigenSparseMatrix<std::complex<double>>;
template class TPZEigenSparseMatrix<std::complex<float>>;
