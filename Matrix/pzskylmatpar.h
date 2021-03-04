/**
 * @file
 * @author Gustavo int64_thin
 * @since 14/12/2000
 * @brief Contains TPZSkylParMatrix class which implements a skyline storage format to parallelized process.
 * @note The Variable PZNTPAR must be defined on a VC++ project, it is referenced on "omnithread.h"
 * @note In order for a PZProject work properly the variable wraps up the following variables: \n
 * __WIN32__,_MBCS,__NT__,OMNITHREAD,_WINSTATIC,_MT,__VC__ \n
 * Also the library LIBCD.lib must be excluded.
 */

#ifndef TSKYLPARMATH
#define TSKYLPARMATH

#include "pzmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzskylmat.h"


#include <stdlib.h>
#include <signal.h>
#include <time.h>
#include <set>

/**
 * @brief Implements a skyline storage format to parallelized process. \ref matrix "Matrix"
 * @ingroup matrix
 */
template<class TVar>
class TPZSkylParMatrix : public TPZSkylMatrix<TVar>
{
public:
	/** @brief Static main for testing */
	static int main();
	static int main_nada();
	
	/** @brief Default constructor */
	TPZSkylParMatrix();
	/** @brief Constructor for given dimension */
	TPZSkylParMatrix(const int64_t dim);
	/** @brief Constructor with number of threads */
	TPZSkylParMatrix(const int64_t dim, const TPZVec<int64_t> &skyline,int NumThreads);
	/** @brief Copy constructor */
	TPZSkylParMatrix(const TPZSkylParMatrix<TVar> &A);
	
	CLONEDEF(TPZSkylParMatrix)
    //    : TPZMatrix(A.Dim(), A.Dim()), fElem(0), fStorage(0) {Copy(A); }
	
	/** @brief Default destructor */
	virtual ~TPZSkylParMatrix();
	
	/**
	 * @name Parallel procedure using pthreads
	 * @brief Implement all data structure used in procedure.
	 * @{
	 */
	
	int Decompose_Cholesky(std::list<int64_t> &singular) override;
	int Decompose_Cholesky() override;
	
	int Decompose_LDLt(std::list<int64_t> &singular) override;
	int Decompose_LDLt() override;
	
	void SetSkyline(const TPZVec<int64_t> &skyline);

	/** @} */
    public:
int ClassId() const override;

private:
	
	static void *ParallelLDLt(void *t);
	static void *ParallelLDLt2(void *t);
	static void *ParallelCholesky(void *t);
	
	/** @brief Determine which column can be decomposed with respect to which column */
	void ColumnToWork(int64_t &lcol, int64_t &lprevcol);
	/** @brief Determine which column has some equations to decompose */
	void ColumnToWork(int64_t &lcol);
	void DecomposeColumnCholesky(int64_t lcol, int64_t lprevcol);
	void DecomposeColumnLDLt(int64_t lcol, int64_t lprevcol);
	void DecomposeColumnLDLt2(int64_t lcol);
	void PrintState();
	
public:
	TPZVec<int> fDec;  
	TPZVec<int64_t> fSkyline;
	
private:
	int64_t fEqDec, fNthreads;
	int * fThreadUsed;
	std::set<int64_t> fColUsed;
	int fNumDecomposed;
	bool fCorrectSingular;
	std::list<int64_t> fSingular;
	
#ifdef OOPARLIB
	
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *CreateInstance(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSkylParMatrix"); }
	virtual int DerivedFrom(const int64_t Classid) const;
	virtual int DerivedFrom(const char *classname) const; // a class with name classname
	
#endif
	
};

#endif
