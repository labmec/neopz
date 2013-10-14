/**
 * @file
 * @author Gustavo Longhin
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

//parallel libraries
#include <pthread.h>

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
	TPZSkylParMatrix(const long dim);
	/** @brief Constructor with number of threads */
	TPZSkylParMatrix(const long dim, const TPZVec<long> &skyline,int NumThreads);
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
	
	int Decompose_Cholesky(std::list<long> &singular);
	int Decompose_Cholesky();
	
	int Decompose_LDLt(std::list<long> &singular);
	int Decompose_LDLt();
	
	void SetSkyline(const TPZVec<long> &skyline);

	/** @} */
	
private:
	
	static void *ParallelLDLt(void *t);
	static void *ParallelLDLt2(void *t);
	static void *ParallelCholesky(void *t);
	
	/** @brief Determine which column can be decomposed with respect to which column */
	void ColumnToWork(long &lcol, long &lprevcol);
	/** @brief Determine which column has some equations to decompose */
	void ColumnToWork(long &lcol);
	void DecomposeColumnCholesky(long lcol, long lprevcol);
	void DecomposeColumnLDLt(long lcol, long lprevcol);
	void DecomposeColumnLDLt2(long lcol);
	void PrintState();
	
public:
	TPZVec<int> fDec;  
	TPZVec<long> fSkyline;
	
private:
	long fEqDec, fNthreads;
	int * fThreadUsed;
	std::set<long> fColUsed;
	int fNumDecomposed;
	bool fCorrectSingular;
	std::list<long> fSingular;
	
#ifdef OOPARLIB
	
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual std::string ClassName() const   { return( "TPZSkylParMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const; // a class with name classname
	
#endif
	
};

#endif
