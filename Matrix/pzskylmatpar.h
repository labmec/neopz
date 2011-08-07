/**
 * @file
 * @brief Contains TPZSkylParMatrix class which implements a skyline storage format to parallelized process.
 */
/**
	Variable settings for PZProject
	by Longhin 14/12/2000

	The Variable PZNTPAR must be defined on a VC++ project,
	In order for a PZProject work properly
	the variable wraps up the following variables :
	__WIN32__,_MBCS,__NT__,OMNITHREAD,_WINSTATIC,_MT,__VC__
	Also the library LIBCD.lib must be excluded.
	It is done putting it on the ignore libraries field under
	Project Settings window.

	PZNTPAR is referenced on "omnithread.h"
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
class TPZSkylParMatrix : public TPZSkylMatrix
{
public:
	static int main();
	static int main_nada();
	TPZSkylParMatrix();
	TPZSkylParMatrix(const int dim);
	TPZSkylParMatrix(const int dim, const TPZVec<int> &skyline,int NumThreads);
	//TPZSkylParMatrix (const int dim,int NumThreads);
	TPZSkylParMatrix(const TPZSkylParMatrix &A);
	
	CLONEDEF(TPZSkylParMatrix)
    //    : TPZMatrix(A.Dim(), A.Dim()), fElem(0), fStorage(0) {Copy(A); }
	
	virtual ~TPZSkylParMatrix();
	///Parallel procedure using pthreads
	///Implement all data structure used in procedure.
	// @{
	int Decompose_Cholesky(std::list<int> &singular);
	int Decompose_Cholesky();
	
	int Decompose_LDLt(std::list<int> &singular);
	int Decompose_LDLt();
	
	void SetSkyline(const TPZVec<int> &skyline);
	// @}
	
private:
	
	static void *ParallelLDLt(void *t);
	static void *ParallelLDLt2(void *t);
	static void *ParallelCholesky(void *t);
	
	/** @brief Determine which column can be decomposed with respect to which column
	 */
	void ColumnToWork(int &lcol, int &lprevcol);
	/// Determine which column has some equations to decompose
	void ColumnToWork(int &lcol);
	void DecomposeColumnCholesky(int lcol, int lprevcol);
	void DecomposeColumnLDLt(int lcol, int lprevcol);
	void DecomposeColumnLDLt2(int lcol);
	void PrintState();
public:
	TPZVec<int> fDec;  
	TPZVec<int> fSkyline;
private:
	int fEqDec, fNthreads;
	int * fThreadUsed;
	std::set<int> fColUsed;
	int fNumDecomposed;
	bool fCorrectSingular;
	std::list<int> fSingular;
	
#ifdef OOPARLIB
	
	//virtual long GetClassID() const    { return TSKYMATRIX_ID; }
	virtual int Unpack( TReceiveStorage *buf );
	static TSaveable *Restore(TReceiveStorage *buf);
	virtual int Pack( TSendStorage *buf ) const;
	virtual char *ClassName() const   { return( "TPZSkylParMatrix"); }
	virtual int DerivedFrom(const long Classid) const;
	virtual int DerivedFrom(const char *classname) const; // a class with name classname
	
#endif
	
};

#endif
