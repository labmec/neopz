
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


#include <pzmatrix.h>
#include <pzvec.h>
#include <pzmanvector.h>
#include <pzskylmat.h>

#ifdef PZNTPAR
	#include <omnithread.h>
#endif

#include <stdlib.h>
#include <signal.h>
//#include <unistd.h>
#include <time.h>

//parallel libraries
#ifndef PZNTPAR
	#include <pthread.h>
#endif
#ifdef OOPARLIB

#include "pzsaveable.h"
#include "pzmatdefs.h"

#endif


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
    //    : TPZMatrix(A.Dim(), A.Dim()), fElem(0), fStorage(0) {Copy(A); }

  virtual ~TPZSkylParMatrix();
  //Parallel procedure using pthreads
  //Implement all data structure used in procedure.
#ifndef PZNTPAR
  int Decompose_Cholesky();//, TPZVec<int> fDec, TPZVec<int> fThreadUsed, int fEqDec, int nthreads, int neq);
#endif
#ifdef OMNITHREAD
  int Decompose_Cholesky();
#endif

  void SetSkyline(const TPZVec<int> &skyline);

 private:

#ifdef OMNITHREAD
  static void *OmniParallelCholesky(void *t);
#endif

#ifndef PZNTPAR 
  static void *ParallelCholesky(void *t);
#endif

  void ColumnToWork(int &lcol, int &lprevcol);
  void DecomposeColumn(int lcol, int lprevcol);
  void PrintState();
 public:
  TPZVec<int> fDec;  
  TPZVec<int> fSkyline;
 private:
  int fEqDec, fNthreads;
  int * fThreadUsed;

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
