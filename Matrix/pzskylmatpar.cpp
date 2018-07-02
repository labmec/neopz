/**
 * @file
 * @brief Contains the implementation of the TPZSkylParMatrix methods.
 */

#include <math.h>
#include <stdlib.h>

#include <stdio.h>  
#include <fstream> 
using namespace std;


#ifdef BLAS
extern "C"{
#include <g2c.h>
#include "fblaswr.h"
};

#endif

#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzskylmatpar.h"

#include "pzlog.h"

#include "pz_pthread.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzskylparmatrix"));
#endif

/** @brief Semaphore */
pthread_mutex_t skymutex = PTHREAD_MUTEX_INITIALIZER;
/** @brief Condition to waiting */
pthread_cond_t condition = PTHREAD_COND_INITIALIZER;

/** @brief Initializing number of iterations to Template SUM */
//const int templatedepth = 10;

//Constructors
template<class TVar>
TPZSkylParMatrix<TVar>::TPZSkylParMatrix(const int64_t dim, const TPZVec<int64_t> &skyline,int NumThreads)

: TPZRegisterClassId(&TPZSkylParMatrix::ClassId),
TPZSkylMatrix<TVar>(dim, skyline),fDec(dim), fCorrectSingular(false)
{
	int i;
	
	fSkyline = skyline;
	
	for(i=0;i<dim;i++) {
		fDec[i]=skyline[i]-1;
		
	}
	fNthreads = NumThreads;
	fEqDec = -1;
}

template<class TVar>
TPZSkylParMatrix<TVar>::TPZSkylParMatrix(const TPZSkylParMatrix<TVar> &copy) : TPZRegisterClassId(&TPZSkylParMatrix::ClassId),
TPZSkylMatrix<TVar>(copy), fDec(copy.fDec),fSkyline(copy.fSkyline),fEqDec(copy.fEqDec),
fNthreads(copy.fNthreads),fThreadUsed(0),fCorrectSingular(false)
{
}

template<class TVar>
TPZSkylParMatrix<TVar>::TPZSkylParMatrix()
: TPZRegisterClassId(&TPZSkylParMatrix::ClassId),
TPZSkylMatrix<TVar>(),fDec(0),fCorrectSingular(false)
{
	//something
}

// TPZSkylParMatrix::TPZSkylParMatrix(const TPZSkylParMatrix &cp):TPZSkylMatrix(cp){
//   fDec = cp.fDec;
//   fSkyline = cp.fSkyline;
//   fNthreads = cp.fNthreads;
//   fEqDec = cp.fEqDec;
// }

template<class TVar>
TPZSkylParMatrix<TVar>::~TPZSkylParMatrix() {
}

template<class TVar>
void TPZSkylParMatrix<TVar>::SetSkyline(const TPZVec<int64_t> &skyline)
{
	fSkyline = skyline;
	int64_t i,neq=this->Dim();
	for(i=0;i<neq;i++) fDec[i]=skyline[i]-1;
	fEqDec=-1;
	fCorrectSingular = false;
	TPZSkylMatrix<TVar>::SetSkyline(skyline);
}

template<class TVar>
void TPZSkylParMatrix<TVar>::PrintState() {
	int i;
	cout << "Decomposed equations " << fEqDec << ' ' << " Thread 'i' working on equation" << ' ';
	for(i=0; i<fNthreads; i++) cout << fThreadUsed[i] << ' ';
	cout << endl;
	cout << "Columns to work " << ' ';
	int64_t neq = this->Dim();
	for(i=0; i<neq; i++) cout << fDec[i] << ' ';
	cout << endl;
	cout << "Singular equations ";
	std::list<int64_t>::iterator it;
	for (it=fSingular.begin(); it!=fSingular.end(); it++) {
		cout << *it << ' ';
	}
	cout.flush();
}

template<class TVar>
void TPZSkylParMatrix<TVar>::ColumnToWork(int64_t &lcol, int64_t &lprevcol) {
	int64_t neq = this->Dim();
	for(lcol=fEqDec+1;lcol<neq;lcol++)
    {
		int i;
		for(i=0;i<fNthreads;i++){
			if(fThreadUsed[i]==lcol){
				//lcol=-1;
				//lprevcol=-1;
				break;
			}
		}
		if(i < fNthreads && fThreadUsed[i] == lcol) continue;
		if(lcol == fEqDec+1 && fDec[lcol] == fEqDec){
			lprevcol=lcol;
			return;
		}else if(fDec[lcol]<fEqDec){
			lprevcol = fDec[lcol] + 1;
			return;
		}
    }
	if(lcol == neq){
		lcol = -1;
	}
}

template<class TVar>
void TPZSkylParMatrix<TVar>::ColumnToWork(int64_t &lcol) {
	int64_t neq = this->Dim();
	int64_t lcolentry = lcol;
#ifdef PZDEBUG
	if(lcolentry < 0 || lcolentry >= neq)
	{
		LOGPZ_ERROR(logger,"ColumnToWork called with wrong lcol argument")
		DebugStop();
	}
#endif
	lcol++;
	lcol = lcol%neq;
	while(lcol != lcolentry)
    {
		if(fColUsed.find(lcol) != fColUsed.end() || fDec[lcol] == lcol) 
		{
			lcol++;
			lcol = lcol%neq;
			continue;
		}
		int decomposeduntil = fDec[lcol]+1;
		if(fDec[decomposeduntil] == decomposeduntil)
		{
#ifdef LOG4CXX
            if (logger->isDebugEnabled()) {
                std::stringstream sout;
                sout << "Will work column " << lcol;
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
			return;
		}
		lcol++;
		if(lcol == neq) lcol = 0;
    }
	if(lcol == lcolentry){
		lcol = -1;
	}
}


template<class TVar>
void * TPZSkylParMatrix<TVar>::ParallelCholesky(void *t) {
	
	TPZSkylParMatrix<TVar> *loc = (TPZSkylParMatrix<TVar> *) t;
	int64_t aux_col = -1;//, k;
	int64_t col, prevcol;
	prevcol=0;
	PZ_PTHREAD_MUTEX_LOCK(&skymutex, "TPZSkylParMatrix::ParallelCholesky()");
	int64_t neq = loc->Dim();
	while (loc->fEqDec < neq-1) {
		loc->ColumnToWork(col, prevcol);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fEqDec > loc->fNthreads){
				PZ_PTHREAD_COND_WAIT(&condition, &skymutex, "TPZSkylParMatrix::ParallelCholesky()");
			}else
			{
				//loc->fNthreads--;
				//cout << "Terminating thread " << pthread_self() << endl;
				break;
			}
		}
		else {
			
			//Registers the working column number
			int i;
			for(i=0;i<loc->fNthreads;i++){
				if(loc->fThreadUsed[i]==-1) {
					loc->fThreadUsed[i]=col;
					aux_col = i;
					break;
				}
			}
			PZ_PTHREAD_COND_SIGNAL(&condition, "TPZSkylParMatrix::ParallelCholesky()");
			PZ_PTHREAD_MUTEX_UNLOCK(&skymutex, "TPZSkylParMatrix::ParallelCholesky()");
			if(loc->fCorrectSingular)
			{
				loc->DecomposeColumn(col, prevcol,loc->fSingular);//loc->DecomposeColumn2(col, prevcol);
			}
			else {
				loc->DecomposeColumn(col, prevcol);
			}
			
			PZ_PTHREAD_MUTEX_LOCK(&skymutex, "TPZSkylParMatrix::ParallelCholesky()");
			loc->fDec[col]=prevcol;
			loc->fThreadUsed[aux_col]=-1;
			if (col==prevcol) {
				if(loc->fEqDec == col-1) {
					loc->fEqDec = col;
					if(!(loc->fEqDec%100) && loc->Dim() > 100) {
						cout << loc->fEqDec << ' ';
						cout.flush();
					}
					if(!(loc->fEqDec%1000)) cout << endl;
					cout.flush();
				}
				else cout << "BIG MISTAKE\n";
			}
			
		}
	}
	PZ_PTHREAD_MUTEX_UNLOCK(&skymutex, "TPZSkylParMatrix::ParallelCholesky()");
	PZ_PTHREAD_COND_BROADCAST(&condition, "TPZSkylParMatrix::ParallelCholesky()");
	cout << endl;
	cout.flush();
	return 0;
	
}

template<class TVar>
void * TPZSkylParMatrix<TVar>::ParallelLDLt(void *t) {
	
	TPZSkylParMatrix<TVar> *loc = (TPZSkylParMatrix<TVar> *) t;
	int64_t aux_col = -1;//, k;
	int64_t col, prevcol;
	prevcol=0;
	PZ_PTHREAD_MUTEX_LOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt()");
	int64_t neq = loc->Dim();
	while (loc->fEqDec < neq-1) {
		loc->ColumnToWork(col, prevcol);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fEqDec > loc->fNthreads){
				PZ_PTHREAD_COND_WAIT(&condition, &skymutex,"TPZSkylParMatrix::ParallelLDLt()");
			}else
			{
				//loc->fNthreads--;
				//cout << "Terminating thread " << pthread_self() << endl;
				break;
			}
		}
		else {
			
			//Registers the working column number
			int i;
			for(i=0;i<loc->fNthreads;i++){
				if(loc->fThreadUsed[i]==-1) {
					loc->fThreadUsed[i]=col;
					aux_col = i;
					break;
				}
			}
			PZ_PTHREAD_COND_SIGNAL(&condition,"TPZSkylParMatrix::ParallelLDLt()");
			PZ_PTHREAD_MUTEX_UNLOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt()");
			loc->DecomposeColumnLDLt(col, prevcol);//loc->DecomposeColumn2(col, prevcol);
			PZ_PTHREAD_MUTEX_LOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt()");
			loc->fDec[col]=prevcol;
			loc->fThreadUsed[aux_col]=-1;
			if (col==prevcol) {
				if(loc->fEqDec == col-1) {
					loc->fEqDec = col;
					if(!(loc->fEqDec%100) && loc->Dim() > 100) {
						cout << loc->fEqDec << ' ';
						cout.flush();
					}
					if(!(loc->fEqDec%1000)) cout << endl;
					cout.flush();
				}
				else cout << "BIG MISTAKE\n";
			}
			
		}
	}
	PZ_PTHREAD_MUTEX_UNLOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt()");
	PZ_PTHREAD_COND_BROADCAST(&condition,"TPZSkylParMatrix::ParallelLDLt()");
	cout << endl;
	cout.flush();
	return 0;
	
}

template<class TVar>
void * TPZSkylParMatrix<TVar>::ParallelLDLt2(void *t) {
	
	TPZSkylParMatrix<TVar> *loc = (TPZSkylParMatrix<TVar> *) t;
	int64_t col = 0, prevcol;
	prevcol=0;
	PZ_PTHREAD_MUTEX_LOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt2()");
	int64_t neq = loc->Dim();
	while (loc->fNumDecomposed < neq) {
		loc->ColumnToWork(col);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fNumDecomposed > loc->fNthreads){
				PZ_PTHREAD_COND_WAIT(&condition, &skymutex,"TPZSkylParMatrix::ParallelLDLt2()");
				col = 0;
			}else
			{
				//loc->fNthreads--;
				break;
			}
		}
		else {
			
			//Registers the working column number
			loc->fColUsed.insert(col);
			PZ_PTHREAD_COND_SIGNAL(&condition,"TPZSkylParMatrix::ParallelLDLt2()");
			PZ_PTHREAD_MUTEX_UNLOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt2()");
			loc->DecomposeColumnLDLt2(col);//loc->DecomposeColumn2(col, prevcol);
			
			PZ_PTHREAD_MUTEX_LOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt2()");
			if(loc->fDec[col] == col)
			{
				loc->fNumDecomposed++;
			}
			loc->fColUsed.erase(col);
		}
	}
	loc->fNthreads--;
	PZ_PTHREAD_COND_BROADCAST(&condition,"TPZSkylParMatrix::ParallelLDLt2()");
#ifdef VC
	cout << "Terminating thread " << pthread_self().x << " numthreads = " << loc->fNthreads << " numdecomposed " << loc->fNumDecomposed <<  endl;
#else
	cout << "Terminating thread " << pthread_self() << " numthreads = " << loc->fNthreads << " numdecomposed " << loc->fNumDecomposed <<  endl;
#endif
	PZ_PTHREAD_MUTEX_UNLOCK(&skymutex,"TPZSkylParMatrix::ParallelLDLt2()");
	cout << endl;
	cout.flush();
	return 0;
	
}

template<class TVar>
int TPZSkylParMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
	bool sing = this->fCorrectSingular;
	this->fCorrectSingular = true;
	int ans = Decompose_Cholesky();
	singular = this->fSingular;
	this->fCorrectSingular = sing;
	return ans;
}

template<class TVar>
int TPZSkylParMatrix<TVar>::Decompose_Cholesky()
{
	
	if (this->fDecomposed == ECholesky) return 1;
	
	int64_t neq = this->Dim();
	
	cout << endl << "TPZSkylParMatrix::Decompose_Cholesky() -> Number of Equations = " << neq << endl;
#ifdef BLAS
	cout << "Using BLAS " << endl;
#endif
	cout.flush();
	
	pthread_t *allthreads = new pthread_t[fNthreads-1];
	int *res = new int[fNthreads];
	
	fThreadUsed = new int[fNthreads];
	int i;
	
	
	//for(i=0; i<neq;i++) fDec[i]=-1;
	//SetSkyline(&skyline);
	
	for(i=0; i<fNthreads;i++) fThreadUsed[i]=-1;
	
	fEqDec=-1;
	
	for(i=0;i<fNthreads-1;i++) {
	  res[i] = PZ_PTHREAD_CREATE(&allthreads[i], NULL, 
				     this->ParallelCholesky, this, __FUNCTION__);
	}

	ParallelCholesky(this);
	for(i=0;i<fNthreads-1;i++) {
	  PZ_PTHREAD_JOIN(allthreads[i], NULL, __FUNCTION__);
	}
	
	delete []allthreads;// fThreadUsed, fDec;
#ifdef LOG4CXX
	if(fSingular.size())
	{
		std::stringstream sout;
		sout << "Singular equations ";
		std::list<int64_t>::iterator it;
		for (it=fSingular.begin(); it!=fSingular.end(); it++) {
			sout << *it << " ";
		}
		LOGPZ_WARN(logger,sout.str())
	}
#endif
	for(i=0;i<this->Dim();i++) {
		fDec[i]=fSkyline[i]-1;
	}
	this->fDecomposed = ECholesky;
	
	fEqDec = -1;
	
	this->fDecomposed = ECholesky;
	
	return 1;
}

template<class TVar>
int TPZSkylParMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
	bool sing = this->fCorrectSingular;
	this->fCorrectSingular = true;
	int ans = Decompose_LDLt();
	singular = this->fSingular;
	this->fCorrectSingular = sing;
	return ans;	
}

template<class TVar>
int TPZSkylParMatrix<TVar>::Decompose_LDLt()
{
	//	return TPZSkylMatrix::Decompose_LDLt();
	if (this->fDecomposed == ELDLt) return 1;
	
	int64_t neq = this->Dim();
	
	cout << endl << "TPZSkylParMatrix::Decompose_LDLt() -> Number of Equations = " << neq << endl;
#ifdef BLAS
	cout << "Using BLAS " << endl;
#endif
	cout.flush();
	int nthreads = fNthreads;
	//	fNthreads=1;
	
	pthread_t *allthreads = new pthread_t[fNthreads-1];
	int *res = new int[fNthreads];
	
	fThreadUsed = new int[fNthreads];
	int i;
	
	
	//for(i=0; i<neq;i++) fDec[i]=-1;
	//SetSkyline(&skyline);
	
	for(i=0; i<fNthreads;i++) fThreadUsed[i]=-1;
	
	fNumDecomposed = 0;
	for(i=0;i<this->Dim();i++) {
		fDec[i]=fSkyline[i]-1;
		if(fDec[i] == i-1) 
		{
			fNumDecomposed++;
			fDec[i] = i;
		}
	}
	fEqDec=-1;
	
	for(i=0;i<nthreads-1;i++) {
	  res[i] = PZ_PTHREAD_CREATE(&allthreads[i], NULL, this->ParallelLDLt2, 
				     this, __FUNCTION__);
	}
	ParallelLDLt2(this);
	for(i=0;i<nthreads-1;i++) {
	  PZ_PTHREAD_JOIN(allthreads[i], NULL, __FUNCTION__);
	}
	
	fNthreads = nthreads;
	delete []allthreads;// fThreadUsed, fDec;
	delete []fThreadUsed;
	fThreadUsed = 0;
	
	for(i=0;i<this->Dim();i++) {
		fDec[i]=fSkyline[i];
	}
	this->fDecomposed = ELDLt;
	
	fEqDec = -1;
	
	this->fDecomposed = ELDLt;
	
	return 1;
}

template<>
void TPZSkylParMatrix<std::complex<double> >::DecomposeColumnCholesky(int64_t col, int64_t prevcol)
{
    DebugStop();
}

template<class TVar>
void TPZSkylParMatrix<TVar>::DecomposeColumnCholesky(int64_t col, int64_t prevcol){
	
	//Cholesky Decomposition
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int64_t skprev, skcol; //prev and col Skyline height respectively
	int64_t minline;
	
	skprev = this->SkyHeight(prevcol);
	skcol = this->SkyHeight(col);
	
	ptrprev = this->Diag(prevcol);
	ptrcol = this->Diag(col);
	
	if((prevcol-skprev) > (col-skcol)){
		minline = prevcol - skprev;
	}else
    {
		minline = col - skcol;
    }
	if(minline > prevcol) {
		cout << "error condition col " << col << " prevcol " << prevcol << "\n";
		cout.flush();
		return;
	}
	TVar *run1 = ptrprev + (prevcol-minline);
	TVar *run2 = ptrcol + (col-minline);
	TVar sum = 0;
	
	//  int64_t n1=1;
	//  int64_t n2=1;  
#ifndef BLAS
#ifdef USETEMPLATE
	while(run1-ptrprev > templatedepth) {
		run1-=templatedepth;
		run2-=templatedepth;
		sum += TemplateSum<templatedepth>(run1--,run2--);
	}
#endif
	
	
	while(run1 != ptrprev) sum += (*run1--)*(*run2--);
#else
	//Blas dot product implementation
	//cout << "Calling BLAS " << endl;
	int64_t n = run1 - ptrprev;
	if (n>0){
		run1-=n-1;
		run2-=n-1;
		sum = ddot_(&n, run1, &n1, run2, &n2);
		
		/*for(int i=0;i<n;i++) {
		 cout << "run1 "<< run1[i] << " run2 " << run2[i] << "\n";
		 }
		 cout << " sum " << sum << " n " << n <<"\n";
		 run1--;
		 run2--;*/
	}
	//cout << "run1-ptrprev " << run1 - ptrprev << "\n";
	
#endif
	*run2-=sum;
	if(run1 != run2){
		*run2 /= *run1;
	}else{
		if (this->fCorrectSingular && IsZero(*run2)) {
			this->fSingular.push_back(col);
			*run2 = 1.;
		} else if (IsZero(*run2)) {
			std::cout << __PRETTY_FUNCTION__ << " Singular pivot at column " << col;
			DebugStop();
		}
		
		*run2=sqrt(*run2);
	}
	
}

template<class TVar>
void TPZSkylParMatrix<TVar>::DecomposeColumnLDLt(int64_t col, int64_t prevcol){
	
	//Cholesky Decomposition
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int64_t skprev, skcol; //prev and col Skyline height respectively
	int64_t minline;
	
	skprev = this->SkyHeight(prevcol);
	skcol = this->SkyHeight(col);
	
	ptrprev = this->Diag(prevcol);
	ptrcol = this->Diag(col);
	
	if((prevcol-skprev) > (col-skcol)){
		minline = prevcol - skprev;
	}else
    {
		minline = col - skcol;
    }
	if(minline > prevcol) {
		cout << "error condition col " << col << " prevcol " << prevcol << "\n";
		cout.flush();
		return;
	}
	int64_t k = minline;
	TVar *run1 = ptrprev + (prevcol-minline);
	TVar *run2 = ptrcol + (col-minline);
	TVar sum = 0;
	
	//  int64_t n1=1;
	//  int64_t n2=1;  
	
	
	while(run1 != ptrprev) sum += (*run1--)*(*run2--)* *this->fElem[k++];
	*run2-=sum;
	if(run1 != run2){
		*run2 /= *run1;
	}else{
#ifdef LOG4CXX
        if (logger->isDebugEnabled())
        {
            std::stringstream sout;
            sout << "col = " << col << " diagonal " << *run2;
            LOGPZ_DEBUG(logger,sout.str())
        }
#endif
		//*run2=sqrt(*run2);
	}
	
}


template<class TVar>
void TPZSkylParMatrix<TVar>::DecomposeColumnLDLt2(int64_t col){
	
	//Cholesky Decomposition
	TVar *ptrprev;     //Pointer to prev column
	TVar *ptrcol;      //Pointer to col column
	int64_t skprev, skcol; //prev and col Skyline height respectively
	int64_t minline;
	int64_t prevcol;
	while((fDec[col] < col && fDec[fDec[col]+1] == fDec[col]+1) || fDec[col] == col-1)
	{
		prevcol = fDec[col]+1;
		skprev = this->SkyHeight(prevcol);
		skcol = this->SkyHeight(col);
		
		ptrprev = this->Diag(prevcol);
		ptrcol = this->Diag(col);
		
		if((prevcol-skprev) > (col-skcol)){
			minline = prevcol - skprev;
		}else
		{
			minline = col - skcol;
		}
		if(minline > prevcol) {
			cout << "error condition col " << col << " prevcol " << prevcol << "\n";
			cout.flush();
			return;
		}
		int64_t k = minline;
		TVar *run1 = ptrprev + (prevcol-minline);
		TVar *run2 = ptrcol + (col-minline);
		TVar sum = 0;
		
		//  int64_t n1=1;
		//  int64_t n2=1;  
		
		
		while(run1 != ptrprev) sum += (*run1--)*(*run2--)* *this->fElem[k++];
		*run2-=sum;
		if(run1 != run2){
			*run2 /= *run1;
		}else{
			if (this->fCorrectSingular &&  fabs(*run2) < 1.e-10 ) {
				this->fSingular.push_back(col);
				*run2 = 1.;
			}
			
#ifdef LOG4CXX
            if (logger->isDebugEnabled())
            {
                std::stringstream sout;
                sout << "col = " << col << " diagonal " << *run2;
                LOGPZ_DEBUG(logger,sout.str())
            }
#endif
			//*run2=sqrt(*run2);
		}
		fDec[col]++;
	}
	
}

template<class TVar>
int TPZSkylParMatrix<TVar>::main_nada()  
{
	
	cout << "Calling constructor"<< endl;
	int64_t tempo1, tempo2, tempo; 
	int64_t dim, nequation, limite, kk, steps;
	kk=0;
	cout << "Entre o numero de equacoes\n";
	cin >> nequation;
	cout << "Entre com o limite\n";
	cin >> limite;
	cout << "Entre com o passo\n";
	cin >> steps;
	char  filename[20];
	int tempo_par[50], tempo_sin[50], passos[50];
	cout << "Entre o nome do Arquivo\n";
	cin >> filename;
	ofstream output(filename,ios::app);
	
	/*  
	 int i;
	 double soma_1=0;
	 double soma_2=0;
	 REAL run1[200];
	 REAL run2[200];
	 
	 const int tempdep = 10;
	 
	 
	 for(i=0;i<tempdep;i++) {
	 run1[i]=i;
	 run2[i]=i;
	 }
	 for(i=0;i<tempdep;i++){
	 soma_1 =soma_1+ run1[i]*run2[i];
	 
	 }
	 soma_2 = soma_2 + TemplateSum<tempdep>(run1,run2);
	 
	 
	 
	 cout << "soma 1 " << soma_1 << endl;
	 cout << "soma 2 " << soma_2 << endl;
	 */
	
	for(dim=nequation;dim<=limite;dim=dim+steps)
    {
		TPZVec<int64_t> skyline(dim);       
		int64_t i;
		for(i=0;i<dim;i++) skyline[i] = 0;//random()%(i+1); //0
		//for(i=0;i<dim;i++) cout << skyline[i] << ' ';
		cout << endl;
		TPZSkylParMatrix<TVar> MySky(dim,skyline,2); 
		int64_t j;
		for(i=0;i<dim;i++) {
			if(i%100==0) {cout << i << ' '; cout.flush();}
			for(j=skyline[i];j<=i;j++) {
				int random = rand();
				double rnd = (random*10.)/RAND_MAX;
				MySky(i,j)=rnd;
				if(i==j) MySky(i,i)=6000.;
			}
		}
		//MySky.SetSkyline(skyline); 
		TPZSkylMatrix<TVar> misael(MySky);
		// MySky.Print("Gustavo");
		//misael.Print("Misael");
		
		//  cout << "Constructed" << endl;
		//cout << "Calling Decompose_Cholesky routine..."<<endl;
		time_t t;
		int nthreads;  
		
		// cout << "Number of threads to work ? "; 
		
		//      cin >> nthreads;
		nthreads=2;
		//cout << endl;
		if(nthreads < 1) {
			cout << "VOCE E BURRO MESMO!!!\n"; 
			nthreads = 1;
		}
		tempo = time(&t);
		MySky.Decompose_LDLt();
		tempo1=time(&t)-tempo;
		MySky.Print("Gustavo",output);
		
		tempo = time(&t);
		misael.Decompose_LDLt();
		tempo2=time(&t)-tempo;
		
		
		misael.Print("misael",output); 
		passos[kk]=dim;
		tempo_par[kk]=tempo1;
		tempo_sin[kk]=tempo2;
		kk++;
    }
	output << "Multi-threaded  ";
	int i;
	for(i=0;i<kk;i++) output << "{"<< passos[i] << "," << tempo_par[i] << "},";
	output << endl;
	output << "Single-Threaded ";
	for(i=0;i<kk;i++) output << "{"<< passos[i] << "," << tempo_sin[i] << "},";
	output << endl;
	output.flush();
	
	return 1;
}

template<class TVar>
int TPZSkylParMatrix<TVar>::ClassId() const{
    return Hash("TPZSkylParMatrix") ^ TPZSkylMatrix<TVar>::ClassId() << 1;
}

template class TPZSkylParMatrix<float>;
template class TPZSkylParMatrix<double>;
template class TPZSkylParMatrix<long double>;

template class TPZSkylParMatrix<std::complex<float> >;
template class TPZSkylParMatrix<std::complex<double> >;
template class TPZSkylParMatrix<std::complex<long double> >;
