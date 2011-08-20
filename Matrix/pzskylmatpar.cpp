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
/*extern "C" {
 #include <cblas.h>
 };*/
#endif

#include "pzfmatrix.h"
#include "pzskylmat.h"
#include "pzskylmatpar.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.matrix.tpzskylparmatrix"));
#endif

/** @brief Semaphore */
pthread_mutex_t skymutex = PTHREAD_MUTEX_INITIALIZER;
/** @brief Condition to waiting */
pthread_cond_t condition = PTHREAD_COND_INITIALIZER;

/** @brief Initializing number of iterations to Template SUM */
const int templatedepth = 10;

//Constructors

/*TPZSkylParMatrix::TPZSkylParMatrix (const int dim,int NumThreads)  : TPZSkylMatrix(dim, dim),fDec(dim)
 {
 fNthreads = NumThreads;
 }
 */
TPZSkylParMatrix::TPZSkylParMatrix(const int dim, const TPZVec<int> &skyline,int NumThreads)

: TPZSkylMatrix(dim, skyline),fDec(dim), fCorrectSingular(false)
{
	int i;
	
	fSkyline = skyline;
	
	for(i=0;i<dim;i++) {
		fDec[i]=skyline[i]-1;
		
	}
	fNthreads = NumThreads;
	fEqDec = -1;
}

TPZSkylParMatrix::TPZSkylParMatrix(const TPZSkylParMatrix &copy) : TPZSkylMatrix(copy), fDec(copy.fDec),fSkyline(copy.fSkyline),fEqDec(copy.fEqDec),
fNthreads(copy.fNthreads),fThreadUsed(0),fCorrectSingular(false)
{
}

TPZSkylParMatrix::TPZSkylParMatrix()
: TPZSkylMatrix(),fDec(0),fCorrectSingular(false)
{
	//something
}

// TPZSkylParMatrix::TPZSkylParMatrix(const TPZSkylParMatrix &cp):TPZSkylMatrix(cp){
//   fDec = cp.fDec;
//   fSkyline = cp.fSkyline;
//   fNthreads = cp.fNthreads;
//   fEqDec = cp.fEqDec;
// }

TPZSkylParMatrix::~TPZSkylParMatrix() {
}

void TPZSkylParMatrix::SetSkyline(const TPZVec<int> &skyline)
{
	fSkyline = skyline;
	int i,neq=Dim();
	for(i=0;i<neq;i++) fDec[i]=skyline[i]-1;
	fEqDec=-1;
	fCorrectSingular = false;
	TPZSkylMatrix::SetSkyline(skyline);
}

void TPZSkylParMatrix::PrintState() {
	int i;
	cout << "Decomposed equations " << fEqDec << ' ' << " Thread 'i' working on equation" << ' ';
	for(i=0; i<fNthreads; i++) cout << fThreadUsed[i] << ' ';
	cout << endl;
	cout << "Columns to work " << ' ';
	int neq = Dim();
	for(i=0; i<neq; i++) cout << fDec[i] << ' ';
	cout << endl;
	cout << "Singular equations ";
	std::list<int>::iterator it;
	for (it=fSingular.begin(); it!=fSingular.end(); it++) {
		cout << *it << ' ';
	}
	cout.flush();
}


void TPZSkylParMatrix::ColumnToWork(int &lcol, int &lprevcol) {
	int neq = Dim();
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

void TPZSkylParMatrix::ColumnToWork(int &lcol) {
	int neq = Dim();
	int lcolentry = lcol;
#ifdef DEBUG
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
			std::stringstream sout;
			sout << "Will work column " << lcol;
			LOGPZ_DEBUG(logger,sout.str())
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



void * TPZSkylParMatrix::ParallelCholesky(void *t) {
	
	TPZSkylParMatrix *loc = (TPZSkylParMatrix *) t;
	int aux_col;//, k;
	int col, prevcol;
	prevcol=0;
	pthread_mutex_lock(&skymutex);
	int neq = loc->Dim();
	while (loc->fEqDec < neq-1) {
		loc->ColumnToWork(col, prevcol);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fEqDec > loc->fNthreads){
				pthread_cond_wait(&condition, &skymutex);
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
			pthread_cond_signal(&condition);
			pthread_mutex_unlock(&skymutex);
			if(loc->fCorrectSingular)
			{
				loc->DecomposeColumn(col, prevcol,loc->fSingular);//loc->DecomposeColumn2(col, prevcol);
			}
			else {
				loc->DecomposeColumn(col, prevcol);
			}
			
			pthread_mutex_lock(&skymutex);
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
	pthread_mutex_unlock(&skymutex);
	pthread_cond_broadcast(&condition);
	cout << endl;
	cout.flush();
	return 0;
	
}

void * TPZSkylParMatrix::ParallelLDLt(void *t) {
	
	TPZSkylParMatrix *loc = (TPZSkylParMatrix *) t;
	int aux_col;//, k;
	int col, prevcol;
	prevcol=0;
	pthread_mutex_lock(&skymutex);
	int neq = loc->Dim();
	while (loc->fEqDec < neq-1) {
		loc->ColumnToWork(col, prevcol);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fEqDec > loc->fNthreads){
				pthread_cond_wait(&condition, &skymutex);
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
			pthread_cond_signal(&condition);
			pthread_mutex_unlock(&skymutex);
			loc->DecomposeColumnLDLt(col, prevcol);//loc->DecomposeColumn2(col, prevcol);
			pthread_mutex_lock(&skymutex);
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
	pthread_mutex_unlock(&skymutex);
	pthread_cond_broadcast(&condition);
	cout << endl;
	cout.flush();
	return 0;
	
}

void * TPZSkylParMatrix::ParallelLDLt2(void *t) {
	
	TPZSkylParMatrix *loc = (TPZSkylParMatrix *) t;
	int col = 0, prevcol;
	prevcol=0;
	pthread_mutex_lock(&skymutex);
	int neq = loc->Dim();
	while (loc->fNumDecomposed < neq) {
		loc->ColumnToWork(col);
		if(col==-1) {
			cout.flush();
			if(neq-loc->fNumDecomposed > loc->fNthreads){
				pthread_cond_wait(&condition, &skymutex);
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
			pthread_cond_signal(&condition);
			pthread_mutex_unlock(&skymutex);
			loc->DecomposeColumnLDLt2(col);//loc->DecomposeColumn2(col, prevcol);
			
			pthread_mutex_lock(&skymutex);
			if(loc->fDec[col] == col)
			{
				loc->fNumDecomposed++;
			}
			loc->fColUsed.erase(col);
		}
	}
	loc->fNthreads--;
	pthread_cond_broadcast(&condition);
	cout << "Terminating thread " << pthread_self() << " numthreads = " << loc->fNthreads << " numdecomposed " << loc->fNumDecomposed <<  endl;
	pthread_mutex_unlock(&skymutex);
	cout << endl;
	cout.flush();
	return 0;
	
}

int TPZSkylParMatrix::Decompose_Cholesky(std::list<int> &singular)
{
	bool sing = this->fCorrectSingular;
	this->fCorrectSingular = true;
	int ans = Decompose_Cholesky();
	singular = this->fSingular;
	this->fCorrectSingular = sing;
	return ans;
}

int TPZSkylParMatrix::Decompose_Cholesky()
{
	
	if (fDecomposed == ECholesky) return 1;
	
	int neq = Dim();
	
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
	
	for(i=0;i<fNthreads-1;i++) res[i] = pthread_create(&allthreads[i],NULL,this->ParallelCholesky, this);
	ParallelCholesky(this);
	for(i=0;i<fNthreads-1;i++) pthread_join(allthreads[i], NULL);
	
	delete []allthreads;// fThreadUsed, fDec;
#ifdef LOG4CXX
	if(fSingular.size())
	{
		std::stringstream sout;
		sout << "Singular equations ";
		std::list<int>::iterator it;
		for (it=fSingular.begin(); it!=fSingular.end(); it++) {
			sout << *it << " ";
		}
		LOGPZ_WARN(logger,sout.str())
	}
#endif
	for(i=0;i<Dim();i++) {
		fDec[i]=fSkyline[i]-1;
	}
	fDecomposed = ECholesky;
	
	fEqDec = -1;
	
	fDecomposed = ECholesky;
	
	return 1;
}

int TPZSkylParMatrix::Decompose_LDLt(std::list<int> &singular)
{
	bool sing = this->fCorrectSingular;
	this->fCorrectSingular = true;
	int ans = Decompose_LDLt();
	singular = this->fSingular;
	this->fCorrectSingular = sing;
	return ans;	
}

int TPZSkylParMatrix::Decompose_LDLt()
{
	//	return TPZSkylMatrix::Decompose_LDLt();
	if (fDecomposed == ELDLt) return 1;
	
	int neq = Dim();
	
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
	for(i=0;i<Dim();i++) {
		fDec[i]=fSkyline[i]-1;
		if(fDec[i] == i-1) 
		{
			fNumDecomposed++;
			fDec[i] = i;
		}
	}
	fEqDec=-1;
	
	for(i=0;i<nthreads-1;i++) res[i] = pthread_create(&allthreads[i],NULL,this->ParallelLDLt2, this);
	ParallelLDLt2(this);
	for(i=0;i<nthreads-1;i++) pthread_join(allthreads[i], NULL);
	
	fNthreads = nthreads;
	delete []allthreads;// fThreadUsed, fDec;
	delete []fThreadUsed;
	fThreadUsed = 0;
	
	for(i=0;i<Dim();i++) {
		fDec[i]=fSkyline[i];
	}
	fDecomposed = ELDLt;
	
	fEqDec = -1;
	
	fDecomposed = ELDLt;
	
	return 1;
}

void TPZSkylParMatrix::DecomposeColumnCholesky(int col, int prevcol){
	
	//Cholesky Decomposition
	REAL *ptrprev;     //Pointer to prev column
	REAL *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	
	skprev = SkyHeight(prevcol);
	skcol = SkyHeight(col);
	
	ptrprev = Diag(prevcol);
	ptrcol = Diag(col);
	
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
	REAL *run1 = ptrprev + (prevcol-minline);
	REAL *run2 = ptrcol + (col-minline);
	REAL sum = 0;
	
	//  long n1=1;
	//  long n2=1;  
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
	long n = run1 - ptrprev;
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
		if (this->fCorrectSingular && *run2 < 1.e-10 ) {
			this->fSingular.push_back(col);
			*run2 = 1.;
		} else if (*run2 < 1.e-25) {
			std::cout << __PRETTY_FUNCTION__ << " Singular pivot at column " << col;
			DebugStop();
		}
		
		*run2=sqrt(*run2);
	}
	
}

void TPZSkylParMatrix::DecomposeColumnLDLt(int col, int prevcol){
	
	//Cholesky Decomposition
	REAL *ptrprev;     //Pointer to prev column
	REAL *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	
	skprev = SkyHeight(prevcol);
	skcol = SkyHeight(col);
	
	ptrprev = Diag(prevcol);
	ptrcol = Diag(col);
	
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
	int k = minline;
	REAL *run1 = ptrprev + (prevcol-minline);
	REAL *run2 = ptrcol + (col-minline);
	REAL sum = 0;
	
	//  long n1=1;
	//  long n2=1;  
	
	
	while(run1 != ptrprev) sum += (*run1--)*(*run2--)* *fElem[k++];
	*run2-=sum;
	if(run1 != run2){
		*run2 /= *run1;
	}else{
#ifdef LOG4CXX
		std::stringstream sout;
		sout << "col = " << col << " diagonal " << *run2;
		LOGPZ_DEBUG(logger,sout.str())
#endif
		//*run2=sqrt(*run2);
	}
	
}



void TPZSkylParMatrix::DecomposeColumnLDLt2(int col){
	
	//Cholesky Decomposition
	REAL *ptrprev;     //Pointer to prev column
	REAL *ptrcol;      //Pointer to col column
	int skprev, skcol; //prev and col Skyline height respectively
	int minline;
	int prevcol;
	while((fDec[col] < col && fDec[fDec[col]+1] == fDec[col]+1) || fDec[col] == col-1)
	{
		prevcol = fDec[col]+1;
		skprev = SkyHeight(prevcol);
		skcol = SkyHeight(col);
		
		ptrprev = Diag(prevcol);
		ptrcol = Diag(col);
		
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
		int k = minline;
		REAL *run1 = ptrprev + (prevcol-minline);
		REAL *run2 = ptrcol + (col-minline);
		REAL sum = 0;
		
		//  long n1=1;
		//  long n2=1;  
		
		
		while(run1 != ptrprev) sum += (*run1--)*(*run2--)* *fElem[k++];
		*run2-=sum;
		if(run1 != run2){
			*run2 /= *run1;
		}else{
			if (this->fCorrectSingular &&  fabs(*run2) < 1.e-10 ) {
				this->fSingular.push_back(col);
				*run2 = 1.;
			}
			
#ifdef LOG4CXX
			std::stringstream sout;
			sout << "col = " << col << " diagonal " << *run2;
			LOGPZ_DEBUG(logger,sout.str())
#endif
			//*run2=sqrt(*run2);
		}
		fDec[col]++;
	}
	
}


int TPZSkylParMatrix::main_nada()  
{
	
	cout << "Calling constructor"<< endl;
	long tempo1, tempo2, tempo; 
	int dim, nequation, limite, kk, steps;
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
		TPZVec<int> skyline(dim);       
		int i;
		for(i=0;i<dim;i++) skyline[i] = 0;//random()%(i+1); //0
		//for(i=0;i<dim;i++) cout << skyline[i] << ' ';
		cout << endl;
		TPZSkylParMatrix MySky(dim,skyline,2); 
		int j; 
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
		TPZSkylMatrix misael(MySky);
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
