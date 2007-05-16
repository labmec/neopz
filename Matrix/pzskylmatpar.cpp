

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



#ifndef PZNTPAR
pthread_mutex_t skymutex = PTHREAD_MUTEX_INITIALIZER;
pthread_cond_t condition = PTHREAD_COND_INITIALIZER;
#endif

#ifdef PZNTPAR
omni_mutex omni_mutexo;
omni_condition omni_condo(&omni_mutexo);
#endif

const int templatedepth = 10;

//Constructors

/*TPZSkylParMatrix::TPZSkylParMatrix (const int dim,int NumThreads)  : TPZSkylMatrix(dim, dim),fDec(dim)
{
	fNthreads = NumThreads;
}
*/
TPZSkylParMatrix::TPZSkylParMatrix(const int dim, const TPZVec<int> &skyline,int NumThreads)

  : TPZSkylMatrix(dim, skyline),fDec(dim)
{
  int i;
  
  fSkyline = skyline;

  for(i=0;i<dim;i++) {
    fDec[i]=skyline[i]-1;

  }
  fNthreads = NumThreads;
  fEqDec = -1;
}

TPZSkylParMatrix::TPZSkylParMatrix()
  : TPZSkylMatrix(),fDec(0)
{
  //something
}

TPZSkylParMatrix::~TPZSkylParMatrix() {
}

void TPZSkylParMatrix::SetSkyline(const TPZVec<int> &skyline)
{
  int i,neq=Dim();
  for(i=0;i<neq;i++) fDec[i]=skyline[i]-1;
  fEqDec=-1;
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


#ifndef PZNTPAR

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
	loc->DecomposeColumn(col, prevcol);//loc->DecomposeColumn2(col, prevcol);
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


int TPZSkylParMatrix::Decompose_Cholesky(std::list<int> &singular)
{
  return Decompose_Cholesky();
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

  for(i=0;i<Dim();i++) {
	fDec[i]=fSkyline[i]-1;
  }
  fDecomposed = ECholesky;

  fEqDec = -1;

  fDecomposed = ECholesky;

  return 1;
}
#endif 

#ifdef OMNITHREAD 
int TPZSkylParMatrix::Decompose_Cholesky(std::list<int> &singular)
{ 
  return Decompose_Cholesky();
}

int TPZSkylParMatrix::Decompose_Cholesky()
{ 

  if (fDecomposed == ECholesky) return 1;

  int neq = Dim();

  cout << "TPZSkylParMatrix::Decompose_Cholesky() -> Number of Equations = " << neq << endl;

		cout.flush();
  
  omni_thread **allthreads = new omni_thread*[fNthreads];
  
  fThreadUsed = new int[fNthreads];
  int i;
  
  
  //for(i=0; i<neq;i++) fDec[i]=-1;
  for(i=0; i<fNthreads;i++) fThreadUsed[i]=-1;
  
  fEqDec=-1; 
 
  
  for(i=0;i<fNthreads;i++) allthreads[i] = omni_thread::create(this->OmniParallelCholesky, this);

  void *return_value;
  for(i=0;i<fNthreads;i++) allthreads[i]->join(&return_value);


  delete[] allthreads;
  delete[] fThreadUsed;
  for(i=0;i<Dim();i++) {
	fDec[i]=fSkyline[i]-1;
  }
  fDecomposed = ECholesky;
  fEqDec = -1;
 
  return 1;
}



void  *TPZSkylParMatrix::OmniParallelCholesky(void *t) {
  //cout << "Omni Parallel..." << omni_thread::self() << "\n";
  //cout.flush();
  TPZSkylParMatrix *loc = (TPZSkylParMatrix *) t;

  int aux_col;

  int col, prevcol;
  prevcol=0;
  omni_mutexo.acquire();
  int neq = loc->Dim();
  while (loc->fEqDec < neq-1) {
    loc->ColumnToWork(col, prevcol);
    if(col==-1) {
      cout.flush();
      if(neq-loc->fEqDec > loc->fNthreads){

		//cout << "Going to cond_wait " << omni_thread::self() << "\n";
		//cout.flush();

		omni_condo.wait();
      }else 
		{

			//loc->fNthreads--;  

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
	omni_condo.signal();
	omni_mutexo.release();
	//loc->DecomposeColumn(col, prevcol);

	loc->TPZSkylMatrix::DecomposeColumn(col, prevcol);

	omni_mutexo.acquire();
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

  omni_mutexo.release();
  omni_condo.broadcast();

  cout << endl;
  cout.flush();
  return NULL;
}
#endif
 
void TPZSkylParMatrix::DecomposeColumn(int col, int prevcol){
 
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
    *run2=sqrt(*run2);
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
      MySky.Decompose_Cholesky();
      tempo1=time(&t)-tempo;
      //MySky.Print("Gustavo");
      
      tempo = time(&t);
      misael.Decompose_Cholesky();
      tempo2=time(&t)-tempo;
 

      //misael.Print("misael"); 
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
