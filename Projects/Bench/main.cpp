#include <vector>
#include <pthread.h>
#include "pzvec.h"
#include "pzfmatrix.h"
#ifdef USING_ATLAS
extern "C"{
     #include <cblas.h>
     };
#endif
using namespace std;
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <time.h>
#include <stdio.h>

#include <CoreServices/CoreServices.h>
#include <vecLib/vDSP.h>
#include <sys/param.h>
#include <sys/sysctl.h>


//#include "main.h"
//#include "JavaMode.h"
#include <StdIO.h>


#ifdef USING_ATLAS
     void cblas_dspr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
     void cblas_dgemv(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TA,
                 const int M, const int N, const double alpha, const double *A,
                 const int lda, const double *X, const int incX,
                 const double beta, double *Y, const int incY);
	void cblas_zaxpy(const int N, const void *alpha, const void *X,
                 const int incX, void *Y, const int incY);

#endif
#ifdef APPLE

#endif

TPZVec<double> vg(9,0.);
TPZVec<double> vr(9,0.);

static void * MultiplySubmat(void * data);
struct SubMatData{
	TPZFMatrix * mat;
	TPZVec<double> vector;
	TPZVec<double> * vr;
	int starting_row;
	pthread_mutex_t * mutex;
};
void * MultiplySubmat(void * data){
	SubMatData * submat = (SubMatData *) data;
	TPZVec<double> vetorb(submat->mat->Rows(),0.);
	int i;
/*	ostringstream filename;
	filename << "file";
	filename << pthread_self();
	string file;
	file = filename.str();
	ofstream out(file.c_str());
	out << "Multiplying on thread\n";
	submat->mat->Print("Matrix 3x3\n", out);
	out << submat->vector.NElements() << endl;;
	
	out << submat->vector[0] << endl;
	out << submat->vector[1] << endl;
	out << submat->vector[2] << endl;
	out.flush();
*/
	#ifdef USING_ATLAS
	CBLAS_ORDER order = CblasColMajor;
	CBLAS_UPLO up_lo = CblasUpper;
	CBLAS_TRANSPOSE TA = CblasNoTrans;
	int sz = submat->mat->Rows();//tamanho;
	long incx = 1;
	double db = -1.;//AuxVec[ilocal];
	cblas_dgemv(order, TA, sz, sz, 1, &submat->mat->operator()(0,0),
			sz, &submat->vector[0], incx, 0, &vetorb[0], incx);
	//cout << "Using ATLAS" << endl;
	#endif

	pthread_mutex_lock(submat->mutex);

	for(i=0;i<submat->vector.NElements();i++){
		//cout << vetorb[i] << " " << endl;
		//cout << submat->starting_row+i << endl;
		(*submat->vr)[submat->starting_row+i]+=vetorb[i];
	}
	cout.flush();
	//contribuir no global
	pthread_mutex_unlock(submat->mutex);
};
        
int TestZAXPY();
int TestZAXPY(){
	int size  = 10;
	TPZFMatrix mat(size,1,1.);
	TPZFMatrix mat_res(size,1,1.);
	double * m = new double[size];
	double * m2 = new double[size];
	int i;
	for(i=0;i<size;i++){
		m[i]=1.;
		m2[i]=1.;
	}
	double alpha = 2.;
	mat.Print("Matriz 1", cout);
	mat_res.Print("Matriz 2", cout);
#ifdef USING_BLAS
        cblas_daxpy(size, alpha, &mat(0,0), 1, &mat_res(0,0), 1);
#endif
	mat.Print("Matriz 1 depois", cout);
	mat_res.Print("Matriz 2 depois", cout);
	cout << "Antes\n";
	for(i=0;i<size;i++){
		cout << m[i] << " ";
		cout << m2[i] << endl;
	}
#ifdef USING_BLAS
	cblas_daxpy(size, alpha, &m[0], 1, &m2[0], 1);
#endif
	cout << "Depois\n";
	for(i=0;i<size;i++){
		cout << m[i] << " ";
		cout << m2[i] << endl;
	}
	delete [] m;
	delete [] m2;
}
#
int main(){
//	TestZAXPY();

 	pthread_mutex_t VectorMutex;
	pthread_mutex_init(&VectorMutex, NULL);
	int i, j, size;
	cin >> size; //Must be a multiple of ndivisions

	TPZFMatrix mat(size,size,0.);
	TPZVec<double> vetor(size);
	TPZVec<double> vetorResult(size, 0.);


	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			mat(i,j)=1+(double) (3000.0*rand()/(RAND_MAX+1.0));
		}
		vetor[i]=1+(double) (30.0*rand()/(RAND_MAX+1.0));
	}
	int ndivisions = 3;
	int division_size = size / ndivisions;
	vector<TPZFMatrix> fMatrices;
	vector<TPZVec <double> > fVectors;
	for(i = 0; i<ndivisions; i++){
		for(j=0;j<ndivisions;j++){
			TPZFMatrix submat;
			mat.GetSub(i*ndivisions,j*ndivisions,division_size, division_size, submat);
			fMatrices.push_back(submat);
		}
	}
	//pegar sub vcs


	//Uncoment this for submatrices visualization
/*
		mat.Print("Matriz inteira", cout);
		vector<TPZFMatrix>::iterator it;
		it = fMatrices.begin();
		for(;it!=fMatrices.end();it++)
			it->Print("Pequenas matrizes", cout);

		for(i=0;i<size;i++)
			cout << vetor[i] << endl;

		cout << "Division size " << division_size << endl;
		*/
		int k;
		for(i = 0; i<ndivisions; i++){
			TPZVec<double> subvec(division_size,0.);
			for(k = 0; k < division_size; k++){
				subvec[k] = vetor[k+(i*division_size)];
				//cout << k << " " << k+(i*division_size) << endl;
			}
			fVectors.push_back(subvec);
		}
/*		double buf;
		double vecout[9];
		int n = 9;

		for (k = 0; k < 9; k++){
			for (i = 0; i < n; i++) {
				buf = 0.;
				for(j = 0; j < n; j++){
				buf += mat(i,j) * vetor[j];
				}
				vecout[i] = buf;
			//      cout << buf << endl;
			}//i
		}
*/
/*		for(i=0;i<size;i++)
			cout << vecout[i] << endl;
		//cin >> k;
*/
	pthread_t * mthreads = new pthread_t[fMatrices.size()];
	int counter=0;
	vector<TPZVec <double> >::iterator vec=fVectors.begin();
	int sr = 0;

	time_t sttime;
	time_t endtime;
	time (& sttime);
	cout << "Start Time "<< sttime << endl;
	cout.flush();
	for(i=0;i<ndivisions*ndivisions;i++){
			
		SubMatData * data = new SubMatData;
		
		data->mat = &fMatrices[i];
		
		data->vector = *vec;
		
		vec++;
		
		data->vr = &vetorResult;
		data->starting_row = sr;//counter * division_size;
		data->mutex = &VectorMutex;
		pthread_create(&mthreads[i], NULL, MultiplySubmat, data);
		counter++;
		if(counter==ndivisions) {
			sr+=division_size;
			counter = 0;
			vec=fVectors.begin();
		}
	}

	for(i=0;i<ndivisions*ndivisions;i++){
		pthread_join(mthreads[i], NULL);
	}

	cout << "Resultado final\n";
	pthread_mutex_lock(&VectorMutex);
	time (& endtime);
	cout << "End Time " << endtime << endl;
	cout.flush();
	int time_elapsed = endtime - sttime;
	cout << "Time Elapsed " << time_elapsed << endl;
	/*for(i=0;i<size;i++){
		cout << vetorResult[i] << endl;
	}*/
	pthread_mutex_unlock(&VectorMutex);
	

}


