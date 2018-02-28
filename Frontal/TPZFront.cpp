/**
 * @file
 * @brief Contains the implementation of the TPZFront methods.
 */

#include "pzsfulmat.h"
#include "TPZFront.h"
#include "pzstack.h"
#include "pzreal.h"
#include <math.h>

#include "tpzeqnarray.h"
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>

using namespace std;

template<class TVar>
TPZFront<TVar>::TPZFront(){
	fExpandRatio = 200;
	fFront = 0;
	fMaxFront=0;
	fWork = 0;
	fNextRigidBodyMode = 0;
	fProductMTData = NULL;
    fDecomposeType = ELU;
}

template<class TVar>
TPZFront<TVar>::TPZFront(int64_t GlobalSize)
{
	fExpandRatio = 200;
	fFront = 0;
	fMaxFront=0;
	fLocal.Resize(GlobalSize);
	int64_t i;
	for(i=0;i<GlobalSize;i++) fLocal[i]=-1;
	fWork = 0;
	fNextRigidBodyMode = GlobalSize;
	fProductMTData = NULL;
    fDecomposeType = ELU;
}

template<class TVar>
TPZFront<TVar>::TPZFront(const TPZFront<TVar> &cp) : fMaxFront(cp.fMaxFront),
fGlobal(cp.fGlobal), fLocal(cp.fLocal),fFront(cp.fFront),fFree(cp.fFree),
fData(cp.fData),
fDecomposeType(cp.fDecomposeType)
{
	fNextRigidBodyMode = cp.fNextRigidBodyMode;
	fExpandRatio = cp.fExpandRatio;
	if (cp.fProductMTData) {
		const int nthreads = cp.fProductMTData->NThreads();
		this->fProductMTData = new STensorProductMTData(nthreads,this);
	}
	else {
		this->fProductMTData = NULL;
	}
}

	
template<class TVar>
TPZFront<TVar>::~TPZFront(){
}

template<class TVar>
void TPZFront<TVar>::PrintGlobal(const char *name, std::ostream& out){
	out << name << endl;
	/*
	 int i, j;
	 for(i=0;i<fLocal.NElements();i++){
	 if(fLocal[i]!=-1) out << i << " ";
	 }
	 out << endl;
	 for(i=0;i<fLocal.NElements();i++){
	 if(fLocal[i]==-1) continue;
	 for(j=0;j<fLocal.NElements();j++){
	 if(fLocal[j]==-1) continue;
	 out << Element(fLocal[i],fLocal[j]) << " ";
	 }
	 out << endl;
	 }
	 out << endl;
	 */
	out << "Not implemented in the abstract Class !" << endl;
	out << "Try one of the subclasses" << endl;
}


template<class TVar>
void TPZFront<TVar>::Print(const char *name, std::ostream& out) const
{
	if(name) out << name << endl;
	/*int i,j,loop_limit;
	 
	 
	 out <<  "Frontal Matrix Size          "<< fFront << endl;
	 out <<  "Maximum Frontal Matrix Size  "<< fMaxFront << endl;
	 
	 out << endl;
	 out <<  "Local Indexation "<< endl;
	 out <<  "Position  "<<" Local index"<< endl;
	 for(i=0;i<fLocal.NElements();i++){
	 out <<   i <<  "         "  << fLocal[i] << endl;
	 }
	 
	 out << endl;
	 out <<  "Global Indexation "<< endl;
	 out <<  "Position  "<<" Global index"<< endl;
	 
	 for(i=0;i<fGlobal.NElements();i++){
	 out <<   i <<  "            "<< fGlobal[i] << endl;
	 }
	 
	 out << endl;
	 out <<  "Local Freed Equatoins  "<< fFree.NElements() << endl;
	 out <<  "position "<<   "Local Equation "<<endl;
	 loop_limit=fFree.NElements();
	 for(i=0;i<loop_limit;i++){
	 out <<  i <<  "             " << fFree[i] << endl;
	 }
	 out <<  "Frontal Matrix "<< endl;
	 if(fData.NElements() > 0) {
	 for(i=0;i<fFront;i++){
	 for(j=0;j<fFront;j++) out << ((i<j) ? Element(i,j) : Element(j,i)) <<  " ";
	 out << endl;
	 }
	 }*/
	
	out << "Not implemented in the abstract Class !" << endl;
	out << "Try one of the subclasses" << endl;
}

template<class TVar>
void TPZFront<TVar>::FreeGlobal(int64_t global)
{
	if(fLocal[global]==-1){
		cout << "TPZFront FreeGlobal was called with wrong parameters !" << endl;
		return;
	}
	int64_t index;
	index=fLocal[global];
	fGlobal[index]=-1;
	fLocal[global]=-1;
	fFree.Push(index);
}

template<class TVar>
int TPZFront<TVar>::Local(int64_t global){
	/*	int index;
	 if(fLocal[global]!=-1) return fLocal[global];
	 if(fFree.NElements()){
	 index=fFree.Pop();
	 }else{
	 index=fFront++;
	 }
	 fLocal[global]=index;
	 // At this point we extend the frontmatrix
	 if(index >= fMaxFront) {
	 cout << "Dynamically expanding the front size to " << index+10 << endl;
	 Expand(index+10);
	 }
	 if(fGlobal.NElements()<=index) fGlobal.Resize(index+1);
	 fGlobal[index]=global;
	 return index;
	 */
	cout << "Not implemented in the abstract Class !" << endl;
	cout << "Try one of the subclasses" << endl;
	return -1;
}

/** Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
template<class TVar>
void TPZFront<TVar>::SymbolicAddKel(TPZVec < int64_t > & destinationindex)
{
	int64_t i, loop_limit, aux;
	loop_limit=destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		aux=destinationindex[i];
		Local(aux);
		fFront = fGlobal.NElements();
	}
	fMaxFront=(fFront<fMaxFront)?fMaxFront:fFront;
	
}

template<class TVar>
void TPZFront<TVar>::SymbolicDecomposeEquations(int64_t mineq, int64_t maxeq)
{
	int64_t i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}



/** Implements tests for TPZFront */
/*
template<class TVar>
void TPZFront<TVar>::main()
{
	int i, j;
	
	 //Populates data structure
	int matsize=6;
	TPZFMatrix<TVar> TestMatrix(matsize,matsize);
	for(i=0;i<matsize;i++) {
		for(j=i;j<matsize;j++) {
			int random = rand();
			float rnd = (float)(random*matsize)/RAND_MAX;
			TestMatrix(i,j)= (TVar)rnd;
			TestMatrix(j,i)=TestMatrix(i,j);
			if(i==j) TestMatrix(i,j)=6000.;
		}
	}
	
	TPZFMatrix<TVar> Prova;
	Prova=TestMatrix;
	
	//	Prova.Decompose_Cholesky();
	Prova.Print("TPZFMatrix<REAL> Cholesky");
	
	TPZFront TestFront(matsize);
	
	
	TPZVec<int64_t> DestIndex(matsize);
	for(i=0;i<matsize;i++) DestIndex[i]=i;
	
	TestFront.SymbolicAddKel(DestIndex);
	TestFront.SymbolicDecomposeEquations(0,matsize-1);
	
	std::string OutFile;
	OutFile = "TPZFrontTest.txt";
	
	ofstream output(OutFile.c_str(),ios::app);
	
	//	TestFront.Compress();
	
	//	TestFront.AllocData();
	
	//	TestFront.AddKel(TestMatrix, DestIndex);
	TPZEqnArray<TVar> Result;
	
//	TestFront.DecomposeEquations(0,0,Result);
//	 
//	 TestFront.Print(OutFile, output);
//	 
//	 ofstream outeqn("TestEQNArray.txt",ios::app);
//	 Result.Print("TestEQNArray.txt",outeqn);
//	 
//	 TestFront.Compress();
//	 
//	 TestFront.Print(OutFile, output);
//	 
	//	TestFront.DecomposeEquations(0,matsize-1,Result);
	ofstream outeqn("TestEQNArray.txt",ios::app);
	
	Result.Print("TestEQNArray.txt",outeqn);
	
	
	TPZFMatrix<TVar> Load(matsize);
	
	for(i=0;i<matsize;i++) {
		int random = rand();
		float rnd = (float)(random*matsize)/RAND_MAX;
		Load(i,0)= (TVar)rnd;
	}
	
	TPZFMatrix<TVar> Load_2(matsize);
	Load_2=Load;
	
	//	Prova.Subst_Forward(&Load);
	//	Prova.Subst_Backward(&Load);
	
	
	DecomposeType decType = ECholesky;
	Prova.SolveDirect(Load, decType);
	
	Load.Print("Load");
	//TestFront.Print(OutFile, output);
	
	Result.EqnForward(Load_2, decType);
	Result.EqnBackward(Load_2, decType);
	
	Load_2.Print("Eqn");
	
}
*/

template<class TVar>
void TPZFront<TVar>::Reset(int64_t GlobalSize)
{
	fData.Resize(0);
	fFree.Resize(0);
	fFront=0;
	fGlobal.Resize(0);
	fLocal.Resize(GlobalSize);
	fLocal.Fill(-1);
	fMaxFront=0;
	fExpandRatio = 200;
	fWork = 0;
	fNextRigidBodyMode = GlobalSize;
}

template<class TVar>
int64_t TPZFront<TVar>::NElements(){
	return fLocal.NElements();
}


template<class TVar>
int64_t TPZFront<TVar>::NFree()
{
	int64_t i;
	int64_t free_eq=0;
	for(i=0;i<fGlobal.NElements();i++)
	{
		if(fGlobal[i]==-1){
			free_eq=free_eq+1;
		}
	}
	return free_eq;
}

template<class TVar>
void TPZFront<TVar>::TensorProductIJ(int ithread,typename TPZFront<TVar>::STensorProductMTData *data)
{
	PZError << "This Method Should only be called in the lower classes" << std::endl;
	DebugStop();
}

template class TPZFront<float>;
template class TPZFront<std::complex<float> >;

template class TPZFront<double>;
template class TPZFront<std::complex<double> >;

template class TPZFront<long double>;
template class TPZFront<std::complex<long double> >;