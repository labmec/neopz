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

TPZFront::TPZFront(){
	fExpandRatio = 200;
	fFront = 0;
	fMaxFront=0;
	fWork = 0;
	fNextRigidBodyMode = 0;
}

TPZFront::TPZFront(int GlobalSize)
{
	fExpandRatio = 200;
	fFront = 0;
	fMaxFront=0;
	fLocal.Resize(GlobalSize);
	int i;
	for(i=0;i<GlobalSize;i++) fLocal[i]=-1;
	fWork = 0;
	fNextRigidBodyMode = GlobalSize;
}

TPZFront::TPZFront(const TPZFront &cp) : fMaxFront(cp.fMaxFront),
fGlobal(cp.fGlobal), fLocal(cp.fLocal),fFront(cp.fFront),fFree(cp.fFree),
fData(cp.fData)
{
	fNextRigidBodyMode = cp.fNextRigidBodyMode;
	fExpandRatio = cp.fExpandRatio;
}

TPZFront::~TPZFront(){
}

void TPZFront::PrintGlobal(const char *name, std::ostream& out){
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


void TPZFront::Print(const char *name, std::ostream& out) const
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

void TPZFront::FreeGlobal(int global)
{
	if(fLocal[global]==-1){
		cout << "TPZFront FreeGlobal was called with wrong parameters !" << endl;
		return;
	}
	int index;
	index=fLocal[global];
	fGlobal[index]=-1;
	fLocal[global]=-1;
	fFree.Push(index);
}
int TPZFront::Local(int global){
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
void TPZFront::SymbolicAddKel(TPZVec < int > & destinationindex)
{
	int i, loop_limit, aux;
	loop_limit=destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		aux=destinationindex[i];
		Local(aux);
		fFront = fGlobal.NElements();
	}
	fMaxFront=(fFront<fMaxFront)?fMaxFront:fFront;
	
}
void TPZFront::SymbolicDecomposeEquations(int mineq, int maxeq)
{
	int i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}



/** Implements tests for TPZFront */
void TPZFront::main()
{
	int i, j;
	/**
	 Populates data structure
	 */
	int matsize=6;
	TPZFMatrix TestMatrix(matsize,matsize);
	for(i=0;i<matsize;i++) {
		for(j=i;j<matsize;j++) {
			int random = rand();
			double rnd = (random*matsize)/RAND_MAX;
			TestMatrix(i,j)=rnd;
			TestMatrix(j,i)=TestMatrix(i,j);
			if(i==j) TestMatrix(i,j)=6000.;
		}
	}
	
	TPZFMatrix Prova;
	Prova=TestMatrix;
	
	//	Prova.Decompose_Cholesky();
	Prova.Print("TPZFMatrix Cholesky");
	
	TPZFront TestFront(matsize);
	
	
	TPZVec<int> DestIndex(matsize);
	for(i=0;i<matsize;i++) DestIndex[i]=i;
	
	TestFront.SymbolicAddKel(DestIndex);
	TestFront.SymbolicDecomposeEquations(0,matsize-1);
	
	std::string OutFile;
	OutFile = "TPZFrontTest.txt";
	
	ofstream output(OutFile.c_str(),ios::app);
	
	//	TestFront.Compress();
	
	//	TestFront.AllocData();
	
	//	TestFront.AddKel(TestMatrix, DestIndex);
	TPZEqnArray Result;
	
	/*TestFront.DecomposeEquations(0,0,Result);
	 
	 TestFront.Print(OutFile, output);
	 
	 ofstream outeqn("TestEQNArray.txt",ios::app);
	 Result.Print("TestEQNArray.txt",outeqn);
	 
	 TestFront.Compress();
	 
	 TestFront.Print(OutFile, output);
	 */
	//	TestFront.DecomposeEquations(0,matsize-1,Result);
	ofstream outeqn("TestEQNArray.txt",ios::app);
	
	Result.Print("TestEQNArray.txt",outeqn);
	
	
	TPZFMatrix Load(matsize);
	
	for(i=0;i<matsize;i++) {
		int random = rand();
		double rnd = (random*matsize)/RAND_MAX;
		Load(i,0)=rnd;
	}
	
	TPZFMatrix Load_2(matsize);
	Load_2=Load;
	
	//	Prova.Subst_Forward(&Load);
	//	Prova.Subst_Backward(&Load);
	
	
	DecomposeType decType = ECholesky;
	Prova.SolveDirect(Load, decType);
	
	Load.Print();
	//TestFront.Print(OutFile, output);
	
	Result.EqnForward(Load_2, decType);
	Result.EqnBackward(Load_2, decType);
	
	Load_2.Print("Eqn");
	
	
	
}
void TPZFront::Reset(int GlobalSize)
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

int TPZFront::NElements(){
	return fLocal.NElements();
}
int TPZFront::NFree()
{
	int i;
	int free_eq=0;
	for(i=0;i<fGlobal.NElements();i++)
	{
		if(fGlobal[i]==-1){
			free_eq=free_eq+1;
		}
	}
	return free_eq;
}
