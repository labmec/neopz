/**
 * @file
 * @brief Contains the implementation of the TPZFrontSym methods.
 */

#include "TPZFrontSym.h"
#include <math.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <fstream>
#include "tpzeqnarray.h"

using namespace std;

#ifdef USING_BLAS
#include "cblas.h"
#endif


#ifdef USING_ATLAS
void cblas_dspr(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                const int N, const double alpha, const double *X,
                const int incX, double *Ap);
#endif
DecomposeType TPZFrontSym::GetDecomposeType() const
{
	return fDecomposeType;
}
void TPZFrontSym::PrintGlobal(const char *name, std::ostream& out){
	int i, j;
	out << name << endl;
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
}
void TPZFrontSym::Print(const char *name, std::ostream& out) const
{
	if(name) out << name << endl;
	int i,j,loop_limit;
	
	
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
	}
}
void TPZFrontSym::AllocData()
{
	fData.Resize(fMaxFront*(fMaxFront+1)/2);
	fGlobal.Fill(-1);
	fData.Fill(0.);
	fFront=0;
	//fLocal.Fill(-1);
}
void TPZFrontSym::Reset(int GlobalSize)
{
	fData.Resize(0);
	fFree.Resize(0);
	fFront=0;
	fGlobal.Resize(0);
	fLocal.Resize(GlobalSize);
	fLocal.Fill(-1);
	fMaxFront=0;
}
int TPZFrontSym::NFree()
{
	return fFree.NElements();
}

int TPZFrontSym::Local(int global){
	int index;
	if(fLocal[global]!=-1) return fLocal[global];
	if(fFree.NElements()){
		index=fFree.Pop();
	}else{
		index=fFront++;
	}
	fLocal[global]=index;
	// At this point we extend the frontmatrix
	if(index >= fMaxFront) {
		//	     cout << endl;
		//		cout << "Dynamically expanding the front size to " << index+fExpandRatio << endl;
		Expand(index+fExpandRatio);
	}
	if(fGlobal.NElements()<=index) fGlobal.Resize(index+1);
	fGlobal[index]=global;
	return index;
}
void TPZFrontSym::FreeGlobal(int global)
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
void TPZFrontSym::DecomposeOneEquation(int ieq, TPZEqnArray &eqnarray)
{
	//eqnarray.SetSymmetric();
    
	int i, ilocal;
	ilocal = Local(ieq);
	TPZVec<REAL> AuxVec(fFront);
	
	for(i=0;i<ilocal;i++) AuxVec[i]=Element(i,ilocal);
	for(i=ilocal;i<fFront;i++) AuxVec[i]=Element(ilocal,i);
	
	if(AuxVec[ilocal]<0) {
		cout << "TPZFront::DecomposeOneEquation AuxVec[ilocal] < 0 " << AuxVec[ilocal] << " ilocal=" << ilocal << " fGlobal=" << fGlobal[ilocal] << endl;
	}
	REAL diag = sqrt(AuxVec[ilocal]);
	for(i=0;i<fFront;i++) AuxVec[i]/=diag;
    
	eqnarray.BeginEquation(ieq);       
	eqnarray.AddTerm(ieq, diag);
	
	//Blas utilizatioin  
#ifdef USING_BLAS     
	CBLAS_ORDER order = CblasColMajor;
	CBLAS_UPLO up_lo = CblasUpper;
	int sz = fFront;
	long incx = 1;
	double db = -1.;//AuxVec[ilocal];
	cblas_dspr(order, up_lo,sz,db,&AuxVec[0],incx,&Element(0,0));
	
#endif
#ifdef USING_ATLAS
	CBLAS_ORDER order = CblasColMajor;
	CBLAS_UPLO up_lo = CblasUpper;
	int sz = fFront;
	long incx = 1;
	double db = -1.;//AuxVec[ilocal];
	cblas_dspr(order, up_lo,sz,db,&AuxVec[0],incx,&Element(0,0));
	//cout << "Using ATLAS" << endl;
#endif
#ifndef USING_BLAS
#ifndef USING_ATLAS
	int j=0;
	for(i=0;i<fFront;i++){
   		for(j=i;j<fFront;j++){
   			Element(i,j)=Element(i,j)-AuxVec[i]*AuxVec[j];
   		}
   	}
#endif
#endif
	
	
    for(i=0;i<fFront;i++) {
		if(i!=ilocal && fGlobal[i]!= -1 && AuxVec[i] != 0.) eqnarray.AddTerm(fGlobal[i],AuxVec[i]);
	}
	eqnarray.EndEquation();
	
	for(i=0;i<ilocal;i++)Element(i,ilocal)=0.;
	for(i=ilocal;i<fFront;i++) Element(ilocal,i)=0.;
	
    FreeGlobal(ieq);
    fDecomposeType=ECholesky;
	//	PrintGlobal("After", output);
}
void TPZFrontSym::AddKel(TPZFMatrix &elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex)
{
	int i, j, ilocal, jlocal, nel;
	nel=sourceindex.NElements();
	for (i = 0; i < nel; i++) {
		// message #1.1.1 to this:TPZFront
		ilocal = this->Local(destinationindex[i]);
		for (j = i; j < nel; j++) {
			// message #1.1.2.1 to this:TPZFront
			jlocal = this->Local(destinationindex[j]);
			
			// message #1.1.2.2 to this:TPZFront
			this->Element(ilocal, jlocal)+=elmat(sourceindex[i],sourceindex[j]);
		}
	}
}
void TPZFrontSym::AddKel(TPZFMatrix &elmat, TPZVec<int> &destinationindex)
{
	int i, j, ilocal, jlocal, nel;
	nel = destinationindex.NElements();
	for(i=0;i<nel;i++){
		ilocal = this->Local(destinationindex[i]);
		for(j=i;j<nel;j++) {
			jlocal=this->Local(destinationindex[j]);
			this->Element(ilocal, jlocal)+=elmat(i,j);
		}
	}
	/*
	 output << "Dest Index " ;
	 for(i=0;i<nel;i++) output << destinationindex[i] << " ";
	 output << endl;
	 elmat.Print("Element matrix ",output);
	 PrintGlobal("After Assemb..." , output);
	 */
}
/*REAL & TPZFrontSym::Element(int i, int j){
 if(i>j){
 int i_temp=i;
 i=j;
 j=i_temp;
 //cout << "Changing row column indexes !"<<endl;
 }
 return fData[(j*(j+1))/2+i];
 }
 */
void TPZFrontSym::Expand(int larger) {
	fData.Resize(larger*(larger+1)/2,0.);
	fMaxFront = larger;
}
void TPZFrontSym::Compress(){
	//	PrintGlobal("Before COmpress",output);
	//	Print("Before Compress", cout);
	TPZStack <int> from;
	int nfound;
	int i, j;
	for(i = 0; i < fFront; i++){
		if(fGlobal[i] != -1) from.Push(i);
	}
	
	/**
	 *First fLocal initialization
	 *Any needed updates is done on next loop
	 */
	nfound = from.NElements();
	for(i=0;i<nfound;i++) {
		fGlobal[i]=fGlobal[from[i]];
		//fGlobal[from[i]] = -1;
		fLocal[fGlobal[i]] = i;
	}
	for(;i<fGlobal.NElements();i++) fGlobal[i] = -1;
	
	if(nfound+fFree.NElements()!=fFront) cout << "TPZFront.Compress inconsistent data structure\n";
	int frontold = fFront;
	fFront = nfound;
	fFree.Resize(0);	
	fGlobal.Resize(fFront);
	if(fData.NElements()==0) return;
	
	for(j = 0; j < nfound; j++){
		for(i = 0; i <= j; i++){
			Element(i,j) = Element(from[i], from[j]);
		}
		//		fGlobal[i] = fGlobal[from[i]];
		//		fLocal[fGlobal[i]] = i;
	}
	for(;j<frontold;j++) {
		for(i=0;i<= j; i++) Element(i,j) = 0.;
	}
	
	//	Print("After Compress", cout);
	//	PrintGlobal("After Compress",output);
}
void TPZFrontSym::SymbolicAddKel(TPZVec < int > & destinationindex)
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
void TPZFrontSym::SymbolicDecomposeEquations(int mineq, int maxeq)
{
	int i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}
void TPZFrontSym::DecomposeEquations(int mineq, int maxeq, TPZEqnArray & eqnarray){
	// message #1.1 to eqnarray:TPZEqnArray
	int ieq;
	eqnarray.Reset();
	eqnarray.SetSymmetric();
	//cout << "Decomposing from " << mineq << " to " << maxeq << "\n";
	//cout.flush();
	for (ieq = mineq; ieq <= maxeq; ieq++) {
		// message #1.2.1 to this:TPZFront
		this->DecomposeOneEquation(ieq, eqnarray);
		//			this->Print("Teste.txt",output);
	}
}
TPZFrontSym::TPZFrontSym(int GlobalSize) : TPZFront(GlobalSize)
{
	fDecomposeType=ECholesky;
}
TPZFrontSym::TPZFrontSym(){
	fDecomposeType=ECholesky;
}


TPZFrontSym::~TPZFrontSym(){}
void TPZFrontSym::main()
{
	int i, j;
	/**
	 * 	Populates data structure
	 */
	int matsize=6;
	TPZFMatrix TestMatrix(matsize,matsize);
	for(i=0;i<matsize;i++) {
		for(j=i;j<matsize;j++) {
			int random = rand();
			double rnd = (random*matsize)/0x7fff;
			TestMatrix(i,j)=rnd;
			TestMatrix(j,i)=TestMatrix(i,j);
			if(i==j) TestMatrix(i,j)=6000.;
		}
	}
	
	TPZFMatrix Prova;
	Prova=TestMatrix;
	
	//	Prova.Decompose_Cholesky();
	Prova.Print("TPZFMatrix Cholesky");
	
	TPZFrontSym TestFront(matsize);
	
	
	TPZVec<int> DestIndex(matsize);
	for(i=0;i<matsize;i++) DestIndex[i]=i;
	
	TestFront.SymbolicAddKel(DestIndex);
	TestFront.SymbolicDecomposeEquations(0,matsize-1);
	
	std::string OutFile;
	OutFile = "TPZFrontSymTest.txt";
	
	ofstream output(OutFile.c_str(),ios::app);
	
	TestFront.Compress();
	
	TestFront.AllocData();
	
	TestFront.AddKel(TestMatrix, DestIndex);
	TPZEqnArray Result;
	
	/*TestFront.DecomposeEquations(0,0,Result);
	 
	 TestFront.Print(OutFile, output);
	 
	 ofstream outeqn("TestEQNArray.txt",ios::app);
	 Result.Print("TestEQNArray.txt",outeqn);
	 
	 TestFront.Compress();
	 
	 TestFront.Print(OutFile, output);
	 */
	TestFront.DecomposeEquations(0,matsize-1,Result);
	ofstream outeqn("TestEQNArray.txt",ios::app);
	
	Result.Print("TestEQNArray.txt",outeqn);
	
	
	TPZFMatrix Load(matsize);
	
	for(i=0;i<matsize;i++) {
		int random = rand();
		double rnd = (random*matsize)/0x7fff;
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

std::string TPZFrontSym::GetMatrixType(){
	return "Symmetric matrix";
}

void TPZFrontSym::ExtractFrontMatrix(TPZFMatrix &front)
{
	int maxeq = fLocal.NElements();
	int mineq = 0;
	for(mineq=0; mineq<maxeq; mineq++) if(fLocal[mineq] != -1) break;
	int numeq = maxeq-mineq;
	front.Redim(numeq,numeq);
	int ieq,jeq;
	for(ieq=mineq;ieq<maxeq;ieq++) {
		if(fLocal[ieq] == -1) continue;
		int il = ieq-mineq;
		for(jeq=ieq;jeq<maxeq;jeq++) {
			if(fLocal[jeq] == -1) continue;
			int jl = jeq-mineq;
			front(il,jl) = this->Element(fLocal[ieq],fLocal[jeq]);
			front(jl,il) = front(il,jl);
		}
	}
	
}
