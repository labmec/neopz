/**
 * @file
 * @brief Contains the implementation of the TPZFrontMatrix methods.
 */
//$Id: TPZFrontMatrix.cpp,v 1.15 2011-05-11 02:10:40 phil Exp $
#include "TPZFrontMatrix.h"
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

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.frontal.frontmatrix"));
static LoggerPtr loggerfw(Logger::getLogger("pz.frontal.frontmatrix.fw"));
#endif

using namespace std;

template<class store, class front>
int TPZFrontMatrix<store, front>::Work(){
	return fFront.Work();
}
template<class store, class front>
void TPZFrontMatrix<store, front>::EquationsToDecompose(TPZVec<int> &destinationindex, int &lower_eq, int &upper_eq)
{
	int i;
	int loop_limit, global;
	loop_limit = destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		global = destinationindex[i];
		fNumElConnected[global]--;
	}
	upper_eq=fLastDecomposed;
	lower_eq=fLastDecomposed+1;
	while(upper_eq < fNumEq-1 && fNumElConnected[upper_eq+1]==0) upper_eq++;
	fLastDecomposed=upper_eq;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Constructor frontmatrix lower_eq "<< lower_eq << " upper_eq " << upper_eq << " fNumElConnected " << fNumElConnected;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/** Initializes the number of elements connected to each equation */
template<class store, class front>
void TPZFrontMatrix<store, front>::SetNumElConnected(TPZVec < int > &numelconnected){
	fNumElConnected.Resize(numelconnected.NElements());
	fNumElConnected=numelconnected;
	fNumElConnectedBackup = fNumElConnected;
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "fNumElConnected " << fNumElConnected;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	//	cout << "Storage Schema -> " << fStorage.GetStorage() << endl; 
	//	cout << "Front Matrix Type -> " << fFront.GetMatrixType() << endl;
#ifdef BLAS
	//     	cout << "Using BLAS" << endl;
#endif
#ifdef USING_ATLAS
	//          cout << "Using ATLAS" << endl;     
#endif
#ifndef USING_BLAS
#ifndef USING_ATLAS
	//     	cout << "Not Using BLAS" << endl;
#endif
#endif
}

/** Add a contribution of a stiffness matrix */
template<class store, class front>
void TPZFrontMatrix<store, front>::AddKel(TPZFMatrix<REAL> & elmat, TPZVec < int > & destinationindex){
	
	// message #1.3 to fFront:TPZFront
	fFront.AddKel(elmat, destinationindex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< fFront.NElements();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
    /*      cout << "destination index" << endl;
	 int i;
	 for(i=0;i<destinationindex.NElements();i++) cout << destinationindex[i] << " ";
	 cout << endl;
	 cout.flush();
	 elmat.Print("Element Matrix");
	 */
	int mineq, maxeq;
	
	EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray AuxEqn;
	if(maxeq >= mineq) {
		
		fFront.DecomposeEquations(mineq,maxeq,AuxEqn);
		CheckCompress();
		fStorage.AddEqnArray(&AuxEqn);
		if(maxeq == this->Rows()-1){
			fStorage.FinishWriting();
			fStorage.ReOpen();
		}
	}
	this->fDecomposed = fFront.GetDecomposeType();
}

/** Add a contribution of a stiffness matrix */
template<class store, class front>
void TPZFrontMatrix<store, front>::AddKel(TPZFMatrix<REAL> & elmat, TPZVec < int > & sourceindex, TPZVec < int > & destinationindex)
{
	fFront.AddKel(elmat, sourceindex, destinationindex);
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
	
	//	EquationsToDecompose(destinationindex);
	//         cout << "AddKel::destination index 2" << endl;
	//     int i;
	//          for(i=0;i<destinationindex.NElements();i++) cout << destinationindex[i] << " ";
	//         cout << endl;
	//         cout.flush();
	//          elmat.Print("AddKel: Element Matrix 2");
	int mineq, maxeq;
	EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray AuxEqn;
	if(maxeq >= mineq) {
		fFront.DecomposeEquations(mineq,maxeq,AuxEqn);
		CheckCompress();
		fStorage.AddEqnArray(&AuxEqn);
		if(maxeq == this->Rows()-1){
			fFront.Reset(0);
			fStorage.FinishWriting();
			fStorage.ReOpen();
		}
		
	}
	this->fDecomposed = fFront.GetDecomposeType();
}


/** Add a contribution of a stiffness matrix using the indexes to compute the frontwidth */
template<class store, class front>
void TPZFrontMatrix<store, front>::SymbolicAddKel(TPZVec < int > & destinationindex)
{
	fFront.SymbolicAddKel(destinationindex);
	int mineq, maxeq;
	EquationsToDecompose(destinationindex, mineq, maxeq);
	
	if(maxeq >= mineq) {
		fFront.SymbolicDecomposeEquations(mineq,maxeq);
		CheckCompress();
	}
	
}


template<class store, class front>
TPZFrontMatrix<store, front>::~TPZFrontMatrix(){
}


template<class store, class front>
TPZFrontMatrix<store, front>::TPZFrontMatrix() : TPZAbstractFrontMatrix()
{
	fFront.Reset();
	fStorage.Reset();
	fNumElConnected.Resize(0);
	fNumElConnectedBackup.Resize(0);
	fLastDecomposed = -1;
	fNumEq=0;
}

template<class store, class front>
TPZFrontMatrix<store, front>::TPZFrontMatrix(int globalsize) : TPZAbstractFrontMatrix(globalsize,globalsize)
{
	fFront.Reset(globalsize);
	fStorage.Reset();
	fNumElConnected.Resize(0);
	fNumElConnectedBackup.Resize(0);
	fLastDecomposed = -1;
	fNumEq=globalsize;
}

template<class store, class front>
void TPZFrontMatrix<store, front>::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	int i;
	out << "Frontal Matrix associated"<< endl;
	fFront.Print(name, out);
	out << "Stored Equations" << endl;
	fStorage.Print(name, out);
	out << "Number of Equations " << fNumEq << endl;
	out << "Number of Elements connected to DF" << endl;
	for(i=0;i<fNumElConnected.NElements();i++){
		out << i << " " << fNumElConnected[i] << endl;
	}
}

template<class store, class front>
void TPZFrontMatrix<store, front>::main()
{
	TPZFMatrix<REAL> KEl1(2,2);
	KEl1(0,0)=4.;
	KEl1(0,1)=6.;
	KEl1(1,0)=6.;
	KEl1(1,1)=12.;
	TPZVec<int> DestInd1(2);
	DestInd1[0]=0;
	DestInd1[1]=1;
	
	
	TPZFMatrix<REAL> KEl2(4,4);
	KEl2(0,0)=4.;
	KEl2(0,1)=-6.;
	KEl2(0,2)=2.;
	KEl2(0,3)=6.;
	
	KEl2(1,1)=12;
	KEl2(1,2)=-6.;
	KEl2(1,3)=-12.;
	
	KEl2(2,2)=4.;
	KEl2(2,3)=6.;
	
	KEl2(3,3)=12.;
	
	int i, j;
	for(i=0;i<4;i++) {
		for(j=i;j<4;j++) KEl2(j,i)=KEl2(i,j);
	}
	
	TPZVec<int> DestInd2(4);
	DestInd2[0]=0;
	DestInd2[1]=1;
	DestInd2[2]=2;
	DestInd2[3]=3;
	
	
	
	TPZFMatrix<REAL> KEl3(2,2);
	KEl3(0,0)=3.;
	KEl3(0,1)=3.;
	KEl1(1,0)=3.;
	KEl3(1,1)=6.;
	
	TPZVec<int> DestInd3(2);
	DestInd3[0]=2;
	DestInd3[1]=3;
	
	
	
	TPZFrontMatrix TestFront(4);
	
	TPZVec<int> NumConnected(4);
	for(i=0;i<4;i++) NumConnected[i]=2;
	
	
	TestFront.SetNumElConnected(NumConnected);
	
	TestFront.SymbolicAddKel(DestInd1);
	TestFront.SymbolicAddKel(DestInd2);
	TestFront.SymbolicAddKel(DestInd3);
	
	TestFront.AllocData();
	TestFront.SetNumElConnected(NumConnected);
	
	TestFront.AddKel(KEl1,DestInd1);
	TestFront.AddKel(KEl2,DestInd2);
	TestFront.AddKel(KEl3,DestInd3);
}

template<class store, class front>
void TPZFrontMatrix<store, front>::CheckCompress()
{
	double nfreerate = ( (double)fFront.NFree() / (double)fFront.FrontSize() ) * 100;
	if(nfreerate>20.) 
	{
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Compressing front nfreerate " << nfreerate << " NFree " << fFront.NFree() << " Front elements " << fFront.FrontSize();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fFront.Compress();
#ifdef LOG4CXX
		{
			std::stringstream sout;
			sout << "Frondwidth after Compress "<< fFront.FrontSize();
			LOGPZ_INFO(loggerfw,sout.str())
		}
#endif
	}
}

template<class store, class front>
void TPZFrontMatrix<store, front>::AllocData()
{
	fFront.AllocData();
	fLastDecomposed=-1;
}

template<class store, class front>
int TPZFrontMatrix<store, front>::Substitution(TPZFMatrix<REAL> *b) const {
	fStorage.Forward(*b, ELU);
	fStorage.Backward(*b, ELU);
	return 1;
}

template<class store, class front>
int TPZFrontMatrix<store, front>::Subst_Forward(TPZFMatrix<REAL> *b) const {
	cout << "Entering Forward Substitution\n";
	cout.flush();
	DecomposeType dec = fFront.GetDecomposeType();
	if(dec != ECholesky) cout << "TPZFrontMatrix::Subst_Forward non matching decomposition\n";
	fStorage.Forward(*b, dec);
	return 1;
}

template<class store, class front>
int TPZFrontMatrix<store, front>::Subst_Backward(TPZFMatrix<REAL> *b) const {
	cout << "Entering Backward Substitution\n";
	cout.flush();
	DecomposeType dec = fFront.GetDecomposeType();
	if(dec != ECholesky) cout << "TPZFrontMatrix::Subst_Forward non matching decomposition\n";
	fStorage.Backward(*b, dec);
	return 1;
}
/*
 template<class store, class front>
 void TPZFrontMatrix<store, front>::SetFileName(const char *name) {
 //     char * bin_template[8]; 
 //     bin_template = "bfXXXXXX"
 //     const char *name = mktemp(template);
 //	fStorage.('w', name);
 }*/

template<class store, class front>
void TPZFrontMatrix<store, front>::SetFileName(char option, const char *name) {
	fStorage.OpenGeneric(option, name);
}


template<class store, class front>
void TPZFrontMatrix<store, front>::FinishWriting() {
	fStorage.FinishWriting();
}
template<class store, class front>
void TPZFrontMatrix<store, front>::ReOpen() {
	fStorage.ReOpen();
}

template<class store, class front>
int TPZFrontMatrix<store, front>::Zero() {
	fStorage.Zero();
	fNumElConnected = fNumElConnectedBackup;
	fLastDecomposed = -1;
	fFront.Reset(fNumEq);
	return 0;
}

class TPZStackEqnStorage;
class TPZFileEqnStorage;
class TPZFrontSym;
class TPZFrontNonSym;

template class TPZFrontMatrix<TPZStackEqnStorage, TPZFrontSym>;
template class TPZFrontMatrix<TPZFileEqnStorage, TPZFrontSym>;
template class TPZFrontMatrix<TPZStackEqnStorage, TPZFrontNonSym>;
template class TPZFrontMatrix<TPZFileEqnStorage, TPZFrontNonSym>;
