/**
 * @file
 * @brief Contains the implementation of the TPZFrontMatrix methods.
 */

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

#ifdef PZ_LOG
static TPZLogger logger("pz.frontal.frontmatrix");
static TPZLogger loggerfw("pz.frontal.frontmatrixfww");
#endif

using namespace std;

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar, store, front>::Work(){
	return fFront.Work();
}
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::EquationsToDecompose(TPZVec<int64_t> &destinationindex, int64_t &lower_eq, int64_t &upper_eq)
{
	int64_t i;
	int64_t loop_limit, global;
	loop_limit = destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		global = destinationindex[i];
		fNumElConnected[global]--;
	}
	upper_eq=fLastDecomposed;
	lower_eq=fLastDecomposed+1;
	while(upper_eq < fNumEq-1 && fNumElConnected[upper_eq+1]==0) upper_eq++;
	fLastDecomposed=upper_eq;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Constructor frontmatrix lower_eq "<< lower_eq << " upper_eq " << upper_eq << " fNumElConnected " << fNumElConnected;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/** Initializes the number of elements connected to each equation */
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::SetNumElConnected(TPZVec < int > &numelconnected){
	fNumElConnected.Resize(numelconnected.NElements());
	fNumElConnected=numelconnected;
	fNumElConnectedBackup = fNumElConnected;
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << "fNumElConnected " << fNumElConnected;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
}

/** Add a contribution of a stiffness matrix */
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & destinationindex){
	
	// message #1.3 to fFront:TPZFront
	fFront.AddKel(elmat, destinationindex);
#ifdef PZ_LOG
    if (loggerfw.isInfoEnabled())
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< fFront.NElements();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
	int64_t mineq, maxeq;
	
	EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray<TVar> AuxEqn;
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
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec < int64_t > & sourceindex, TPZVec < int64_t > & destinationindex)
{
	fFront.AddKel(elmat, sourceindex, destinationindex);
#ifdef PZ_LOG
    if (loggerfw.isInfoEnabled())
	{
		std::stringstream sout;
		sout << "Frondwidth after AddKel "<< fFront.FrontSize();
		LOGPZ_INFO(loggerfw,sout.str())
	}
#endif
	
	int64_t mineq, maxeq;
	EquationsToDecompose(destinationindex, mineq, maxeq);
	TPZEqnArray<TVar> AuxEqn;
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
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::SymbolicAddKel(TPZVec < int64_t > & destinationindex)
{
	fFront.SymbolicAddKel(destinationindex);
	int64_t mineq, maxeq;
	EquationsToDecompose(destinationindex, mineq, maxeq);
	
	if(maxeq >= mineq) {
		fFront.SymbolicDecomposeEquations(mineq,maxeq);
		CheckCompress();
	}
	
}


template<class TVar, class store, class front>
TPZFrontMatrix<TVar,store, front>::~TPZFrontMatrix(){
}


template<class TVar, class store, class front>
TPZFrontMatrix<TVar,store, front>::TPZFrontMatrix() : TPZRegisterClassId(&TPZFrontMatrix::ClassId),
TPZAbstractFrontMatrix<TVar>()
{
	fFront.Reset();
	fStorage.Reset();
	fNumElConnected.Resize(0);
	fNumElConnectedBackup.Resize(0);
	fLastDecomposed = -1;
	fNumEq=0;
}

template<class TVar, class store, class front>
TPZFrontMatrix<TVar,store, front>::TPZFrontMatrix(int64_t globalsize) : TPZRegisterClassId(&TPZFrontMatrix::ClassId),
TPZAbstractFrontMatrix<TVar>(globalsize,globalsize)
{
	fFront.Reset(globalsize);
	fStorage.Reset();
	fNumElConnected.Resize(0);
	fNumElConnectedBackup.Resize(0);
	fLastDecomposed = -1;
	fNumEq=globalsize;
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::ClassId() const{
    return Hash("TPZFrontMatrix") ^ TPZAbstractFrontMatrix<TVar>::ClassId() << 1 ^ store().ClassId() << 2 ^ front().ClassId() << 3;
}

template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::Print(const char *name, std::ostream& out,const MatrixOutputFormat form) const
{
	int64_t i;
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


template<class TVar, class store, class front>
int64_t TPZFrontMatrix<TVar,store,front>::Size() const
{
	PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return -1;
}
template<class TVar, class store, class front>
TVar* &TPZFrontMatrix<TVar,store,front>::Elem()
{
	PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	static TVar* t{nullptr};
  return t;
}
template<class TVar, class store, class front>
const TVar* TPZFrontMatrix<TVar,store,front>::Elem() const
{
	PZError<<__PRETTY_FUNCTION__;
  PZError<<"ERROR:Should not be called\n.Aborting...\n";
  DebugStop();
	return nullptr;
}

template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::main()
{
	TPZFMatrix<TVar> KEl1(2,2);
	KEl1(0,0)=4.;
	KEl1(0,1)=6.;
	KEl1(1,0)=6.;
	KEl1(1,1)=12.;
	TPZVec<int64_t> DestInd1(2);
	DestInd1[0]=0;
	DestInd1[1]=1;
	
	
	TPZFMatrix<TVar> KEl2(4,4);
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
	
	TPZVec<int64_t> DestInd2(4);
	DestInd2[0]=0;
	DestInd2[1]=1;
	DestInd2[2]=2;
	DestInd2[3]=3;

	TPZFMatrix<TVar> KEl3(2,2);
	KEl3(0,0)=3.;
	KEl3(0,1)=3.;
	KEl1(1,0)=3.;
	KEl3(1,1)=6.;
	
	TPZVec<int64_t> DestInd3(2);
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

template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::CheckCompress()
{
	double nfreerate = ( (double)fFront.NFree() / (double)fFront.FrontSize() ) * 100;
	if(nfreerate>20.) 
	{
#ifdef PZ_LOG
        if (logger.isDebugEnabled())
		{
			std::stringstream sout;
			sout << "Compressing front nfreerate " << nfreerate << " NFree " << fFront.NFree() << " Front elements " << fFront.FrontSize();
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
		fFront.Compress();
#ifdef PZ_LOG
        if (loggerfw.isInfoEnabled())
		{
			std::stringstream sout;
			sout << "Frondwidth after Compress "<< fFront.FrontSize();
			LOGPZ_INFO(loggerfw,sout.str())
		}
#endif
	}
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::SolveDirect( TPZFMatrix<TVar> &B , DecomposeType dt, std::list<int64_t> &singular) {
    if (fFront.GetDecomposeType() != dt) {
        DebugStop();
    }
#ifdef PZ_LOG2
    if (logger.isDebugEnabled())
    {
        
        std::stringstream sout;
        B.Print("On input " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Subst_Forward(&B);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B.Print("After forward and diagonal " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    Subst_Backward(&B);
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        B.Print("Final result " , sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    return 1;
}


template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::AllocData()
{
	fFront.AllocData();
	fLastDecomposed=-1;
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::Substitution(TPZFMatrix<TVar> *b) const {
	fStorage.Forward(*b, ELU);
	fStorage.Backward(*b, ELU);
	return 1;
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::Subst_Forward(TPZFMatrix<TVar> *b) const {
	cout << "Entering Forward Substitution\n";
	cout.flush();
	DecomposeType dec = fFront.GetDecomposeType();
//	if(dec != ECholesky) cout << "TPZFrontMatrix::Subst_Forward non matching decomposition\n";
	fStorage.Forward(*b, dec);
	return 1;
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::Subst_Backward(TPZFMatrix<TVar> *b) const {
	cout << "Entering Backward Substitution\n";
	cout.flush();
	DecomposeType dec = fFront.GetDecomposeType();
//	if(dec != ECholesky) cout << "TPZFrontMatrix::Subst_Forward non matching decomposition\n";
	fStorage.Backward(*b, dec);
	return 1;
}

template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::SetFileName(char option, const char *name) {
	fStorage.OpenGeneric(option, name);
}

template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::FinishWriting() {
	fStorage.FinishWriting();
}
template<class TVar, class store, class front>
void TPZFrontMatrix<TVar,store, front>::ReOpen() {
	fStorage.ReOpen();
}

template<class TVar, class store, class front>
int TPZFrontMatrix<TVar,store, front>::Zero() {
	fStorage.Zero();
	fNumElConnected = fNumElConnectedBackup;
	fLastDecomposed = -1;
	fFront.Reset(fNumEq);
	return 0;
}

template<class TVar>
class TPZStackEqnStorage;
template<class TVar>
class TPZFileEqnStorage;
template<class TVar>
class TPZFrontSym;
template<class TVar>
class TPZFrontNonSym;


template class TPZFrontMatrix<float,TPZStackEqnStorage<float>, TPZFrontSym<float> >;
template class TPZFrontMatrix<float,TPZFileEqnStorage<float>, TPZFrontSym<float> >;
template class TPZFrontMatrix<float,TPZStackEqnStorage<float>, TPZFrontNonSym<float> >;
template class TPZFrontMatrix<float,TPZFileEqnStorage<float>, TPZFrontNonSym<float> >;

template class TPZFrontMatrix<double,TPZStackEqnStorage<double>, TPZFrontSym<double> >;
template class TPZFrontMatrix<double,TPZFileEqnStorage<double>, TPZFrontSym<double> >;
template class TPZFrontMatrix<double,TPZStackEqnStorage<double>, TPZFrontNonSym<double> >;
template class TPZFrontMatrix<double,TPZFileEqnStorage<double>, TPZFrontNonSym<double> >;

template class TPZFrontMatrix<long double,TPZStackEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZFrontMatrix<long double,TPZFileEqnStorage<long double>, TPZFrontSym<long double> >;
template class TPZFrontMatrix<long double,TPZStackEqnStorage<long double>, TPZFrontNonSym<long double> >;
template class TPZFrontMatrix<long double,TPZFileEqnStorage<long double>, TPZFrontNonSym<long double> >;

template class TPZFrontMatrix<std::complex<float> ,TPZStackEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZFrontMatrix<std::complex<float> ,TPZFileEqnStorage<std::complex<float> >, TPZFrontSym<std::complex<float> > >;
template class TPZFrontMatrix<std::complex<float> ,TPZStackEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;
template class TPZFrontMatrix<std::complex<float> ,TPZFileEqnStorage<std::complex<float> >, TPZFrontNonSym<std::complex<float> > >;

template class TPZFrontMatrix<std::complex<double> ,TPZStackEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZFrontMatrix<std::complex<double> ,TPZFileEqnStorage<std::complex<double> >, TPZFrontSym<std::complex<double> > >;
template class TPZFrontMatrix<std::complex<double> ,TPZStackEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;
template class TPZFrontMatrix<std::complex<double> ,TPZFileEqnStorage<std::complex<double> >, TPZFrontNonSym<std::complex<double> > >;

template class TPZFrontMatrix<std::complex<long double> ,TPZStackEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZFrontMatrix<std::complex<long double> ,TPZFileEqnStorage<std::complex<long double> >, TPZFrontSym<std::complex<long double> > >;
template class TPZFrontMatrix<std::complex<long double> ,TPZStackEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;
template class TPZFrontMatrix<std::complex<long double> ,TPZFileEqnStorage<std::complex<long double> >, TPZFrontNonSym<std::complex<long double> > >;
