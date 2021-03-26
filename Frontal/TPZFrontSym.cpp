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

template<class TVar>
DecomposeType TPZFrontSym<TVar>::GetDecomposeType() const
{
	return this->fDecomposeType;
}
template<class TVar>
void TPZFrontSym<TVar>::PrintGlobal(const char *name, std::ostream& out){
	int64_t i, j;
	out << name << endl;
	for(i=0;i<this->fLocal.NElements();i++){
		if(this->fLocal[i]!=-1) out << i << " ";
	}
	out << endl;
	for(i=0;i<this->fLocal.NElements();i++){
		if(this->fLocal[i]==-1) continue;
		for(j=0;j<this->fLocal.NElements();j++){
			if(this->fLocal[j]==-1) continue;
			out << Element(this->fLocal[i],this->fLocal[j]) << " ";
		}
		out << endl;
	}
	out << endl;
}
template<class TVar>
void TPZFrontSym<TVar>::Print(const char *name, std::ostream& out) const
{
	if(name) out << name << endl;
	int64_t i,j,loop_limit;
	
	
	out <<  "Frontal Matrix Size          "<< this->fFront << endl;
	out <<  "Maximum Frontal Matrix Size  "<< this->fMaxFront << endl;
	
	out << endl;
	out <<  "Local Indexation "<< endl;
	out <<  "Position  "<<" Local index"<< endl;
	for(i=0;i<this->fLocal.NElements();i++){
		out <<   i <<  "         "  << this->fLocal[i] << endl;
	}
	
	out << endl;
	out <<  "Global Indexation "<< endl;
	out <<  "Position  "<<" Global index"<< endl;
	
	for(i=0;i<this->fGlobal.NElements();i++){
		out <<   i <<  "            "<< this->fGlobal[i] << endl;
	}
	
	out << endl;
	out <<  "Local Freed Equatoins  "<< this->fFree.NElements() << endl;
	out <<  "position "<<   "Local Equation "<<endl;
	loop_limit=this->fFree.NElements();
	for(i=0;i<loop_limit;i++){
		out <<  i <<  "             " << this->fFree[i] << endl;
	}
	out <<  "Frontal Matrix "<< endl;
	if(this->fData.NElements() > 0) {
		for(i=0;i<this->fFront;i++){
			for(j=0;j<this->fFront;j++) out << ((i<j) ? Element(i,j) : Element(j,i)) <<  " ";
			out << endl;
		}
	}
}
template<class TVar>
void TPZFrontSym<TVar>::AllocData()
{
	this->fData.Resize(this->fMaxFront*(this->fMaxFront+1)/2);
	this->fGlobal.Fill(-1);
	this->fData.Fill(0.);
	this->fFront=0;
	//this->fLocal.Fill(-1);
}
template<class TVar>
void TPZFrontSym<TVar>::Reset(int64_t GlobalSize)
{
	this->fData.Resize(0);
	this->fFree.Resize(0);
	this->fFront=0;
	this->fGlobal.Resize(0);
	this->fLocal.Resize(GlobalSize);
	this->fLocal.Fill(-1);
	this->fMaxFront=0;
}
template<class TVar>
int64_t TPZFrontSym<TVar>::NFree()
{
	return this->fFree.NElements();
}

template<class TVar>
int TPZFrontSym<TVar>::Local(int64_t global){
	int64_t index;
	if(this->fLocal[global]!=-1) return this->fLocal[global];
	if(this->fFree.NElements()){
		index=this->fFree.Pop();
	}else{
		index=this->fFront++;
	}
	this->fLocal[global]=index;
	// At this point we extend the frontmatrix
	if(index >= this->fMaxFront) {
		//	     cout << endl;
		//		cout << "Dynamically expanding the front size to " << index+fExpandRatio << endl;
		Expand(index+this->fExpandRatio);
	}
	if(this->fGlobal.NElements()<=index) this->fGlobal.Resize(index+1);
	this->fGlobal[index]=global;
	return index;
}
template<class TVar>
void TPZFrontSym<TVar>::FreeGlobal(int64_t global)
{
	if(this->fLocal[global]==-1){
		cout << "TPZFront FreeGlobal was called with wrong parameters !" << endl;
		return;
	}
	int64_t index;
	index=this->fLocal[global];
	this->fGlobal[index]=-1;
	this->fLocal[global]=-1;
	this->fFree.Push(index);
}

template<>
void TPZFrontSym<std::complex<float> >::DecomposeOneEquation(int64_t ieq, TPZEqnArray<std::complex<float> > &eqnarray)
{
	DebugStop();
}
template<>
void TPZFrontSym<std::complex<double> >::DecomposeOneEquation(int64_t ieq, TPZEqnArray<std::complex<double> > &eqnarray)
{
	DebugStop();
}
template<>
void TPZFrontSym<std::complex<long double> >::DecomposeOneEquation(int64_t ieq, TPZEqnArray<std::complex<long double> > &eqnarray)
{
	DebugStop();
}

template<class TVar>
void TPZFrontSym<TVar>::DecomposeOneEquation(int64_t ieq, TPZEqnArray<TVar> &eqnarray)
{
	//eqnarray.SetSymmetric();
	
	int64_t i, ilocal;
	ilocal = Local(ieq);
	TPZManVector<TVar> AuxVec(this->fFront);
	
	for(i=0;i<ilocal;i++) AuxVec[i]=Element(i,ilocal);
	for(i=ilocal;i<this->fFront;i++) AuxVec[i]=Element(ilocal,i);
	
	if(AuxVec[ilocal]<0 && this->fDecomposeType == ECholesky) {
		cout << "TPZFront::DecomposeOneEquation AuxVec[ilocal] < 0 " << AuxVec[ilocal] << " ilocal=" << ilocal << " fGlobal=" << this->fGlobal[ilocal] << endl;
	}
    TVar diag;
    if (this->fDecomposeType == ECholesky) {
        diag = sqrt(AuxVec[ilocal]);
        for(i=0;i<this->fFront;i++) AuxVec[i]/=diag;
    }
    else
    {
        diag = AuxVec[ilocal];
        for(i=0;i<this->fFront;i++) AuxVec[i]/=diag;
    }
	
	eqnarray.BeginEquation(ieq);       
	eqnarray.AddTerm(ieq, diag);
	
	
    if (this->fDecomposeType == ECholesky)
    {
        if(this->fProductMTData){
            this->ProductTensorMT( AuxVec, AuxVec );
        }
        else {
            for(int64_t j=0;j<this->fFront;j++){
                for(int64_t i=0;i<=j;i++){
                    Element(i,j)-=AuxVec[i]*AuxVec[j];
                }
            }		
        }

    }
    else
    {
        if(this->fProductMTData){
            this->fProductMTData->fDiagonal = diag;
            this->ProductTensorMT( AuxVec, AuxVec );
        }
        else {
            for(int64_t j=0;j<this->fFront;j++){
                for(int64_t i=0;i<=j;i++){
                    Element(i,j)-=AuxVec[i]*AuxVec[j]*diag;
                }
            }
        }
        
    }
		// #endif
	// #endif
	
	for(i=0;i<this->fFront;i++) {
		if(i!=ilocal && this->fGlobal[i]!= -1 && AuxVec[i] != 0.) eqnarray.AddTerm(this->fGlobal[i],AuxVec[i]);
	}
	eqnarray.EndEquation();
	
	for(i=0;i<ilocal;i++)Element(i,ilocal)=0.;
	for(i=ilocal;i<this->fFront;i++) Element(ilocal,i)=0.;
	
	FreeGlobal(ieq);
	//	PrintGlobal("After", output);
}


template <class TVar>
void TPZFrontSym<TVar>::TensorProductIJ(int ithread,typename TPZFront<TVar>::STensorProductMTData *data){
  if(!data) DebugStop();
#ifdef PZDEBUG
  TPZFrontSym<TVar> * matrix = dynamic_cast<TPZFrontSym<TVar> * > (data->fMat);
  if(matrix != this) DebugStop();
#endif
  while(data->fRunning){
    data->fWorkSem[ithread].Wait();
    if(!data->fRunning) break;
    const int n = data->fAuxVecCol->NElements();
    const int Nthreads = data->NThreads();
      
      if (data->fMat->GetDecomposeType() == ELDLt) {
          for(int j = 0+ithread; j < n; j += Nthreads){
              int i = 0;
              const TVar RowVal = data->fAuxVecRow->operator[](j)*data->fDiagonal;
              TVar * ColValPtr = &(data->fAuxVecCol->operator[](i));
              TVar * elemPtr = &this->Element4JGreatEqualI(i,j);
              for( ; i <= j; i++, ColValPtr++, elemPtr++ ){
                  (*elemPtr) -= (*ColValPtr) * RowVal;
              }///i
          }///j

      }
		else if(data->fMat->GetDecomposeType() == ECholesky)
        {
    for(int j = 0+ithread; j < n; j += Nthreads){
      int i = 0;
      const TVar RowVal = data->fAuxVecRow->operator[](j);
      TVar * ColValPtr = &(data->fAuxVecCol->operator[](i));
      TVar * elemPtr = &this->Element4JGreatEqualI(i,j);
      for( ; i <= j; i++, ColValPtr++, elemPtr++ ){
        (*elemPtr) -= (*ColValPtr) * RowVal;
      }///i
    }///j
        }
    data->WorkDone();
  }///while
}///void


template<class TVar>
void TPZFrontSym<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex)
{
	int64_t i, j, ilocal, jlocal, nel;
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
template<class TVar>
void TPZFrontSym<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int64_t> &destinationindex)
{
	int64_t i, j, ilocal, jlocal, nel;
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
/*TVar & TPZFrontSym<TVar>::Element(int i, int j){
 if(i>j){
 int i_temp=i;
 i=j;
 j=i_temp;
 //cout << "Changing row column indexes !"<<endl;
 }
 return fData[(j*(j+1))/2+i];
 }
 */
template<class TVar>
void TPZFrontSym<TVar>::Expand(int larger) {
	this->fData.Resize(larger*(larger+1)/2,0.);
	this->fMaxFront = larger;
}
template<class TVar>
void TPZFrontSym<TVar>::Compress(){
	//	PrintGlobal("Before COmpress",output);
	//	Print("Before Compress", cout);
	TPZStack <int64_t> from;
	int64_t nfound;
	int64_t i, j;
	for(i = 0; i < this->fFront; i++){
		if(this->fGlobal[i] != -1) from.Push(i);
	}
	
	/**
	 *First this->fLocal initialization
	 *Any needed updates is done on next loop
	 */
	nfound = from.NElements();
	for(i=0;i<nfound;i++) {
		this->fGlobal[i]=this->fGlobal[from[i]];
		//fGlobal[from[i]] = -1;
		this->fLocal[this->fGlobal[i]] = i;
	}
	for(;i<this->fGlobal.NElements();i++) this->fGlobal[i] = -1;
	
	if(nfound+this->fFree.NElements()!=this->fFront) cout << "TPZFront.Compress inconsistent data structure\n";
	int frontold = this->fFront;
	this->fFront = nfound;
	this->fFree.Resize(0);	
	this->fGlobal.Resize(this->fFront);
	if(this->fData.NElements()==0) return;
	
	for(j = 0; j < nfound; j++){
		for(i = 0; i <= j; i++){
			Element(i,j) = Element(from[i], from[j]);
		}
		//		fGlobal[i] = fGlobal[from[i]];
		//		this->fLocal[fGlobal[i]] = i;
	}
	for(;j<frontold;j++) {
		for(i=0;i<= j; i++) Element(i,j) = 0.;
	}
	
	//	Print("After Compress", cout);
	//	PrintGlobal("After Compress",output);
}
template<class TVar>
void TPZFrontSym<TVar>::SymbolicAddKel(TPZVec < int64_t > & destinationindex)
{
	int64_t i, loop_limit, aux;
	loop_limit=destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		aux=destinationindex[i];
		Local(aux);
		this->fFront = this->fGlobal.NElements();
	}
	this->fMaxFront=(this->fFront<this->fMaxFront)?this->fMaxFront:this->fFront;
	
}
template<class TVar>
void TPZFrontSym<TVar>::SymbolicDecomposeEquations(int64_t mineq, int64_t maxeq)
{
	int64_t i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}
template<class TVar>
void TPZFrontSym<TVar>::DecomposeEquations(int64_t mineq, int64_t maxeq, TPZEqnArray<TVar> & eqnarray){
	// message #1.1 to eqnarray:TPZEqnArray
	int64_t ieq;
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
template<class TVar>
TPZFrontSym<TVar>::TPZFrontSym(int64_t GlobalSize) : TPZRegisterClassId(&TPZFrontSym<TVar>::ClassId),
TPZFront<TVar>(GlobalSize)
{
	this->fDecomposeType=ECholesky;
}
template<class TVar>
TPZFrontSym<TVar>::TPZFrontSym(): TPZRegisterClassId(&TPZFrontSym<TVar>::ClassId){
	this->fDecomposeType=ECholesky;
}


template<class TVar>
TPZFrontSym<TVar>::~TPZFrontSym(){}

template<class TVar>
void TPZFrontSym<TVar>::main()
{
	int i, j;
	/**
	 * 	Populates data structure
	 */
	int matsize=6;
	TPZFMatrix<TVar> TestMatrix(matsize,matsize);
	for(i=0;i<matsize;i++) {
		for(j=i;j<matsize;j++) {
			int random = rand();
			double rnd = (random*matsize)/0x7fff;
			TestMatrix(i,j)=rnd;
			TestMatrix(j,i)=TestMatrix(i,j);
			if(i==j) TestMatrix(i,j)=6000.;
		}
	}
	
	TPZFMatrix<TVar> Prova;
	Prova=TestMatrix;
	
	//	Prova.Decompose_Cholesky();
	Prova.Print("TPZFMatrix<TVar> Cholesky");
	
	TPZFrontSym TestFront(matsize);
	
	
	TPZVec<int64_t> DestIndex(matsize);
	for(i=0;i<matsize;i++) DestIndex[i]=i;
	
	TestFront.SymbolicAddKel(DestIndex);
	TestFront.SymbolicDecomposeEquations(0,matsize-1);
	
	std::string OutFile;
	OutFile = "TPZFrontSymTest.txt";
	
	ofstream output(OutFile.c_str(),ios::app);
	
	TestFront.Compress();
	
	TestFront.AllocData();
	
	TestFront.AddKel(TestMatrix, DestIndex);
	TPZEqnArray<TVar> Result;
	
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
	
	
	TPZFMatrix<TVar> Load(matsize);
	
	for(i=0;i<matsize;i++) {
		int random = rand();
		double rnd = (random*matsize)/0x7fff;
		Load(i,0)=rnd;
	}
	
	TPZFMatrix<TVar> Load_2(matsize);
	Load_2=Load;
	
	//	Prova.Subst_Forward(&Load);
	//	Prova.Subst_Backward(&Load);
	
	
	DecomposeType decType = ECholesky;
	Prova.SolveDirect(Load, decType);
	
	Load.Print("load");
	//TestFront.Print(OutFile, output);
	
	Result.EqnForward(Load_2, decType);
	Result.EqnBackward(Load_2, decType);
	
	Load_2.Print("Eqn");
	
	
	
}

template<class TVar>
std::string TPZFrontSym<TVar>::GetMatrixType(){
	return "Symmetric matrix";
}

template<class TVar>
void TPZFrontSym<TVar>::ExtractFrontMatrix(TPZFMatrix<TVar> &front)
{
	int64_t maxeq = this->fLocal.NElements();
	int64_t mineq = 0;
	for(mineq=0; mineq<maxeq; mineq++) if(this->fLocal[mineq] != -1) break;
	int64_t numeq = maxeq-mineq;
	front.Redim(numeq,numeq);
	int64_t ieq,jeq;
	for(ieq=mineq;ieq<maxeq;ieq++) {
		if(this->fLocal[ieq] == -1) continue;
		int64_t il = ieq-mineq;
		for(jeq=ieq;jeq<maxeq;jeq++) {
			if(this->fLocal[jeq] == -1) continue;
			int64_t jl = jeq-mineq;
			front(il,jl) = this->Element(this->fLocal[ieq],this->fLocal[jeq]);
			front(jl,il) = front(il,jl);
		}
	}
	
}

template<class TVar>
int TPZFrontSym<TVar>::ClassId() const{
    return Hash("TPZFrontSym") ^ TPZFront<TVar>::ClassId() << 1;
}

template class TPZFrontSym<float>;
template class TPZFrontSym<std::complex<float> >;

template class TPZFrontSym<double>;
template class TPZFrontSym<std::complex<double> >;

template class TPZFrontSym<long double>;
template class TPZFrontSym<std::complex<long double> >;
