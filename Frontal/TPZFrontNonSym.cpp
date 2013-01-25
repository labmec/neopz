/**
 * @file
 * @brief Contains the implementation of the TPZFrontNonSym methods.
 */

#include "TPZFrontNonSym.h"
#include "tpzeqnarray.h"

#ifdef USING_BLAS
void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double  *X, const int incX,
                const double  *Y, const int incY, double  *A, const int lda);
#endif
#ifdef USING_ATLAS
void cblas_dger(const enum CBLAS_ORDER order, const int M, const int N,
                const double alpha, const double  *X, const int incX,
                const double  *Y, const int incY, double  *A, const int lda);
void cblas_dscal(const int N, const double alpha, double *X, const int incX);
void cblas_dcopy(const int N, const double *X, const int incX,
                 double *Y, const int incY);
#endif

/** @brief Initializing tolerance for current implementations */
const REAL TOL=1.e-10;

using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("pz.frontstrmatrix.frontnonsym"));
#endif

template<class TVar>
DecomposeType TPZFrontNonSym<TVar>::GetDecomposeType() const{
	return fDecomposeType;
}

template<class TVar>
void TPZFrontNonSym<TVar>::PrintGlobal(const char *name, std::ostream& out){
	int i, j;
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
	out.flush();
}

template<class TVar>
void TPZFrontNonSym<TVar>::Print(const char *name, std::ostream& out) const
{
	if(name) out << name << endl;
	int i,j,loop_limit;
	
	
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
	out <<  "Local Freed Equations  "<< this->fFree.NElements() << endl;
	out <<  "position "<<   "Local Equation "<<endl;
	loop_limit=this->fFree.NElements();
	for(i=0;i<loop_limit;i++){
		out <<  i <<  "             " << this->fFree[i] << endl;
	}
	out <<  "Frontal Matrix "<< endl;
	if(this->fData.NElements() > 0) {
		for(i=0;i<this->fFront;i++){
			for(j=0;j<this->fFront;j++) out << Element(i,j) <<  " ";
			out << endl;
		}
	}
}

template<class TVar>
void TPZFrontNonSym<TVar>::AllocData()
{
	this->fData.Resize(this->fMaxFront*this->fMaxFront);
	this->fGlobal.Fill(-1);
	this->fData.Fill(0.);
	this->fFront=0;
	//this->fLocal.Fill(-1);
}

template<class TVar>
void TPZFrontNonSym<TVar>::Reset(int GlobalSize)
{
	this->fData.Resize(0);
	this->fFree.Resize(0);
	this->fFront=0;
	this->fGlobal.Resize(0);
	this->fLocal.Resize(GlobalSize);
	this->fLocal.Fill(-1);
	this->fMaxFront=0;
	this->fWork=0;
	this->fNextRigidBodyMode = GlobalSize;
}

template<class TVar>
int TPZFrontNonSym<TVar>::NFree()
{
	return this->fFree.NElements();
}

template<class TVar>
int TPZFrontNonSym<TVar>::Local(int global){
	int index;
	
	if(this->fLocal[global]!=-1) return this->fLocal[global];
	if(this->fFree.NElements()){
		index=this->fFree.Pop();
	}else{
		index=this->fFront;
	}
	// At this point we extend the frontmatrix
	if(index >= this->fMaxFront) {
		//	     cout << endl;
		//		cout << "Dynamically expanding the front size to " << index+this->fExpandRatio << endl;
		Expand(index+this->fExpandRatio);
	}
	this->fLocal[global]=index;
	if(index==this->fFront) this->fFront++;
	if(this->fGlobal.NElements()<=index) this->fGlobal.Resize(index+1);
	this->fGlobal[index]=global;
	return index;
}

template<class TVar>
void TPZFrontNonSym<TVar>::FreeGlobal(int global)
{
	if(this->fLocal[global]==-1){
		cout << "TPZFront FreeGlobal was called with wrong parameters !" << endl;
		return;
	}
	int index;
	index=this->fLocal[global];
	this->fGlobal[index]=-1;
	this->fLocal[global]=-1;
	this->fFree.Push(index);
}

template<>
void TPZFrontNonSym<std::complex<float> >::DecomposeOneEquation(int ieq, TPZEqnArray<std::complex<float> > &eqnarray)
{
    DebugStop();
}
template<>
void TPZFrontNonSym<std::complex<double> >::DecomposeOneEquation(int ieq, TPZEqnArray<std::complex<double> > &eqnarray)
{
    DebugStop();
}
template<>
void TPZFrontNonSym<std::complex<long double> >::DecomposeOneEquation(int ieq, TPZEqnArray<std::complex<long double> > &eqnarray)
{
    DebugStop();
}
template<class TVar>
void TPZFrontNonSym<TVar>::DecomposeOneEquation(int ieq, TPZEqnArray<TVar> &eqnarray)
{
	//std::cout<<" fNextRigidBodyMode AQQQQ "<<this->fNextRigidBodyMode<<endl;
	
	//	if (this->fNextRigidBodyMode > this->fLocal.NElements()|| this->fNextRigidBodyMode == this->fLocal.NElements()) {
	//	DebugStop();
	//	}
	//eqnarray.SetNonSymmetric();
	int i, ilocal;
	ilocal = Local(ieq);

#ifdef LOG4CXX
	{
       	double diagonal=fabs(Element(ilocal,ilocal));
		std::stringstream sout;
		sout<<" Valor da Diagonal ( " << ilocal<< ", "<< ilocal<< " ) = "<< diagonal<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
	
	if(fabs(Element(ilocal,ilocal)) <TOL) 
	{
		Local(this->fNextRigidBodyMode);
	}
	
	TPZVec<TVar> AuxVecRow(this->fFront);
	TPZVec<TVar> AuxVecCol(this->fFront);
	
#ifdef USING_ATLAS
	cblas_dcopy(this->fFront, &Element(0, ilocal), 1, &AuxVecCol[0], 1);
	cblas_dcopy(this->fFront, &Element(ilocal, 0), this->fMaxFront, &AuxVecRow[0], 1);
#else
	for(i=0;i<this->fFront;i++){
		AuxVecCol[i]=Element(i,ilocal);
		AuxVecRow[i]=Element(ilocal,i);
	}
#endif
	
	TVar diag = AuxVecRow[ilocal];
	
	if (fabs(diag) < TOL ) {
		if (this->fNextRigidBodyMode < this->fLocal.NElements()) {
			int jlocal = Local(this->fNextRigidBodyMode);
			
			
			AuxVecRow[ilocal] = 1.;
			AuxVecRow[jlocal] = -1.;
			AuxVecCol[jlocal] = -1.;
			diag = 1.;
			Element(ilocal, jlocal) = -1.;
			Element(jlocal, ilocal) = -1.;
			Element(jlocal, jlocal) =  1.;
			this->fNextRigidBodyMode++;
		}
	}
	
#ifdef USING_ATLAS
	
	cblas_dscal(this->fFront, (1/diag), &AuxVecCol[0], 1);
	
#else
	for(i=0;i<this->fFront;i++){
		AuxVecCol[i]/=diag;
	}
#endif
	this->fWork+=this->fFront*this->fFront;
#ifdef USING_ATLAS
	//Blas utilizatioin
	CBLAS_ORDER order = CblasColMajor;
	//     CBLAS_UPLO up_lo = CblasUpper;
	long sz = this->fFront;
	long incx = 1;
	long stride = this->fMaxFront;
	double db = -1.;//AuxVec[ilocal];
	//resultado=cblas_dger(sz,&sz,&db,(double)&AuxVecCol[0],&incx,&AuxVecRow[0],&incx,&Element(0,0),&stride);
	cblas_dger(order, sz, sz, db,
			   &AuxVecCol[0], incx,
			   &AuxVecRow[0], incx, &Element(0,0), stride);
#endif
#ifdef USING_BLAS
	//Blas utilizatioin  
	CBLAS_ORDER order = CblasColMajor;
	//     CBLAS_UPLO up_lo = CblasUpper;
	long sz = this->fFront;
	long incx = 1;
	long stride = this->fMaxFront;
	double db = -1.;//AuxVec[ilocal];
	//resultado=cblas_dger(sz,&sz,&db,(double)&AuxVecCol[0],&incx,&AuxVecRow[0],&incx,&Element(0,0),&stride);
	cblas_dger(order, sz, sz, db,
			   &AuxVecCol[0], incx,
			   &AuxVecRow[0], incx, &Element(0,0), stride);
	
#endif
#ifndef	USING_ATLAS
#ifndef USING_BLAS
	int j;
	for(i=0;i<this->fFront;i++){
		for(j=0;j<this->fFront;j++) Element(i,j)-=AuxVecCol[i]*AuxVecRow[j];
	}
	
	/*     Print("After correct elimination",cout);
	 */
#endif
#endif
    AuxVecRow[ilocal]=1.;
    eqnarray.BeginEquation(ieq);
    eqnarray.AddTerm(ieq,diag);
    for(i=0;i<this->fFront;i++) {
		if(i!=ilocal && this->fGlobal[i]!= -1 && AuxVecRow[i] != 0.) eqnarray.AddTerm(this->fGlobal[i],AuxVecRow[i]);
	}
    eqnarray.EndEquation();
	
	eqnarray.BeginEquation(ieq);//, diag);
	eqnarray.AddTerm(ieq,1.);
	for(i=0;i<this->fFront;i++) {
        if(i!=ilocal && this->fGlobal[i]!= -1 && AuxVecCol[i] != 0.) eqnarray.AddTerm(this->fGlobal[i],AuxVecCol[i]);
	}
    eqnarray.EndEquation();
	
#ifdef USING_ATLAS
	TPZVec<double> zero(this->fFront, 0.);
	cblas_dcopy(this->fFront, &Element(0, ilocal), 1, &zero[0], 1);
	cblas_dcopy(this->fFront, &Element(ilocal, 0), this->fMaxFront, &zero[0], 1);
#else
	for(i=0;i<this->fFront;i++){
        Element(i,ilocal)=0.;
        Element(ilocal,i)=0.;
    }
#endif
	
    FreeGlobal(ieq);
    fDecomposeType=ELU;
	//	PrintGlobal("After", output);
}

template<class TVar>
void TPZFrontNonSym<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex)
{
	int i, j, ilocal, jlocal, nel;
	nel=sourceindex.NElements();
	for (i = 0; i < nel; i++) {
		// message #1.1.1 to this:TPZFront
		ilocal = this->Local(destinationindex[i]);
		for (j = 0; j < nel; j++) {
			// message #1.1.2.1 to this:TPZFront
			jlocal = this->Local(destinationindex[j]);
			
			// message #1.1.2.2 to this:TPZFront
			this->Element(ilocal, jlocal)+=elmat(sourceindex[i],sourceindex[j]);
		}
	}
}

template<class TVar>
void TPZFrontNonSym<TVar>::AddKel(TPZFMatrix<TVar> &elmat, TPZVec<int> &destinationindex)
{
    int i, j, ilocal, jlocal, nel;
    nel = destinationindex.NElements();
    for(i=0;i<nel;i++){
        ilocal = this->Local(destinationindex[i]);
		for(j=0;j<nel;j++) {
            jlocal=this->Local(destinationindex[j]);
            this->Element(ilocal, jlocal)+=elmat(i,j);
        }
    }
}

template<class TVar>
void TPZFrontNonSym<TVar>::Expand(int larger) {
	//	PrintGlobal("Antes do Expande");
	this->fData.Resize(larger*larger,0.);
    int i,j;
	for(j=this->fFront-1;j>=0;j--){
		for(i=this->fFront-1;i>=0;i--){
			this->fData[j*larger + i]=this->fData[j*this->fMaxFront + i];
			if(j) this->fData[j*this->fMaxFront+i] = 0.;
		}
	}
	this->fMaxFront = larger;
}

template<class TVar>
void TPZFrontNonSym<TVar>::Compress(){
	//	PrintGlobal("Before COmpress");
	//	Print("Before Compress", cout);
	TPZStack <int> from;
	int nfound;
	int i, j;
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
		//this->fGlobal[from[i]] = -1;
		this->fLocal[this->fGlobal[i]] = i;
	}
	for(;i<this->fGlobal.NElements();i++) this->fGlobal[i] = -1;
	
	if(nfound+this->fFree.NElements()!=this->fFront) cout << "TPZFront.Compress inconsistent data structure\n";
	
	this->fFront = nfound;
	this->fFree.Resize(0);	
	this->fGlobal.Resize(this->fFront);
	if(this->fData.NElements()==0) return;
	
	for(j = 0; j < nfound; j++){
		for(i = 0; i < nfound; i++){
			Element(i,j) = Element(from[i], from[j]);
			if(from[i]!=i || from[j]!=j) Element(from[i],from[j])=0.;
		}
	}
}

template<class TVar>
void TPZFrontNonSym<TVar>::SymbolicAddKel(TPZVec < int > & destinationindex)
{
	int i, loop_limit, aux;
	loop_limit=destinationindex.NElements();
	for(i=0;i<loop_limit;i++){
		aux=destinationindex[i];
		Local(aux);
		this->fFront = this->fGlobal.NElements();
	}
	this->fMaxFront=(this->fFront<this->fMaxFront)?this->fMaxFront:this->fFront;
	
}

template<class TVar>
void TPZFrontNonSym<TVar>::SymbolicDecomposeEquations(int mineq, int maxeq)
{
	int i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}

template<class TVar>
void TPZFrontNonSym<TVar>::DecomposeEquations(int mineq, int maxeq, TPZEqnArray<TVar> & eqnarray){
	// message #1.1 to eqnarray:TPZEqnArray
	int ieq;
	
	eqnarray.Reset();
	eqnarray.SetNonSymmetric();
	
	for (ieq = mineq; ieq <= maxeq; ieq++) {
		// message #1.2.1 to this:TPZFront
		this->DecomposeOneEquation(ieq, eqnarray);
		//			this->Print("Teste.txt",output);
	}
}

template<class TVar>
TPZFrontNonSym<TVar>::TPZFrontNonSym(int GlobalSize) : TPZFront<TVar>(GlobalSize)
{
	fDecomposeType=ELU;
	this->fWork=0;
}

template<class TVar>
TPZFrontNonSym<TVar>::TPZFrontNonSym() : TPZFront<TVar>() {
	fDecomposeType=ELU;
	this->fWork=0;
}

template<class TVar>
TPZFrontNonSym<TVar>::TPZFrontNonSym(const TPZFrontNonSym &cp) : TPZFront<TVar>(cp) , fDecomposeType(cp.fDecomposeType) {
}

template<class TVar>
TPZFrontNonSym<TVar>::~TPZFrontNonSym(){}



template<class TVar>
void TPZFrontNonSym<TVar>::main()
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
	
	TPZFrontNonSym TestFront(matsize);
	
	
	TPZVec<int> DestIndex(matsize);
	for(i=0;i<matsize;i++) DestIndex[i]=i;
	
	TestFront.SymbolicAddKel(DestIndex);
	TestFront.SymbolicDecomposeEquations(0,matsize-1); 
	
	std::string OutFile;
	OutFile = "TPZFrontNonSymTest.txt";
	
	ofstream output(OutFile.c_str(),ios::app);
	
	TestFront.Compress();
	
	TestFront.AllocData();
	
	TestFront.AddKel(TestMatrix, DestIndex);
	TPZEqnArray<TVar> Result;

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
	
	DecomposeType decType = ECholesky;
	Prova.SolveDirect(Load, decType);
	
	Load.Print();
	//TestFront.Print(OutFile, output);
	
	Result.EqnForward(Load_2, decType);
	Result.EqnBackward(Load_2, decType);
	
	Load_2.Print("Eqn");
}

template<class TVar>
std::string TPZFrontNonSym<TVar>::GetMatrixType(){
	return "Non symmetric matrix";
}

template<class TVar>
void TPZFrontNonSym<TVar>::ExtractFrontMatrix(TPZFMatrix<TVar> &front)
{
	// Extend the front with the non initialized rigid body modes
	int ieq;
	int maxeq = this->fLocal.NElements();
	for (ieq = this->fNextRigidBodyMode; ieq< maxeq; ieq++) {
		int ilocal = Local(ieq);
		Element(ilocal, ilocal) = 1.;
	}
	
	int mineq = 0;
	for(mineq=0; mineq<maxeq; mineq++) if(this->fLocal[mineq] != -1) break;
	int numeq = maxeq-mineq;
	front.Redim(numeq,numeq);
	int jeq;
	for(ieq=mineq;ieq<maxeq;ieq++) {
		if(this->fLocal[ieq] == -1) continue;
		int il = ieq-mineq;
		for(jeq=0;jeq<maxeq;jeq++) {
			if(this->fLocal[jeq] == -1) continue;
			int jl = jeq-mineq;
			front(il,jl) = this->Element(this->fLocal[ieq],this->fLocal[jeq]);
		}
	}
	
}

template class TPZFrontNonSym<float>;
template class TPZFrontNonSym<double>;
template class TPZFrontNonSym<long double>;

template class TPZFrontNonSym<std::complex<float> >;
template class TPZFrontNonSym<std::complex<double> >;
template class TPZFrontNonSym<std::complex<long double> >;
