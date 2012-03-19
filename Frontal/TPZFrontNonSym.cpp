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

DecomposeType TPZFrontNonSym::GetDecomposeType() const{
	return fDecomposeType;
}
void TPZFrontNonSym::PrintGlobal(const char *name, std::ostream& out = cout){
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
	out.flush();
}
void TPZFrontNonSym::Print(const char *name, std::ostream& out = cout) const
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
	out <<  "Local Freed Equations  "<< fFree.NElements() << endl;
	out <<  "position "<<   "Local Equation "<<endl;
	loop_limit=fFree.NElements();
	for(i=0;i<loop_limit;i++){
		out <<  i <<  "             " << fFree[i] << endl;
	}
	out <<  "Frontal Matrix "<< endl;
	if(fData.NElements() > 0) {
		for(i=0;i<fFront;i++){
			for(j=0;j<fFront;j++) out << Element(i,j) <<  " ";
			out << endl;
		}
	}
}
void TPZFrontNonSym::AllocData()
{
	fData.Resize(fMaxFront*fMaxFront);
	fGlobal.Fill(-1);
	fData.Fill(0.);
	fFront=0;
	//fLocal.Fill(-1);
}
void TPZFrontNonSym::Reset(int GlobalSize)
{
	fData.Resize(0);
	fFree.Resize(0);
	fFront=0;
	fGlobal.Resize(0);
	fLocal.Resize(GlobalSize);
	fLocal.Fill(-1);
	fMaxFront=0;
	fWork=0;
	fNextRigidBodyMode = GlobalSize;
}

int TPZFrontNonSym::NFree()
{
	return fFree.NElements();
}

int TPZFrontNonSym::Local(int global){
	int index;
	
	if(fLocal[global]!=-1) return fLocal[global];
	if(fFree.NElements()){
		index=fFree.Pop();
	}else{
		index=fFront;
	}
	// At this point we extend the frontmatrix
	if(index >= fMaxFront) {
		//	     cout << endl;
		//		cout << "Dynamically expanding the front size to " << index+fExpandRatio << endl;
		Expand(index+fExpandRatio);
	}
	fLocal[global]=index;
	if(index==fFront) fFront++;
	if(fGlobal.NElements()<=index) fGlobal.Resize(index+1);
	fGlobal[index]=global;
	return index;
}
void TPZFrontNonSym::FreeGlobal(int global)
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
void TPZFrontNonSym::DecomposeOneEquation(int ieq, TPZEqnArray &eqnarray)
{
	std::cout<<" fNextRigidBodyMode AQQQQ "<<fNextRigidBodyMode<<endl;
	
	//	if (fNextRigidBodyMode > fLocal.NElements()|| fNextRigidBodyMode == fLocal.NElements()) {
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
		Local(fNextRigidBodyMode);
	}
	
	TPZVec<REAL> AuxVecRow(fFront);
	TPZVec<REAL> AuxVecCol(fFront);
	
#ifdef USING_ATLAS
	cblas_dcopy(fFront, &Element(0, ilocal), 1, &AuxVecCol[0], 1);
	cblas_dcopy(fFront, &Element(ilocal, 0), fMaxFront, &AuxVecRow[0], 1);
#else
	for(i=0;i<fFront;i++){
		AuxVecCol[i]=Element(i,ilocal);
		AuxVecRow[i]=Element(ilocal,i);
	}
#endif
	//memcpy(&AuxVecCol[0], &Element(0,ilocal), fFront);
	//for(i=0;i<fFront;i++) AuxVecRow[i]=Element(ilocal,i);
	//memcpy(&AuxVecRow[0], &Element(ilocal, 0), fFront);
	
	
	REAL diag = AuxVecRow[ilocal];
	
	if (fabs(diag) < TOL ) {
		if (fNextRigidBodyMode < fLocal.NElements()) {
			int jlocal = Local(fNextRigidBodyMode);
			
			
			AuxVecRow[ilocal] = 1.;
			AuxVecRow[jlocal] = -1.;
			AuxVecCol[jlocal] = -1.;
			diag = 1.;
			Element(ilocal, jlocal) = -1.;
			Element(jlocal, ilocal) = -1.;
			Element(jlocal, jlocal) =  1.;
			fNextRigidBodyMode++;
		}
	}
	
#ifdef USING_ATLAS
	
	cblas_dscal(fFront, (1/diag), &AuxVecCol[0], 1);
	
#else
	for(i=0;i<fFront;i++){
		AuxVecCol[i]/=diag;
	}
#endif
	fWork+=fFront*fFront;
#ifdef USING_ATLAS
	//Blas utilizatioin
	CBLAS_ORDER order = CblasColMajor;
	//     CBLAS_UPLO up_lo = CblasUpper;
	long sz = fFront;
	long incx = 1;
	long stride = fMaxFront;
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
	long sz = fFront;
	long incx = 1;
	long stride = fMaxFront;
	double db = -1.;//AuxVec[ilocal];
	//resultado=cblas_dger(sz,&sz,&db,(double)&AuxVecCol[0],&incx,&AuxVecRow[0],&incx,&Element(0,0),&stride);
	cblas_dger(order, sz, sz, db,
			   &AuxVecCol[0], incx,
			   &AuxVecRow[0], incx, &Element(0,0), stride);
	
#endif
#ifndef	USING_ATLAS
#ifndef USING_BLAS
	int j;
	for(i=0;i<fFront;i++){
		for(j=0;j<fFront;j++) Element(i,j)-=AuxVecCol[i]*AuxVecRow[j];
	}
	
	/*     Print("After correct elimination",cout);
	 */
#endif
#endif
    AuxVecRow[ilocal]=1.;
    eqnarray.BeginEquation(ieq);
    eqnarray.AddTerm(ieq,diag);
    for(i=0;i<fFront;i++) {
		if(i!=ilocal && fGlobal[i]!= -1 && AuxVecRow[i] != 0.) eqnarray.AddTerm(fGlobal[i],AuxVecRow[i]);
	}
    eqnarray.EndEquation();
	
	eqnarray.BeginEquation(ieq);//, diag);
	eqnarray.AddTerm(ieq,1.);
	for(i=0;i<fFront;i++) {
        if(i!=ilocal && fGlobal[i]!= -1 && AuxVecCol[i] != 0.) eqnarray.AddTerm(fGlobal[i],AuxVecCol[i]);
	}
    eqnarray.EndEquation();
	
#ifdef USING_ATLAS
	TPZVec<double> zero(fFront, 0.);
	cblas_dcopy(fFront, &Element(0, ilocal), 1, &zero[0], 1);
	cblas_dcopy(fFront, &Element(ilocal, 0), fMaxFront, &zero[0], 1);
#else
	for(i=0;i<fFront;i++){
        Element(i,ilocal)=0.;
        Element(ilocal,i)=0.;
    }
#endif
	
    FreeGlobal(ieq);
    fDecomposeType=ELU;
	//	PrintGlobal("After", output);
}
void TPZFrontNonSym::AddKel(TPZFMatrix<REAL> &elmat, TPZVec<int> &sourceindex,  TPZVec<int> &destinationindex)
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
void TPZFrontNonSym::AddKel(TPZFMatrix<REAL> &elmat, TPZVec<int> &destinationindex)
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
	/*
	 output << "Dest Index " ;
	 for(i=0;i<nel;i++) output << destinationindex[i] << " ";
	 output << endl;
	 elmat.Print("Element matrix ",output);
	 PrintGlobal("After Assemb..." , output);
	 */
}

void TPZFrontNonSym::Expand(int larger) {
	//	PrintGlobal("Antes do Expande");
	fData.Resize(larger*larger,0.);
    int i,j;
	for(j=fFront-1;j>=0;j--){
		for(i=fFront-1;i>=0;i--){
			fData[j*larger + i]=fData[j*fMaxFront + i];
			if(j) fData[j*fMaxFront+i] = 0.;
		}
	}
	fMaxFront = larger;
	//	PrintGlobal("Depois do Expande");
	
}
void TPZFrontNonSym::Compress(){
	//	PrintGlobal("Before COmpress");
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
	
	fFront = nfound;
	fFree.Resize(0);	
	fGlobal.Resize(fFront);
	if(fData.NElements()==0) return;
	
	for(j = 0; j < nfound; j++){
		for(i = 0; i < nfound; i++){
			Element(i,j) = Element(from[i], from[j]);
			if(from[i]!=i || from[j]!=j) Element(from[i],from[j])=0.;
		}
		//		fGlobal[i] = fGlobal[from[i]];
		//		fLocal[fGlobal[i]] = i;
	}
	
	//	Print("After Compress", cout);
	//	PrintGlobal("After Compress",output);
}
void TPZFrontNonSym::SymbolicAddKel(TPZVec < int > & destinationindex)
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
void TPZFrontNonSym::SymbolicDecomposeEquations(int mineq, int maxeq)
{
	int i;
	for(i=mineq;i<=maxeq;i++) FreeGlobal(i);
}
void TPZFrontNonSym::DecomposeEquations(int mineq, int maxeq, TPZEqnArray & eqnarray){
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
TPZFrontNonSym::TPZFrontNonSym(int GlobalSize) : TPZFront(GlobalSize)
{
	fDecomposeType=ELU;
	fWork=0;
}

TPZFrontNonSym::TPZFrontNonSym() : TPZFront() {
	fDecomposeType=ELU;
	fWork=0;
}

TPZFrontNonSym::TPZFrontNonSym(const TPZFrontNonSym &cp) : TPZFront(cp) , fDecomposeType(cp.fDecomposeType) {
}

TPZFrontNonSym::~TPZFrontNonSym(){}



void TPZFrontNonSym::main()
{
	int i, j;
	/**
	 * 	Populates data structure
	 */
	int matsize=6;
	TPZFMatrix<REAL> TestMatrix(matsize,matsize);
	for(i=0;i<matsize;i++) {
		for(j=i;j<matsize;j++) {
			int random = rand();
			double rnd = (random*matsize)/0x7fff;
			TestMatrix(i,j)=rnd;
			TestMatrix(j,i)=TestMatrix(i,j);
			if(i==j) TestMatrix(i,j)=6000.;
		}
	}
	
	TPZFMatrix<REAL> Prova;
	Prova=TestMatrix;
	
	//	Prova.Decompose_Cholesky();
	Prova.Print("TPZFMatrix<REAL> Cholesky");
	
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
	
	
	TPZFMatrix<REAL> Load(matsize);
	
	for(i=0;i<matsize;i++) {
		int random = rand();
		double rnd = (random*matsize)/0x7fff;
		Load(i,0)=rnd;
	}
	
	TPZFMatrix<REAL> Load_2(matsize);
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

std::string TPZFrontNonSym::GetMatrixType(){
	return "Non symmetric matrix";
}

void TPZFrontNonSym::ExtractFrontMatrix(TPZFMatrix<REAL> &front)
{
	// Extend the front with the non initialized rigid body modes
	int ieq;
	int maxeq = fLocal.NElements();
	for (ieq = fNextRigidBodyMode; ieq< maxeq; ieq++) {
		int ilocal = Local(ieq);
		Element(ilocal, ilocal) = 1.;
	}
	
	int mineq = 0;
	for(mineq=0; mineq<maxeq; mineq++) if(fLocal[mineq] != -1) break;
	int numeq = maxeq-mineq;
	front.Redim(numeq,numeq);
	int jeq;
	for(ieq=mineq;ieq<maxeq;ieq++) {
		if(fLocal[ieq] == -1) continue;
		int il = ieq-mineq;
		for(jeq=0;jeq<maxeq;jeq++) {
			if(fLocal[jeq] == -1) continue;
			int jl = jeq-mineq;
			front(il,jl) = this->Element(fLocal[ieq],fLocal[jeq]);
		}
	}
	
}
