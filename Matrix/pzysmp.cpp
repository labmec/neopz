/**
 * @file
 * @brief Contains the implementation of the TPZFYsmpMatrix methods.
 */

#include "pzysmp.h"
#include "pzfmatrix.h"
#include "pzvec.h"

#include <memory.h>
#include <string>
#include <map>
#include <pthread.h>

#include "pz_pthread.h"

// #ifdef USING_BLAS
// 	double cblas_ddoti(const int N, const double *X, const int *indx,
//                    const double *Y);
// #endif

using namespace std;

// ****************************************************************************
//
// Constructors and the destructor
//
// ****************************************************************************
template<class TVar>
void TPZFYsmpMatrix<TVar>::InitializeData(){}

template<class TVar>
void TPZFYsmpMatrix<TVar>::MultiplyDummy(TPZFYsmpMatrix<TVar> & B, TPZFYsmpMatrix<TVar> & Res){
    long i,j,k;
    if (B.Rows()!=this->Rows()) return;
    long rows = this->Rows();
    REAL aux=0.;
    for(i=0;i<rows;i++){
        for(j=0;j<rows;j++){
            for(k=0;k<rows;k++){
				// C[i][j] += A[i][k]*B[k][j];
				aux+=GetVal(i,k)*B.GetVal(i,k);
			}
			Res.PutVal(i,j,aux);
			aux=0.;
		}
    }
}

// ****************************************************************************
//
// Constructor
//
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix(const TPZVerySparseMatrix<TVar> &cp) : TPZMatrix<TVar>
()
{
	*this = cp;
}

template<class TVar>
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix() : TPZMatrix<TVar>(), fIA(),fJA(),fA(),fDiag()
{
}

template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator=(const TPZFYsmpMatrix<TVar> &cp) {
	// Deletes everything associated with a TPZFYsmpMatrix
	TPZMatrix<TVar>::operator=(cp);
    fIA = cp.fIA;
    fA = new TPZVec<TVar>((TPZVec<TVar>&)cp.fA);
    fJA = cp.fJA;
    fDiag = new TPZVec<TVar>((TPZVec<TVar>&)cp.fDiag);
	return *this;
}

template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator=(const TPZVerySparseMatrix<TVar> &cp)
{
	// Deletes everything associated with a TPZFYsmpMatrix
	long nrows = cp.Rows();
	
	long count = 0, c = 0, r = 0;
	
	count = cp.fExtraSparseData.size();
    fJA = new TPZVec<long>(count);
    fA = new TPZVec<TVar>(count);
    fIA = new TPZVec<long>(nrows+1);
	(fIA)->operator[](0) = 0;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;
	
	typename map< pair<long,long>, TVar>::const_iterator it;
	c = 0;
	r = 0;
	for(it=cp.fExtraSparseData.begin(); it!= cp.fExtraSparseData.end(); it++)
	{
		long row = it->first.first;
		if(r != row)
		{
			r++;
			while(r < row) 
			{
				IA[r] = c;
				r++;
			}
			IA[row] = c;
			r = row;
		}
		long col = it->first.second;
		JA[c] = col;
		TVar val = it->second;
		A[c] = val;
		c++;
	}
	r++;
	while(r<=nrows)
	{
		IA[r] = c;
		r++;
	}
	return *this;
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::PutVal(const long row, const long col, const TVar &Value){
    long k;
    int flag=0;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;

    for(k=IA[row];k<IA[row+1];k++){
		if(JA[k]==col){
			flag=1;
			A[k]=Value;
			break;
		}
    }
    if(!flag) 
    {
		cout << "TPZFYsmpMatrix::PutVal: Non existing position on sparse matrix: line = " << row << " column " << col << endl;
		DebugStop();
		return 0;
    }
    else
    {
		return 1;
    }
}
template<class TVar>
void TPZFYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<long> & destinationindex){
    long i,j,k = 0;
    TVar value=0.;
    long ipos,jpos;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;

    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(i,j);
            //cout << "j= " << j << endl;
            if(value != 0.){
                //cout << "fIA[ipos] " << fIA[ipos] << "     fIA[ipos+1] " << fIA[ipos+1] << endl;
                int flag = 0;
				k++;
				if(k >= IA[ipos] && k < IA[ipos+1] && JA[k]==jpos)
				{ // OK -> elements in sequence
					A[k]+=value;
					flag = 1;
				}else
				{
					for(k=IA[ipos];k<IA[ipos+1];k++){
						if(JA[k]==jpos || JA[k]==-1){
							//cout << "fJA[k] " << fJA[k] << " jpos "<< jpos << "   " << value << endl;
							//cout << "k " << k << "   "<< jpos << "   " << value << endl;
							flag=1;
							if(JA[k]==-1){
								JA[k]=jpos;
								A[k]=value;
								// cout << jpos << "   " << value << endl;
								break;
							}else{
								A[k]+=value;
								break;
							}
						}
					}
				}
                if(!flag) cout << "TPZFYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << endl;         }
        }
    }
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<long> & sourceindex, TPZVec<long> & destinationindex){
	long i,j,k = 0;
	TVar value=0.;
	long ipos,jpos;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;

	for(i=0;i<sourceindex.NElements();i++){
		for(j=0;j<sourceindex.NElements();j++){
			ipos=destinationindex[i];
			jpos=destinationindex[j];
			value=elmat.GetVal(sourceindex[i],sourceindex[j]);
            //cout << "j= " << j << endl;
			if(value != 0.){
                //cout << "fIA[ipos] " << fIA[ipos] << "     fIA[ipos+1] " << fIA[ipos+1] << endl;
				int flag = 0;
				k++;
				if(k >= IA[ipos] && k < IA[ipos+1] && JA[k]==jpos)
				{ // OK -> elements in sequence
					A[k]+=value;
					flag = 1;
				}else
				{
					for(k=IA[ipos];k<IA[ipos+1];k++){
						if(JA[k]==jpos || JA[k]==-1){
							//cout << "fJA[k] " << fJA[k] << " jpos "<< jpos << "   " << value << endl;
							//cout << "k " << k << "   "<< jpos << "   " << value << endl;
							flag=1;
							if(JA[k]==-1){
								JA[k]=jpos;
								A[k]=value;
								// cout << jpos << "   " << value << endl;
								break;
							}else{
								A[k]+=value;
								break;
							}
						}
					}
				}
				if(!flag) cout << "TPZFYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << endl;         }
		}
	}
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::AddKelOld(TPZFMatrix<TVar> & elmat, TPZVec < int > & destinationindex){
	long i=0;
	long j=0;
	long ilocal=0;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;

	//  int jlocal=0;
	long nel = destinationindex.NElements();
	std::multimap<long,long> mapindex;
	std::multimap<long,long>::iterator hint = mapindex.begin();
	for(i=0;i<nel;i++){
		ilocal = destinationindex[i];
		hint = mapindex.insert(hint,std::make_pair(ilocal,i));
		//    mapindex[ilocal] = i;
	}
	for(i=0;i<nel;i++){
		ilocal = destinationindex[i];
		long jfirst = IA[ilocal];
		long jlast = IA[ilocal+1];
		long *Aptr = &JA[jfirst];
		long *AptrLast = &JA[jlast];
		j=0;
		std::multimap<long,long>::iterator itelmat = mapindex.begin();
		while(j<nel && Aptr != AptrLast)
		{
			if(*Aptr == (*itelmat).first)
			{
				long jel = (*itelmat).second;
				A[jfirst] += elmat(i,jel);
				itelmat++;
				if(itelmat != mapindex.end() && (*itelmat).second != jel)
				{
					Aptr++;
					jfirst++;
				}
				j++;
			}
			else if(*Aptr < (*itelmat).first)
			{
				Aptr++;
				jfirst++;
			}
			else if(*Aptr > (*itelmat).second)
			{
				std::cout << __PRETTY_FUNCTION__ << " inconsistent\n";
				long *iptr = &JA[jfirst];
				while(iptr < AptrLast) 
				{
					cout << *iptr << " ";
					iptr++;
				}
				cout << endl;
				std::multimap<long,long>::iterator itelmat2 = mapindex.begin();
				for(;itelmat2 != mapindex.end(); itelmat2++)
				{
					cout << (*itelmat2).first << "/" << (*itelmat2).second << " ";
				}
				cout << endl;
				
			}
		}
		if(j!= nel)
		{
			std::cout << __PRETTY_FUNCTION__ << " inconsistent2 j = " << j << " nel " << nel << "\n";
			long *iptr = &JA[jfirst];
			while(iptr < AptrLast) 
			{
				cout << *iptr << " ";
				iptr++;
			}
			cout << endl;
			std::multimap<long,long>::iterator itelmat2 = mapindex.begin();
			for(;itelmat2 != mapindex.end(); itelmat2++)
			{
				cout << (*itelmat2).first << "/" << (*itelmat2).second << " ";
			}
			cout << endl;
		}
	}
	
}

template<class TVar>
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix(const long rows,const long cols ) : TPZMatrix<TVar>(rows,cols) {
	// Constructs an empty TPZFYsmpMatrix
	//    fSolver = -1;
	fSymmetric = 0;
	//    fMaxIterations = 4;
	//    fSORRelaxation = 1.;
	fDiag = new TPZVec<TVar>;
	fA = new TPZVec<TVar>;
	fIA = new TPZVec<long>(1,0);
	fJA = new TPZVec<long>;
#ifdef CONSTRUCTOR
	cerr << "TPZFYsmpMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZFYsmpMatrix<TVar>::~TPZFYsmpMatrix() {
	// Deletes everything associated with a TPZFYsmpMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZFYsmpMatrix()\n";
#endif
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

template<class TVar>
const TVar & TPZFYsmpMatrix<TVar>::GetVal(const long row,const long col ) const {
	// Get the matrix entry at (row,col) without bound checking
	
	// Now look through the requested row and see if there is anything
	// in column col
	/*  int loccol = col+1;
	 for(int ic=fIA[row]-1 ; ic < fIA[row+1]-1; ic++ ) {
	 if ( fJA[ic] == loccol ) return fA[ic];
	 }*/
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	long loccol = col;
	for(long ic=IA[row] ; ic < IA[row+1]; ic++ ) {
		if ( JA[ic] == loccol && JA[ic] != -1 ) return A[ic];
	}
	return this->gZero;
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZFYsmpMatrix<TVar>::MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
									 TPZFMatrix<TVar> &z,
									 const TVar alpha,const TVar beta,const int opt,const int stride )  {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
	{
		cout << "\nERROR! in TPZVerySparseMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
		return;
	}
	
	long  ir, ic, icol, xcols;
	xcols = x.Cols();
	TVar sum;
	long  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	for(ic=0; ic<xcols; ic++) {
		TVar *zp = &(z(0,ic));
		if(beta != 0) {
			const TVar *yp = &(y.g(0,0));
			TVar *zlast = zp+r*stride;
			if(beta != 1. || (&z != &y && stride != 1)) {
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			}
			else if(&z != &y) {
				memcpy(zp,yp,r*sizeof(TVar));
			}
		} else {
			TVar *zp = &(z(0,0)), *zlast = zp+r*stride;
			while(zp != zlast) {
				*zp = 0.;
				zp += stride;
			}
		}
	}
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;

	// Compute alpha * A * x
	if(xcols == 1 && stride == 1 && opt == 0)
	{
		if(this->Cols() != x.Rows()*stride || this->Rows() != y.Rows()*stride)
		{
			cout << "\nERROR! in TPZFYsmpMatrix::MultiplyAddMT: incompatible dimensions in opt=false\n";
			return;
		} 
		for(ir=0; ir<r; ir++) {
			long icolmin = IA[ir];
			long icolmax = IA[ir+1];
			const TVar *xptr = &(x.g(0,0));
			TVar *Aptr = &A[0];
			long *JAptr = &JA[0];
			for(sum = 0.0, icol=icolmin; icol<icolmax; icol++ ) {
				sum += Aptr[icol] * xptr[JAptr[icol]];
			}
			z(ir,0) += alpha * sum;
		}
	}
	else 
	{
		for(ic=0; ic<xcols; ic++) {
			if(opt == 0) {
				
				for(ir=0; ir<this->Rows(); ir++) {
					for(sum = 0.0, icol=IA[ir]; icol<IA[ir+1]; icol++ ) {
						sum += A[icol] * x.g((JA[icol])*stride,ic);
					}
					z(ir*stride,ic) += alpha * sum;
				}
			}
			
			// Compute alpha * A^T * x
			else 
			{
				if (this->Rows() != x.Rows()*stride || this->Cols() != y.Rows()*stride)
				{
					cout << "\nERROR! in TPZFYsmpMatrix::MultiplyAddMT: incompatible dimensions in opt=true\n";
					return; 
				}
				long jc;
				long icol;
				for(ir=0; ir<this->Rows(); ir++) {
					for(icol=IA[ir]; icol<IA[ir+1]; icol++ ) {
						if(JA[icol]==-1) break; //Checa a exist�cia de dado ou n�
						jc = JA[icol];
						TVar aval = A[icol];
						//cout << "FA["<<icol<<"] = "<<aval<< " * x["<< ir<<"] ="<< x.Get(ir,ic)<< endl;
						z(jc*stride,ic) += alpha * aval * x.g(ir*stride,ic);
					}
				}
			}
		}
	}
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
void TPZFYsmpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
								   TPZFMatrix<TVar> &z,
								   const TVar alpha,const TVar beta,const int opt,const int stride ) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
	{
		cout << "\nERROR! em TPZFYsmpMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
		return;
	}
	
	long  ic, xcols;
	xcols = x.Cols();
	long  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	for(ic=0; ic<xcols; ic++) {
		TVar *zp = &(z(0,ic));
		if(beta != 0) {
			const TVar *yp = &(y.g(0,0));
			TVar *zlast = zp+r*stride;
			if(beta != 1. || (&z != &y && stride != 1)) {
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp += stride;
					yp += stride;
				}
			}
			else if(&z != &y) {
				memcpy(zp,yp,r*sizeof(REAL));
			}
		} else {
			TVar *zp = &(z(0,0)), *zlast = zp+r*stride;
			while(zp != zlast) {
				*zp = 0.;
				zp += stride;
			}
		}
	}
	/*
	 TPZFYsmpMatrix *target;
	 int fFirsteq;
	 int fLasteq;
	 TPZFMatrix<>*fX;
	 TPZFMatrix<>*fZ;
	 REAL fAlpha;
	 int fOpt;
	 int fStride;
	 */
    const int numthreads = 2;
	pthread_t allthreads[numthreads];
	TPZMThread alldata[numthreads];
	int res[numthreads];
	int i;
	int eqperthread = r/numthreads;
	int firsteq = 0;
	for(i=0;i<numthreads;i++) 
	{
		alldata[i].target = this;
		alldata[i].fFirsteq = firsteq;
		alldata[i].fLasteq = firsteq+eqperthread;
		firsteq += eqperthread;
		if(i==numthreads-1) alldata[i].fLasteq = this->Rows();
		alldata[i].fX = &x;
		alldata[i].fZ = &z;
		alldata[i].fAlpha = alpha;
		alldata[i].fOpt = opt;
		alldata[i].fStride = stride;
		res[i] = PZ_PTHREAD_CREATE(&allthreads[i], NULL, 
					   ExecuteMT, &alldata[i], __FUNCTION__);
	}
	for(i=0;i<numthreads;i++) {
	  PZ_PTHREAD_JOIN(allthreads[i], NULL, __FUNCTION__);
	}
	
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZFYsmpMatrix<TVar>::Print(const char *title, ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	if(form != EInputFormat) {
		out << "\nTFYsmpMatrix Print: " << title << '\n'
		<< "\tRows    = " << this->Rows()  << '\n'
		<< "\tColumns = " << this->Cols() << '\n';
		long i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<this->Rows(); i++) {
			out << i      << '\t'
			<< IA[i] << '\t'
			<< JA[i] << '\t'
			<< A[i]  << '\n';
		}
		for(i=this->Rows()+1; i<IA[this->Rows()]; i++) {
			out << i      << "\t\t"
			<< JA[i] << '\t'
			<< A[i]  << '\n';
		}
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZFYsmpMatrix<TVar>::ComputeDiagonal() {
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;

	if(Diag.size()) return;
	long rows = this->Rows();
	Diag.resize(rows);
	for(long ir=0; ir<rows; ir++) {
		Diag[ir] = GetVal(ir,ir);
	}
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::SolveSOR( long &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
							  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &/*scratch*/,
							  const REAL overrelax, REAL &tol,
							  const int FromCurrent,const int direction ) {
	
	//    if(!fDiag) ComputeDiagonal();
	long irStart = 0,irLast = this->Rows(),irInc = 1;
    TPZVec<long> &IA = fIA;
    TPZVec<long> &JA = fJA;
    TPZVec<TVar> &A = fA;
    TPZVec<TVar> &Diag = fDiag;

	if(direction < 0) {
		irStart = this->Rows()-1;
		irLast = -1;
		irInc = -1;
	}
	if(!FromCurrent) x.Zero();
	TVar eqres = 2.*tol;
	long iteration;
	for(iteration=0; iteration<numiterations && eqres >= tol; iteration++) {
		eqres = 0.;
		long ir=irStart;
		while(ir != irLast) {
			TVar xnewval=rhs.g(ir,0);
			for(long ic=IA[ir]; ic<IA[ir+1]; ic++) {
				xnewval -= A[ic] * x(JA[ic],0);
			}
			eqres += xnewval*xnewval;
			x(ir,0) += overrelax*(xnewval/Diag[ir]);
			ir += irInc;
		}
		eqres = sqrt(eqres);
	}
	tol = eqres;
	numiterations = iteration;
	if(residual) this->Residual(x,rhs,*residual);
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::Zero()
{
    fA->Fill(TVar(0.));
    fDiag->Fill(TVar(0.));
	return 1;
}


/**
 * Solves the linear system using Jacobi method. \n
 * @param numiterations The number of interations for the process.
 * @param F The right hand side of the system.
 * @param result The solution.
 * @param residual Returns F - A*U which is the solution residual.
 * @param scratch Available manipulation area on memory.
 * @param tol The tolerance value.
 * @param FromCurrent It starts the solution based on FromCurrent. Obtaining solution FromCurrent + 1.
 */
template<class TVar>
void TPZFYsmpMatrix<TVar>::SolveJacobi(long & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result, TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent)
{
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	if(!Diag.size()) {
		cout << "TPZSYsmpMatrix::Jacobi cannot be called without diagonal\n";
		numiterations = 0;
		if(residual) {
			this->Residual(result,F,*residual);
			tol = sqrt(Norm(*residual));
		}
		return;
	}
	long c = F.Cols();
	long r = this->Rows();
	long it=0;
	if(FromCurrent) {
		this->Residual(result,F,scratch);
		for(long ic=0; ic<c; ic++) {
			for(long i=0; i<r; i++) {
				result(i,ic) += scratch(i,ic)/(Diag)[i];
			}
		}
	} else 
	{
		for(long ic=0; ic<c; ic++) {
			for(long i=0; i<r; i++) {
				result(i,ic) = F.GetVal(i,ic)/(Diag)[i];
			}
		}
	}
	if(it<numiterations)
	{
		this->Residual(result,F,scratch);
		TVar res = Norm(scratch);
		for(long it=1; it<numiterations && res > tol; it++) {
			for(long ic=0; ic<c; ic++) {
				for(long i=0; i<r; i++) {
					result(i,ic) += (scratch)(i,ic)/(Diag)[i];
				}
			}
			this->Residual(result,F,scratch);
			res = Norm(scratch);
		}
	}
	if(residual) *residual = scratch;
}

template<class TVar>
void *TPZFYsmpMatrix<TVar>::ExecuteMT(void *entrydata)
{
    //DebugStop();
	TPZMThread *data = (TPZMThread *) entrydata;
	const TPZFYsmpMatrix *mat = data->target;
	TVar sum;
	long xcols = data->fX->Cols();
    long ic,ir,icol;
    TPZVec<long> &IA = *(mat->fIA.operator->());
    TPZVec<long> &JA = *(mat->fJA.operator->());
    TPZVec<TVar> &A = *(mat->fA.operator->());

	// Compute alpha * A * x
	if(xcols == 1 && data->fStride == 1 && data->fOpt == 0)
	{
		for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
			long icolmin = IA[ir];
			long icolmax = IA[ir+1];
			const TVar *xptr = &(data->fX->g(0,0));
			TVar *Aptr = &A[0];
			long *JAptr = &JA[0];
			for(sum = 0.0, icol=icolmin; icol<icolmax; icol++ ) {
				sum += Aptr[icol] * xptr[JAptr[icol]];
			}
			data->fZ->operator()(ir,0) += data->fAlpha * sum;
		}
	}
	else 
	{
		for(ic=0; ic<xcols; ic++) {
			if(data->fOpt == 0) {
				
				for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
					for(sum = 0.0, icol=IA[ir]; icol<IA[ir+1]; icol++ ) {
						sum += A[icol] * data->fX->g((JA[icol])*data->fStride,ic);
					}
					data->fZ->operator()(ir*data->fStride,ic) += data->fAlpha * sum;
				}
			}
			
			// Compute alpha * A^T * x
			else 
			{
				long jc;
				long icol;
				for(ir=data->fFirsteq; ir<data->fLasteq; ir++) 
				{
					for(icol=IA[ir]; icol<IA[ir+1]; icol++ )
					{
						if(JA[icol]==-1) break; //Checa a exist�cia de dado ou n�
						jc = JA[icol];
						data->fZ->operator()(jc*data->fStride,ic) += data->fAlpha * A[icol] * data->fX->g(ir,ic);
					}
				}
				
			}
		}
	}
	return 0;
}
static int  FindCol(long *colf, long *coll, long col)
{
//    if(col == *colf && col == *coll) {
//        return 0;
//    }
	if(col == *colf || col == *coll) return col;
	long *begin = colf;
	long *end = coll;
	while (begin != end)
	{
		long dist = (end-begin)/2;
		long *mid = begin+dist;
		if(*mid == col) return (mid-colf);
		else if(*mid > col) end=mid;
		else begin = mid;
	}
	return -1;
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::GetSub(const long sRow,const long sCol,const long rowSize,
								 const long colSize, TPZFMatrix<TVar> & A ) const {
	long ir;
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &Av = fA;
    const TPZVec<TVar> &Diag = fDiag;

	for(ir=0; ir<rowSize; ir++)
	{
		long row = sRow+ir;
		long colfirst = IA[row];
		long collast = IA[row+1];
		long iacol = FindCol(&JA[0]+colfirst,&JA[0]+collast-1,sCol);
		long ic;
		for(ic=0; ic<colSize; ic++) A(ir,ic) = Av[iacol+colfirst];
	}
	return 0;
}


template<class TVar>
void TPZFYsmpMatrix<TVar>::GetSub(const TPZVec<long> &indices,TPZFMatrix<TVar> &block) const
{
	std::map<long,long> indord;
	long i,size = indices.NElements();
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	for(i=0; i<size; i++)
	{
		indord[indices[i]] = i;
	}
	std::map<long,long>::iterator itset,jtset;
	for(itset = indord.begin(); itset != indord.end(); itset++)
	{
		long *jfirst = &JA[0]+IA[(*itset).first];
		long *jlast = &JA[0]+IA[(*itset).first+1]-1;
		//    int row = (*itset).first;
		for(jtset = indord.begin(); jtset != indord.end(); jtset++)
		{
			long col = FindCol(jfirst,jlast,(*jtset).first);
			long dist = jfirst+col-&JA[0];
			block((*itset).second,(*jtset).second) = A[dist];
			jfirst += col+1;
		}
	}
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::GetSub(const TPZVec<long> &i_indices, const TPZVec<long> &j_indices,TPZFMatrix<TVar> &block) const
{
    std::map<long,long> indord,jndord;
    long i,j;
    long i_size = i_indices.NElements();
    long j_size = j_indices.NElements();
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;
    
    for(i=0; i<i_size; i++)
    {
        indord[i_indices[i]] = i;
    }
    for(j=0; j<j_size; j++)
    {
        jndord[j_indices[j]] = j;
    }
    
    std::map<long,long>::iterator itset,jtset;
    for(itset = indord.begin(); itset != indord.end(); itset++)
    {
        long *jfirst = &JA[0]+IA[(*itset).first];
        long *jlast = &JA[0]+IA[(*itset).first+1]-1;
        //    int row = (*itset).first;
        for(jtset = jndord.begin(); jtset != jndord.end(); jtset++)
        {
            int last = (*itset).first;

            long col = FindCol(jfirst,jlast,(*jtset).first);
            long dist = jfirst+col-&JA[0];
            block((*itset).second,(*jtset).second) = A[dist];
            jfirst += col+1;
        }
    }
}

/*
 * Perform row update of the sparse matrix
 */
template<class TVar>
void TPZFYsmpMatrix<TVar>::RowLUUpdate(long sourcerow, long destrow)
{
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	long *sourcefirst = &JA[0]+IA[sourcerow];
	long *sourcelast = &JA[0]+IA[sourcerow+1]-1;
	long sourcecol = FindCol(sourcefirst,sourcelast,sourcerow);
	if(sourcecol < 0)
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " source not found\n";
		return;
	}
	long sourcedist = sourcefirst+sourcecol-&JA[0];
	long *destfirst = &JA[0]+IA[destrow];
	long *destlast = &JA[0]+IA[destrow+1]-1;
	long destcol = FindCol(destfirst,destlast,destrow);
	long destdist = destfirst+destcol-&JA[0];
	if(destcol < 0)
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " destrow not found\n";
		return;
	}
	if(A[sourcedist] < 1.e-15)
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " small pivot " << A[sourcedist] << "\n";
		return;
	}
	TVar mult = A[destdist]/A[sourcedist];
	if(mult == 0.) return;
	destdist++;
	sourcedist++;
	while(destdist < IA[destrow+1] && sourcedist < IA[sourcerow+1])
	{
		if(JA[destdist] == JA[sourcedist])
		{
			A[destdist] -= A[sourcedist]*mult;
			destdist++;
			sourcedist++;
		}
		else if(JA[destdist] < JA[sourcedist])
		{
			destdist++;
		}
		else
		{
			sourcedist++;
		}
	}
	
}

/**
 * Decomposes the current matrix using LU decomposition.
 */
template<class TVar>
int TPZFYsmpMatrix<TVar>::Decompose_LU(std::list<long> &singular)
{
	return Decompose_LU();
}
template<class TVar>
int TPZFYsmpMatrix<TVar>::Decompose_LU()
{
	long row;
	long neq = this->Rows();
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	for(row=1; row<neq; row++)
	{
		//    int firstcol = fIA[row];
		long lastcol = IA[row+1];
		long colind = 0;
		if(JA[lastcol-1] < row) continue;
		while(JA[colind] < row)
		{
			RowLUUpdate(JA[colind],row);
			colind++;
		}
	}
	this->fDecomposed=1;
	return 1;
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
	long row;
	long bcol = B->Cols();
	long col;
	long neq = this->Rows();
    const TPZVec<long> &IA = fIA;
    const TPZVec<long> &JA = fJA;
    const TPZVec<TVar> &A = fA;
    const TPZVec<TVar> &Diag = fDiag;

	
	// forward substitution
	for(row=0; row<neq; row++)
	{
		long firstrow = IA[row];
		long lastrow = IA[row+1];
		if(JA[firstrow] > row || JA[lastrow-1] < row)
		{
			cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " inconsistent column information for row " << row << endl;
			continue;
		}
		long rowcounter = firstrow;
		while(JA[rowcounter] < row)
		{
			for(col=0; col<bcol; col++)
			{
				(*B)(row,col) -= A[rowcounter]*(*B)(JA[rowcounter],col);
			}
		}
		for(col=0; col<bcol; col++)
		{
			(*B)(row,col) /= A[rowcounter];
		}
	}
	// backward substitution
	for(row = neq-1; row >= 0; row--)
	{
		long firstrow = IA[row];
		long lastrow = IA[row+1];
		long col = FindCol(&JA[0]+firstrow,&JA[0]+lastrow-1,row);
		if(col < 0)
		{
			cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " inconsistent column information for row " << row << endl;
			continue;
		}
		long coldist = firstrow+col+1;
		while(coldist < lastrow)
		{
			for(col=0; col<bcol; col++)
			{
				(*B)(row,col) -= A[coldist]*(*B)(JA[coldist],col);
			}
		}
	}
	return 1;
}

template class TPZFYsmpMatrix<long double>;
template class TPZFYsmpMatrix<double>;
template class TPZFYsmpMatrix<float>;
template class TPZFYsmpMatrix<int>;
