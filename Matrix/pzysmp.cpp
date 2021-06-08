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
#include <thread>
#include <vector>
#include "tpzverysparsematrix.h"
#include "pzstack.h"

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
    int64_t i,j,k;
    if (B.Rows()!=this->Rows()) return;
    int64_t rows = this->Rows();
    TVar aux=0.;
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
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix(const TPZVerySparseMatrix<TVar> &cp) : TPZRegisterClassId(&TPZFYsmpMatrix::ClassId),
TPZMatrix<TVar>
()
{
	*this = cp;
}

template<class TVar>
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix() : TPZRegisterClassId(&TPZFYsmpMatrix::ClassId),
TPZMatrix<TVar>(), fIA(1,0),fJA(),fA(),fDiag()
{
}

template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator=(const TPZVerySparseMatrix<TVar> &cp)
{
	// Deletes everything associated with a TPZFYsmpMatrix
	int64_t nrows = cp.Rows();
	
	int64_t count = 0, c = 0, r = 0;
	
	count = cp.fExtraSparseData.size();
	fJA.resize(count);
    fA.resize(count);
    fIA.resize(nrows+1);
	fIA[0] = 0;
	
	typename map< pair<int64_t,int64_t>, TVar>::const_iterator it;
	c = 0;
	r = 0;
	for(it=cp.fExtraSparseData.begin(); it!= cp.fExtraSparseData.end(); it++)
	{
		int64_t row = it->first.first;
		if(r != row)
		{
			r++;
			while(r < row) 
			{
				fIA[r] = c;
				r++;
			}
			fIA[row] = c;
			r = row;
		}
		int64_t col = it->first.second;
		fJA[c] = col;
		TVar val = it->second;
		fA[c] = val;
		c++;
	}
	r++;
	while(r<=nrows)
	{
		fIA[r] = c;
		r++;
	}
	return *this;
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::PutVal(const int64_t row, const int64_t col, const TVar &Value){
    int64_t k;
    int flag=0;
    for(k=fIA[row];k<fIA[row+1];k++){
		if(fJA[k]==col){
			flag=1;
			fA[k]=Value;
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
void TPZFYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex){
    int64_t i,j,k = 0;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            value=elmat.GetVal(i,j);
            //cout << "j= " << j << endl;
            if(!IsZero(value)){
                //cout << "fIA[ipos] " << fIA[ipos] << "     fIA[ipos+1] " << fIA[ipos+1] << endl;
                int flag = 0;
				k++;
				if(k >= fIA[ipos] && k < fIA[ipos+1] && fJA[k]==jpos)
				{ // OK -> elements in sequence
					fA[k]+=value;
					flag = 1;
				}else
				{
					for(k=fIA[ipos];k<fIA[ipos+1];k++){
						if(fJA[k]==jpos || fJA[k]==-1){
							//cout << "fJA[k] " << fJA[k] << " jpos "<< jpos << "   " << value << endl;
							//cout << "k " << k << "   "<< jpos << "   " << value << endl;
							flag=1;
							if(fJA[k]==-1){
								fJA[k]=jpos;
								fA[k]=value;
								// cout << jpos << "   " << value << endl;
								break;
							}else{
								fA[k]+=value;
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
void TPZFYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
	int64_t i,j,k = 0;
	TVar value=0.;
	int64_t ipos,jpos;
	for(i=0;i<sourceindex.NElements();i++){
		for(j=0;j<sourceindex.NElements();j++){
			ipos=destinationindex[i];
			jpos=destinationindex[j];
			value=elmat.GetVal(sourceindex[i],sourceindex[j]);
            //cout << "j= " << j << endl;
			if(!IsZero(value)){
                //cout << "fIA[ipos] " << fIA[ipos] << "     fIA[ipos+1] " << fIA[ipos+1] << endl;
				int flag = 0;
				k++;
				if(k >= fIA[ipos] && k < fIA[ipos+1] && fJA[k]==jpos)
				{ // OK -> elements in sequence
					fA[k]+=value;
					flag = 1;
				}else
				{
					for(k=fIA[ipos];k<fIA[ipos+1];k++){
						if(fJA[k]==jpos || fJA[k]==-1){
							//cout << "fJA[k] " << fJA[k] << " jpos "<< jpos << "   " << value << endl;
							//cout << "k " << k << "   "<< jpos << "   " << value << endl;
							flag=1;
							if(fJA[k]==-1){
								fJA[k]=jpos;
								fA[k]=value;
								// cout << jpos << "   " << value << endl;
								break;
							}else{
								fA[k]+=value;
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
	int64_t i=0;
	int64_t j=0;
	int64_t ilocal=0;
	//  int jlocal=0;
	int64_t nel = destinationindex.NElements();
	std::multimap<int64_t,int64_t> mapindex;
	std::multimap<int64_t,int64_t>::iterator hint = mapindex.begin();
	for(i=0;i<nel;i++){
		ilocal = destinationindex[i];
		hint = mapindex.insert(hint,std::make_pair(ilocal,i));
		//    mapindex[ilocal] = i;
	}
	for(i=0;i<nel;i++){
		ilocal = destinationindex[i];
		int64_t jfirst = fIA[ilocal];
		int64_t jlast = fIA[ilocal+1];
		int64_t *Aptr = &fJA[jfirst];
		int64_t *AptrLast = &fJA[jlast];
		j=0;
		std::multimap<int64_t,int64_t>::iterator itelmat = mapindex.begin();
		while(j<nel && Aptr != AptrLast)
		{
			if(*Aptr == (*itelmat).first)
			{
				int64_t jel = (*itelmat).second;
				fA[jfirst] += elmat(i,jel);
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
				int64_t *iptr = &fJA[jfirst];
				while(iptr < AptrLast) 
				{
					cout << *iptr << " ";
					iptr++;
				}
				cout << endl;
				std::multimap<int64_t,int64_t>::iterator itelmat2 = mapindex.begin();
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
			int64_t *iptr = &fJA[jfirst];
			while(iptr < AptrLast) 
			{
				cout << *iptr << " ";
				iptr++;
			}
			cout << endl;
			std::multimap<int64_t,int64_t>::iterator itelmat2 = mapindex.begin();
			for(;itelmat2 != mapindex.end(); itelmat2++)
			{
				cout << (*itelmat2).first << "/" << (*itelmat2).second << " ";
			}
			cout << endl;
		}
	}
	
}

template<class TVar>
TPZFYsmpMatrix<TVar>::TPZFYsmpMatrix(const int64_t rows,const int64_t cols ) :
TPZRegisterClassId(&TPZFYsmpMatrix::ClassId),TPZMatrix<TVar>(rows,cols) {
	// Constructs an empty TPZFYsmpMatrix
	//    fSolver = -1;
	fSymmetric = 0;
	//    fMaxIterations = 4;
	//    fSORRelaxation = 1.;
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
const TVar TPZFYsmpMatrix<TVar>::GetVal(const int64_t row,const int64_t col ) const {
	// Get the matrix entry at (row,col) without bound checking
	
	// Now look through the requested row and see if there is anything
	// in column col
	/*  int loccol = col+1;
	 for(int ic=fIA[row]-1 ; ic < fIA[row+1]-1; ic++ ) {
	 if ( fJA[ic] == loccol ) return fA[ic];
	 }*/
	int64_t loccol = col;
	for(int64_t ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
		if ( fJA[ic] == loccol && fJA[ic] != -1 ){
			return fA[ic];
		}
	}
	return (TVar) 0;
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::CheckTypeCompatibility(const TPZMatrix<TVar>*A,
																									const TPZMatrix<TVar>*B)const
{
  auto incompatSparse = [](){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: incompatible matrices\n.Aborting...\n";
    DebugStop();
  };
	auto aPtr = dynamic_cast<const TPZFYsmpMatrix<TVar>*>(A);
  auto bPtr = dynamic_cast<const TPZFYsmpMatrix<TVar>*>(B);
  if(!aPtr || !bPtr){
    incompatSparse();
  }
	bool check{false};
	const auto nIA = aPtr->fIA.size();
	for(auto i = 0; i < nIA; i++){
		check = check || aPtr->fIA[i] != bPtr->fIA[i];
	}

	const auto nJA = aPtr->fJA.size();
	for(auto i = 0; i < nJA; i++){
		check = check || aPtr->fJA[i] != bPtr->fJA[i];
	}
	if(check) incompatSparse();
}

// ****************************************************************************
//
// Multiply and Multiply-Add
//
// ****************************************************************************

template<class TVar>
TPZFYsmpMatrix<TVar> TPZFYsmpMatrix<TVar>::operator+(const TPZFYsmpMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
  const auto sizeA = res.fA.size();
  for(auto i = 0; i < sizeA; i++) res.fA[i] += mat.fA[i];
	return res;
}
template<class TVar>
TPZFYsmpMatrix<TVar> TPZFYsmpMatrix<TVar>::operator-(const TPZFYsmpMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
  const auto sizeA = res.fA.size();
  for(auto i = 0; i < sizeA; i++) res.fA[i] -= mat.fA[i];
	return res;
}

template<class TVar>
TPZFYsmpMatrix<TVar> TPZFYsmpMatrix<TVar>::operator*(const TVar alpha) const
{
	auto res(*this);
	for(auto &el : res.fA) el *= alpha;
	return res;
}

template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator+=(const TPZFYsmpMatrix<TVar> &A )
{
	TPZFYsmpMatrix<TVar> res((*this)+A);
	*this = res;
	return *this;
}
template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator-=(const TPZFYsmpMatrix<TVar> &A )
{
	TPZFYsmpMatrix<TVar> res((*this)-A);
	*this = res;
	return *this;
}
template<class TVar>
TPZFYsmpMatrix<TVar> &TPZFYsmpMatrix<TVar>::operator*=(const TVar val)
{
	TPZFYsmpMatrix<TVar> res((*this)*val);
	*this = res;
	return *this;
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::MultAddMT(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
									 TPZFMatrix<TVar> &z,
									 const TVar alpha,const TVar beta,const int opt )  {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	if(x.Cols() != y.Cols() || x.Cols() != z.Cols() || y.Rows() != z.Rows() )
	{
		cout << "\nERROR! in TPZVerySparseMatrix::MultiplyAdd : incompatible dimensions in x, y or z\n";
		return;
	}
	
	int64_t  ir, ic, icol, xcols;
	xcols = x.Cols();
	TVar sum;
	int64_t  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	for(ic=0; ic<xcols; ic++) {
		TVar *zp = &(z(0,ic));
		if(!IsZero(beta)){
			const TVar *yp = &(y.g(0,0));
			TVar *zlast = zp+r;

            if(&z != &y) {
				memcpy(zp,yp,r*sizeof(TVar));
			}
		} else {
			TVar *zp = &(z(0,0)), *zlast = zp+r;
			while(zp != zlast) {
				*zp = 0.;
				zp ++;
			}
		}
	}
	// Compute alpha * A * x
	if(xcols == 1 && opt == 0)
	{
		if(this->Cols() != x.Rows() || this->Rows() != y.Rows())
		{
			cout << "\nERROR! in TPZFYsmpMatrix::MultiplyAddMT: incompatible dimensions in opt=false\n";
			return;
		} 
		for(ir=0; ir<r; ir++) {
			int64_t icolmin = fIA[ir];
			int64_t icolmax = fIA[ir+1];
			const TVar *xptr = &(x.g(0,0));
			TVar *Aptr = &fA[0];
			int64_t *JAptr = &fJA[0];
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
					for(sum = 0.0, icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
						sum += fA[icol] * x.g((fJA[icol]),ic);
					}
					z(ir,ic) += alpha * sum;
				}
			}
			
			// Compute alpha * A^T * x
			else 
			{
				if (this->Rows() != x.Rows() || this->Cols() != y.Rows())
				{
					cout << "\nERROR! in TPZFYsmpMatrix::MultiplyAddMT: incompatible dimensions in opt=true\n";
					return; 
				}
				int64_t jc;
				int64_t icol;
				for(ir=0; ir<this->Rows(); ir++) {
					for(icol=fIA[ir]; icol<fIA[ir+1]; icol++ ) {
						if(fJA[icol]==-1) break; //Checa a exist�cia de dado ou n�
						jc = fJA[icol];
						TVar aval = fA[icol];
						//cout << "FA["<<icol<<"] = "<<aval<< " * x["<< ir<<"] ="<< x.Get(ir,ic)<< endl;
						z(jc,ic) += alpha * aval * x.g(ir,ic);
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
								   const TVar alpha,const TVar beta,const int opt) const {
	// computes z = beta * y + alpha * opt(this)*x
	//          z and x cannot share storage
	
#ifdef PZDEBUG
    if ((!opt && this->Cols() != x.Rows()) || (opt && this->Rows() != x.Rows())) {
        std::cout << "TPZFMatrix::MultAdd matrix x with incompatible dimensions>" ;
        return;
    }
    if(!IsZero(beta) && ((!opt && this->Rows() != y.Rows()) || (opt && this->Cols() != y.Rows()) || y.Cols() != x.Cols())) {
        std::cout << "TPZFMatrix::MultAdd matrix y with incompatible dimensions>";
        return;
    }
#endif
    if(!opt) {
        if(z.Cols() != x.Cols() || z.Rows() != this->Rows()) {
            z.Redim(this->Rows(),x.Cols());
        }
    } else {
        if(z.Cols() != x.Cols() || z.Rows() != this->Cols()) {
            z.Redim(this->Cols(),x.Cols());
        }
    }
    if(this->Cols() == 0) {
        z.Zero();
        return;
    }

    int64_t  ic, xcols;
	xcols = x.Cols();
	int64_t  r = (opt) ? this->Cols() : this->Rows();
	
	// Determine how to initialize z
	for(ic=0; ic<xcols; ic++) {
		TVar *zp = &(z(0,ic));
		if(!IsZero(beta)){
			const TVar *yp = &(y.g(0,0));
			TVar *zlast = zp+r;
			if(!IsZero(beta-(TVar)1.)){
				while(zp < zlast) {
					*zp = beta * (*yp);
					zp ++;
					yp ++;
				}
			}
			else if(&z != &y) {
				memcpy(zp,yp,r*sizeof(REAL));
			}
		} else {
			TVar *zp = &(z(0,ic)), *zlast = zp+r;
			while(zp != zlast) {
				*zp = 0.;
				zp ++;
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
  */
  int numthreads = 2;
  if(opt) numthreads = 1;
  std::vector<std::thread> allthreads;
	TPZVec<TPZMThread> alldata(numthreads);
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
    allthreads.push_back(std::thread(ExecuteMT, &alldata[i]));
	}
	for(i=0;i<numthreads;i++) {
    allthreads[i].join();
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
	if(form == EInputFormat) {
		out << "\nTFYsmpMatrix Print: " << title << '\n'
		<< "\tRows    = " << this->Rows()  << '\n'
		<< "\tColumns = " << this->Cols() << '\n';
		int64_t i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<this->Rows(); i++) {
			out << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(i=this->Rows()+1; i<fIA[this->Rows()]; i++) {
			out << i      << "\t\t"
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
	} else {
		TPZMatrix<TVar>::Print(title,out,form);
	}
}


// ****************************************************************************
//
// Various solvers
//
// ****************************************************************************

template<class TVar>
void TPZFYsmpMatrix<TVar>::ComputeDiagonal() {
	if(fDiag.size()) return;
	int64_t rows = this->Rows();
	fDiag.resize(rows);
	for(int64_t ir=0; ir<rows; ir++) {
		fDiag[ir] = GetVal(ir,ir);
	}
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::SolveSOR( int64_t &numiterations, const TPZFMatrix<TVar> &rhs, TPZFMatrix<TVar> &x,
							  TPZFMatrix<TVar> *residual, TPZFMatrix<TVar> &/*scratch*/,
							  const REAL overrelax, REAL &tol,
							  const int FromCurrent,const int direction ) {
	
	//    if(!fDiag) ComputeDiagonal();
	int64_t irStart = 0,irLast = this->Rows(),irInc = 1;
	if(direction < 0) {
		irStart = this->Rows()-1;
		irLast = -1;
		irInc = -1;
	}
	if(!FromCurrent) x.Zero();
	RTVar eqres = 2.*tol;
	int64_t iteration;
	for(iteration=0; iteration<numiterations && eqres >= tol; iteration++) {
		TVar local_eqres = 0.;
		int64_t ir=irStart;
		while(ir != irLast) {
			TVar xnewval=rhs.g(ir,0);
			for(int64_t ic=fIA[ir]; ic<fIA[ir+1]; ic++) {
				xnewval -= fA[ic] * x(fJA[ic],0);
			}
			local_eqres += xnewval*xnewval;
			x(ir,0) += (TVar)overrelax*(xnewval/fDiag[ir]);
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
    fA.Fill(TVar(0.));
    fDiag.Fill(TVar(0.));
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
void TPZFYsmpMatrix<TVar>::SolveJacobi(int64_t & numiterations, const TPZFMatrix<TVar> & F, TPZFMatrix<TVar> & result, TPZFMatrix<TVar> * residual, TPZFMatrix<TVar> & scratch, REAL & tol, const int FromCurrent)
{
	if(!fDiag.size()) {
		cout << "TPZSYsmpMatrix::Jacobi cannot be called without diagonal\n";
		numiterations = 0;
		if(residual) {
			this->Residual(result,F,*residual);
			tol = sqrt(Norm(*residual));
		}
		return;
	}
	int64_t c = F.Cols();
	int64_t r = this->Rows();
	int64_t it=0;
	if(FromCurrent) {
		this->Residual(result,F,scratch);
		for(int64_t ic=0; ic<c; ic++) {
			for(int64_t i=0; i<r; i++) {
				result(i,ic) += scratch(i,ic)/(fDiag)[i];
			}
		}
	} else 
	{
		for(int64_t ic=0; ic<c; ic++) {
			for(int64_t i=0; i<r; i++) {
				result(i,ic) = F.GetVal(i,ic)/(fDiag)[i];
			}
		}
	}
	if(it<numiterations)
	{
		this->Residual(result,F,scratch);
		TVar res = Norm(scratch);
		bool cond = [res,tol](){
			if constexpr (is_complex<TVar>::value){
				return abs(res)>tol;
			}else{
				return res>tol;
			}
		}();
	
		for(int64_t it=1; it<numiterations && cond; it++) {
			for(int64_t ic=0; ic<c; ic++) {
				for(int64_t i=0; i<r; i++) {
					result(i,ic) += (scratch)(i,ic)/(fDiag)[i];
				}
			}
			this->Residual(result,F,scratch);
			res = Norm(scratch);
			if constexpr (is_complex<TVar>::value) {
				cond = abs(res) > tol;
			} else {
				cond = res > tol;
			}
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
	int64_t xcols = data->fX->Cols();
	int64_t ic,ir,icol;
	// Compute alpha * A * x
	if(xcols == 1 && data->fOpt == 0)
	{
		for(ir=data->fFirsteq; ir<data->fLasteq; ir++) {
			int64_t icolmin = mat->fIA[ir];
			int64_t icolmax = mat->fIA[ir+1];
			const TVar *xptr = &(data->fX->g(0,0));
			TVar *Aptr = &mat->fA[0];
			int64_t *JAptr = &mat->fJA[0];
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
					for(sum = 0.0, icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ ) {
						sum += mat->fA[icol] * data->fX->g((mat->fJA[icol]),ic);
					}
					data->fZ->operator()(ir,ic) += data->fAlpha * sum;
				}
			}
			
			// Compute alpha * A^T * x
			else 
			{
				int64_t jc;
				int64_t icol;
				for(ir=data->fFirsteq; ir<data->fLasteq; ir++) 
				{
					for(icol=mat->fIA[ir]; icol<mat->fIA[ir+1]; icol++ ) 
					{
						if(mat->fJA[icol]==-1) break; //Checa a exist�cia de dado ou n�
						jc = mat->fJA[icol];
						data->fZ->operator()(jc,ic) += data->fAlpha * mat->fA[icol] * data->fX->g(ir,ic);
					}
				}
				
			}
		}
	}
	return 0;
}
static int  FindCol(int64_t *colf, int64_t *coll, int64_t col)
{
	if(col == *colf) return 0;
	int64_t *begin = colf;
	int64_t *end = coll;
	while (begin != end)
	{
		int64_t dist = (end-begin)/2;
		int64_t *mid = begin+dist;
		if(*mid == col) return (mid-colf);
		else if(*mid > col) end=mid;
		else begin = mid;
	}
	return -1;
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::GetSub(const int64_t sRow,const int64_t sCol,const int64_t rowSize,
								 const int64_t colSize, TPZFMatrix<TVar> & A ) const {
	int64_t ir;
	for(ir=0; ir<rowSize; ir++)
	{
		int64_t row = sRow+ir;
		int64_t colfirst = fIA[row];
		int64_t collast = fIA[row+1];
		int64_t iacol = FindCol(&fJA[0]+colfirst,&fJA[0]+collast-1,sCol);
		int64_t ic;
		for(ic=0; ic<colSize; ic++) A(ir,ic) = fA[iacol+colfirst];
	}
	return 0;
}


template<class TVar>
void TPZFYsmpMatrix<TVar>::GetSub(const TPZVec<int64_t> &indices,TPZFMatrix<TVar> &block) const
{
	std::map<int64_t,int64_t> indord;
	int64_t i,size = indices.NElements();
	for(i=0; i<size; i++)
	{
		indord[indices[i]] = i;
	}
	std::map<int64_t,int64_t>::iterator itset,jtset;
	for(itset = indord.begin(); itset != indord.end(); itset++)
	{
		int64_t *jfirst = &fJA[0]+fIA[(*itset).first];
		int64_t *jlast = &fJA[0]+fIA[(*itset).first+1]-1;
		//    int row = (*itset).first;
		for(jtset = indord.begin(); jtset != indord.end(); jtset++)
		{
			int64_t col = FindCol(jfirst,jlast,(*jtset).first);
			int64_t dist = jfirst+col-&fJA[0];
			block((*itset).second,(*jtset).second) = fA[dist];
			jfirst += col+1;
		}
	}
}

/*
 * Perform row update of the sparse matrix
 */
template<class TVar>
void TPZFYsmpMatrix<TVar>::RowLUUpdate(int64_t sourcerow, int64_t destrow)
{
	int64_t *sourcefirst = &fJA[0]+fIA[sourcerow];
	int64_t *sourcelast = &fJA[0]+fIA[sourcerow+1]-1;
	int64_t sourcecol = FindCol(sourcefirst,sourcelast,sourcerow);
	if(sourcecol < 0)
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " source not found\n";
		return;
	}
	int64_t sourcedist = sourcefirst+sourcecol-&fJA[0];
	int64_t *destfirst = &fJA[0]+fIA[destrow];
	int64_t *destlast = &fJA[0]+fIA[destrow+1]-1;
	int64_t destcol = FindCol(destfirst,destlast,destrow);
	int64_t destdist = destfirst+destcol-&fJA[0];
	if(destcol < 0)
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " destrow not found\n";
		return;
	}

	if(IsZero(fA[sourcedist]))
	{
		cout << __PRETTY_FUNCTION__ << " at line " << __LINE__ << " small pivot " << fA[sourcedist] << "\n";
		return;
	}
	TVar mult = fA[destdist]/fA[sourcedist];
	if(IsZero(mult)) return;
	destdist++;
	sourcedist++;
	while(destdist < fIA[destrow+1] && sourcedist < fIA[sourcerow+1])
	{
		if(fJA[destdist] == fJA[sourcedist])
		{
			fA[destdist] -= fA[sourcedist]*mult;
			destdist++;
			sourcedist++;
		}
		else if(fJA[destdist] < fJA[sourcedist])
		{
			destdist++;
		}
		else
		{
			sourcedist++;
		}
	}
	
}
/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZFYsmpMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (symmetric && nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    orig.AutoFill(nrow,ncol,symmetric);
    
    TPZVec<int64_t> IA(nrow+1);
    TPZStack<int64_t> JA;
    TPZStack<TVar> A;
    IA[0] = 0;
    TPZVec<std::set<int64_t> > eqs(nrow);
    for (int64_t row=0; row<nrow; row++) {
        if(nrow == ncol) eqs[row].insert(row);
        for (int64_t col = 0; col<ncol; col++) {
            REAL test = rand()*1./RAND_MAX;
            if (test > 0.5) {
                eqs[row].insert(col);
                if (symmetric) {
                    eqs[col].insert(row);
                }
            }
        }
    }
    int64_t pos=0;
    for (int64_t row=0; row< nrow; row++) {
        for (std::set<int64_t>::iterator col = eqs[row].begin(); col != eqs[row].end(); col++) {
            JA.Push(*col);
            A.Push(orig(row,*col));
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}


/**
 * Decomposes the current matrix using LU decomposition.
 */
template<class TVar>
int TPZFYsmpMatrix<TVar>::Decompose_LU(std::list<int64_t> &singular)
{
	return Decompose_LU();
}
template<class TVar>
int TPZFYsmpMatrix<TVar>::Decompose_LU()
{
    
#ifdef USING_MKL
    if(this->IsDecomposed() == ELU) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
	typename TPZPardisoSolver<TVar>::MStructure str =
        this->IsSymmetric() ?
		TPZPardisoSolver<TVar>::MStructure::ESymmetric:
		TPZPardisoSolver<TVar>::MStructure::ENonSymmetric;
	typename TPZPardisoSolver<TVar>::MSystemType sysType =
		this->IsSymmetric() ?
		TPZPardisoSolver<TVar>::MSystemType::ESymmetric:
		TPZPardisoSolver<TVar>::MSystemType::ENonSymmetric;
	typename TPZPardisoSolver<TVar>::MProperty prop =
		this->IsDefPositive() ?
		TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
		TPZPardisoSolver<TVar>::MProperty::EIndefinite;
	fPardisoControl.SetStructure(str);
	fPardisoControl.SetMatrixType(sysType,prop);
    fPardisoControl.Decompose(this);
    this->SetIsDecomposed(ELU);
    return 1;
#endif
    
	int64_t row;
	int64_t neq = this->Rows();
	for(row=1; row<neq; row++)
	{
		//    int firstcol = fIA[row];
		int64_t lastcol = fIA[row+1];
		int64_t colind = 0;
		if(fJA[lastcol-1] < row) continue;
		while(fJA[colind] < row)
		{
			RowLUUpdate(fJA[colind],row);
			colind++;
		}
	}
	this->fDecomposed=1;
	return 1;
}

template<class TVar>
void TPZFYsmpMatrix<TVar>::SetIsDecomposed(int val)
{
	if(val)
		fPardisoControl.fDecomposed = true;
	TPZBaseMatrix::SetIsDecomposed(val);
}

template<class TVar>
int TPZFYsmpMatrix<TVar>::Substitution( TPZFMatrix<TVar> *B ) const
{
    
#ifdef USING_MKL
    TPZFMatrix<TVar> x(*B);
    fPardisoControl.Solve(this,*B,x);
    *B = x;
    return 1;
#endif
    
	int64_t row;
	int64_t bcol = B->Cols();
	int64_t col;
	int64_t neq = this->Rows();
	
	// forward substitution
	for(row=0; row<neq; row++)
	{
		int64_t firstrow = fIA[row];
		int64_t lastrow = fIA[row+1];
		if(fJA[firstrow] > row || fJA[lastrow-1] < row)
		{
			cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " inconsistent column information for row " << row << endl;
			continue;
		}
		int64_t rowcounter = firstrow;
		while(fJA[rowcounter] < row)
		{
			for(col=0; col<bcol; col++)
			{
				(*B)(row,col) -= fA[rowcounter]*(*B)(fJA[rowcounter],col);
			}
		}
		for(col=0; col<bcol; col++)
		{
			(*B)(row,col) /= fA[rowcounter];
		}
	}
	// backward substitution
	for(row = neq-1; row >= 0; row--)
	{
		int64_t firstrow = fIA[row];
		int64_t lastrow = fIA[row+1];
		int64_t col = FindCol(&fJA[0]+firstrow,&fJA[0]+lastrow-1,row);
		if(col < 0)
		{
			cout << __PRETTY_FUNCTION__ << " " << __LINE__ << " inconsistent column information for row " << row << endl;
			continue;
		}
		int64_t coldist = firstrow+col+1;
		while(coldist < lastrow)
		{
			for(col=0; col<bcol; col++)
			{
				(*B)(row,col) -= fA[coldist]*(*B)(fJA[coldist],col);
			}
		}
	}
	return 1;
}


template<class TVar>
int TPZFYsmpMatrix<TVar>::ClassId() const{
    return Hash("TPZFYsmpMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template class TPZFYsmpMatrix<long double>;
template class TPZFYsmpMatrix<double>;
template class TPZFYsmpMatrix<float>;
template class TPZFYsmpMatrix<std::complex<long double>>;
template class TPZFYsmpMatrix<std::complex<double>>;
template class TPZFYsmpMatrix<std::complex<float>>;
