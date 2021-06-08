/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */

#include <memory.h>

#include "pzsysmp.h"
#include "pzfmatrix.h"
#include "pzstack.h"
// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix() : TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
TPZMatrix<TVar>() {

#ifdef CONSTRUCTOR
    cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}


template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix(const int64_t rows,const int64_t cols ) : TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
TPZMatrix<TVar>(rows,cols) {

#ifdef CONSTRUCTOR
	cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}

template<class TVar>
TPZSYsmpMatrix<TVar>::~TPZSYsmpMatrix() {
	// Deletes everything associated with a TPZSYsmpMatrix
#ifdef DESTRUCTOR
	cerr << "~TPZSYsmpMatrix()\n";
#endif
}

template <class TVar> int TPZSYsmpMatrix<TVar>::Zero() {
  fA.Fill(0.);
  fDiag.Fill(0.);
#ifndef USING_MKL
  TPZMatrix<TVar>::fDecomposed = ENoDecompose;
#endif
  return 0;
}

// ****************************************************************************
//
// Find the element of the matrix at (row,col) in the stencil matrix
//
// ****************************************************************************

template<class TVar>
const TVar TPZSYsmpMatrix<TVar>::GetVal(const int64_t row,const int64_t col ) const {
    if (row > col) {
        for(int ic=fIA[col] ; ic < fIA[col+1]; ic++ ) {
            if ( fJA[ic] == row ) {
                if constexpr (is_complex<TVar>::value){
                    return std::conj(fA[ic]);
                }else{
                    return fA[ic];
                }
            }
        }
        return (TVar)0;
    }
	for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
		if ( fJA[ic] == col ) return fA[ic];
	}
	return (TVar)0;
}

/** @brief Put values without bounds checking \n
 *  This method is faster than "Put" if DEBUG is defined.
 */
template<class TVar>
int TPZSYsmpMatrix<TVar>::PutVal(const int64_t r,const int64_t c,const TVar & val )
{
    // Get the matrix entry at (row,col) without bound checking
    int64_t row(r),col(c);
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
    }
    for(int ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fA[ic] = val;
            return 0;
        }
    }
    if (val != (TVar(0.))) {
        DebugStop();
    }
    return 0;
    
}



template<class TVar>
void TPZSYsmpMatrix<TVar>::CheckTypeCompatibility(const TPZMatrix<TVar>*A,
                                                  const TPZMatrix<TVar>*B)const
{
  auto incompatSparse = [](){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nERROR: incompatible matrices\n.Aborting...\n";
    DebugStop();
  };
  auto aPtr = dynamic_cast<const TPZSYsmpMatrix<TVar>*>(A);
  auto bPtr = dynamic_cast<const TPZSYsmpMatrix<TVar>*>(B);
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
TPZSYsmpMatrix<TVar> TPZSYsmpMatrix<TVar>::operator+(const TPZSYsmpMatrix<TVar>&mat) const
{
  CheckTypeCompatibility(this,&mat);
	auto res(*this);
  const auto sizeA = res.fA.size();
  for(auto i = 0; i < sizeA; i++) res.fA[i] += mat.fA[i];
	return res;
}

template<class TVar>
TPZSYsmpMatrix<TVar> TPZSYsmpMatrix<TVar>::operator-(const TPZSYsmpMatrix<TVar>&mat) const
{
	CheckTypeCompatibility(this,&mat);
	auto res(*this);
  const auto sizeA = res.fA.size();
  for(auto i = 0; i < sizeA; i++) res.fA[i] -= mat.fA[i];
	return res;
}

template<class TVar>
TPZSYsmpMatrix<TVar> TPZSYsmpMatrix<TVar>::operator*(const TVar alpha) const
{
	auto res(*this);
	for(auto &el : res.fA) el *= alpha;
	return res;
}

template<class TVar>
TPZSYsmpMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator+=(const TPZSYsmpMatrix<TVar> &A )
{
	TPZSYsmpMatrix<TVar> res((*this)+A);
	*this = res;
	return *this;
}
template<class TVar>
TPZSYsmpMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator-=(const TPZSYsmpMatrix<TVar> &A )
{
	TPZSYsmpMatrix<TVar> res((*this)-A);
	*this = res;
	return *this;
}
template<class TVar>
TPZSYsmpMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator*=(const TVar val)
{
	TPZSYsmpMatrix<TVar> res((*this)*val);
	*this = res;
	return *this;
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
							 TPZFMatrix<TVar> &z,
							 const TVar alpha,const TVar beta,const int opt) const {	
	// Determine how to initialize z
    this->PrepareZ(y,z,beta,opt);
	
	// Compute alpha * A * x
    const int64_t ncols = x.Cols();
    const int64_t nrows = this->Rows();

    
    if constexpr (is_complex<TVar>::value){
        auto GetMyVal = [](const int64_t ir, const int64_t ic,
                           const bool opt, const TVar val){
            if((ir <= ic && !opt)||(ir >= ic && opt)) return val;
            else return std::conj(val);
        };
        for (int64_t col=0; col<ncols; col++){
            for(int64_t row=0; row<nrows; row++) {
                for(int64_t iv=fIA[row]; iv<fIA[row+1]; iv++) {
                    const int64_t ic = fJA[iv];
                    const TVar val = GetMyVal(row,ic,opt,fA[iv]);
                    z(row,col) += alpha * val * x.GetVal(ic,col);
                    if(row != ic){
                        z(ic,col) += alpha* std::conj(val) * x.GetVal(row,col);
                    }
                }
            }
        }
    }else{
        for (int64_t col=0; col<ncols; col++){
            for(int64_t row=0; row<nrows; row++) {
                for(int64_t iv=fIA[row]; iv<fIA[row+1]; iv++) {
                    const int64_t ic = fJA[iv];
                    z(row,col) += alpha*fA[iv] * x.GetVal(ic,col);
                    if(row != ic){
                        z(ic,col) += alpha*fA[iv] * x.GetVal(row,col);
                    }
                }
            }
        }
    }
}

// ****************************************************************************
//
// Print the matrix
//
// ****************************************************************************

template<class TVar>
void TPZSYsmpMatrix<TVar>::Print(const char *title, std::ostream &out ,const MatrixOutputFormat form) const {
	// Print the matrix along with a identification title
	if(form == EInputFormat) {
		out << "\nTSYsmpMatrix Print: " << title << '\n'
        << "\tNon zero elements    = " << fA.size()  << '\n'
		<< "\tRows    = " << this->Rows()  << '\n'
		<< "\tColumns = " << this->Cols() << '\n';
		int i;
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(i=0; i<=this->Rows(); i++) {
			out << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(i=this->Rows()+1; i<fIA[this->Rows()]-1; i++) {
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
void TPZSYsmpMatrix<TVar>::ComputeDiagonal() {
	if(!fDiag.size()) fDiag.resize(this->Rows());
	for(int ir=0; ir<this->Rows(); ir++) {
		fDiag[ir] = GetVal(ir,ir);
	}
}

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, int symmetric)
{
    if (!symmetric || nrow != ncol) {
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
        eqs[row].insert(row);
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
            if(*col >= row)
            {
                JA.Push(*col);
                A.Push(orig(row,*col));
            }
        }
        IA[row+1] = JA.size();
    }
    TPZMatrix<TVar>::Resize(nrow,ncol);
    SetData(IA, JA, A);
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & destinationindex){
    int64_t i,j,k = 0;
    TVar value=0.;
    int64_t ipos,jpos;
    for(i=0;i<elmat.Rows();i++){
        for(j=0;j<elmat.Rows();j++){
            ipos=destinationindex[i];
            jpos=destinationindex[j];
            if(jpos<ipos) continue;
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
                if(!flag) std::cout << "TPZSYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;         }
        }
    }
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::AddKel(TPZFMatrix<TVar> & elmat, TPZVec<int64_t> & sourceindex, TPZVec<int64_t> & destinationindex){
	int64_t i,j,k = 0;
	TVar value=0.;
	int64_t ipos,jpos;
	for(i=0;i<sourceindex.NElements();i++){
		for(j=0;j<sourceindex.NElements();j++){
			ipos=destinationindex[i];
			jpos=destinationindex[j];
            if(jpos<ipos) continue;
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
				if(!flag) std::cout << "TPZSYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;         }
		}
	}
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::SetIsDecomposed(int val)
{
	if(val)
		fPardisoControl.fDecomposed = true;
	TPZBaseMatrix::SetIsDecomposed(val);
}

#ifdef USING_MKL

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{
    Decompose_LDLt();
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt()
{
    if(this->IsDecomposed() == ELDLt) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    typename TPZPardisoSolver<TVar>::MStructure str =
        TPZPardisoSolver<TVar>::MStructure::ESymmetric;
    typename TPZPardisoSolver<TVar>::MSystemType sysType =
		TPZPardisoSolver<TVar>::MSystemType::ESymmetric;
	typename TPZPardisoSolver<TVar>::MProperty prop =
		this->IsDefPositive() ?
		TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
		TPZPardisoSolver<TVar>::MProperty::EIndefinite;
    fPardisoControl.SetStructure(str);
	fPardisoControl.SetMatrixType(sysType,prop);
    fPardisoControl.Decompose(this);
    this->SetIsDecomposed(ELDLt);
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky()
{
    if(this->IsDecomposed() == ECholesky) return 1;
    if (this->IsDecomposed() != ENoDecompose) {
        DebugStop();
    }
    this->fDefPositive = true;
    typename TPZPardisoSolver<TVar>::MStructure str =
        TPZPardisoSolver<TVar>::MStructure::ESymmetric;
    typename TPZPardisoSolver<TVar>::MSystemType sysType =
		TPZPardisoSolver<TVar>::MSystemType::ESymmetric;
	typename TPZPardisoSolver<TVar>::MProperty prop =
		this->IsDefPositive() ?
		TPZPardisoSolver<TVar>::MProperty::EPositiveDefinite:
		TPZPardisoSolver<TVar>::MProperty::EIndefinite;
    fPardisoControl.SetStructure(str);
	fPardisoControl.SetMatrixType(sysType,prop);
    fPardisoControl.Decompose(this);
    this->SetIsDecomposed(ELDLt);
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    return Decompose_Cholesky();
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    fPardisoControl.Solve(this,*b,x);
    *b = x;
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    TPZFMatrix<TVar> x(*b);
    fPardisoControl.Solve(this,*b,x);
    *b = x;
    return 1;
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    return 1;
}

#else
//perhaps we could default to a less eficient implementation for
//solving. Perhaps removing the DebugStop() on PutVal() would be enough?
#define NOMKL \
    PZError<<__PRETTY_FUNCTION__<<" is not available if NeoPZ ";\
    PZError<<"was not configured with MKL. Aborting..."<<std::endl;\
    DebugStop();\
    return -1;


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt(std::list<int64_t> &singular)
{    
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_LDLt()
{
    NOMKL
    
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky()
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Decompose_Cholesky(std::list<int64_t> &singular)
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LForward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_LBackward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Diag( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Forward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}


template<class TVar>
int TPZSYsmpMatrix<TVar>::Subst_Backward( TPZFMatrix<TVar>* b ) const
{
    NOMKL
}
#endif


template<class TVar>
int TPZSYsmpMatrix<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}
template class TPZSYsmpMatrix<double>;
template class TPZSYsmpMatrix<float>;
template class TPZSYsmpMatrix<long double>;
template class TPZSYsmpMatrix<std::complex<float>>;
template class TPZSYsmpMatrix<std::complex<double>>;
template class TPZSYsmpMatrix<std::complex<long double>>;
