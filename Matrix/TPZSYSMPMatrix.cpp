/**
 * @file
 * @brief Contains the implementation of the TPZSYsmpMatrix methods.
 */

#include <memory.h>

#include "TPZSYSMPMatrix.h"
#include "pzfmatrix.h"
#include "pzstack.h"
#include "TPZParallelUtils.h"
// ****************************************************************************
// 
// Constructors and the destructor
// 
// ****************************************************************************

template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix() : TPZRegisterClassId(&TPZSYsmpMatrix::ClassId),
                                         TPZMatrix<TVar>() {
  this->fSymProp = SymProp::Herm;
#ifdef CONSTRUCTOR
    cerr << "TPZSYsmpMatrix(int rows,int cols)\n";
#endif
}


template<class TVar>
TPZSYsmpMatrix<TVar>::TPZSYsmpMatrix(const int64_t rows,const int64_t cols) :
    TPZRegisterClassId(&TPZSYsmpMatrix::ClassId), TPZMatrix<TVar>(rows,cols) {
  this->fSymProp = SymProp::Herm;
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
        for(int64_t ic=fIA[col] ; ic < fIA[col+1]; ic++ ) {
            if ( fJA[ic] == row ) {
                if constexpr (is_complex<TVar>::value){
                    return this->fSymProp == SymProp::Herm ? std::conj(fA[ic]) : fA[ic];
                }else{
                    return fA[ic];
                }
            }
        }
        return (TVar)0;
    }
	for(int64_t ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
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
    TVar valcp{val};
    if (r > c) {
        int64_t temp = r;
        row = col;
        col = temp;
        if constexpr (is_complex<TVar>::value){
            valcp = this->fSymProp == SymProp::Herm? std::conj(val) : val;
        }
    }
    for(int64_t ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
        if ( fJA[ic] == col )
        {
            fA[ic] = valcp;
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
#ifdef PZDEBUG
	CheckTypeCompatibility(this, &A);
#endif
	const int nnzero = this->fA.size();
	for(int i = 0; i < nnzero; i++){
		this->fA[i] += A.fA[i];
	}
	return *this;
}
template<class TVar>
TPZSYsmpMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator-=(const TPZSYsmpMatrix<TVar> &A )
{
#ifdef PZDEBUG
	CheckTypeCompatibility(this, &A);
#endif
	const int nnzero = this->fA.size();
	for(int i = 0; i < nnzero; i++){
		this->fA[i] -= A.fA[i];
	}
	return *this;
}
template<class TVar>
TPZMatrix<TVar> &TPZSYsmpMatrix<TVar>::operator*=(const TVar val)
{
	for(auto &el : this->fA) el *= val;
	return *this;
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::MultAdd(const TPZFMatrix<TVar> &x,const TPZFMatrix<TVar> &y,
							 TPZFMatrix<TVar> &z,
							 const TVar alpha,const TVar beta,const int opt) const {	

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

    // Determine how to initialize z
    this->PrepareZ(y,z,beta,opt);
	
	// Compute alpha * A * x
    const int64_t ncols = x.Cols();
    const int64_t nrows = this->Rows();

    const bool must_conj =
      is_complex<TVar>::value && this->GetSymmetry() == SymProp::Herm;


    if(must_conj){
      if constexpr (is_complex<TVar>::value){
        auto GetMyVal = [](const int64_t ir, const int64_t ic,
                           const bool opt, const TVar val){
          if((ir <= ic && !opt)||(ir > ic && opt)) return val;
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
      }//no need for else
    }
    else{
      for (int64_t col=0; col<ncols; col++){
        for(int64_t row=0; row<nrows; row++) {
          for(int64_t iv=fIA[row]; iv<fIA[row+1]; iv++) {
            const int64_t ic = fJA[iv];
            z(row,col) += alpha * fA[iv] * x.GetVal(ic,col);
            if(row != ic){
              z(ic,col) += alpha * fA[iv] * x.GetVal(row,col);
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
		
		out << "\tIA\tJA\tA\n"
		<< "\t--\t--\t-\n";
		for(int64_t i=0; i<=this->Rows(); i++) {
			out << i      << '\t'
			<< fIA[i] << '\t'
			<< fJA[i] << '\t'
			<< fA[i]  << '\n';
		}
		for(int64_t i=this->Rows()+1; i<fIA[this->Rows()]-1; i++) {
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

template<class TVar>
void TPZSYsmpMatrix<TVar>::SetSymmetry (SymProp sp){
    if(sp == SymProp::NonSym){
        PZError<<__PRETTY_FUNCTION__
               <<"\nTrying to set matrix with symmetric storage as non symmetric\n"
               <<"Aborting..."<<std::endl;
        DebugStop();
    }
    TPZBaseMatrix::SetSymmetry(sp);
}

/** @brief Fill matrix storage with randomic values */
/** This method use GetVal and PutVal which are implemented by each type matrices */
template<class TVar>
void TPZSYsmpMatrix<TVar>::AutoFill(int64_t nrow, int64_t ncol, SymProp sym)
{
    if (sym == SymProp::NonSym || nrow != ncol) {
        DebugStop();
    }
    TPZFMatrix<TVar> orig;
    orig.AutoFill(nrow,ncol,sym);
    
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
                eqs[col].insert(row);
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
            if(!flag) std::cout << "TPZSYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;         
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
            if(!flag) std::cout << "TPZSYsmpMatrix::AddKel: Non existing position on sparse matrix: line =" << ipos << " column =" << jpos << std::endl;
		}
	}
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::AddKelAtomic(TPZFMatrix<TVar>&elmat, TPZVec<int64_t> &sourceindex,  TPZVec<int64_t> &destinationindex){

    int64_t nelem = sourceindex.NElements();
    int64_t icoef,jcoef,ieq,jeq,ieqs,jeqs;
    double prevval;
    int64_t row, col;

    for(icoef=0; icoef<nelem; icoef++) {
        ieq = destinationindex[icoef];
        ieqs = sourceindex[icoef];
        for(jcoef=icoef; jcoef<nelem; jcoef++) {
            jeq = destinationindex[jcoef];
            jeqs = sourceindex[jcoef];
            {
                row = ieq; col = jeq;
                if (row > col) {
                    int64_t temp = row;
                    row = col;
                    col = temp;
                }
                int64_t ic;
                for(ic=fIA[row] ; ic < fIA[row+1]; ic++ ) {
                    if ( fJA[ic] == col )
                    {
                        pzutils::AtomicAdd(fA[ic],elmat(ieqs,jeqs));
                        break;
                    }
                }
                if (ic == fIA[row+1]) {
                    if (elmat(ieqs,jeqs) != (TVar(0.))){
                        DebugStop();

                    }

                }
            }

        }

    }
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::ComputeFillIn(const int64_t resolution, TPZFMatrix<REAL> &fillin) const{
    int64_t nequations = this->Rows();
    int64_t divider = nequations/resolution;
    if(divider*resolution != nequations) divider++;
    REAL factor = 1./(divider*divider);
    fillin.Redim(resolution,resolution);
    
    for(int64_t i = 0; i < nequations ; i++) {
        int64_t indexfirst = fIA[i];
        int64_t indexlast = fIA[i+1];
        for(int64_t c = indexfirst; c < indexlast ; c++) {
            int64_t j = fJA[c];
            fillin(i/divider,j/divider) += factor;
            fillin(j/divider,i/divider) += factor;
        }
    }
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::ComputeFillIn(const int64_t resolution, const TPZVec<long long>& perm, TPZFMatrix<REAL> &fillin) const{
    int64_t nequations = this->Rows();
    if(perm.size() != nequations) DebugStop();
    int64_t divider = nequations/resolution;
    if(divider*resolution != nequations) divider++;
    REAL factor = 1./(divider*divider);
    fillin.Redim(resolution,resolution);
    
    TPZVec<long long> invperm(perm.size());
    for (int64_t i = 0; i < nequations; i++) {
        invperm[perm[i]] = i;
    }
    
    
    for(int64_t i = 0; i < nequations ; i++) {
        int64_t indexfirst = fIA[i];
        int64_t indexlast = fIA[i+1];
        for(int64_t c = indexfirst; c < indexlast ; c++) {
            int64_t j = fJA[c];
            fillin(invperm[i]/divider,invperm[j]/divider) += factor;
            fillin(invperm[j]/divider,invperm[i]/divider) += factor;
        }
    }
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::Permute(const TPZVec<long long>& perm) {
            
    int64_t nequations = this->Rows();
    if(perm.size() != nequations) DebugStop();
    TPZVec<int64_t> tempIA(fIA.size()), tempJA(fJA.size(),-1);
    TPZVec<TVar> tempA(fA.size());
    TPZVec<int64_t> nterms(nequations);
    
    TPZVec<long long> invperm(perm.size());
    for (int64_t i = 0; i < nequations; i++) {
        invperm[perm[i]] = i;
    }
    
    for(int64_t i = 0 ; i < nequations ; i++) {
        int64_t indexfirst = fIA[i];
        int64_t indexlast = fIA[i+1];
        const int64_t linesize = indexlast - indexfirst;
        nterms[invperm[i]] = linesize;
    }
        
    int64_t JAsize = 0;
    tempIA[0] = 0;
    for(int64_t i = 0 ; i < nequations ; i++) {
        tempIA[i+1] = tempIA[i] + nterms[i];
    }
    
    for(int64_t i = 0 ; i < nequations ; i++) {
        int64_t indexfirst = fIA[i];
        int64_t indexlast = fIA[i+1];
        
        for(int64_t c = indexfirst; c < indexlast ; c++) {
            int64_t j = fJA[c];
            int64_t inew = invperm[i], jnew = invperm[j];
            int64_t indexfirstnew = tempIA[inew];
            int64_t indexlastnew = tempIA[inew+1];
            for(int64_t cnew = indexfirstnew; cnew < indexlastnew ; cnew++) {
                if(tempJA[cnew] == -1){
                    tempJA[cnew] = jnew;
                    tempA[cnew] = fA[c];
                    break;
                }
                if(tempJA[cnew] == jnew) DebugStop();
            }
        }
    }
    
    fIA = tempIA;
    fJA = tempJA;
    fA = tempA;
}

template<class TVar>
int TPZSYsmpMatrix<TVar>::ClassId() const{
    return Hash("TPZSYsmpMatrix") ^ TPZMatrix<TVar>::ClassId() << 1;
}


template<class TVar>
void TPZSYsmpMatrix<TVar>::Read(TPZStream &buf, void *context){
	TPZMatrix<TVar>::Read(buf,context);
	buf.Read(fIA);
	buf.Read(fJA);
	buf.Read(fA);
	buf.Read(fDiag);	
}

template<class TVar>
void TPZSYsmpMatrix<TVar>::Write(TPZStream &buf, int withclassid) const{
	TPZMatrix<TVar>::Write(buf,withclassid);
	buf.Write(fIA);
	buf.Write(fJA);
	buf.Write(fA);
	buf.Write(fDiag);
}

#define TEMPL_INST(T) \
  template class TPZRestoreClass<TPZSYsmpMatrix<T>>;\
  template class TPZSYsmpMatrix<T>;

TEMPL_INST(double)
TEMPL_INST(float)
TEMPL_INST(long double)
TEMPL_INST(std::complex<float>)
TEMPL_INST(std::complex<double>)
TEMPL_INST(std::complex<long double>)

#undef TEMPL_INST