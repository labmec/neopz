#include "TPZJacobiPrecond.h"
#include "TPZFMatrixRef.h"

template<class TVar>
TPZJacobiPrecond<TVar>::TPZJacobiPrecond(TPZAutoPointer<TPZMatrix<TVar> >  refmat)
{
  if(!refmat){
    PZError<<__PRETTY_FUNCTION__
           <<"\nInvalid matrix as input! nullptr"<<std::endl;
    DebugStop();
  }
  this->fReferenceMatrix = refmat;
  this->SetMatrix(new TPZFMatrix<TVar>());
}


template<class TVar>
void TPZJacobiPrecond<TVar>::UpdateFrom(TPZAutoPointer<TPZBaseMatrix> mat)
{
  if(mat && mat == this->fReferenceMatrix){
    const auto sz = mat->Rows();
    if(sz!= this->Matrix()->Rows()){this->Matrix()->Resize(sz,1);}
    auto ms = this->Matrix()->Storage().Elem();
    auto matval = TPZAutoPointerDynamicCast<TPZMatrix<TVar>>(mat);
    for(int i = 0; i < sz; i++){
      ms[i] = 1./(matval->GetVal(i,i));
    }
  }
}

template<class TVar>
void TPZJacobiPrecond<TVar>::Solve(const TPZFMatrix<TVar> &F, TPZFMatrix<TVar> &result, TPZFMatrix<TVar> *residual)
{
  const auto nr = F.Rows();

#ifdef PZDEBUG
  if(nr != this->Matrix()->Rows()){
    PZError<<__PRETTY_FUNCTION__
           <<"\nIncompatible sizes!\n"
           <<"\tmatrix diagonal size: "<<this->Matrix()->Rows()<<'\n'
           <<"\trhs size: "<<nr<<'\n'
           <<"Aborting..."<<std::endl;
    DebugStop();
  }

#endif
  const auto nc = F.Cols();
  this->fScratch = F;
  const TPZFMatrix<TVar> fcp = F;//F and result might be the sme

  result.Redim(nr,nc);

  const int ncolors = this->fColorVec.size();
  if(ncolors > 0){
    for(int ic = 0; ic < nc; ic++){
      for(int ico = 0; ico < ncolors; ico++){
        auto &current_nodes = this->fColorVec[ico];
        for(auto ir : current_nodes){
          //we opt for bound checking in this loop
          const auto res =
            result(ir,ic) + this->fScratch(ir,ic)*this->Matrix()->operator()(ir,0);
          result.PutVal(ir,ic, res);
        }
        this->Matrix()->Residual(result, fcp, this->fScratch);
      }
    }  
  }else{
    //we performed bound checking at the beginning of this function
    auto ms = this->Matrix()->Storage().Elem();
    auto ss = this->fScratch.Storage().Elem();
    auto rs = result.Storage().Elem();
    for(int ic = 0; ic < nc; ic++){
      for(auto ir = 0; ir < nr; ir++){
        rs[nr*ic+ir] += ss[nr*ic+ir] * ms[ir];
      }
    }
  }
  if(residual) *residual = this->fScratch;
}

template<class TVar>
void TPZJacobiPrecond<TVar>::SetColoring(const TPZVec<int64_t> &colors,
                                         const int numcolors)
{
  fColorVec.Resize(numcolors);
  TPZVec<int64_t> colorcount(numcolors,0);
  for(auto c : colors){
    colorcount[c]++;
  }
  for(auto c = 0; c < numcolors; c++){
    fColorVec[c].Resize(colorcount[c]);
  }
  colorcount.Fill(0);
  const int nbl = colors.size();
  for(auto ibl = 0; ibl < nbl; ibl++){
    auto c = colors[ibl];
    auto pos = (colorcount[c])++;
    fColorVec[c][pos] = ibl;
  }
}

template class TPZJacobiPrecond<STATE>;
template class TPZJacobiPrecond<CSTATE>;