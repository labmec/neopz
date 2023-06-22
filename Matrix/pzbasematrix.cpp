#include "pzbasematrix.h"
#include "TPZStream.h"
#include "Hash/TPZHash.h"


void TPZBaseMatrix::GetRowIndices(const int64_t i, TPZVec<int64_t> &indices) const{
  indices.Resize(fCol);
  for(int i = 0; i < fCol; i++){indices[i] = i;}
}

void TPZBaseMatrix::SetSymmetry(SymProp sp){
  if(fRow!=fCol && sp != SymProp::NonSym){
    PZError<<__PRETTY_FUNCTION__
           <<"\nTrying to set a non-square matrix as symmetric/hermitian\n"
           <<"Aborting..."<<std::endl;
    DebugStop();
  }
  fSymProp = sp;
}
void TPZBaseMatrix::Read(TPZStream &buf, void *context){
  buf.Read(&fRow);
  buf.Read(&fCol);
  int temp;
  buf.Read(&temp);
  fDecomposed = (DecomposeType) temp;
  buf.Read(&temp);
  fDefPositive = temp == 1 ? 1 : 0;
}

void TPZBaseMatrix::Write(TPZStream &buf, int withclassid) const{
  buf.Write(&fRow);
  buf.Write(&fCol);
  int temp = (int) fDecomposed;
  buf.Write(&temp);
  temp = fDefPositive ? 1 : 0;
  buf.Write(&temp);
}


int TPZBaseMatrix::ClassId() const{
  return Hash("TPZBaseMatrix");
}
