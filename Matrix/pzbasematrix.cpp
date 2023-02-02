#include "pzbasematrix.h"
#include "TPZStream.h"
#include "Hash/TPZHash.h"


void TPZBaseMatrix::Read(TPZStream &buf, void *context){
  buf.Read(&fRow);
  buf.Read(&fCol);
    int temp;
  buf.Read(&temp);
    fDecomposed = (DecomposeType) temp;
  buf.Read(&fDefPositive);
}

void TPZBaseMatrix::Write(TPZStream &buf, int withclassid) const{
  buf.Write(&fRow);
  buf.Write(&fCol);
    int temp;
  buf.Write(&temp);
  buf.Write(&fDefPositive);
}


int TPZBaseMatrix::ClassId() const{
  return Hash("TPZBaseMatrix");
}
