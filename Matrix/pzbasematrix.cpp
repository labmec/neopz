#include "pzbasematrix.h"
#include "TPZStream.h"
#include "Hash/TPZHash.h"


void TPZBaseMatrix::Read(TPZStream &buf, void *context){
  buf.Read(&fRow);
  buf.Read(&fCol);
  buf.Read(&fDecomposed);
  buf.Read(&fDefPositive);
}

void TPZBaseMatrix::Write(TPZStream &buf, int withclassid) const{
  buf.Write(&fRow);
  buf.Write(&fCol);
  buf.Write(&fDecomposed);
  buf.Write(&fDefPositive);
}


int TPZBaseMatrix::ClassId() const{
  return Hash("TPZBaseMatrix");
}