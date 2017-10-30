
#include "TPZPlasticState.h"

template <>
void TPZPlasticState<STATE>::Read(TPZStream& buf, void* context){
    fEpsT.Read(buf, context);
    fEpsP.Read(buf, context);
    buf.Read(&fAlpha);
}
   
template <>
void TPZPlasticState<STATE>::Write(TPZStream& buf, int withclassid) const{
    fEpsT.Write(buf, withclassid);
    fEpsP.Write(buf, withclassid);
    buf.Write(&fAlpha);
}