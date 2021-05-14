#include "TPZBndCond.h"
#include "TPZMaterial.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
    
TPZMaterial * TPZBndCond::Material(){
    return fMaterial;
}

void TPZBndCond::Print(std::ostream & out) const {
	out << "BC type: " << fType << '\n';
}

int TPZBndCond::ClassId() const{
    return Hash("TPZBndCond");
}

void TPZBndCond::Read(TPZStream& buf, void* context){
    //let us NOT read anything from the material here
    buf.Read(&fType);
}

void TPZBndCond::Write(TPZStream& buf, int withclassid) const{
    //let us NOT write anything from the material here
    buf.Write(&fType);
}