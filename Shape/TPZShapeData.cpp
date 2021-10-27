#include "TPZShapeData.h"
#include <string>

TPZShapeData::TPZShapeData() :
    TPZRegisterClassId(&TPZShapeData::ClassId)
{

}

std::string TPZShapeData::ShapeFunctionType() const
{
    switch(fShapeType){
        case EEmpty:
            return "NotInitialized";
        case EScalarShape:
            return "Scalar";
        case EVecandShape:
            return "Vector combined with Scalar";
        case EVecShape:
            return "Vector shape";
        default:
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Invalid value\n";
            DebugStop();
    }
    return "All Wrong!\n";
}


void TPZShapeData::Print(std::ostream &out) const
{
    out << "Shape function type " << ShapeFunctionType() << std::endl;
    fPhi.Print("phi",out);
    fDPhi.Print("dphi",out);
    out << "fMasterDirections" << fMasterDirections << std::endl;
    if (fSDVecShapeIndex.size()) {
        out << "VecShapeIndex: ";
        for (int64_t i = 0; i < fSDVecShapeIndex.size(); i++) {
            out << fSDVecShapeIndex[i].first << '/' << fSDVecShapeIndex[i].second << ' ';
        }
        out << '\n';
    }
}

/** Print the data in a format suitable for Mathematica */
void TPZShapeData::PrintMathematica(std::ostream &out) const
{
    fPhi.Print("phi = ",out,EMathematicaInput);
    fDPhi.Print("dphi = ",out,EMathematicaInput);
    out << "fMasterDirections" << fMasterDirections << std::endl;
}


int TPZShapeData::ClassId() const
{
    return Hash("TPZShapeData");
}
void TPZShapeData::Write(TPZStream &buf, int withclassid) const
{

}
    
void TPZShapeData::Read(TPZStream &buf, void *context)
{
}
