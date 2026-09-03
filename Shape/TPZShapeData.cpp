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
    out << "H1 data" << std::endl;
    fH1.fPhi.Print("phi",out);
    fH1.fDPhi.Print("dphi",out);
    out << "HDiv data" << std::endl;
    out << "fMasterDirections" << fHDiv.fMasterDirections << std::endl;
    if (fHDiv.fSDVecShapeIndex.size()) {
        out << "VecShapeIndex: ";
        for (int64_t i = 0; i < fHDiv.fSDVecShapeIndex.size(); i++) {
            out << fHDiv.fSDVecShapeIndex[i].first << '/' << fHDiv.fSDVecShapeIndex[i].second << ' ';
        }
        out << '\n';
    }
    out << "HCurl data" << std::endl;
    out << "fMasterDirections" << fHCurl.fMasterDirections << std::endl;
    if (fHCurl.fSDVecShapeIndex.size()) {
        out << "VecShapeIndex: ";
        for (int64_t i = 0; i < fHCurl.fSDVecShapeIndex.size(); i++) {
            out << fHCurl.fSDVecShapeIndex[i].first << '/' << fHCurl.fSDVecShapeIndex[i].second << ' ';
        }
        out << '\n';
    }
}

/** Print the data in a format suitable for Mathematica */
void TPZShapeData::PrintMathematica(std::ostream &out) const
{
    fH1.fPhi.Print("H1phi = ",out,EMathematicaInput);
    fH1.fDPhi.Print("H1dphi = ",out,EMathematicaInput);
    out << "HDivfMasterDirections = " << fHDiv.fMasterDirections << std::endl;
    out << "HCurlfMasterDirections = " << fHCurl.fMasterDirections << std::endl;
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
