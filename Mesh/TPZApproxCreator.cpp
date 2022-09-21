//
// Created by victor on 21/09/2022.
//

#include "TPZApproxCreator.h"


TPZApproxCreator::TPZApproxCreator(TPZGeoMesh *gmesh) : fGeoMesh(gmesh)
{
}

/**Insert a material object in the datastructure*/
int TPZApproxCreator::InsertMaterialObject(int matid, TPZMaterial * mat) {
    if(!mat) DebugStop();
    std::pair<int, TPZMaterial*> matinfo;
    matinfo.first = matid;
    matinfo.second = mat;
    fMaterialVec.Push(matinfo);
    return fMaterialVec.size();
}

int TPZApproxCreator::InsertBCMaterialObject(int matid, TPZBndCond *mat)  {
    if(!mat) DebugStop();
    std::pair<int, TPZBndCond*> matinfo;
    matinfo.first = matid;
    matinfo.second = mat;
    fBCMaterialVec.Push(matinfo);
    return fBCMaterialVec.size();
}



void TPZApproxCreator::ComputePeriferalMaterialIds(int base)
{
    if(base < 2) base = 2;
    int max_matid = 0;
    int64_t nel = fGeoMesh->NElements();
    for (int64_t el = 0; el < nel; el++) {
        TPZGeoEl *gel = fGeoMesh->Element(el);
        if(!gel) continue;
        max_matid = std::max(max_matid,gel->MaterialId());
    }
    int remain = max_matid % base;
    int matid_base = max_matid-remain + base;
    fHybridizationData.fWrapMatId = matid_base;
    fHybridizationData.fInterfaceMatId = matid_base+base;
    fHybridizationData.fLagrangeMatId = matid_base + 2*base;
    fHybridizationData.fSecondInterfaceMatId = matid_base + 3*base;
    fHybridizationData.fSecondLagrangeMatId = matid_base + 4*base;
}

void TPZApproxCreator::AddHybridizationGeoElements(){
    std::cout << "\n\nImplement me\n\n";
    DebugStop();
}

void TPZApproxCreator::Print(std::ostream &ofs){
    std::stringstream ss;
    ss << "\n";
    int numStars = 70;
    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    std::string hybridization = [](HybridizationType hybtype){
        switch (hybtype){
            case HybridizationType::ENone:
                return  "ENone";
            case HybridizationType::EStandard:
                return "EStandard";
            case HybridizationType::EStandardSquared:
                return "EStandardSquared";
            case HybridizationType::ESemi:
                return "ESemi";
            default:
                DebugStop();
        }
    }(this->fHybridType);

    std::string problem = [](ProblemType probtype){
        switch (probtype){
            case ProblemType::ENone:
                return  "ENone";
            case ProblemType::EElastic:
                return "EElastic";
            case ProblemType::EDarcy:
                return "EDarcy";
            case ProblemType::EStokes:
                return "EStokes";
            default:
                DebugStop();
        }
    }(this->fProbType);

    ss << __PRETTY_FUNCTION__ << "\nApproximation space creation set up:\n\n";
    ss << "fHybridType:\t" << hybridization << "fProbType:\t" << problem;
    ss <<"\nExtra Internal Pressure Order: " << fExtraInternalPOrder;
    ss << "\nDefault Pressure Order:\t" << fDefaultPOrder;
    if(this->fHybridType != HybridizationType::ENone) {
        ss << "\nfHybridizeBCLevel: " << fHybridizationData.fHybridizeBCLevel;
    }
    ss << "\nfShouldCondense:\t" << fShouldCondense;
    ss << "\nfIsEnhancedSpaces" << fIsEnhancedSpaces;
    ss << "\n\nmatids:\n\t{";
    for(auto it : fMaterialVec){
        ss << it.first <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    ss << "bcmatids:\n\t{";
    for(auto it : fBCMaterialVec){
        ss << it.first <<", ";
    }
    ss.seekp(-2,ss.cur); //return stringstream head by two positions, as consequence, this data will be erased.
    ss << "}\n";

    if(this->fHybridType == HybridizationType::EStandard || this->fHybridType == HybridizationType::EStandardSquared) {
        ss << "peripheralMatids:\n";
        ss << "\tfWrapMatId = " << fHybridizationData.fWrapMatId;
        ss << "\n\tfInterfaceMatId = " << fHybridizationData.fInterfaceMatId;
        ss << "\n\tfLagrangeMatId = " << fHybridizationData.fLagrangeMatId;
        if( this->fHybridType == HybridizationType::EStandardSquared) {
            ss << "\n\tfSecondInterfaceMatId = " << fHybridizationData.fSecondInterfaceMatId;
            ss << "\n\tfSecondLagrangeMatId = " << fHybridizationData.fSecondLagrangeMatId;
        }
        ss << "\n";
    }

    for(int istar = 0; istar < numStars ; istar++)
        ss <<"*";
    ss <<"\n";

    ofs << ss.str();
}