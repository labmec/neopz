#include "Electromagnetics/TPZWaveguideModalAnalysisPML.h"
#include "TPZMaterialDataT.h"

using namespace std::complex_literals;

TPZWaveguideModalAnalysisPML::TPZWaveguideModalAnalysisPML(
    const int id,const TPZWaveguideModalAnalysis &mat,
    const bool &attX, REAL &pmlBeginX,
    const bool &attY, REAL &pmlBeginY,
    const REAL &alphaMax, const REAL &d) :
    TPZWaveguideModalAnalysis(mat),
    fAttX(attX), fAttY(attY), fPmlBeginX(pmlBeginX),fPmlBeginY(pmlBeginY),
    fAlphaMax(alphaMax), fD(d)
{
    this->SetId(id);
    if(fAlphaMax < 0){//for the attenuation to happen this value must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"alpha max is invalid : "<<alphaMax<<std::endl;
        DebugStop();
    }
     
    if(fD < 0){ // pml width must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"PML width is invalid : "<<d<<std::endl;
        DebugStop();
    }
    if(!fAttX && !fAttY){//a pml should attenuate at some direction
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"PML must attenuate at at least one direction"<<std::endl;
        DebugStop();
    }
}

TPZWaveguideModalAnalysisPML * TPZWaveguideModalAnalysisPML::NewMaterial() const
{
    return new TPZWaveguideModalAnalysisPML;
}

void TPZWaveguideModalAnalysisPML::ComputeSParameters(const TPZVec<REAL> &x,
                                                      CSTATE &sx,
                                                      CSTATE &sy)
{
    /*****************CALCULATE S PML PARAMETERS*************************
     * In the current application, the waveguide's cross section is always
     * in the xy-plane. Therefore, sz will always be unity, and omitted for
     * the folllowing calculations. The same principle applies, for instance,
     * for the z-component of the hcurl functions, the x and y components of
     * their curl and so on.
     */

    sx = 1;
    sy = 1;
    if(fAttX){
        sx = 1. - 1i * fAlphaMax * ((x[0]-fPmlBeginX) / fD ) * ((x[0]-fPmlBeginX) / fD );
    }
    if(fAttY){
        sy = 1. - 1i * fAlphaMax * ((x[1]-fPmlBeginY) / fD ) * ((x[1]-fPmlBeginY) / fD );
    }
}
void TPZWaveguideModalAnalysisPML::Contribute(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec, REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
    
    TPZManVector<CSTATE,3> oldUr(fUr);
    TPZManVector<CSTATE,3> oldEr(fEr);
    CSTATE sx{1}, sy{1};
    ComputeSParameters(datavec[0].x,sx,sy);
    fUr={
        oldUr[0] * sy / sx,
        oldUr[1] * sx / sy,
        oldUr[2] * sy * sx
    };
    fEr={
        oldEr[0] * sy / sx,
        oldEr[1] * sx / sy,
        oldEr[2] * sy * sx
    };
    TPZWaveguideModalAnalysis::Contribute(datavec,weight,ek,ef);
    fUr=oldUr;
    fEr=oldEr;
}

void TPZWaveguideModalAnalysisPML::Solution(
    const TPZVec<TPZMaterialDataT<CSTATE>> &datavec,
    int var, TPZVec<CSTATE> &solout)
{
    TPZManVector<CSTATE,3> oldUr(fUr);
    TPZManVector<CSTATE,3> oldEr(fEr);
    CSTATE sx{1}, sy{1};
    ComputeSParameters(datavec[0].x,sx,sy);
    fUr={
        oldUr[0] * sy / sx,
        oldUr[1] * sx / sy,
        oldUr[2] * sy * sx
    };
    fEr={
        oldEr[0] * sy / sx,
        oldEr[1] * sx / sy,
        oldEr[2] * sy * sx
    };
    TPZWaveguideModalAnalysis::Solution(datavec,var,solout);
    fUr=oldUr;
    fEr=oldEr;
}

int TPZWaveguideModalAnalysisPML::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const
{
    int pmax = 0;
    for (int ip=0;  ip<elPMaxOrder.size(); ip++)
    {
        if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }

    const int integrationorder = 4+2*pmax;

    return  integrationorder;
}