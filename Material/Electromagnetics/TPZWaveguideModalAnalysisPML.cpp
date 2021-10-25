#include "Electromagnetics/TPZWaveguideModalAnalysisPML.h"
#include "TPZMaterialDataT.h"

using namespace std::complex_literals;

TPZWaveguideModalAnalysisPML::TPZWaveguideModalAnalysisPML(
    const int id, const TPZWaveguideModalAnalysis &mat):
    TPZWaveguideModalAnalysis(mat)
{
    this->SetId(id);
}

void TPZWaveguideModalAnalysisPML::SetAttX(const REAL pmlBegin,
                                           const STATE alpha,
                                           const REAL d)
{
    if(d < 0){ // pml width must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"PML width is invalid : "<<d<<std::endl;
        DebugStop();
    }
    if(alpha < 0){//for the attenuation to happen this value must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"alpha max is invalid : "<<alpha<<std::endl;
        DebugStop();
    }
    fAttX = true;
    fPmlBeginX = pmlBegin;
    fAlphaMaxX = alpha;
    fDX = d;
}

void TPZWaveguideModalAnalysisPML::SetAttY(const REAL pmlBegin,
                                           const STATE alpha,
                                           REAL d)
{
    if(d < 0){ // pml width must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"PML width is invalid : "<<d<<std::endl;
        DebugStop();
    }
    if(alpha < 0){//for the attenuation to happen this value must be positive
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"alpha max is invalid : "<<alpha<<std::endl;
        DebugStop();
    }
    fAttY = true;
    fPmlBeginY = pmlBegin;
    fAlphaMaxY = alpha;
    fDY = d;
}

TPZWaveguideModalAnalysisPML * TPZWaveguideModalAnalysisPML::NewMaterial() const
{
    return new TPZWaveguideModalAnalysisPML;
}

void TPZWaveguideModalAnalysisPML::ComputeSParameters(const TPZVec<REAL> &x,
                                                      CSTATE &sx,
                                                      CSTATE &sy) const
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
        sx = 1. - 1i * fAlphaMaxX * ((x[0]-fPmlBeginX) / fDX )
            * ((x[0]-fPmlBeginX) / fDX );
    }
    if(fAttY){
        sy = 1. - 1i * fAlphaMaxY * ((x[1]-fPmlBeginY) / fDY ) *
            ((x[1]-fPmlBeginY) / fDY );
    }
}

void TPZWaveguideModalAnalysisPML::GetPermittivity(
  const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const
{
    TPZWaveguideModalAnalysis::GetPermittivity(x,er);
    CSTATE sx{1}, sy{1};
    ComputeSParameters(x,sx,sy);

    er[0] *= sy / sx;
    er[1] *= sx / sy;
    er[2] *= sy * sx;
}

void TPZWaveguideModalAnalysisPML::GetPermeability(
  const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const
{
    TPZWaveguideModalAnalysis::GetPermeability(x,ur);
    CSTATE sx{1}, sy{1};
    ComputeSParameters(x,sx,sy);
    
    ur[0] *= sy / sx;
    ur[1] *= sx / sy;
    ur[2] *= sy * sx;
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