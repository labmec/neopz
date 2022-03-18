#include "TPZPlanarWGScatteringPML.h"
#include "TPZMaterialDataT.h"

using namespace std::complex_literals;

TPZPlanarWGScatteringPML::TPZPlanarWGScatteringPML(
    const int id, const TPZPlanarWGScattering &mat):
    TPZPlanarWGScattering(mat)
{
    this->SetId(id);
}

void TPZPlanarWGScatteringPML::SetAttX(const REAL pmlBegin,
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

void TPZPlanarWGScatteringPML::SetAttY(const REAL pmlBegin,
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

TPZPlanarWGScatteringPML * TPZPlanarWGScatteringPML::NewMaterial() const
{
    return new TPZPlanarWGScatteringPML;
}

void TPZPlanarWGScatteringPML::ComputeSParameters(const TPZVec<REAL> &x,
                                                  CSTATE &sx,
                                                  CSTATE &sy,
                                                  CSTATE &sz) const
{
    /*****************CALCULATE S PML PARAMETERS*************************
     * In the current application, the waveguide's cross section is always
     * in the xy-plane. Therefore, in the propagation problem of the planar
     * waveguide, the domain is the yz-plane (assuming symmetry in the 
     * x-direction).
     */
    if(fAttX){
        const auto distx = (x[0]-fPmlBeginX) / fDX;
        sx = 1. - 1i * fAlphaMaxX * distx * distx;
    }
    if(fAttY){
        const auto disty = (x[1]-fPmlBeginY) / fDY;
        sy = 1. - 1i * fAlphaMaxY * disty * disty;
    }
}

TPZVec<CSTATE> TPZPlanarWGScatteringPML::GetPermittivity(
  const TPZVec<REAL> &x) const
{
    auto er = TPZPlanarWGScattering::GetPermittivity(x);
    CSTATE sx{1}, sy{1}, sz{1};
    ComputeSParameters(x,sx,sy,sz);

    er[0] *= sy * sz / sx;
    er[1] *= sz * sx / sy;
    er[2] *= sx * sy / sz;
    return er;
}

TPZVec<CSTATE> TPZPlanarWGScatteringPML::GetPermeability(
  const TPZVec<REAL> &x) const
{
    auto ur = TPZPlanarWGScattering::GetPermeability(x);
    CSTATE sx{1}, sy{1}, sz{1};
    ComputeSParameters(x,sx,sy,sz);
    
    ur[0] *= sy * sz / sx;
    ur[1] *= sz * sx / sy;
    ur[2] *= sx * sy / sz;
    return ur;
}

int TPZPlanarWGScatteringPML::IntegrationRuleOrder(const int elPMaxOrder) const
{
    const int integrationorder = 2+2*elPMaxOrder;

    return  integrationorder;
}