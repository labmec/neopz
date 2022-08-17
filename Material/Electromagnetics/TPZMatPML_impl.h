#ifndef _TPZMATPML_IMPL_H_
#define _TPZMATPML_IMPL_H_
#include "TPZMatPML.h"
#include "TPZMaterialDataT.h"

template<class TMAT>
TPZMatPML<TMAT>::TPZMatPML(
  const int id, const TMAT &mat):
  TMAT(mat)
{
  this->SetId(id);
}

template<class TMAT>
void TPZMatPML<TMAT>::SetAttX(const REAL pmlBegin,
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

template<class TMAT>
void TPZMatPML<TMAT>::SetAttY(const REAL pmlBegin,
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

template<class TMAT>
void TPZMatPML<TMAT>::SetAttZ(const REAL pmlBegin,
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
  fAttZ = true;
  fPmlBeginZ = pmlBegin;
  fAlphaMaxZ = alpha;
  fDZ = d;
}

template<class TMAT>
TPZMatPML<TMAT> * TPZMatPML<TMAT>::NewMaterial() const
{
  return new TPZMatPML<TMAT>(*this);
}

template<class TMAT>
void TPZMatPML<TMAT>::ComputeSParameters(const TPZVec<REAL> &x,
                                         CSTATE &sx,
                                         CSTATE &sy,
                                         CSTATE &sz) const
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
  sz = 1;
  static constexpr CSTATE imag{0,1};
  if(fAttX){
    sx = 1. - imag * fAlphaMaxX * ((x[0]-fPmlBeginX) / fDX )
      * ((x[0]-fPmlBeginX) / fDX );
  }
  if(fAttY){
    sy = 1. - imag * fAlphaMaxY * ((x[1]-fPmlBeginY) / fDY ) *
      ((x[1]-fPmlBeginY) / fDY );
  }
  if(fAttZ){
    sz = 1. - imag * fAlphaMaxZ * ((x[2]-fPmlBeginZ) / fDZ ) *
      ((x[2]-fPmlBeginZ) / fDZ );
  }
}

template<class TMAT>
void TPZMatPML<TMAT>::GetPermittivity(
  const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const
{
  TMAT::GetPermittivity(x,er);
  CSTATE sx{1}, sy{1}, sz{1};
  ComputeSParameters(x,sx,sy,sz);

  er[0] *= (sz*sy) / sx;
  er[1] *= (sx*sz) / sy;
  er[2] *= (sy*sx) / sz;
}

template<class TMAT>
void TPZMatPML<TMAT>::GetPermeability(
  const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const
{
  TMAT::GetPermeability(x,ur);
  CSTATE sx{1}, sy{1}, sz{1};
  ComputeSParameters(x,sx,sy,sz);
    
  ur[0] *= (sz*sy) / sx;
  ur[1] *= (sx*sz) / sy;
  ur[2] *= (sy*sx) / sz;
}

template<class TMAT>
int TPZMatPML<TMAT>::ClassId() const
{
  return Hash("TPZMatPML") ^
    TMAT::ClassId() << 1;
}


template<class TMAT>
int TPZSingleSpacePML<TMAT>::IntegrationRuleOrder(const int elPMaxOrder) const
{
  const int integrationorder = 2+2*elPMaxOrder;

  return  integrationorder;
}


template<class TMAT>
int TPZCombinedSpacesPML<TMAT>::IntegrationRuleOrder(const TPZVec<int> &elPMaxOrder) const
{
  int pmax = 0;
  for (int ip=0;  ip<elPMaxOrder.size(); ip++)
    {
      if(elPMaxOrder[ip] > pmax) pmax = elPMaxOrder[ip];
    }

  const int integrationorder = 2+2*pmax;

  return  integrationorder;
}
#endif /* _TPZMATPML_IMPL_H_ */

