#ifndef _TPZCYLINDRICALPML_IMPL_H_
#define _TPZCYLINDRICALPML_IMPL_H_
#include "TPZCylindricalPML.h"
#include "TPZMaterialDataT.h"

template<class TMAT>
void TPZCylindricalPML<TMAT>::SetAttR(const REAL pmlBegin,
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
  fAttR = true;
  fPmlBeginR = pmlBegin;
  fAlphaMaxR = alpha;
  fDR = d;
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::SetAttZ(const REAL pmlBegin,
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
TPZCylindricalPML<TMAT> * TPZCylindricalPML<TMAT>::NewMaterial() const
{
  return new TPZCylindricalPML<TMAT>(*this);
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::ComputeSParameters(const REAL &r,
                                                 const REAL &z,
                                                 CSTATE &sr,
                                                 CSTATE &sz) const
{
  /*************CALCULATE S PML PARAMETERS*************/

  static constexpr CSTATE imag{0,1};
  if(fAttR){
    const auto dr = (r-fPmlBeginR);
    sr = 1. - imag * fAlphaMaxR * dr * dr /(fDR*fDR);
  }
  if(fAttZ){
    const auto dz = ((z-fPmlBeginZ) / fDZ );
    sz = 1. - imag * fAlphaMaxZ * dz * dz;
      
  }
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::GetPermittivity(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const
{
  TMAT::GetPermittivity(x,er);
  CSTATE sr{1}, sz{1};
  const auto r = sqrt(x[0]*x[0]+x[1]*x[1]);
  const auto z = x[2];
  ComputeSParameters(r,z,sr,sz);
  const auto imagsr = sr.imag();
  const auto sx = CSTATE(1. + 1i*imagsr*x[0]/r);
  const auto sy = CSTATE(1. + 1i*imagsr*x[1]/r);
  const auto dets = sx*sy*sz;
  TPZFNMatrix<9,CSTATE> smat(3,3,0.), tmp(3,3,0.);
  smat.PutVal(0,0,sx);
  smat.PutVal(1,1,sy);
  smat.PutVal(2,2,sz);
  smat.Multiply(er,tmp);
  tmp.Multiply(smat,er);
  er *= dets;
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::GetPermeability(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const
{
  TMAT::GetPermeability(x,ur);
  CSTATE sr{1}, sz{1};
  const auto r = sqrt(x[0]*x[0]+x[1]*x[1]);
  const auto z = x[2];
  ComputeSParameters(r,z,sr,sz);
  const auto imagsr = sr.imag();
  const auto sx = CSTATE(1. + 1i*imagsr*x[0]/r);
  const auto sy = CSTATE(1. + 1i*imagsr*x[1]/r);
  const auto dets = sx*sy*sz;
  TPZFNMatrix<9,CSTATE> smat(3,3,0.), tmp(3,3,0.);
  smat.PutVal(0,0,sx);
  smat.PutVal(1,1,sy);
  smat.PutVal(2,2,sz);
  smat.Multiply(ur,tmp);
  tmp.Multiply(smat,ur);
  ur *= dets;
}

template<class TMAT>
int TPZCylindricalPML<TMAT>::ClassId() const
{
  return Hash("TPZCylindricalPML") ^
    TMAT::ClassId() << 1;
}

#endif /* _TPZCYLINDRICALPML_IMPL_H_ */

