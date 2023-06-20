#ifndef _TPZCARTESIANPML_IMPL_H_
#define _TPZCARTESIANPML_IMPL_H_
#include "TPZCartesianPML.h"
#include "TPZMaterialDataT.h"

template<class TMAT>
void TPZCartesianPML<TMAT>::SetAttX(const REAL pmlBegin,
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
void TPZCartesianPML<TMAT>::SetAttY(const REAL pmlBegin,
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
void TPZCartesianPML<TMAT>::SetAttZ(const REAL pmlBegin,
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
TPZCartesianPML<TMAT> * TPZCartesianPML<TMAT>::NewMaterial() const
{
  return new TPZCartesianPML<TMAT>(*this);
}

template<class TMAT>
void TPZCartesianPML<TMAT>::ComputeSParameters(const TPZVec<REAL> &x,
                                         CSTATE &sx,
                                         CSTATE &sy,
                                         CSTATE &sz) const
{
  /*************CALCULATE S PML PARAMETERS*************/

  static constexpr CSTATE imag{0,1};
  if(fAttX){
    const auto dx = ((x[0]-fPmlBeginX) / fDX );
    sx = 1. - imag * fAlphaMaxX * dx * dx;
  }
  if(fAttY){
    const auto dy = ((x[1]-fPmlBeginY) / fDY );
    sy = 1. - imag * fAlphaMaxY * dy * dy;
      
  }
  if(fAttZ){
    const auto dz = ((x[2]-fPmlBeginZ) / fDZ );
    sz = 1. - imag * fAlphaMaxZ * dz * dz;
      
  }
}

template<class TMAT>
void TPZCartesianPML<TMAT>::GetPermittivity(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const
{
  TMAT::GetPermittivity(x,er);
  CSTATE sx{1}, sy{1}, sz{1};
  ComputeSParameters(x,sx,sy,sz);
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
void TPZCartesianPML<TMAT>::GetPermeability(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const
{
  TMAT::GetPermeability(x,ur);
  CSTATE sx{1}, sy{1}, sz{1};
  ComputeSParameters(x,sx,sy,sz);
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
int TPZCartesianPML<TMAT>::ClassId() const
{
  return Hash("TPZCartesianPML") ^
    TMAT::ClassId() << 1;
}

#endif /* _TPZCARTESIANPML_IMPL_H_ */

