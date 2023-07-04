#ifndef _TPZCYLINDRICALPML_IMPL_H_
#define _TPZCYLINDRICALPML_IMPL_H_
#include "TPZCylindricalPML.h"
#include "TPZMaterialDataT.h"
using namespace std::complex_literals;

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
void TPZCylindricalPML<TMAT>::ComputeTransformMat(TPZFMatrix<CSTATE> &mat,
                           const REAL r,
                           const REAL phi,
                           const REAL z) const
{
  /*************CALCULATE S PML PARAMETERS*************/
  CSTATE sr{1};
  CSTATE sphi{1};
  CSTATE sz{1};
  constexpr CSTATE imag = 1i;
  if(fAttR){
    const CSTATE ar = this->fAlphaMaxR*(1.+0.00i);
    const auto &rmin = this->fPmlBeginR;
    const auto &dr = this->fDR;
    const CSTATE diffr = r-rmin;
    sr = 1.- imag*ar * diffr*diffr/(dr*dr);
    const CSTATE rt = r - imag*ar*diffr*diffr*diffr/(3*dr*dr);
    sphi = rt/r;
  }
  if(fAttZ){
    const auto dz = ((z-fPmlBeginZ) / fDZ );
    sz = 1. - imag * fAlphaMaxZ * dz * dz;
  }
  mat.Redim(3,3);
  mat.Put(0,0,1./sr);
  mat.Put(1,1,1./sphi);
  mat.Put(2,2,1./sz);
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::GetPermittivity(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &er) const
{
  const auto rho = sqrt(x[0]*x[0]+x[1]*x[1]);
  const auto phi = std::atan2(x[1],x[0]);
  const auto z = x[2];
  TPZFNMatrix<9,CSTATE> smat,t1,t2;
  //this is actually J^-1, not J (just to avoid inverting it)
  this->ComputeTransformMat(smat,rho,phi,z);
  
  const CSTATE detmatinv = smat.Get(0,0)*smat.Get(1,1)*smat.Get(2,2);
  TMAT::GetPermittivity(x,er);
  smat.Multiply(er,t1);
  t1.Multiply(smat,er);
  er *= 1./detmatinv;
  //now we convert back to cartesian coords
  auto RotationMatrix = [](TPZFMatrix<CSTATE> &mat, STATE theta){
    mat.Redim(3,3);
    mat.Put(0,0, std::cos(theta));
    mat.Put(0,1,-std::sin(theta));
    mat.Put(1,0, std::sin(theta));
    mat.Put(1,1, std::cos(theta));
    mat.Put(2,2,1);
  };
  RotationMatrix(t1,phi);
  t1.Multiply(er,t2);
  RotationMatrix(t1,-phi);
  t2.Multiply(t1,er);
}

template<class TMAT>
void TPZCylindricalPML<TMAT>::GetPermeability(
  const TPZVec<REAL> &x,TPZFMatrix<CSTATE> &ur) const
{
  const auto rho = sqrt(x[0]*x[0]+x[1]*x[1]);
  const auto phi = std::atan2(x[1],x[0]);
  const auto z = x[2];
  TPZFNMatrix<9,CSTATE> smat,t1,t2;
  //this is actually J^-1, not J (just to avoid inverting it)
  this->ComputeTransformMat(smat,rho,phi,z);
  
  const CSTATE detmatinv = smat.Get(0,0)*smat.Get(1,1)*smat.Get(2,2);
  TMAT::GetPermeability(x,ur);
  smat.Multiply(ur,t1);
  t1.Multiply(smat,ur);
  ur *= 1./detmatinv;
  //now we convert back to cartesian coords
  auto RotationMatrix = [](TPZFMatrix<CSTATE> &mat, STATE theta){
    mat.Redim(3,3);
    mat.Put(0,0, std::cos(theta));
    mat.Put(0,1,-std::sin(theta));
    mat.Put(1,0, std::sin(theta));
    mat.Put(1,1, std::cos(theta));
    mat.Put(2,2,1);
  };
  RotationMatrix(t1,phi);
  t1.Multiply(ur,t2);
  RotationMatrix(t1,-phi);
  t2.Multiply(t1,ur);
}

template<class TMAT>
int TPZCylindricalPML<TMAT>::ClassId() const
{
  return Hash("TPZCylindricalPML") ^
    TMAT::ClassId() << 1;
}

#endif /* _TPZCYLINDRICALPML_IMPL_H_ */

