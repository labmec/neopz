#include "TPZPeriodicWgma.h"

#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;

TPZPeriodicWgma *TPZPeriodicWgma::NewMaterial() const {
  return new TPZPeriodicWgma(*this);
}

void TPZPeriodicWgma::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                          REAL weight, TPZFMatrix<CSTATE> &ek,
                                          TPZFMatrix<CSTATE> &ef) {
  TPZFNMatrix<9,CSTATE> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  CSTATE cGx{0}, cGy{0}, cS{0};
  switch (fMode) {
  case ModeType::TE:
    cGx = 1. / ur.GetVal(1,1);
    cGy = 1. / ur.GetVal(0,0);
    cS = er.GetVal(2,2);
    break;
  case ModeType::TM:
    cGx = 1. / er.GetVal(1,1);
    cGy = 1. / er.GetVal(0,0);
    cS =  ur.GetVal(2,2);
    break;
  }
  switch(this->fAssembling){
  case EWhichMatrix::A:
    ContributeInternalA(cGx, cGy, cS, data, weight, ek, ef);
    return;
  case EWhichMatrix::B:
    ContributeInternalB(cGx, cGy, cS, data, weight, ek, ef);
    return;
  case EWhichMatrix::NDefined:
    DebugStop();
    return;
  }
}

void TPZPeriodicWgma::ContributeInternalA(
    const CSTATE cGradX, const CSTATE cGradY, const CSTATE cScal,
    const TPZMaterialDataT<CSTATE> &data, REAL weight, TPZFMatrix<CSTATE> &ek,
    TPZFMatrix<CSTATE> &ef) {
  const auto &phi = data.phi;
  const int nshape = phi.Rows();

  TPZFNMatrix<3, REAL> dphix(3, phi.Rows(), 0.);
  {
    const TPZFMatrix<REAL> &dphidaxes = data.dphix;
    TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, dphix, data.axes);
  }
  const REAL k0 = fScaleFactor * 2*M_PI/fLambda;
  for (int i = 0; i < nshape; i++) {
    for (int j = 0; j < nshape; j++) {
      const STATE gradX = dphix.GetVal(0, i) * dphix.GetVal(0, j);
      const STATE gradY = dphix.GetVal(1, i) * dphix.GetVal(1, j);
      const STATE cross =
        phi.GetVal(i,0) * dphix.GetVal(0, j) - 
        phi.GetVal(j,0) * dphix.GetVal(0, i);
      const STATE phiIphiJ = phi.GetVal(i, 0) * phi.GetVal(j, 0);
      CSTATE stiff = 0;
      stiff -= gradX * cGradX + gradY * cGradY;
      stiff -= 1i *fBeta * cGradX * cross;
      stiff += k0 * k0 * cScal * phiIphiJ;
      ek(i, j) += weight * stiff;
    } // for j
  }   // for i
}

void TPZPeriodicWgma::ContributeInternalB(
    const CSTATE cGradX, const CSTATE cGradY, const CSTATE cScal,
    const TPZMaterialDataT<CSTATE> &data, REAL weight, TPZFMatrix<CSTATE> &ek,
    TPZFMatrix<CSTATE> &ef) {
  const auto &phi = data.phi;
  const int nshape = phi.Rows();

  for (int i = 0; i < nshape; i++) {
    for (int j = 0; j < nshape; j++) {
      const STATE phiIphiJ = phi.GetVal(i, 0) * phi.GetVal(j, 0);
      CSTATE stiff = cGradX * phiIphiJ;
      ek(i, j) += weight * stiff;
    } // for j
  }   // for i
}

void TPZPeriodicWgma::ContributeBC(
    const TPZMaterialDataT<CSTATE> &data, REAL weight, TPZFMatrix<CSTATE> &ek,
    TPZFMatrix<CSTATE> &ef, TPZBndCondT<CSTATE> &bc) {
  switch(this->fAssembling){
  case EWhichMatrix::A:
    ContributeBCInternalA(data, weight, ek, ef, bc);
    return;
  case EWhichMatrix::B://nothing to do here
  case EWhichMatrix::NDefined:
    return;
  }
}

void TPZPeriodicWgma::ContributeBCInternalA(
    const TPZMaterialDataT<CSTATE> &data, REAL weight,
    TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef, TPZBndCondT<CSTATE> &bc) {
  const auto &phi = data.phi;
  const auto &BIG = TPZMaterial::fBigNumber;

  const CSTATE v1 = bc.Val1()(0, 0);
  const CSTATE v2 = bc.Val2()[0];
  constexpr STATE tol = std::numeric_limits<STATE>::epsilon();
  if (std::abs(v2) > tol) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "\nThis method supports only homogeneous boundary conditions.\n";
    std::cout << "Stopping now..." << std::endl;
    DebugStop();
  }
  switch (bc.Type()) {
  case 0: {
    const int nshape = phi.Rows();
    for (int i = 0; i < nshape; i++) {
      ef(i, 0) += weight * BIG * v2 * phi(i, 0);
      for (int j = 0; j < nshape; j++) {
        const STATE stiff = phi(i, 0) * phi(j, 0) * BIG;
        ek(i, j) += stiff * weight;
      }
    }
    break;
  }
  case 1:
    /// PMC condition just adds zero to both matrices. nothing to do here....
    break;
  case 2:
    /// periodic conditions are treated at a mesh level
    break;
  default:
    PZError << __PRETTY_FUNCTION__;
    PZError << "\nThis module supports only dirichlet and neumann boundary "
               "conditions.\n";
    PZError << "Stopping now..." << std::endl;
    DebugStop();
    break;
  }
}

int TPZPeriodicWgma::ClassId() const
{
  return Hash("TPZPeriodicWgma") ^
    TPZScalarField::ClassId() << 1;
}

void TPZPeriodicWgma::Read(TPZStream& buf, void* context)
{
  TPZScalarField::Read(buf,context);
  buf.Read(&fBeta);
}

void TPZPeriodicWgma::Write(TPZStream& buf, int withclassid) const
{
  TPZScalarField::Write(buf,withclassid);
  buf.Write(&fBeta);
}