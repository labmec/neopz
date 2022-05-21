#include "TPZPeriodicWgma.h"

#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;

TPZPeriodicWgma::TPZPeriodicWgma() : TBase() {}

TPZPeriodicWgma::TPZPeriodicWgma(
    int id, int dim, const CSTATE er, const CSTATE ur, const STATE lambda,
    const ModeType mode, const REAL &scale)
  : TBase(id), fScaleFactor(scale), fMode(mode), fDim(dim) {
  SetWavelength(lambda);
  SetPermittivity(er);
  SetPermeability(ur);
}

TPZPeriodicWgma *TPZPeriodicWgma::NewMaterial() const {
  return new TPZPeriodicWgma(*this);
}

void TPZPeriodicWgma::SetWavelength(STATE lambda) {
  if (lambda < 0) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "Setting negative wavelength. Aborting..\n";
    DebugStop();
  }
  fLambda = lambda;
}

void TPZPeriodicWgma::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const {
  ur = {fUr, fUr, fUr};
}

void TPZPeriodicWgma::SetPermeability(const CSTATE ur) {
  if (std::real(ur) < 0) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "Setting negative permeability. Aborting..\n";
    DebugStop();
  }
  fUr = ur;
}

void TPZPeriodicWgma::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE>&er) const {
  er = {fEr, fEr, fEr};
}

void TPZPeriodicWgma::SetPermittivity(const CSTATE er) {
  if (std::real(er) < 0) {
    PZError << __PRETTY_FUNCTION__;
    PZError << "Setting negative permeability. Aborting..\n";
    DebugStop();
  }
  fEr = er;
}

void TPZPeriodicWgma::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                          REAL weight, TPZFMatrix<CSTATE> &ek,
                                          TPZFMatrix<CSTATE> &ef) {
  TPZManVector<CSTATE,3> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  CSTATE cGx{0}, cGy{0}, cS{0};
  switch (fMode) {
  case ModeType::TE:
    cGx = 1. / ur[1];
    cGy = 1. / ur[0];
    cS = er[2];
    break;
  case ModeType::TM:
    cGx = er[1];
    cGy = er[0];
    cS = 1. / ur[2];
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
      const STATE cross = dphix.GetVal(0, i) * phi.GetVal(j,0) -
        dphix.GetVal(0, j) * phi.GetVal(i,0);
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
  CSTATE cGradX{-1};
  TPZManVector<CSTATE,3> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  switch (fMode) {
  case ModeType::TE:
    cGradX = 1. / ur[1];
    break;
  case ModeType::TM:
    cGradX = er[1];
    break;
  }
  switch(this->fAssembling){
  case EWhichMatrix::A:
    ContributeBCInternalA(cGradX, data, weight, ek, ef, bc);
    return;
  case EWhichMatrix::B://nothing to do here
  case EWhichMatrix::NDefined:
    return;
  }
}

void TPZPeriodicWgma::ContributeBCInternalA(
    const CSTATE coeffGradX, const TPZMaterialDataT<CSTATE> &data, REAL weight,
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
  case 2:// periodic conditions are treated at a mesh level
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
//! Variable index of a given solution
int TPZPeriodicWgma::VariableIndex(const std::string &name) const {
  if (strcmp(name.c_str(), "Field_re") == 0)
    return 0;
  if (strcmp(name.c_str(), "Field_abs") == 0)
    return 1;
  if (strcmp(name.c_str(), "Deriv_re") == 0)
    return 2;
  if (strcmp(name.c_str(), "Deriv_abs") == 0)
    return 3;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int TPZPeriodicWgma::NSolutionVariables(int var) const {
  switch(var){
  case 0: //field (real part)
  case 1: //field (abs val)
    return 1;
  case 2://deriv (real part)
  case 3://deriv (abs val)
    return 2;
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void TPZPeriodicWgma::Solution(const TPZMaterialDataT<CSTATE> &data,
                                        int var, TPZVec<CSTATE> &solout) {
  switch (var) {
  case 0:
    solout[0] = std::real(data.sol[0][0]);
    break;
  case 1:
    solout[0] = std::abs(data.sol[0][0]);
    break;
  case 2:
    solout[0] = std::real(data.dsol[0][0]);
    solout[1] = std::real(data.dsol[0][1]);
    break;
  case 3:
    solout[0] = std::abs(data.dsol[0][0]);
    solout[1] = std::abs(data.dsol[0][1]);
    break;
  default:
    TBase::Solution(data, var, solout);
  }
}