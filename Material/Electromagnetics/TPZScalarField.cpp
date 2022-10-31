#include "TPZScalarField.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;


TPZScalarField::TPZScalarField(int id, const CSTATE er,const CSTATE ur,
                                             const STATE lambda, const ModeType mode,
                                             const REAL scale) :
  TBase(id), fScaleFactor(scale), fMode(mode)
{
  SetWavelength(lambda);
  SetPermittivity(er);
  SetPermeability(ur);
}

void TPZScalarField::SetWavelength(STATE lambda)
{
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    fLambda = lambda;
}

void TPZScalarField::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const
{
  ur = {fUr,fUr,fUr};
}


 void TPZScalarField::SetPermeability(const CSTATE ur)
{
    if (std::real(ur) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fUr = ur;
}

void TPZScalarField::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const
{
  er = {fEr,fEr,fEr};
}

void TPZScalarField::SetPermittivity(const CSTATE er)
{
    if (std::real(er) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fEr = er;
}


int TPZScalarField::ClassId() const
{
  return Hash("TPZScalarField")
    ^ TBase::ClassId() << 1;
}


static constexpr int intTE{0};
static constexpr int intTM{1};

void TPZScalarField::Read(TPZStream& buf, void* context)
{
  buf.Read(&fUr);
  buf.Read(&fEr);
  buf.Read(&fLambda);
  buf.Read(&fScaleFactor);
  int i = 0;
  buf.Read(&i);
  if(i == intTE){
    fMode = ModeType::TE;
  }else if (i == intTM){
    fMode = ModeType::TM;
  }else{
    DebugStop();
  }
}

void TPZScalarField::Write(TPZStream& buf, int withclassid) const
{
  buf.Write(&fUr);
  buf.Write(&fEr);
  buf.Write(&fLambda);
  buf.Write(&fScaleFactor);

  switch(fMode){
  case ModeType::TE:
    buf.Write(&intTE);
  case ModeType::TM:
    buf.Write(&intTM);
  }
}


//! Variable index of a given solution
int TPZScalarField::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "Field_real") == 0) return 0;
  if( strcmp(name.c_str(), "Field_imag") == 0) return 1;
  if( strcmp(name.c_str(), "Field_abs") == 0) return 2;
  if( strcmp(name.c_str(), "Field_phase") == 0) return 3;
  if( strcmp(name.c_str(), "Deriv_real") == 0) return 4;
  if( strcmp(name.c_str(), "Deriv_imag") == 0) return 5;
  if( strcmp(name.c_str(), "Deriv_abs") == 0) return 6;
  if( strcmp(name.c_str(), "Deriv_phase") == 0) return 7;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int TPZScalarField::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
  case 1: //field (imag val)
  case 2: //field (abs val)
  case 3: //field (phase)
    return 1;
  case 4://deriv (real part)
  case 5://deriv (imag val)
  case 6://deriv (abs val)
  case 7://deriv (phase)
    return this->Dimension();
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void TPZScalarField::Solution(const TPZMaterialDataT<CSTATE> &data,
              int var, TPZVec<CSTATE> &solout)
{
  TPZFNMatrix<3, CSTATE> dsolx(3, 1, 0.);
  const auto dim = this->Dimension();//1 or 2
  if(dim == 2){
    const TPZFMatrix<CSTATE> &dsoldaxes = data.dsol[0];
    TPZAxesTools<CSTATE>::Axes2XYZ(dsoldaxes, dsolx, data.axes);
  }else{
    dsolx = data.dsol[0];
  }
  
  switch (var) {
  case 0:
    solout[0] = std::real(data.sol[0][0]);
    break;
  case 1:
    solout[0] = std::imag(data.sol[0][0]);
    break;
  case 2:
    solout[0] = std::abs(data.sol[0][0]);
    break;
  case 3:
    solout[0] = std::arg(data.sol[0][0]);
    break;
  case 4:
    for(int ix = 0; ix < dim; ix++){
      solout[ix] = std::real(dsolx.GetVal(ix,0));
    }
    break;
  case 5:
    for(int ix = 0; ix < dim; ix++){
      solout[ix] = std::imag(dsolx.GetVal(ix,0));
    }
    break;
  case 6:
    for(int ix = 0; ix < dim; ix++){
      solout[ix] = std::abs(dsolx.GetVal(ix,0));
    }
    break;
  case 7:
    for(int ix = 0; ix < dim; ix++){
      solout[ix] = std::arg(dsolx.GetVal(ix,0));
    }
    break;
  default:
    TBase::Solution(data, var, solout);
  }
}
