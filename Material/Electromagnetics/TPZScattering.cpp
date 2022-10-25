#include "TPZScattering.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;


TPZScattering::TPZScattering(int id, const CSTATE er,
                             const CSTATE ur, const STATE lambda,
                             const REAL scale)
  : TBase(id),fLambda(lambda), fScaleFactor(scale)
{
  SetPermeability(ur);
  SetPermittivity(er);
}
TPZScattering::TPZScattering(int id, const TPZVec<CSTATE>& er,
                             const TPZVec<CSTATE> &ur,
                             const STATE lambda,
                             const REAL scale)
  : TBase(id),fLambda(lambda), fScaleFactor(scale)
{
  SetPermeability(ur);
  SetPermittivity(er);
}

TPZScattering * TPZScattering::NewMaterial() const
{
  return new TPZScattering(*this);
}


void TPZScattering::GetPermittivity(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &er) const
{
    er = fEr;
}

void TPZScattering::GetPermeability(
  [[maybe_unused]] const TPZVec<REAL> &x,TPZVec<CSTATE> &ur) const
{
    ur = fUr;
}

void TPZScattering::SetPermeability(CSTATE ur)
{
    if (std::real(ur) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fUr = {ur,ur,ur};
}

void TPZScattering::SetPermeability(const TPZVec<CSTATE>& ur)
{
    if(ur.size()!=3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of ur != 3. Aborting...\n";
        DebugStop();
    }
    for(const auto &iur : ur){
        if (std::real(iur) <0){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Setting negative permeability. Aborting..\n";
            DebugStop();
        }
    }
    fUr = ur;
}
void TPZScattering::SetPermittivity(CSTATE er)
{
    if (std::real(er) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fEr = {er,er,er};
}

void TPZScattering::SetPermittivity(const TPZVec<CSTATE>&er)
{
    if(er.size()!=3){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"\nSize of er != 3. Aborting...\n";
        DebugStop();
    }
    for(const auto &ier : er){
        if (std::real(ier) <0){
            PZError<<__PRETTY_FUNCTION__;
            PZError<<"Setting negative permitivitty. Aborting..\n";
            DebugStop();
        }
    }
    fEr = er;
}

void TPZScattering::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                       REAL weight,
                                       TPZFMatrix<CSTATE> &ek,
                                       TPZFMatrix<CSTATE> &ef)
{
  TPZManVector<CSTATE,3> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  const int nshape = data.phi.Rows();
  const auto &phi_real = data.phi;
  const auto &curl_phi_real = data.curlphi;
  
  const STATE k0 = fScaleFactor * 2*M_PI/fLambda;
  TPZFNMatrix<9,CSTATE> er_mat(3,3,0.);
  TPZFNMatrix<9,CSTATE> ur_inv_mat(3,3,0.);
  for(int x = 0; x < 3; x++){
    er_mat(x,x) = er[x] * k0 * k0;
    ur_inv_mat(x,x) = 1./ur[x];
  }
  
  //making complex version of phi
  TPZFNMatrix<3000,CSTATE> phi(nshape,3);
  TPZFNMatrix<3000,CSTATE> curl_phi(3,nshape);
  for(int i = 0; i < nshape; i++){
    for(int x = 0; x < 3; x++){
      phi.PutVal(i,x,phi_real.GetVal(i,x));
      curl_phi.PutVal(x,i,curl_phi_real.GetVal(x,i));
    }
  }
  //initializing with correct dimensions
  TPZFNMatrix<3000,CSTATE> phi_t(3, nshape);
  TPZFNMatrix<3000,CSTATE> curl_phi_t(nshape, 3);

  phi.Transpose(&phi_t);
  curl_phi.Transpose(&curl_phi_t);
    
  ek += ((curl_phi_t * (ur_inv_mat * curl_phi)) -  phi*(er_mat*phi_t))*weight;
}

void TPZScattering::ContributeBC(const TPZMaterialDataT<CSTATE> &data,
                                         REAL weight,
                                         TPZFMatrix<CSTATE> &ek,
                                         TPZFMatrix<CSTATE> &ef,
                                         TPZBndCondT<CSTATE> &bc)
{
  const auto &phi = data.phi;
  const auto& BIG = TPZMaterial::fBigNumber;
    
  const CSTATE v1 = bc.Val1()(0,0);
  const CSTATE v2 = bc.Val2()[0];
  constexpr STATE tol = std::numeric_limits<STATE>::epsilon();
  if(std::abs(v2) > tol){
    PZError<<__PRETTY_FUNCTION__;
    PZError<<"\nThis method supports only homogeneous boundary conditions.\n";
    std::cout<<"Stopping now..."<<std::endl;
    DebugStop();
  }
  switch ( bc.Type() )
    {
    case 0:{
      const int nshape=phi.Rows();
      for(int i = 0 ; i<nshape ; i++){
        CSTATE load{0};
        for(int x = 0; x < 3; x++){
          load += weight * BIG * v2 * phi(i,x);
        }
        ef(i,0) += load;
        for(int j=0;j<nshape;j++){
          STATE stiff{0};
          for(int x = 0; x < 3; x++){
            stiff += phi(i,x) * phi(j,x) * BIG ;
          }
          ek(i,j) += stiff*weight;
        }
      }
      break;
    }
    case 1:
      ///PMC condition just adds zero to both matrices. nothing to do here....
      break;
    default:
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
      PZError<<"Stopping now..."<<std::endl;
      DebugStop();
      break;
    }
}

int TPZScattering::ClassId() const {
  return Hash("TPZScattering") ^
    TBase::ClassId() << 1;


}

//! Variable index of a given solution
int TPZScattering::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "Field_real") == 0) return 0;
  if( strcmp(name.c_str(), "Field_imag") == 0) return 1;
  if( strcmp(name.c_str(), "Field_abs") == 0) return 2;
  if( strcmp(name.c_str(), "Deriv_real") == 0) return 3;
  if( strcmp(name.c_str(), "Deriv_imag") == 0) return 4;
  if( strcmp(name.c_str(), "Deriv_abs") == 0) return 5;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int TPZScattering::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
  case 1: //field (imag val)
  case 2: //field (abs val)
  case 3://deriv (real part)
  case 4://deriv (imag val)
  case 5://deriv (abs val)
    return this->Dimension();
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void TPZScattering::Solution(const TPZMaterialDataT<CSTATE> &data,
              int var, TPZVec<CSTATE> &solout)
{

  const auto &sol = data.sol[0];
  const auto &curlsol = data.curlsol[0];

  const auto op = [var](auto &val){
    switch(var){
    case 0:
    case 3:
      return std::real(val);
    case 1:
    case 4:
      return std::imag(val);
    case 2:
    case 5:
      return std::abs(val);
    default:
      DebugStop();
      return std::real(val);
    }
  };

  const auto &val = var < 3 ? sol : curlsol;

  for(auto x = 0; x < 3; x++){
    solout[x] = op(val[x]);
  }
}

#include "TPZMatPML.h"
template class TPZSingleSpacePML<TPZScattering>;
