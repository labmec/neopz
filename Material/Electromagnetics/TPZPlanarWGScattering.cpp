#include "TPZPlanarWGScattering.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;


TPZPlanarWGScattering::TPZPlanarWGScattering() : TBase()
{
  
}

TPZPlanarWGScattering::TPZPlanarWGScattering(int id, const CSTATE ur,const CSTATE er,
                                             const STATE lambda, const ModeType mode,
                                             const REAL &scale) :
  TBase(id), fScaleFactor(scale), fMode(mode)
{
  SetWavelength(lambda);
  SetPermeability(ur);
  SetPermittivity(er);
}

  

TPZPlanarWGScattering * TPZPlanarWGScattering::NewMaterial() const
{
  return new TPZPlanarWGScattering(*this);
}

void TPZPlanarWGScattering::SetWavelength(STATE lambda)
{
    if (lambda <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative wavelength. Aborting..\n";
        DebugStop();
    }
    fLambda = lambda;
}


TPZVec<CSTATE>
TPZPlanarWGScattering::GetPermeability([[maybe_unused]] const TPZVec<REAL> &x) const
{ return {fUr,fUr,fUr};}

 void TPZPlanarWGScattering::SetPermeability(const CSTATE ur)
{
    if (std::real(ur) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fUr = ur;
}

TPZVec<CSTATE>
TPZPlanarWGScattering::GetPermittivity([[maybe_unused]] const TPZVec<REAL> &x) const
{ return {fEr,fEr,fEr};}

void TPZPlanarWGScattering::SetPermittivity(const CSTATE er)
{
    if (std::real(er) <0){
        PZError<<__PRETTY_FUNCTION__;
        PZError<<"Setting negative permeability. Aborting..\n";
        DebugStop();
    }
    fEr = er;
}

void TPZPlanarWGScattering::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                       REAL weight,
                                       TPZFMatrix<CSTATE> &ek,
                                       TPZFMatrix<CSTATE> &ef)
{
  auto er = GetPermittivity(data.x);
  auto ur = GetPermeability(data.x);
  CSTATE cGx{0}, cGy{0}, cS{0};
  switch(fMode){
  case ModeType::TE:
    cGx = 1./ur[1];
    cGy = 1./ur[0];
    cS = er[2];
    break;
  case ModeType::TM:
    cGx = er[1];
    cGy = er[0];
    cS = 1./ur[2];
    break;
  }
  ContributeInternal(cGx, cGy, cS, data, weight, ek, ef);
}

void
TPZPlanarWGScattering::ContributeInternal(const CSTATE cGradX,
                                          const CSTATE cGradY,
                                          const CSTATE cScal,
                                          const TPZMaterialDataT<CSTATE> &data,
                                          REAL weight,
                                          TPZFMatrix<CSTATE> &ek,
                                          TPZFMatrix<CSTATE> &ef)
{
  const int nshape = data.phi.Rows();
  const auto &phi = data.phi;

  //antigo
  // const auto &dphix = data.dphix;
  //novo
  TPZFNMatrix<3,REAL> dphix(3, phi.Rows(), 0.);
  {
    const TPZFMatrix<REAL> &dphidaxes = data.dphix;
    TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, dphix, data.axes);
  }
  
  const STATE k0 = fScaleFactor * 2*M_PI/fLambda;
  const auto er = GetPermittivity(data.x);
  const auto ur = GetPermeability(data.x);
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
      const STATE gradX = dphix.GetVal(0,i) * dphix.GetVal(0,j);
      const STATE gradY = dphix.GetVal(1,i) * dphix.GetVal(1,j);
      const STATE phiIphiJ = phi.GetVal(i,0) * phi.GetVal(j,0);
      CSTATE stiff = 0;
      stiff += gradX * cGradX + gradY * cGradY;
      stiff -= k0*k0*cScal*phiIphiJ;
      ek(i, j) += weight * stiff;
		}//for j
	}//for i

}

void TPZPlanarWGScattering::ContributeBC(const TPZMaterialDataT<CSTATE> &data,
                                         REAL weight,
                                         TPZFMatrix<CSTATE> &ek,
                                         TPZFMatrix<CSTATE> &ef,
                                         TPZBndCondT<CSTATE> &bc)
{
  CSTATE cGradX{-1};
  switch(fMode){
  case ModeType::TE:
    cGradX = 1./GetPermeability(data.x)[1];
    break;
  case ModeType::TM:
    cGradX = GetPermittivity(data.x)[1];
    break;
  }
   ContributeBCInternal(cGradX, data, weight, ek, ef, bc);
}


void
TPZPlanarWGScattering::ContributeBCInternal(const CSTATE coeffGradX,
                                            const TPZMaterialDataT<CSTATE> &data,
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
        ef(i,0) += weight * BIG * v2 * phi(i,0);
        for(int j=0;j<nshape;j++){
          const STATE stiff = phi(i,0) * phi(j,0) * BIG ;
          ek(i,j) += stiff*weight;
        }
      }
      break;
    }
    case 1:
      ///PMC condition just adds zero to both matrices. nothing to do here....
      break;
    case 2:{
      ///Source term
      TPZFNMatrix<9,CSTATE> dummy;
      TPZManVector<CSTATE,2> res(2);
      bc.ForcingFunctionBC()(data.x,res,dummy);
      const auto Am = res[0];//amplitude
      const auto beta = res[1];//beta
      const int nshape=phi.Rows();
      for(int i = 0 ; i<nshape ; i++){
        const CSTATE load = (phi(i,0) * coeffGradX) * -2i * beta * Am;
        ef(i,0) += weight * load;
      }
      break;
    }
    default:
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
      PZError<<"Stopping now..."<<std::endl;
      DebugStop();
      break;
    }
}
//! Variable index of a given solution
int TPZPlanarWGScattering::VariableIndex(const std::string &name) const
{
  if( strcmp(name.c_str(), "Field_re") == 0) return 0;
  if( strcmp(name.c_str(), "Field_abs") == 0) return 1;
  return TPZMaterial::VariableIndex(name);
}
//! Number of variables associated with a given solution
int TPZPlanarWGScattering::NSolutionVariables(int var) const
{
  switch(var){
  case 0: //field (real part)
    return 1;
  case 1: //field (abs val)
    return 1;
  default:
    return TPZMaterial::NSolutionVariables(var);
  }
}
//! Computes the solution at an integration point
void TPZPlanarWGScattering::Solution(const TPZMaterialDataT<CSTATE> &data,
              int var, TPZVec<CSTATE> &solout)
{
  switch(var){
  case 0:
    solout[0] = std::real(data.sol[0][0]);
    break;
  case 1:
    solout[0] = std::abs(data.sol[0][0]);
    break;
  default:
    TBase::Solution(data,var,solout);
  } 
}