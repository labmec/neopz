#include "TPZPlanarWgScatt.h"
#include "TPZMaterialDataT.h"
#include <pzaxestools.h>
using namespace std::complex_literals;

TPZPlanarWgScatt * TPZPlanarWgScatt::NewMaterial() const
{
  return new TPZPlanarWgScatt(*this);
}

void TPZPlanarWgScatt::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                       REAL weight,
                                       TPZFMatrix<CSTATE> &ek,
                                       TPZFMatrix<CSTATE> &ef)
{
  TPZManVector<CSTATE,3> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  CSTATE cGx{0}, cGy{0}, cS{0};
  switch(fMode){
  case ModeType::TE:
    cGx = 1./ur[1];
    cGy = 1./ur[0];
    cS = er[2];
    break;
  case ModeType::TM:
    cGx = 1./er[1];
    cGy = 1./er[0];
    cS = ur[2];
    break;
  }
  const int nshape = data.phi.Rows();
  const auto &phi = data.phi;

  TPZFNMatrix<3,REAL> dphix(3, phi.Rows(), 0.);
  {
    const TPZFMatrix<REAL> &dphidaxes = data.dphix;
    TPZAxesTools<REAL>::Axes2XYZ(dphidaxes, dphix, data.axes);
  }
  
  const STATE k0 = fScaleFactor * 2*M_PI/fLambda;
  
	for(int i = 0; i < nshape; i++){
		for(int j = 0; j < nshape; j++){
      const STATE gradX = dphix.GetVal(0,i) * dphix.GetVal(0,j);
      const STATE gradY = dphix.GetVal(1,i) * dphix.GetVal(1,j);
      const STATE phiIphiJ = phi.GetVal(i,0) * phi.GetVal(j,0);
      CSTATE stiff = 0;
      stiff += gradX * cGx + gradY * cGy;
      stiff -= k0*k0*cS*phiIphiJ;
      ek(i, j) += weight * stiff;
		}//for j
	}//for i
}

void TPZPlanarWgScatt::ContributeBC(const TPZMaterialDataT<CSTATE> &data,
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
    default:
      PZError<<__PRETTY_FUNCTION__;
      PZError<<"\nThis module supports only dirichlet and neumann boundary conditions.\n";
      PZError<<"Stopping now..."<<std::endl;
      DebugStop();
      break;
    }
}

int TPZPlanarWgScatt::ClassId() const {
  return Hash("TPZPlanarWgScatt") ^
    TPZScalarField::ClassId() << 1;


}

#include "TPZMatPML.h"
template class TPZSingleSpacePML<TPZPlanarWgScatt>;
