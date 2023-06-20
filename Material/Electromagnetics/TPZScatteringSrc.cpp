#include "TPZScatteringSrc.h"
#include <TPZMaterialDataT.h>

using namespace std::complex_literals;

//! Unique identifier for serialization purposes
int TPZScatteredSol3D::ClassId() const
{
  return Hash("TPZScatteredSol3D");
}

void TPZScatteredSol3D::Write(TPZStream &buf, int withclassid) const
{
  buf.Write(sol);
}
  //! Read from stream(serialization method)
void TPZScatteredSol3D::Read(TPZStream &buf, void *context)
{
  buf.Read(sol);
}

void TPZScatteredSol3D::Print(std::ostream &out) const
{
  out << "sol ";
  for(auto s : sol){
    out << s<<" ";
  }
}


void TPZScatteringSrc::FillDataRequirements(
  TPZMaterialData &data) const
{
  data.fNeedsNormal = true;
}

//! Contribution to the integration point
void TPZScatteringSrc::Contribute(const TPZMaterialDataT<CSTATE> &data,
                                     REAL weight, TPZFMatrix<CSTATE> &ef)
{
  //index of integration point
  const int gp_index = data.intGlobPtIndex;
  const auto &mem_item = this->MemItem(gp_index);
  
  const auto &phi_real = data.phi;
  const auto beta = this->Beta();
  const int nshape=phi_real.Rows();


  //making complex version of phi
  TPZFNMatrix<3000,CSTATE> phi(3,nshape);
  for(int i = 0; i < nshape; i++){
    for(int x = 0; x < 3; x++){
      phi.PutVal(x,i,phi_real.GetVal(i,x));
    }
  }
  

  TPZFNMatrix<9,CSTATE> ur_inv;
  GetPermeability(data.x,ur_inv);
  ur_inv.Decompose(ELU);
  
  ///Source term
  const auto src_cte = 2. * 1i * beta * fScaleFactor;
  
  const auto &sol = mem_item.sol;
  TPZFNMatrix<3,CSTATE> src_mat(3,1,0.);
  src_mat(0,0) = -1.*src_cte * sol[0];
  src_mat(1,0) =  1.*src_cte * sol[1];

  //Contribution
  ur_inv.Substitution(&src_mat);
  ef.AddContribution(0, 0, phi, true, src_mat, false, weight);
}

TPZScatteringSrc * TPZScatteringSrc::NewMaterial() const
{
  return new TPZScatteringSrc(*this);
}
//! Unique identifier for serialization purposes
int TPZScatteringSrc::ClassId() const
{
  return
    Hash("TPZScatteringSrc")
    ^
    TPZScattering::ClassId() << 1
    ^
    TPZMatWithMem<TPZScatteredSol3D>::ClassId() << 2;
}
//! Write to stream(serialization method)
void TPZScatteringSrc::Write(TPZStream &buf, int withclassid) const
{
  TPZScattering::Write(buf,withclassid);
  TPZMatWithMem<TPZScatteredSol3D>::Write(buf, withclassid);
  buf.Write(&fBeta);
}
//! Read from stream(serialization method)
void TPZScatteringSrc::Read(TPZStream &buf, void *context)
{
  TPZScattering::Read(buf,context);
  TPZMatWithMem<TPZScatteredSol3D>::Read(buf, context);
  buf.Read(&fBeta);
}

#include "TPZCartesianPML.h"
template class TPZSingleSpaceCartesianPML<TPZScatteringSrc>;