#include "TPZPlanarWgScattSrc.h"
#include <TPZMaterialDataT.h>

using namespace std::complex_literals;

//! Unique identifier for serialization purposes
int TPZScatteredSol2D::ClassId() const
{
  return Hash("TPZScatteredSol2D");
}

void TPZScatteredSol2D::Write(TPZStream &buf, int withclassid) const
{
  buf.Write(&sol);
  buf.Write(dsol);
}
  //! Read from stream(serialization method)
void TPZScatteredSol2D::Read(TPZStream &buf, void *context)
{
  buf.Read(&sol);
  buf.Read(dsol);
}

void TPZScatteredSol2D::Print(std::ostream &out) const
{
  out << "sol "<<sol<<'\n';
  out<< "dsol ";
  for(auto d : dsol){
    out << d<<" ";
  }
}


void TPZPlanarWgScattSrc::FillDataRequirements(
  TPZMaterialData &data) const
{
  data.fNeedsNormal = true;
}

//! Contribution to the integration point
void TPZPlanarWgScattSrc::Contribute(
  const TPZMaterialDataT<CSTATE> &data, REAL weight,
  TPZFMatrix<CSTATE> &ek, TPZFMatrix<CSTATE> &ef)
{
  //index of integration point
  const int gp_index = data.intGlobPtIndex;
  const auto &mem_item = this->MemItem(gp_index);
  
  const auto &phi = data.phi;
  const auto beta = this->Beta();
  const int nshape=phi.Rows();
  
  const auto sol = mem_item.sol;
  const auto &dsol = mem_item.dsol;
  const auto &n = data.normal;

  TPZManVector<CSTATE,3> er,ur;
  GetPermittivity(data.x,er);
  GetPermeability(data.x,ur);
  CSTATE cGx{0};
  switch(fMode){
  case ModeType::TE:
    cGx = 1./ur[1];
    break;
  case ModeType::TM:
    cGx = 1./er[1];
    break;
  }

  //dsol \cdot normal vec //FOR NOW LET US RESTRICT TO SOURCES ALIGNED WITH Y AXIS
  // CSTATE dsol_n = 0;
  // for(int ix = 0; ix < 3; ix++) { dsol_n += dsol[ix] * n[ix];}
  // std::cout<<"pt"<<std::endl;
  // std::cout<<"\tnormal"<<n<<std::endl;
  // std::cout<<"\tdata.x"<<data.x<<std::endl;
  // std::cout<<"\t sol.x"<<mem_item.x<<std::endl;
  // std::cout<<"\t   sol"<<mem_item.sol<<std::endl;
  // std::cout<<"\tdsol.n"<<dsol_n<<std::endl;
  
  
  ///Source term
  const auto src = 2. *( -1i * beta * sol + dsol[0]);

  //Contribution
  for(int i = 0 ; i<nshape ; i++){
    const CSTATE load = src * cGx * phi(i,0);
    ef(i,0) += weight * load;
  }
}

TPZPlanarWgScattSrc * TPZPlanarWgScattSrc::NewMaterial() const
{
  return new TPZPlanarWgScattSrc(*this);
}
//! Unique identifier for serialization purposes
int TPZPlanarWgScattSrc::ClassId() const
{
  return
    Hash("TPZPlanarWgScattSrc")
    ^
    TPZPlanarWgScatt::ClassId() << 1
    ^
    TPZMatWithMem<TPZScatteredSol2D>::ClassId() << 2;
}
//! Write to stream(serialization method)
void TPZPlanarWgScattSrc::Write(TPZStream &buf, int withclassid) const
{
  TPZPlanarWgScatt::Write(buf,withclassid);
  TPZMatWithMem<TPZScatteredSol2D>::Write(buf, withclassid);
  buf.Write(&fBeta);
}
//! Read from stream(serialization method)
void TPZPlanarWgScattSrc::Read(TPZStream &buf, void *context)
{
  TPZPlanarWgScatt::Read(buf,context);
  TPZMatWithMem<TPZScatteredSol2D>::Read(buf, context);
  buf.Read(&fBeta);
}