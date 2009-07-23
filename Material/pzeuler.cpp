//$Id: pzeuler.cpp,v 1.1 2009-07-23 20:36:41 fortiago Exp $

#include "pzeuler.h"

#include "pzartdiff.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzreal.h"
#include <math.h>
#include "pzstring.h"
#include "pzsave.h"
#include "pzerror.h"


TPZEulerEquation::~TPZEulerEquation(){

}

TPZEulerEquation::TPZEulerEquation(int nummat, REAL gamma) : 
                  TPZDiscontinuousGalerkin(nummat),fAUSMFlux(gamma){
  fGamma = gamma;
}

TPZEulerEquation::TPZEulerEquation():TPZDiscontinuousGalerkin(),fAUSMFlux(-1.){
  fGamma = -1.;
}

TPZEulerEquation::TPZEulerEquation(const TPZEulerEquation &cp) : 
                TPZDiscontinuousGalerkin(cp),fAUSMFlux(cp.fAUSMFlux){
  this->fGamma = cp.fGamma;
}

TPZAutoPointer<TPZMaterial> TPZEulerEquation::NewMaterial(){
  return new TPZEulerEquation(*this);
}

int TPZEulerEquation::NStateVariables() {
  return 5;///U = (rho, rhou, rhov, rhow, rhoe)
}

int TPZEulerEquation::Dimension(){
  return 3;
}

void TPZEulerEquation::Print(std::ostream &out) {
  TPZDiscontinuousGalerkin::Print(out);
  out << "fGamma = " << fGamma << "\n";
}

int TPZEulerEquation::VariableIndex(const std::string &name) {
  if( !strcmp(name.c_str(),"density")  )     return 1;//rho
  if( !strcmp(name.c_str(),"velocity") )     return 2;//(u,v,w)
  if( !strcmp(name.c_str(),"energy")   )     return 3;//E
  if( !strcmp(name.c_str(),"pressure") )     return 4;//p
  if( !strcmp(name.c_str(),"solution") )     return 5;//(ro,u,v,w,E)
  if( !strcmp(name.c_str(),"normvelocity") ) return 6;//sqrt(u+v+w)
  if( !strcmp(name.c_str(),"Mach") )         return 7;//sqrt(u+v+w)/c
  cout << "TPZEulerEquation::VariableIndex not defined\n";
  return TPZMaterial::VariableIndex(name);
}

int TPZEulerEquation::NSolutionVariables(int var){

  if(var == 1 || var == 3 || var == 4 || var == 6 || var == 7) return 1;
  if(var == 2) return Dimension();
  if(var == 5) return NStateVariables();

  cout << "TPZEulerEquation::NSolutionVariables not defined\n";
  return 0;
}

void TPZEulerEquation::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

  if(fabs(Sol[0]) < 1.e-10){
    PZError << "\nTPZEulerEquation::Solution: Density almost null\n"
            << "Density = " << Sol[0] << endl;
  }

  if(var == 1) {
    Solout.Resize(1);
    Solout[0] = Sol[0];//density
    return;
  } else if(var == 2) {
    int dim = Dimension();
    Solout.Resize(dim);
    for(int i=0;i<dim;i++) Solout[i] = Sol[i+1]/Sol[0];//velocity vector
    return;
  } else if(var == 3) {
    Solout.Resize(1);
    int pos = Dimension() + 1;
    Solout[0] = Sol[pos];//energy
    return;
  } else if(var == 4) {
    Solout.Resize(1);
    Solout[0] = Pressure(Sol);//pressure
    return;
  } else if(var == 5) {
    int nstate = NStateVariables();
    Solout.Resize(nstate);
    for(int i=0;i<nstate;i++) Solout[i] = Sol[i];//(ro,ro*u,ro*v,ro*w,E)
    return;
  } else if(var == 6) {
    int nstate = NStateVariables();
    Solout.Resize(1);
    REAL ro2 = Sol[0]*Sol[0];
    REAL veloc = 0.0;
    for(int i=1;i<nstate-1;i++) veloc += Sol[i]*Sol[i];//velocity vector
    Solout[0] = sqrt(veloc/ro2);
    return;
  } else if(var == 7) {
//    int nstate = NStateVariables();
    Solout.Resize(1);
    const REAL cspeed = this->cSpeed(Sol);
    const REAL us = this->uRes(Sol);
    Solout[0] = us / cspeed;
    return;
  } else {
    //cout << "TPZEulerEquation::Solution variable in the base class\n";
    TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
  }
}

void TPZEulerEquation::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
///nothing to be done here
  cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}

void TPZEulerEquation::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){
///nothing to be done here
  cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}

void TPZEulerEquation::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
  this->ContributeInterface(data,weight,ef);
  cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}


void TPZEulerEquation::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ef){
  TPZFNMatrix<15> FluxLeft(5,3), FluxRight(5,3), Flux(5,3), NormalFlux(5,1);;
  this->ComputeFlux(data.soll,FluxLeft);
  this->ComputeFlux(data.solr,FluxRight);

  TPZFNMatrix<3> normalVec(3,1);
  for(int i = 0; i < 3; i++) normalVec(i,0) = data.normal[i];

  #warning fluxo que nao funciona, apenas para desenvolvimento
  Flux = FluxLeft;
  Flux += FluxRight;
  Flux *= 0.5;

  Flux.Multiply(normalVec,NormalFlux);

  for(int i = 0; i < 5; i++) ef(i,0)   += -1. * weight*Flux(i,0);
  for(int i = 0; i < 5; i++) ef(i+5,0) += +1. * weight*Flux(i,0);

}///void

void TPZEulerEquation::ContributeBC(TPZMaterialData &data,
                                    REAL weight,
                                    TPZFMatrix &ek, TPZFMatrix &ef,
                                    TPZBndCond &bc){
  #warning implement it
}

void TPZEulerEquation::ContributeBCInterface(TPZMaterialData &data,
                                              REAL weight,
                                              TPZFMatrix &ek,TPZFMatrix &ef,
                                              TPZBndCond &bc){
#warning implement it
}

void TPZEulerEquation::ContributeBCInterface(TPZMaterialData &data,
                                              REAL weight,
                                              TPZFMatrix &ef,
                                              TPZBndCond &bc){
#warning implement it
}

REAL TPZEulerEquation::Pressure(TPZVec<REAL> &U){

  if(fabs(U[0]) < 1.e-6){
    PZError << "TPZEulerEquation::Pressure - Negative or too small density"
            << U[0] << std::endl;
    DebugStop();
  }

  REAL press = 0.0;

  //U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro w , ro e)
  REAL rho_velocity = ( U[1]*U[1] + U[2]*U[2] + U[3]*U[3] )/U[0];
  press = ((this->fGamma-1.)*( U[4] - REAL(0.5) * rho_velocity ));

  if(press < 0){
    REAL temp = (this->fGamma-1.)*U[4];
    PZError << "TPZEulerEquation::Pressure Negative pressure: " << press << " (gama-1)*E = " << temp << std::endl;
    DebugStop();
  }
  return press;

}///method


REAL TPZEulerEquation::cSpeed(TPZVec<REAL> & sol){

  if(sol[0] < REAL(1e-10)){
      PZError << "TPZEulerEquation(::cSpeed Too small or negative density\n";
      DebugStop();
   }

   const REAL press = this->Pressure(sol);
   const REAL temp = this->fGamma * press;

   if(temp < REAL(1e-10)) // too low or negative
   {
      PZError << "TPZEulerEquation::cSpeed Too low or negative numerator\n";
   }
   const REAL c = sqrt(this->fGamma * press/ sol[0]);
   return c;

}///method

REAL TPZEulerEquation::uRes(TPZVec<REAL> & sol){
  const REAL temp = sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3];
  if(temp < REAL(1e-40)){
    PZError << "TPZEulerEquation::uRes Zero Velocity\n";
    DebugStop();
  }
  const REAL us = sqrt(temp)/sol[0];
  return us;
}

void TPZEulerEquation::ComputeFlux(TPZVec<REAL> &sol, TPZFMatrix & F){
  const double rho = sol[0];
  const double rhoU = sol[1];
  const double rhoV = sol[2];
  const double rhoW = sol[3];
  const double rhoE = sol[4];
  const double u = rhoU/rho;
  const double v = rhoV/rho;
  const double w = rhoW/rho;
  const double p = this->Pressure(sol);
  F.Resize(5,3);
  F(0,0) = rhoU;        F(0,1) = rhoV;        F(0,2) = rhoW;
  F(1,0) = rhoU*u+p;    F(1,1) = rhoU*v;      F(1,2) = rhoU*w;
  F(2,0) = rhoV*u;      F(2,1) = rhoV*v+p;    F(2,2) = rhoV*w;
  F(3,0) = rhoW*u;      F(3,1) = rhoW*v;      F(3,2) = rhoW*w+p;
  F(4,0) = (rhoE+p)*u;  F(4,1) = (rhoE+p)*v;  F(4,2) = (rhoE+p)*w;
}///void
