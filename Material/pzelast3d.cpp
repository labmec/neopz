
#include "pzelast3d.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pztempmat.h"
#include "pzmanvector.h"
#include <math.h>
#include <fstream>

TPZElasticity3D::TPZElasticity3D(int nummat, REAL E, REAL poisson, TPZVec<REAL> &force) : TPZMaterial(nummat){
  this->fE = E;
  this->fPoisson = poisson;
  this->fForce = force;
}//method

TPZElasticity3D::~TPZElasticity3D(){}

void TPZElasticity3D::Print(std::ostream & out){
  out << "\nTPZElasticity3D material:\n";
  out << "\tfE       = " << this->fE << std::endl;
  out << "\tfPoisson = " << this->fPoisson << std::endl;
  out << "\tfForce   = " << this->fForce << std::endl;
  out << "\tBase class print\n";
  TPZMaterial::Print(out);
  out << "End of TPZElasticity3D::Print\n";
}

void TPZElasticity3D::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,REAL weight,
                        TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef){
  
  const int phr = phi.Rows();
  if(this->fForcingFunction){
    this->fForcingFunction(x,this->fForce);
  }
  
  //this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
                          //{dvdy*dudx, dvdy*dudy, dvdy*dudz},
                          //{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
  TPZFNMatrix<9> Deriv(3,3);
  const REAL E  = this->fE;
  const REAL nu = this->fPoisson;
  const REAL C1 = E / (4.+4.*nu);
  const REAL C2 = E * nu / (-1. + nu + 2.*nu*nu);
  const REAL C3 = E * (nu - 1.) / (-1. + nu +2. * nu * nu);
  
  int in;
  for(in = 0; in < phr; in++) {
    int kd;
    for(kd = 0; kd < 3; kd++){
      ef(in*3+kd, 0) += weight* fForce[kd] * phi(in,0);
    }//kd
    REAL val;
    for( int jn = 0; jn < phr; jn++ ) {
    
      //Compute Deriv matrix
      for(int ud = 0; ud < 3; ud++){
        for(int vd = 0; vd < 3; vd++){
          Deriv(vd,ud) = dphi(vd,in)*dphi(ud,jn);
        }//ud
      }//vd
      
      //First equation Dot[Sigma1, gradV1]
      val = ( Deriv(1,1) + Deriv(2,2) ) * C1 + Deriv(0,0) * C3;
      ek(in*3+0,jn*3+0) += weight * val;
      
      val = Deriv(1,0) * C1 - Deriv(0,1) * C2;
      ek(in*3+0,jn*3+1) += weight * val;
      
      val = Deriv(2,0) * C1 - Deriv(0,2) * C2;
      ek(in*3+0,jn*3+2) += weight * val;
           
      //Second equation Dot[Sigma2, gradV2]
      val = Deriv(0,1) * C1 - Deriv(1,0) * C2;
      ek(in*3+1,jn*3+0) += weight * val;
      
      val = ( Deriv(0,0) + Deriv(2,2) ) * C1 + Deriv(1,1) * C3;
      ek(in*3+1,jn*3+1) += weight * val;
      
      val = Deriv(2,1) * C1 - Deriv(1,2) * C2;
      ek(in*3+1,jn*3+2) += weight * val;
      
      //Third equation Dot[Sigma3, gradV3]
      val = Deriv(0,2) * C1 - Deriv(2,0) * C2;
      ek(in*3+2,jn*3+0) += weight * val;
      
      val = Deriv(1,2) * C1 - Deriv(2,1) * C2;
      ek(in*3+2,jn*3+1) += weight * val;
      
      val = ( Deriv(0,0) + Deriv(1,1) ) * C1 + Deriv(2,2) * C3;
      ek(in*3+2,jn*3+2) += weight * val;
      
    }//jn
  }//in
  
#ifdef DEBUG   
   if ( !ek.VerifySymmetry( 1.e-8 ) ) PZError << __PRETTY_FUNCTION__ << "\nERROR - NON SYMMETRIC MATRIX" << std::endl;
#endif
}//method

void TPZElasticity3D::ContributeBC(TPZVec<REAL> &x,TPZVec<REAL> &sol,REAL weight,
                          TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
  const REAL BIGNUMBER  = 1.e18;

  const int phr = phi.Rows();
  int in,jn;
  REAL v2[3];
  v2[0] = bc.Val2()(0,0);
  v2[1] = bc.Val2()(1,0);
  v2[2] = bc.Val2()(2,0);

  switch (bc.Type()) {
  case 0: // Dirichlet condition
    for(in = 0 ; in < phr; in++) {
      ef(3*in+0,0) += BIGNUMBER * v2[0] * phi(in,0) * weight;
      ef(3*in+1,0) += BIGNUMBER * v2[1] * phi(in,0) * weight;        
      ef(3*in+2,0) += BIGNUMBER * v2[2] * phi(in,0) * weight;        
      
      for (jn = 0 ; jn < phr; jn++) {
        ek(3*in+0,3*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
        ek(3*in+1,3*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
        ek(3*in+2,3*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
      }//jn
    }//in
    break;

  case 1: // Neumann condition
    for(in = 0 ; in < phi.Rows(); in++) {
      ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
      ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
      ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
    }//in
    break;
    
  default:
    PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
  }//switch
                          
}//method

int TPZElasticity3D::VariableIndex(char *name){

}

int TPZElasticity3D::NSolutionVariables(int var){

}

void TPZElasticity3D::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

}

void TPZElasticity3D::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudx, 
                    TPZFMatrix &axes, TPZVec<REAL> &flux,
                    TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values){
                    
}
