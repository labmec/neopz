//$Id: pzelast3d.cpp,v 1.2 2005-11-28 13:42:06 tiago Exp $

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
#ifdef DEBUG
  if (force.NElements() != 3) PZError << __PRETTY_FUNCTION__ << " - error!" << std::endl;
#endif  
  int i;  
  this->fForce.Resize(3);
  for(i = 0; i < 3; i++) this->fForce[i] = force[i];
  //Default directions is {1,0,0}
  this->fPostProcessDirection.Resize(3);
  this->fPostProcessDirection.Fill(0.);
  this->fPostProcessDirection[0] = 1.;
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
  if(!strcmp("Displacement",name))  return TPZElasticity3D::EDisplacement;
  if(!strcmp("DisplacementX",name))  return TPZElasticity3D::EDisplacementX;
  if(!strcmp("DisplacementY",name))  return TPZElasticity3D::EDisplacementY;
  if(!strcmp("DisplacementZ",name))  return TPZElasticity3D::EDisplacementZ;
  if(!strcmp("PrincipalStress", name))  return TPZElasticity3D::EPrincipalStress;
  if(!strcmp("PrincipalStrain", name))  return TPZElasticity3D::EPrincipalStrain;
  if(!strcmp("VonMises",    name))  return TPZElasticity3D::EVonMisesStress;
  if(!strcmp("Stress",     name))  return TPZElasticity3D::EStress;
  if(!strcmp("Strain",     name))  return TPZElasticity3D::EStrain;
  PZError << "TPZElasticity3D::VariableIndex Error\n";
  return -1;
}

int TPZElasticity3D::NSolutionVariables(int var){
  if(var == TPZElasticity3D::EDisplacement)    return 3;
  if(var == TPZElasticity3D::EDisplacementX)   return 1;
  if(var == TPZElasticity3D::EDisplacementY)   return 1;
  if(var == TPZElasticity3D::EDisplacementZ)   return 1;
  if(var == TPZElasticity3D::EPrincipalStress) return 3;
  if(var == TPZElasticity3D::EPrincipalStrain) return 3;
  if(var == TPZElasticity3D::EVonMisesStress)  return 1;  
  if(var == TPZElasticity3D::EStress)          return 3;  
  if(var == TPZElasticity3D::EStrain)          return 3;  
  PZError << "TPZElasticity3D::NSolutionVariables Error\n";
  return -1;
}

void TPZElasticity3D::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

  if(var == TPZElasticity3D::EDisplacement){
    int i;
    for(i = 0; i < 3; i++){
      Solout[i] = Sol[i];
    }//for
    return;
  }//TPZElasticity3D::EDisplacement
  
  if(var == TPZElasticity3D::EDisplacementX){
    int i;
    Solout[0] = Sol[0];
    return;
  }//TPZElasticity3D::EDisplacementX
  
  if(var == TPZElasticity3D::EDisplacementY){
    int i;
    Solout[0] = Sol[1];
    return;
  }//TPZElasticity3D::EDisplacementY  
  
  if(var == TPZElasticity3D::EDisplacementZ){
    int i;
    Solout[0] = Sol[2];
    return;
  }//TPZElasticity3D::EDisplacementZ  
  
  if(var == TPZElasticity3D::EPrincipalStress){
   
  }//TPZElasticity3D::EPrincipalStress

  
  if(var == TPZElasticity3D::EPrincipalStrain){
  
  }//TPZElasticity3D::EPrincipalStrain

  
  if(var == TPZElasticity3D::EVonMisesStress){
    
  }//TPZElasticity3D::EVonMisesStress
  
  if(var == TPZElasticity3D::EStress){
    TPZFNMatrix<6> Stress(6,1);
    this->ComputeStressVector(Stress, DSol);
    this->ApplyDirection(Stress, Solout);
    return;
  }//TPZElasticity3D::EStress
  
  if(var == TPZElasticity3D::EStrain){
    TPZFNMatrix<6> Strain(6,1);
    this->ComputeStrainVector(Strain, DSol);
    this->ApplyDirection(Strain, Solout);
    return;
  }//TPZElasticity3D::EStrain
  
}//Solution

void TPZElasticity3D::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudx, 
                    TPZFMatrix &axes, TPZVec<REAL> &flux,
                    TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values){
  #warning Implement me, Tiago                    
}

void TPZElasticity3D::ComputeStrainTensor(TPZFMatrix &Strain, TPZFMatrix &DSol){
  Strain.Redim(3,3);
  Strain(0,0) = DSol(0,0);
  Strain(0,1) = 0.5 * ( DSol(1,0) + DSol(0,1) );
  Strain(0,2) = 0.5 * ( DSol(2,0) + DSol(0,2) );
  Strain(1,0) = Strain(0,1);
  Strain(1,1) = DSol(1,1);
  Strain(1,2) = 0.5 * ( DSol(2,1) + DSol(1,2) );
  Strain(2,0) = Strain(0,2);
  Strain(2,1) = Strain(1,2);
  Strain(2,2) = DSol(2,2);
}

void TPZElasticity3D::ComputeStrainVector(TPZFMatrix &Strain, TPZFMatrix &DSol){
  Strain.Redim(6,1);
  Strain(0,0) = DSol(0,0);
  Strain(1,0) = DSol(1,1);
  Strain(2,0) = DSol(2,2);
  Strain(3,0) = 0.5 * ( DSol(1,0) + DSol(0,1) );
  Strain(4,0) = 0.5 * ( DSol(2,0) + DSol(0,2) );
  Strain(5,0) = 0.5 * ( DSol(2,1) + DSol(1,2) );
}

void TPZElasticity3D::ComputeStressVector(TPZFMatrix &Stress, TPZFMatrix &DSol){
  REAL const1 = -1. + this->fPoisson;
  REAL const2 = -1. + this->fPoisson + 2. * fPoisson * fPoisson;
  const REAL E = this->fE;
  const REAL ni = this->fPoisson;
  Stress.Redim(6,1);
  Stress(0,0) = E * ( DSol(0,0) * const1 - ( DSol(1,1) + DSol(2,2) ) * ni ) / const2;
  Stress(1,0) = E * ( DSol(1,1) * const1 - ( DSol(0,0) + DSol(2,2) ) * ni ) / const2;
  Stress(2,0) = E * ( DSol(2,2) * const1 - ( DSol(0,0) + DSol(1,1) ) * ni ) / const2;
  
  REAL const3 = 4. * ( 1. + ni );
  Stress(3,0) = E * ( DSol(1,0) + DSol(0,1) ) / const3;
  Stress(4,0) = E * ( DSol(2,0) + DSol(0,2) ) / const3;
  Stress(5,0) = E * ( DSol(3,0) + DSol(0,3) ) / const3;
}

void TPZElasticity3D::ApplyDirection(TPZFMatrix &StrVec, TPZVec<REAL> &Out){
  Out.Resize(3);
  TPZVec<REAL> &Dir = this->fPostProcessDirection;
  Out[0] = Dir[0] * StrVec(0,0) + Dir[1] * StrVec(3,0) + Dir[2] * StrVec(4,0);
  Out[1] = Dir[1] * StrVec(1,0) + Dir[0] * StrVec(3,0) + Dir[2] * StrVec(5,0);
  Out[2] = Dir[2] * StrVec(2,0) + Dir[0] * StrVec(4,0) + Dir[1] * StrVec(5,0);
}

