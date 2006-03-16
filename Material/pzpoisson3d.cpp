// -*- c++ -*-

//$Id: pzpoisson3d.cpp,v 1.20 2006-03-16 01:45:54 tiago Exp $

#include "pzpoisson3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

using namespace std;
REAL TPZMatPoisson3d::gAlfa = 0.5;

TPZMatPoisson3d::TPZMatPoisson3d(int nummat, int dim) : TPZDiscontinuousGalerkin(nummat), fXf(0.), fDim(dim), fSD(1.) {
  fK = 1.;
  fC = 0.;
  fConvDir[0] = 1.;
  fConvDir[1] = 0.;
  fConvDir[2] = 0.;
  fPenaltyConstant = 1.;
  this->SetNonSymmetric();
  this->SetRightK(fK);
}

void TPZMatPoisson3d::SetParameters(REAL diff, REAL conv, TPZVec<REAL> &convdir) {
  fK = diff;
  this->SetRightK(fK);
  fC = conv;
  int d;
  for(d=0; d<fDim; d++) fConvDir[d] = convdir[d];
}

void TPZMatPoisson3d::GetParameters(REAL &diff, REAL &conv, TPZVec<REAL> &convdir) {
  diff = fK;
  conv = fC;
  int d;
  for(d=0; d<fDim; d++) convdir[d] = fConvDir[d];
}

TPZMatPoisson3d::~TPZMatPoisson3d() {
}

int TPZMatPoisson3d::NStateVariables() {
  return 1;
}

void TPZMatPoisson3d::Print(ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "Laplace operator multiplier fK "<< fK << endl;
  out << "Laplace operator multiplier fK of right neighbour " << this->fRightK << endl;
  out << "Convection coeficient fC " << fC << endl;
  out << "Convection direction " << fConvDir[0] << ' ' << fConvDir[1] << ' ' <<  fConvDir[2] << endl;
  out << "Forcing vector fXf " << fXf << endl;
  out << "Penalty constant " << fPenaltyConstant << endl;
  out << "Base Class properties : \n";
  TPZMaterial::Print(out);
}

void TPZMatPoisson3d::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &  dsol ,REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {

  int phr = phi.Rows();

  if(fForcingFunction) {            // phi(in, 0) = phi_in
    TPZManVector<REAL> res(1);
    fForcingFunction(x,res);       // dphi(i,j) = dphi_j/dxi
    fXf = res[0];
  }
  REAL delx = 0.;
  REAL ConvDirAx[3] = {0.};
  if(fC != 0.0) {
    int di,dj;
    delx = 0.;
    for(di=0; di<fDim; di++) {
      for(dj=0; dj<fDim; dj++) {
        delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
      }
    }
    delx = 2./delx;
      
    
    switch(fDim) {
      case 1:
//        delx = jacinv(0,0);
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        break;
      case 2:
//        delx = jacinv(0,0)*jacinv(1,1)-jacinv(1,0)*jacinv(0,1);
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        break;
      case 3:
//        delx = jacinv(0,0)*jacinv(1,1)*jacinv(2,2)+
//          jacinv(0,1)*jacinv(1,2)*jacinv(2,0)+
//          jacinv(1,0)*jacinv(2,1)*jacinv(0,2)-
//          jacinv(2,0)*jacinv(1,1)*jacinv(0,2)-
//          jacinv(1,0)*jacinv(0,1)*jacinv(2,2)-
//          jacinv(2,1)*jacinv(2,1)*jacinv(0,0);
        ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
        ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
        ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
        break;
      default:
        PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
    }
  }
    
  //Equacao de Poisson
  for( int in = 0; in < phr; in++ ) {
    int kd;
    REAL dphiic = 0;
    for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*dphi(kd,in);
    ef(in, 0) += - weight * ( fXf*phi(in,0) 
//			      + 0.5*fSD*delx*fC*dphiic*fXf
			      );
    for( int jn = 0; jn < phr; jn++ ) {
      for(kd=0; kd<fDim; kd++) {
        ek(in,jn) += weight * (
          +fK * ( dphi(kd,in) * dphi(kd,jn) ) 
          -fC * ( ConvDirAx[kd]* dphi(kd,in) * phi(jn) )
//	  +0.5 * fSD * delx * fC * dphiic * dphi(kd,jn)* ConvDirAx[kd]
          );
      }
    }
  }
    
    if (this->IsSymetric()){    
      if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
    
}


void TPZMatPoisson3d::ContributeBC(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*sol*/,REAL weight,
				     TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {

  int phr = phi.Rows();
  short in,jn;
  REAL v2[1];
  v2[0] = bc.Val2()(0,0);

  switch (bc.Type()) {
  case 0 :			// Dirichlet condition
    for(in = 0 ; in < phr; in++) {
      ef(in,0) += gBigNumber * v2[0] * phi(in,0) * weight;
      for (jn = 0 ; jn < phr; jn++) {
	ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
      }
    }
    break;
  case 1 :			// Neumann condition
    for(in = 0 ; in < phi.Rows(); in++) {
      ef(in,0) += v2[0] * phi(in,0) * weight;
    }
    break;
  case 2 :		// condi�o mista
    for(in = 0 ; in < phi.Rows(); in++) {
      ef(in, 0) += v2[0] * phi(in, 0) * weight;
      for (jn = 0 ; jn < phi.Rows(); jn++) {
	ek(in,jn) += bc.Val1()(0,0) * phi(in,0) *
	  phi(jn,0) * weight;     // peso de contorno => integral de contorno
      }
    }
  }
  
  if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
    if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
  }
}

/** returns the variable index associated with the name*/
int TPZMatPoisson3d::VariableIndex(char *name){
  if(!strcmp("Displacement6",name))   return  0;
  if(!strcmp("Solution",name))        return  1;
  if(!strcmp("Derivate",name))        return  2;
  if(!strcmp("KDuDx",name))           return  3;
  if(!strcmp("KDuDy",name))           return  4;
  if(!strcmp("KDuDz",name))           return  5;
  if(!strcmp("NormKDu",name))         return  6;
  if(!strcmp("MinusKGradU",name))     return  7;
  if(!strcmp("POrder",name))          return 10;
  cout << "TPZMatPoisson3d::VariableIndex Error\n";
  return -1;
}

int TPZMatPoisson3d::NSolutionVariables(int var){

  if(var == 0) return 6;
  if(var == 1) return 1;
  if(var == 2) return fDim;
  if ((var == 3) || (var == 4) || (var == 5) || (var == 6)) return 1;
  if (var == 7) return fDim;
  if(var == 10) return 1;
  return TPZMaterial::NSolutionVariables(var);
  //  cout << "TPZMatPoisson3d::NSolutionVariables Error\n";
  //  return 0;
}

void TPZMatPoisson3d::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){

  if(var == 0 || var == 1) Solout[0] = Sol[0];//function
  if(var == 2) {
    int id;
    for(id=0 ; id<fDim; id++) {
      Solout[id] = DSol(id,0);//derivate
    }
  }//var == 2
  if (var == 3){ //KDuDx
    Solout[0] = DSol(0,0) * this->fK;
  }//var ==3
  if (var == 4){ //KDuDy
    Solout[0] = DSol(1,0) * this->fK;
  }//var == 4 
  if (var == 5){ //KDuDz
    Solout[0] = DSol(2,0) * this->fK;
  }//var == 5
  if (var == 6){ //NormKDu
    int id;
    REAL val = 0.;
    for(id=0 ; id<fDim; id++) {
      val += (DSol(id,0) * this->fK) * (DSol(id,0) * this->fK);
    }
    Solout[0] = sqrt(val);
  }//var == 6
  if (var == 7){ //MinusKGradU
    int id;
    REAL val = 0.;
    for(id=0 ; id<fDim; id++) {
      Solout[id] = -1. * this->fK * DSol(id,0);
    }
    Solout[0] = sqrt(val);
  }//var == 7  
}//method

void TPZMatPoisson3d::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
  //Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

//ofstream ErroFile("erro.txt");
void TPZMatPoisson3d::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

  TPZManVector<REAL> sol(1),dsol(3,0.);
  Solution(u,dudx,axes,1,sol);
  Solution(u,dudx,axes,2,dsol);
  REAL dx[3] = {0.};
  int id,jd;
  for(id=0; id<fDim; id++) for(jd=0; jd<3; jd++) dx[jd] += dsol[id]*axes(id,jd);
  //values[1] : eror em norma L2
  values[1]  = (sol[0] - u_exact[0])*(sol[0] - u_exact[0]);
  //values[2] : erro em semi norma H1
  values[2] = 0.;
  for(id=0; id<fDim; id++) {
    values[2]  += fK*(dx[id] - du_exact(id,0))*(dx[id] - du_exact(id,0));
  }
  //values[0] : erro em norma H1 <=> norma Energia
  values[0]  = values[1]+values[2];
  
//  ErroFile << x[0] << "  " << x[1] << "  " << u_exact[0] << "  " << sol[0] << "  " <<  du_exact(0,0) << "  " << dsol[0] << "  " << du_exact(1,0) << "  " << dsol[1] << endl;
}


void TPZMatPoisson3d::InterfaceJumps(TPZVec<REAL> &x, TPZVec<REAL> &leftu,  TPZVec<REAL> &leftNormalDeriv,
                    TPZVec<REAL> &rightu, TPZVec<REAL> &rightNormalDeriv,
                    TPZVec<REAL> &values){
  values.Resize(3);
  values[1]  = (leftu[0] - rightu[0]) * (leftu[0] - rightu[0]);

  values[2] = (leftNormalDeriv[0] - rightNormalDeriv[0]) * (leftNormalDeriv[0] - rightNormalDeriv[0]);

  values[0]  = values[1]+values[2];
}//method

void TPZMatPoisson3d::BCInterfaceJumps(TPZVec<REAL> &leftu,
                                       TPZBndCond &bc,
                                       TPZVec<REAL> &values){
  if(bc.Type() == 0){ //DIRICHLET
    REAL f = bc.Val2()(0,0);
    values.Resize(3);
    values[1]  = (leftu[0] - f) * (leftu[0] - f);  
    values[2] = 0.0; //it is zero because it is not possible 
                     //to be computed since only u is restrained 
                     //by Dirichlet BC in that formulation/implementation
    values[0]  = values[1] + values[2];    
  }//if                 
}//method

#ifdef _AUTODIFF
void TPZMatPoisson3d::ContributeEnergy(TPZVec<REAL> &x,
			      TPZVec<FADFADREAL> &sol,
			      TPZVec<FADFADREAL> &dsol,
			      FADFADREAL &U,
			      REAL weight)
{
      int dim = dsol.NElements()/sol.NElements();

      //Equa�o de Poisson

      int i, eqs = dsol.NElements()/dim;
      if(sol.NElements() != 1) PZError << "";

//cout << "FADREAL init : \n" << FADREAL(weight * fXf(0,0));

      //FADFADREAL Buff;

      U+= sol[0] * FADREAL(weight * fXf);

      switch(dim)
      {
      case 1:
             U+=fK*(dsol[0] * dsol[0])*FADREAL(weight/2.); // U=((du/dx)^2)/2

	 break;
      case 2:
             U+=fK*(dsol[0] * dsol[0] +
	         dsol[1] * dsol[1])*(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2
             /*Buff  = dsol[0] * dsol[0];
             Buff += dsol[1] * dsol[1];
	     U += Buff * FADREAL(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2*/
	 break;
      case 3:
             U+=fK*(dsol[0] * dsol[0] +
                 dsol[1] * dsol[1] +
	         dsol[2] * dsol[2])*(weight/2.); // U=((du/dx)^2+(du/dy)^2+(du/dz)^2)/2*/
             /*Buff  = dsol[0] * dsol[0];
             Buff += dsol[1] * dsol[1];
             Buff += dsol[2] * dsol[2];
	     U += Buff * FADREAL(weight/2.); //  U=((du/dx)^2+(du/dy)^2+(du/dz)^2)/2*/
	 break;
      }
      //cout << "\nCalcEnergy\n" << U;

}

void TPZMatPoisson3d::ContributeBCEnergy(TPZVec<REAL> & x,TPZVec<FADFADREAL> & sol, FADFADREAL &U, REAL weight, TPZBndCond &bc)
{
  //int i, phr=sol[0].size();

  FADFADREAL solMinBC = sol[0] - FADREAL( bc.Val2()(0,0) );

//cout << "\nsolution " << sol[0];

  switch (bc.Type()) {
  case 0 :	// Dirichlet condition
    // U += 1/2* Big * weight * Integral((u - u0)^2 dOmega)
    U += (solMinBC * solMinBC) * FADREAL(weight * gBigNumber / 2.);
    break;
  case 1 :	// Neumann condition
    // U -= weight * Integral([g].u dOmega)
    U -= sol[0] * FADREAL( bc.Val2()(0,0) * weight);
    break;
  case 2 :	// condi�o mista
    // U += 1/2 * weight * Integral(<(u-u0), [g].(u-u0)> dOmega)
    U += ( solMinBC * /*scalar*/ FADREAL(bc.Val1()(0,0)) * /*matrix oprt*/ solMinBC ) * FADREAL(weight / 2.);
    break;

  }
}

#endif

/** Termos de penalidade. */
void TPZMatPoisson3d::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                          TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                          TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                          TPZFMatrix &ek,TPZFMatrix &ef, int LeftPOrder, int RightPOrder, REAL faceSize){

  this->ContributeInterface( x, solL, solR, dsolL, dsolR, weight, normal, phiL, phiR, dphiL, dphiR, ek, ef );
  
  REAL leftK, rightK;
  leftK  = this->fK;
  rightK = this->GetRightK();
  
  int nrowl = phiL.Rows();
  int nrowr = phiR.Rows();
  int il,jl,ir,jr;
//penalty = <A p^2>/h 
  const REAL penalty = fPenaltyConstant * (0.5 * (leftK*LeftPOrder*LeftPOrder + rightK*RightPOrder*RightPOrder)) / faceSize;
  
  // 1) left i / left j
  for(il=0; il<nrowl; il++) {
    for(jl=0; jl<nrowl; jl++) {
      ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
    }
  }
  
  // 2) right i / right j
  for(ir=0; ir<nrowr; ir++) {
    for(jr=0; jr<nrowr; jr++) {
      ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
    }
  }
  
  // 3) left i / right j
  for(il=0; il<nrowl; il++) {
    for(jr=0; jr<nrowr; jr++) {
      ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
    }
  }
  
  // 4) right i / left j
  for(ir=0; ir<nrowr; ir++) {
    for(jl=0; jl<nrowl; jl++) {
      ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
    }
  }
  
}

/** Termos de penalidade. */
void TPZMatPoisson3d::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
                                            TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc, int POrder, REAL faceSize){
  
  this->ContributeBCInterface( x, solL, dsolL, weight, normal, phiL, dphiL, ek, ef, bc);

  int il,jl,nrowl;
  nrowl = phiL.Rows(); 
  const REAL penalty = fPenaltyConstant * fK * POrder * POrder / faceSize; //Ap^2/h
  REAL outflow = 0.;
  for(il=0; il<fDim; il++) outflow += fConvDir[il]*normal[il];
   
  switch(bc.Type()) {
  case 0: // DIRICHLET  
    for(il=0; il<nrowl; il++) {
      ef(il,0) += weight * penalty * phiL(il,0) * bc.Val2()(0,0);
      for(jl=0; jl<nrowl; jl++) {
        ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
      }
    }
     
    break;
  case 1: // Neumann
    if(outflow > 0.)
    {
      for(il=0; il<nrowl; il++)
      {
        for(jl=0; jl<nrowl; jl++)
        {
          ek(il,jl) += weight * outflow * fC * phiL(il,0) * phiL(jl,0);
        }
      }
    }
    //nothing to be done
    break;
  default:
    PZError << "TPZMatPoisson3d::Wrong boundary condition type\n";
    break;
  }    
  
}

void TPZMatPoisson3d::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
                                          TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
                                          TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
                                          TPZFMatrix &ek,TPZFMatrix &ef){

  int nrowl = phiL.Rows();
  int nrowr = phiR.Rows();
  int il,jl,ir,jr,id;
  
  //Convection term
  REAL ConvNormal = 0.;
  for(id=0; id<fDim; id++) ConvNormal += fConvDir[id]*normal[id];
  if(ConvNormal > 0.) {
    for(il=0; il<nrowl; il++) {
      for(jl=0; jl<nrowl; jl++) {
        ek(il,jl) += weight * fC * ConvNormal * phiL(il)*phiL(jl);
      }
    }
    for(ir=0; ir<nrowr; ir++) {
      for(jl=0; jl<nrowl; jl++) {
        ek(ir+nrowl,jl) -= weight * fC * ConvNormal * phiR(ir) * phiL(jl);
      }
    }
  } else {
    for(ir=0; ir<nrowr; ir++) {
      for(jr=0; jr<nrowr; jr++) {
        ek(ir+nrowl,jr+nrowl) -= weight * fC * ConvNormal * phiR(ir) * phiR(jr);
      }
    }
    for(il=0; il<nrowl; il++) {
      for(jr=0; jr<nrowr; jr++) {
        ek(il,jr+nrowl) += weight * fC * ConvNormal * phiL(il) * phiR(jr);
      }
    }
  }
  
  
  //diffusion term
  REAL leftK, rightK;
  leftK  = this->fK;
  rightK = this->GetRightK();
  
  // 1) phi_I_left, phi_J_left
  for(il=0; il<nrowl; il++) {
    REAL dphiLinormal = 0.;
    for(id=0; id<fDim; id++) {
      dphiLinormal += dphiL(id,il)*normal[id];
    }
    for(jl=0; jl<nrowl; jl++) {
      REAL dphiLjnormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiLjnormal += dphiL(id,jl)*normal[id];
      }
      ek(il,jl) += weight * leftK * (
        this->fSymmetry * 0.5*dphiLinormal*phiL(jl,0)-0.5*dphiLjnormal*phiL(il,0)
      );
    }
  }
  
  // 2) phi_I_right, phi_J_right
  for(ir=0; ir<nrowr; ir++) {
    REAL dphiRinormal = 0.;
    for(id=0; id<fDim; id++) {
      dphiRinormal += dphiR(id,ir)*normal[id];
    }
    for(jr=0; jr<nrowr; jr++) {
      REAL dphiRjnormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiRjnormal += dphiR(id,jr)*normal[id];
      }
      ek(ir+nrowl,jr+nrowl) += weight * rightK * (
        this->fSymmetry * (-0.5 * dphiRinormal * phiR(jr) ) + 0.5 * dphiRjnormal * phiR(ir)
      );
    }
  }
  
  // 3) phi_I_left, phi_J_right
  for(il=0; il<nrowl; il++) {
    REAL dphiLinormal = 0.;
    for(id=0; id<fDim; id++) {
      dphiLinormal += dphiL(id,il)*normal[id];
    }
    for(jr=0; jr<nrowr; jr++) {
      REAL dphiRjnormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiRjnormal += dphiR(id,jr)*normal[id];
      }
      ek(il,jr+nrowl) += weight * (
        this->fSymmetry * (-0.5 * dphiLinormal * leftK * phiR(jr) ) - 0.5 * dphiRjnormal * rightK * phiL(il)
      );
    }
  }
  
  // 4) phi_I_right, phi_J_left
  for(ir=0; ir<nrowr; ir++) {
    REAL dphiRinormal = 0.;
    for(id=0; id<fDim; id++) {
      dphiRinormal += dphiR(id,ir)*normal[id];
    }
    for(jl=0; jl<nrowl; jl++) {
      REAL dphiLjnormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiLjnormal += dphiL(id,jl)*normal[id];
      }
      ek(ir+nrowl,jl) += weight * (
        this->fSymmetry * 0.5 * dphiRinormal * rightK * phiL(jl) + 0.5 * dphiLjnormal * leftK * phiR(ir)
      );
    }
  }
  
  if (this->IsSymetric()){
    if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
  }
}

void TPZMatPoisson3d::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
					    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {

  //  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
  int il,jl,nrowl,id;
  nrowl = phiL.Rows();
  REAL ConvNormal = 0.;
  for(id=0; id<fDim; id++) ConvNormal += fConvDir[id]*normal[id];
  switch(bc.Type()) {
  case 0: // DIRICHLET
    
    //Diffusion
    for(il=0; il<nrowl; il++) {
      REAL dphiLinormal = 0.;
      for(id=0; id<fDim; id++) {
        dphiLinormal += dphiL(id,il)*normal[id];
      }
      ef(il,0) += weight*fK*dphiLinormal*bc.Val2()(0,0) * this->fSymmetry;
      for(jl=0; jl<nrowl; jl++) {
        REAL dphiLjnormal = 0.;
        for(id=0; id<fDim; id++) {
          dphiLjnormal += dphiL(id,jl)*normal[id];
        }
        ek(il,jl) += weight*fK*(this->fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0));
      }
    }
    
    //Convection
    if(ConvNormal > 0.) {
      for(il=0; il<nrowl; il++) {
        for(jl=0; jl<nrowl; jl++) {
          ek(il,jl) += weight * fC * ConvNormal * phiL(il)*phiL(jl);
        }
      }
    } else {
      for(il=0; il<nrowl; il++) {
        ef(il,0) -= weight * fC * ConvNormal * bc.Val2()(0,0) * phiL(il);
      }
    }
      
    break;
  case 1: // Neumann
    for(il=0; il<nrowl; il++) {
      ef(il,0) += weight*phiL(il,0)*bc.Val2()(0,0);
    }
    break;
  default:
    PZError << "TPZMatPoisson3d::Wrong boundary condition type\n";
    break;
  }
    if (this->IsSymetric()){
      if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
}

void TPZMatPoisson3d::InterfaceErrors(TPZVec<REAL> &/*x*/,
                                      TPZVec<REAL> &leftu, TPZFMatrix &leftdudx, /* TPZFMatrix &leftaxes,*/ 
				      TPZVec<REAL> &rightu, TPZFMatrix &rightdudx, /* TPZFMatrix &rightaxes,*/ 
                                      TPZVec<REAL> &/*flux*/,
				      TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values, 
				      TPZVec<REAL> normal, REAL elsize) {
#ifndef DOS
#warning Metodo nao funcional
#endif
  TPZManVector<REAL,3> Lsol(1), Ldsol(3,0.), Rsol(1), Rdsol(3,0.);

  TPZFMatrix fake_axes(fDim,fDim,0.);  

  Solution(leftu,leftdudx,fake_axes,1,Lsol);
  Solution(leftu,leftdudx,fake_axes,2,Ldsol);

  Solution(rightu,rightdudx,fake_axes,1,Rsol);
  Solution(rightu,rightdudx,fake_axes,2,Rdsol);

#ifdef DEBUG
  if ( (leftdudx.Rows() != rightdudx.Rows()) || (leftdudx.Rows() != du_exact.Rows()) ){
    PZError << "TPZMatPoisson3d::InterfaceErrors - Left and right matrices should have" 
	    << endl 
	    << "same sizes in internal boundaries." 
	    << endl;
    exit (-1);
  }
#endif

  REAL Ldsolnormal = 0., Rdsolnormal = 0., ExactDNormal = 0.;
  for(int id = 0; id < fDim; id++) {
    Ldsolnormal  += Ldsol[id] * normal[id];
    Rdsolnormal  += Rdsol[id] * normal[id];
    ExactDNormal += du_exact(id, 0) * normal[id];
  }
 
  values.Resize(3);
  REAL aux;

  //values[1] : eror em norma L2

  //Jump aprox. solution - jump of exact solution i.e. zero
  aux = (Lsol - Rsol);

  //*= h ^ -gAlfa
  aux *= pow(elsize, -1.0 * gAlfa);
  values[1] = aux * aux;

  //values[2] : erro em semi norma H1
  values[2] = 0.;

  for(int id=0; id<fDim; id++) {
    //Normal gradient average <grad V> = 0.5 * (grad_left.n + grad_right.n)
    aux = 0.5 * (Ldsolnormal + Rdsolnormal);
    //<grad V> - <grad exact> = <grad V> - grad exact
    aux = aux - ExactDNormal;
    //*= h ^ gAlfa
    aux *= pow(elsize, gAlfa);
    values[2]  += aux * aux;
  }
  //values[0] : erro em norma H1 <=> norma Energia
  values[0]  = values[1]+values[2];
}


