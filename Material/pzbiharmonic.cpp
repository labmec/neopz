//$Id: pzbiharmonic.cpp,v 1.1 2003-11-27 16:02:39 igor Exp $

#include "pzbiharmonic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

REAL TPZBiharmonic::gLambda1 = -1.0;
REAL TPZBiharmonic::gLambda2 = -1.0; 
REAL TPZBiharmonic::gSigma      =  1.0;
REAL TPZBiharmonic::gL_alpha    =  6.0;
REAL TPZBiharmonic::gM_alpha   =  3.0;
REAL TPZBiharmonic::gL_betta     =  4.0;
REAL TPZBiharmonic::gM_betta    =  1.0;

TPZBiharmonic::TPZBiharmonic(int nummat, int dim) : TPZDiscontinuousGalerkin(nummat), fXf(1,1,0.), fXk(1,1,0.), fDim(dim) {
}

TPZBiharmonic::~TPZBiharmonic() {
}

void TPZBiharmonic::Print(ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "properties : \n";
  TPZMaterial::Print(out);

  fXf.Print ("Matrix fxf", out ) ;

  fXk.Print ("Matrix fxk", out ) ;         
}

void TPZBiharmonic::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &  dsol ,REAL weight,TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {

  int phr = phi.Rows();

  if(fForcingFunction) {            // phi(in, 0) = phi_in
    TPZManVector<REAL> res(1);
    fForcingFunction(x,res);       // dphi(i,j) = dphi_j/dxi
    fXf(0,0) = res[0];
  }
  //Equação de Poisson
  for( int in = 0; in < phr; in++ ) {
    ef(in, 0) +=  weight * fXf(0,0)*phi(in,0);
    for( int jn = 0; jn < phr; jn++ ) { 
      ek(in,jn) +=fXk(0,0)* weight * ( dphi(fDim,in) * dphi(fDim,jn) );
    }
  }
}


void TPZBiharmonic::ContributeBC(TPZVec<REAL> &/*x*/,TPZVec<REAL> &/*sol*/,REAL weight,
				     TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {
  //<!> usa ??????
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
  case 2 :		// condiçao mista
    for(in = 0 ; in < phi.Rows(); in++) {
      ef(in, 0) += v2[0] * phi(in, 0) * weight;
      for (jn = 0 ; jn < phi.Rows(); jn++) {
	ek(in,jn) += bc.Val1()(0,0) * phi(in,0) *
	  phi(jn,0) * weight;     // peso de contorno => integral de contorno
      }
    }
  }
}

/** returns the variable index associated with the name*/
int TPZBiharmonic::VariableIndex(char *name){
  if(!strcmp("Displacement6",name))   return  0;
  if(!strcmp("Solution",name))        return  1;
  if(!strcmp("Derivate",name))        return  2;
  if(!strcmp("POrder",name))          return 10;
  cout << "TPZBiharmonic::VariableIndex Error\n";
  return -1;
}

int TPZBiharmonic::NSolutionVariables(int var){

  if(var == 0) return 6;
  if(var == 1) return 1;
  if(var == 2) return fDim;
  if(var == 10) return 1;
  return TPZMaterial::NSolutionVariables(var);
  //  cout << "TPZBiharmonic::NSolutionVariables Error\n";
  //  return 0;
}

void TPZBiharmonic::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){

  if(var == 0 || var == 1) Solout[0] = Sol[0];//function
  if(var == 2) {
    int id;
    for(id=0 ; id<fDim; id++) {
      Solout[id] = DSol(id,0);//derivate
    }
  }
}

void TPZBiharmonic::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
  //Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZBiharmonic::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

  TPZVec<REAL> sol(1),dsol(3,0.);
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
    values[2]  += (dx[id] - du_exact(id,0))*(dx[id] - du_exact(id,0));
  }
  //values[0] : erro em norma H1 <=> norma Energia
  values[0]  = values[1]+values[2];
}


void TPZBiharmonic::ContributeInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL,TPZVec<REAL> &solR,TPZFMatrix &dsolL,
				   TPZFMatrix &dsolR,REAL weight,TPZVec<REAL> &normal,TPZFMatrix &phiL,
				   TPZFMatrix &phiR,TPZFMatrix &dphiL,TPZFMatrix &dphiR,
					  TPZFMatrix &ek,TPZFMatrix &ef){

  int nrowl = phiL.Rows();
  int nrowr = phiR.Rows();
  int il,jl,ir,jr,id;
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
      ek(il,jl) += weight*(
			   0.5*dphiLinormal*phiL(jl,0)-0.5*dphiLjnormal*phiL(il,0)
			   );
    }
  }
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
      ek(ir+nrowl,jr+nrowl) += weight*(
			   -0.5*dphiRinormal*phiR(jr)+0.5*dphiRjnormal*phiR(ir)
			   );
    }
  }
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
      ek(il,jr+nrowl) += weight*(
			   -0.5*dphiLinormal*phiR(jr)-0.5*dphiRjnormal*phiL(il)
			   );
    }
  }
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
      ek(ir+nrowl,jl) += weight*(
			   +0.5*dphiRinormal*phiL(jl)+0.5*dphiLjnormal*phiR(ir)
			   );
    }
  }
}

void TPZBiharmonic::ContributeBCInterface(TPZVec<REAL> &x,TPZVec<REAL> &solL, TPZFMatrix &dsolL, REAL weight, TPZVec<REAL> &normal,
					    TPZFMatrix &phiL,TPZFMatrix &dphiL, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc) {

  //  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
  int il,jl,nrowl,id;
  nrowl = phiL.Rows();
  switch(bc.Type()) {
  case 0: // DIRICHLET
    for(il=0; il<nrowl; il++) {
      REAL dphiLinormal = 0.;
      for(id=0; id<fDim; id++) {
	dphiLinormal += dphiL(id,il)*normal[id];
      }
      ef(il,0) += weight*dphiLinormal*bc.Val2()(0,0);
      for(jl=0; jl<nrowl; jl++) {
	REAL dphiLjnormal = 0.;
	for(id=0; id<fDim; id++) {
	  dphiLjnormal += dphiL(id,jl)*normal[id];
	}
	ek(il,jl) += weight*(dphiLinormal*phiL(jl,0)-dphiLjnormal*phiL(il,0));
      }
    }
    break;
  case 1: // Neumann
    for(il=0; il<nrowl; il++) {
      ef(il,0) += weight*phiL(il,0)*bc.Val2()(0,0);
    }
    break;
  default:
    PZError << "TPZBiharmonic::Wrong boundary condition type\n";
    break;
  }
}
