#include "pzpoisson3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

int TPZMatPoisson3d::problema = 0;

TPZMatPoisson3d::TPZMatPoisson3d(int nummat, int dim) : TPZMaterial(nummat), fXf(1,1,0.), fDim(dim) {
}

TPZMatPoisson3d::~TPZMatPoisson3d() {
}

int TPZMatPoisson3d::NStateVariables() {
  return 1;
}

void TPZMatPoisson3d::Print(ostream &out) {
  out << "name of material : " << Name() << "\n";
  out << "properties : \n";           
  TPZMaterial::Print(out);
}

void TPZMatPoisson3d::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &  dsol ,REAL weight,TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &dphi,TPZFMatrix &ek,TPZFMatrix &ef) {

  int phr = phi.Rows();

  if(fForcingFunction) {            // phi(in, 0) = phi_in
    TPZManVector<REAL> res(1);
    fForcingFunction(x,res);       // dphi(i,j) = dphi_j/dxi
    fXf(0,0) = res[0];
  }
  int dim = dphi.Rows();

  if(problema==1) {
    //projeção L2 da carga fXf : fForcingFunction
    for( int in = 0; in < phr; in++ ) {
      ef(in, 0) += weight * fXf(0,0) * phi(in, 0);
      for( int jn = 0; jn < phr; jn++ ) {
	ek(in,jn) += weight * phi(in,0)* phi(jn,0);
      }
    }
  } else
    if(problema==2) {
      //Equação de Poisson
      for( int in = 0; in < phr; in++ ) {
        int kd;
        for(kd=0; kd<dim; kd++) {
        	ef(in, 0) += - weight * dsol(kd, 0) * dphi(kd, in);
        }
	ef(in, 0) += weight * fXf(0,0) * phi(in, 0);
	for( int jn = 0; jn < phr; jn++ ) {
	  for(kd=0; kd<dim; kd++) {
	    ek(in,jn) += weight * ( dphi(kd,in) * dphi(kd,jn) );
	  }
	}
      }
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
int TPZMatPoisson3d::VariableIndex(char *name){
  if(!strcmp("Displacement6",name))   return  0;
  if(!strcmp("Solution",name))        return  1;
  if(!strcmp("Derivate",name))        return  2;
  if(!strcmp("POrder",name))          return 10;
  cout << "TPZMatPoisson3d::VariableIndex Error\n";
  return -1;
}

int TPZMatPoisson3d::NSolutionVariables(int var){

  if(var == 0 || var == 1 || var == 2 || var == 10) return 1;
  cout << "TPZMatPoisson3d::NSolutionVariables Error\n";
  return 0;
}

void TPZMatPoisson3d::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &/*axes*/,int var,TPZVec<REAL> &Solout){

  if(var == 0 || var == 1) Solout[0] = Sol[0];//function
  if(var == 2) {
    Solout[0] = DSol(0,0);//derivate
    Solout[1] = DSol(1,0);//derivate
    Solout[2] = DSol(2,0);//derivate
  }
}

void TPZMatPoisson3d::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
  //Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZMatPoisson3d::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
			       TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
			       TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

  TPZVec<REAL> sol(1),dsol(3);
  Solution(u,dudx,axes,1,sol);
  Solution(u,dudx,axes,2,dsol);
  REAL dx = dsol[0]*axes(0,0)+dsol[1]*axes(1,0)+dsol[2]*axes(2,0);
  REAL dy = dsol[0]*axes(0,1)+dsol[1]*axes(1,1)+dsol[2]*axes(2,1);
  REAL dz = dsol[0]*axes(0,2)+dsol[1]*axes(1,2)+dsol[2]*axes(2,2);
  //values[1] : eror em norma L2
  values[1]  = pow(sol[0] - u_exact[0],2.0);
  //values[2] : erro em semi norma H1
  values[2]  = pow(dx - du_exact(0,0),2.0);
  values[2] += pow(dy - du_exact(1,0),2.0);
  values[2] += pow(dz - du_exact(2,0),2.0);
  //values[0] : erro em norma H1 <=> norma Energia
  values[0]  = values[1]+values[2];
}


#ifdef _AUTODIFF
void TPZMatPoisson3d::ContributeEnergy(TPZVec<REAL> &x,
			      TPZVec<FADFADREAL> &sol,
			      TPZVec<FADFADREAL> &dsol,
			      FADFADREAL &U,
			      REAL weight)
{
      int dim = dsol.NElements()/sol.NElements();

      //Equação de Poisson

      int i, eqs = dsol.NElements()/dim;
      if(sol.NElements() != 1) PZError << "";

//cout << "FADREAL init : \n" << FADREAL(weight * fXf(0,0));

      //FADFADREAL Buff;

      U+= sol[0] * FADREAL(weight * fXf(0,0));

      switch(dim)
      {
      case 1:
             U+=(dsol[0] * dsol[0])*FADREAL(weight/2.); // U=((du/dx)^2)/2

	 break;
      case 2:
             U+=(dsol[0] * dsol[0] +
	         dsol[1] * dsol[1])*(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2
             /*Buff  = dsol[0] * dsol[0];
             Buff += dsol[1] * dsol[1];
	     U += Buff * FADREAL(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2*/
	 break;
      case 3:
             U+=(dsol[0] * dsol[0] +
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
  case 2 :	// condiçao mista
    // U += 1/2 * weight * Integral(<(u-u0), [g].(u-u0)> dOmega)
    U += ( solMinBC * /*scalar*/ FADREAL(bc.Val1()(0,0)) * /*matrix oprt*/ solMinBC ) * FADREAL(weight / 2.);
    break;

  }
}

#endif
