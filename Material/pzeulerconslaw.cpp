//$Id: pzeulerconslaw.cpp,v 1.11 2003-11-20 21:39:06 erick Exp $

#include "pzeulerconslaw.h"
//#include "TPZDiffusionConsLaw.h"
#include "pzartdiff.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzreal.h"
#include <math.h>
#include "pzstring.h"

TPZEulerConsLaw2::~TPZEulerConsLaw2(){

}

TPZEulerConsLaw2::TPZEulerConsLaw2(int nummat,REAL timeStep,
			REAL gamma,int dim,
			TPZArtDiffType artdiff) :
			TPZConservationLaw2(nummat,timeStep,dim),
			fArtDiff(artdiff, gamma),
			fDiff(Explicit_TD),
			fConvVol(Explicit_TD),
			fConvFace(Explicit_TD)
{
  fGamma = gamma;
}

TPZEulerConsLaw2::TPZEulerConsLaw2() :
			TPZConservationLaw2(-1, 0, 3),
			fArtDiff(LeastSquares_AD, 1.4),
			fDiff(),
			fConvVol(),
			fConvFace()
{
  fGamma = 1.4;
}

void TPZEulerConsLaw2::SetTimeDiscr(TPZTimeDiscr Diff, TPZTimeDiscr ConvVol, TPZTimeDiscr ConvFace)
{
   fDiff = Diff;
   fConvVol = ConvVol;
   fConvFace = ConvFace;
}

REAL TPZEulerConsLaw2::OptimalCFL(int degree)
{
   return 1./((2.0*(REAL)degree) + 1.0);
}

void TPZEulerConsLaw2::SetTimeStep(REAL maxveloc,REAL deltax,int degree)
{
  REAL CFL = fCFL;
  // Notice that fCFL remains 0, so that optimal CFL will
  // be computed unless CFL is redefined.
  if(CFL == 0.0) CFL = OptimalCFL(degree);

  REAL deltaT = CFL*deltax/maxveloc;
  //cout << "TPZCompMesh::Delta Time : " << deltaT << endl;
  TPZConservationLaw2::SetTimeStep(deltaT);
}

int TPZEulerConsLaw2::NStateVariables() {
  return (2 + Dimension());//U = (rho, rhou, rhov, rhow, rhoe)
}

REAL TPZEulerConsLaw2::Pressure(TPZVec<REAL> &U)
{
   REAL press;
   Pressure(press, U);
   return press;
}

void TPZEulerConsLaw2::Print(ostream &out) {

  TPZMaterial::Print(out);

  TPZConservationLaw2::Print(out);
  out << "Artificial Diffusion: " <<
  fArtDiff.DiffusionName().Str() << endl;
  out << "Number of State Variables: " << NStateVariables() << endl;
  out << "Number of Fluxes: " << NFluxes() << endl;

  switch(fDiff)
  {
     case(Explicit_TD):
        cout << "Explicit Diffusive term\n";
	break;
     case(ApproxImplicit_TD):
        cout << "ApproxImplicit Diffusive term\n";
	break;
     case(Implicit_TD):
        cout << "Implicit Diffusive term\n";
	break;
     default:
        cout << "No Diffusive term\n";
  }

  switch(fConvVol)
  {
     case(Explicit_TD):
        cout << "Explicit Volume Convective term\n";
	break;
     case(Implicit_TD):
        cout << "Implicit Volume Convective term\n";
	break;
     default:
        cout << "No Volume Convective term\n";
  }

  switch(fConvFace)
  {
     case(Explicit_TD):
        cout << "Explicit Face Convective term\n";
	break;
     case(Implicit_TD):
        cout << "Implicit Face Convective term\n";
	break;
     default:
        cout << "No Face Convective term\n";
  }


}

int TPZEulerConsLaw2::VariableIndex(char *name) {
  if( !strcmp(name,"density")  )     return 1;//rho
  if( !strcmp(name,"velocity") )     return 2;//(u,v,w)
  if( !strcmp(name,"energy")   )     return 3;//E
  if( !strcmp(name,"pressure") )     return 4;//p
  if( !strcmp(name,"solution") )     return 5;//(ro,u,v,w,E)
  if( !strcmp(name,"normvelocity") ) return 6;//sqrt(u²+v²+w²)
  cout << "TPZEulerConsLaw2::VariableIndex not defined\n";
  return TPZMaterial::VariableIndex(name);
}

int TPZEulerConsLaw2::NSolutionVariables(int var){

  if(var == 1 || var == 3 || var == 4 || var == 6) return 1;
  if(var == 2) return Dimension();
  if(var == 5) return NStateVariables();

  cout << "TPZEulerConsLaw2::NSolutionVariables not defined\n";
  return 0;
}

int TPZEulerConsLaw2::NFluxes()
{
  return Dimension();
}

//-----------------Solutions


void TPZEulerConsLaw2::ComputeSolLeft(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft){

  if(!bcleft){
    PZError << "TPZEulerConsLaw2::ComputeSolLeft null bundary condition return\n";
    return;
  }
  int i,nstate = NStateVariables();
  REAL vpn=0.;
  switch (bcleft->Type()){
  case 0://Dirichlet
  case 1://Neumann
  case 2://Mista
    PZError << "TPZEulerConsLaw2::ComputeSolLeft boundary condition error\n";
    break;
  case 3://Dirichlet: nada a fazer a CC é a correta
    for(i=0;i<nstate; i++) soll[i] = bcleft->Val2().operator()(i,0);//<!>erick
    break;
  case 4://recuperar valor da solu¢ão MEF direita: saida livre
    for(i=0;i<nstate; i++) soll[i] = solr[i];
    break;
  case 5://parede
    for(i=1;i<nstate-1;i++) vpn += solr[i]*normal[i-1];//v.n
    for(i=1;i<nstate-1;i++) soll[i] = solr[i] - 2.0*vpn*normal[i-1];
    soll[0] = solr[0];
    soll[nstate-1] = solr[nstate-1];
    break;
  case 6://não refletivas (campo distante)
    for(i=0;i<nstate;i++) soll[i] = solr[i];
    break;
  default:
    PZError << "TPZEulerConsLaw2::ContributeInterface Boundary Condition Type  does Not Exist\n";
  }
}


void TPZEulerConsLaw2::ComputeSolRight(TPZVec<REAL> &solr,TPZVec<REAL> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright){

  if(!bcright){
    PZError << "TPZEulerConsLaw2::ComputeSolLeft null bundary condition return\n";
    return;
  }
  int i,nstate = NStateVariables();
  REAL vpn=0.;
  switch (bcright->Type()){
  case 0://Dirichlet
  case 1://Neumann
  case 2://Mista
    PZError << "TPZEulerConsLaw2::ComputeSolLeft boundary condition error\n";
    break;
  case 3://Dirichlet: nada a fazer a CC é a correta
    for(i=0;i<nstate; i++) solr[i] = bcright->Val2().operator()(i,0);//<!>erick
    break;
  case 4://recuperar valor da solu¢ão MEF esquerda: saida livre
    for(i=0;i<nstate; i++) solr[i] = soll[i];
    break;
  case 5://condi¢ão de parede
    for(i=1;i<nstate-1;i++) vpn += soll[i]*normal[i-1];//v.n
    for(i=1;i<nstate-1;i++) solr[i] = soll[i] - 2.0*vpn*normal[i-1];
    solr[0] = soll[0];
    solr[nstate-1] = soll[nstate-1];
    break;
  case 6://não refletivas (campo distante)
    for(i=0;i<nstate;i++) solr[i] = soll[i];
    break;
  default:
    PZError << "TPZEulerConsLaw2::ContributeInterface Boundary Condition Type does Not Exist\n";
  }
}


void TPZEulerConsLaw2::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

  if(fabs(Sol[0]) < 1.e-10) {
    cout << "\nTPZEulerConsLaw2::Solution: Density almost null\n";
    cout << "Density = " << Sol[0] << endl;
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
  } else {
    //cout << "TPZEulerConsLaw2::Solution variable in the base class\n";
    TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
  }
}


void TPZEulerConsLaw2::Flux(TPZVec<REAL> &x,TPZVec<REAL> &Sol,TPZFMatrix &DSol,
			   TPZFMatrix &axes,TPZVec<REAL> &flux) {
  TPZVec<REAL> Fx,Fy,Fz;
  Flux(Sol,Fx,Fy,Fz);
  int cap = Sol.NElements();
  int nstate = NStateVariables(),i;
  if(cap != nstate){
    PZError << "\nTPZEulerConsLaw2::Flux data size error\n";
    flux.Resize(0);
    return;
  }
  if(nstate == 3){
    flux.Resize(3);
    for(i=0;i<3;i++) flux[i] = Fx[i];
    return;
  } else
  if(nstate == 4){
    flux.Resize(8);
    for(i=0;i<4;i++) flux[i] = Fx[i];
    for(i=4;i<8;i++) flux[i] = Fy[i - 4];
    return;
  } else
  if(nstate == 5){
    flux.Resize(15);
    for(i=00;i<05;i++) flux[i] = Fx[i];
    for(i=05;i<10;i++) flux[i] = Fy[i-5];
    for(i=10;i<15;i++) flux[i] = Fz[i-10];
  }
}


void TPZEulerConsLaw2::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
                               TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
                               TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {

  cout << "\nTPZEulerConsLaw2::Errors not implemented yet, program exit\n";
  exit(-1);

//   TPZVec<REAL> sol(1),dsol(3);
//   Solution(u,dudx,axes,1,sol);
//   Solution(u,dudx,axes,2,dsol);
//   REAL dx = dsol[0]*axes(0,0)+dsol[1]*axes(1,0)+dsol[2]*axes(2,0);
//   REAL dy = dsol[0]*axes(0,1)+dsol[1]*axes(1,1)+dsol[2]*axes(2,1);
//   REAL dz = dsol[0]*axes(0,2)+dsol[1]*axes(1,2)+dsol[2]*axes(2,2);
//   //values[1] : eror em norma L2
//   values[1]  = pow(sol[0] - u_exact[0],2.0);
//   //values[2] : erro em semi norma H1
//   values[2]  = pow(dx - du_exact(0,0),2.0);
//   values[2] += pow(dy - du_exact(1,0),2.0);
//   values[2] += pow(dz - du_exact(2,0),2.0);
//   //values[0] : erro em norma H1 <=> norma Energia
//   values[0]  = values[1]+values[2];
}

void TPZEulerConsLaw2::SetDelta(REAL delta)
{
   fArtDiff.SetDelta(delta);
}

/*

void TPZEulerConsLaw2::Flux(TPZVec<REAL> &U,TPZVec<REAL> &Fx,TPZVec<REAL> &Fy,TPZVec<REAL> &Fz)
{
   TFlux(U, Fx, Fy, Fz);
}
*/


//------------------Differentiable variables setup

#ifdef _AUTODIFF

void TPZEulerConsLaw2::PrepareFAD(
		TPZVec<REAL> & sol, TPZFMatrix & dsol,
		TPZFMatrix & phi, TPZFMatrix & dphi,
		TPZVec<FADREAL> & FADsol,
		TPZVec<FADREAL> & FADdsol)
{
   int nState = NStateVariables();
   int nShape = phi.Rows();
   int i_state, i_shape, k;
   int nDer = nState * nShape;

   // initializing the differentiable variables
   FADREAL defaultFAD(nDer, 0., 0.);
   if(defaultFAD.dx(0)==1.)PZError << "\nError: FAD doesn't have default constructor for parameters: (number of derivatives, default value, default derivative value) !";
   FADsol.Resize(nState);
   FADsol.Fill(defaultFAD);

   FADdsol.Resize(nState * fDim);
   FADdsol.Fill(defaultFAD);

   // copying the solution and spatial derivative values
   for(i_state = 0; i_state < nState; i_state++)
   {
      FADsol[i_state].val() = sol[i_state];
      for(k = 0; k < fDim; k ++)
         FADdsol[i_state * fDim + k].val() = dsol(k,i_state);
   }

   // preparing the coefficient derivatives
   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_state + i_shape * nState;
         FADsol[i_state].fastAccessDx(index)=phi(i_shape,0);
         for(k = 0; k < fDim; k++)
            FADdsol[i_state * fDim + k].fastAccessDx(index)=dphi(k,i_shape);
      }
}

void TPZEulerConsLaw2::PrepareInterfaceFAD(
		TPZVec<REAL> &solL,TPZVec<REAL> &solR,
		TPZFMatrix &phiL,TPZFMatrix &phiR,
		TPZVec<FADREAL> & FADsolL,
		TPZVec<FADREAL> & FADsolR)
{
   int nState = NStateVariables();
   int nShapeL = phiL.Rows();
   int nShapeR = phiR.Rows();
   int i_state, i_shape;
   int nDerL = nState * nShapeL;
   int nDerR = nState * nShapeR;

   // initializing the differentiable variables
   FADREAL defaultFAD(nDerL + nDerR, 0., 0.);
   if(defaultFAD.dx(0)==1.)PZError << "\nError: FAD doesn't have default constructor for parameters: (number of derivatives, default value, default derivative value) !";
   FADsolL.Resize(nState);
   FADsolL.Fill(defaultFAD);

   FADsolR.Resize(nState);
   FADsolR.Fill(defaultFAD);

   // copying the solution and spatial derivatives values
   for(i_state = 0; i_state < nState; i_state++)
   {
      FADsolL[i_state].val() = solL[i_state];
      FADsolR[i_state].val() = solR[i_state];
   }

   // preparing the coefficient derivatives
   for(i_shape = 0; i_shape < nShapeL; i_shape++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_state + i_shape * nState;
         FADsolL[i_state].fastAccessDx(index)=phiL(i_shape,0);
         //FADsolR[i_state].fastAccessDx(index)=phiL(i_shape,0);
      }
   for(i_shape = 0/*nShapeL*/; i_shape < /*nShapeL +*/ nShapeR; i_shape++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_state + (i_shape + nShapeL) * nState;
         //FADsolL[i_state].fastAccessDx(index)=phiR(i_shape,0);
         FADsolR[i_state].fastAccessDx(index)=phiR(i_shape,0);
      }
}

#endif

//----------------Contributions

void TPZEulerConsLaw2::Contribute(TPZVec<REAL> &x,TPZFMatrix &jacinv,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight, TPZFMatrix &axes,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
//return;

   // initial guesses for sol
   // fForcingFunction is null at iterations > 0
   if(fForcingFunction)
   {
      TPZVec<REAL> res;
      int i, nState = NStateVariables();
      fForcingFunction(x, res);
      for(i = 0; i < nState; i++)
         sol[i] = res[i];
   }

   if(fContributionTime == Last_CT)
   {
       ContributeLast(x, jacinv,
			sol, dsol,
			weight,
			phi, dphi,
			ef);
       return;
   }

   if(fContributionTime == Advanced_CT)
   {
       ContributeAdv(x, jacinv,
			sol, dsol,
			weight,
			phi, dphi,
			ek, ef);
       return;
   }

   PZError << "TPZEulerConsLaw2::Contribute> Unhandled Contribution Time";
}

void TPZEulerConsLaw2::ContributeLast(TPZVec<REAL> &x,TPZFMatrix &jacinv,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
// contributing the explicit parcell of the residual to the
// rhs.

   // the parcell T2 is always explicit.
   ContributeExplT2(x,sol,weight,phi,ef);

   // contributing volume-based quantities
   // diffusive term
   if (fDiff == Explicit_TD)
         ContributeExplDiff(x,sol,dsol,weight,phi,dphi,ef);

   // Volume Convective term
   if (fConvVol == Explicit_TD)
         ContributeExplConvVol(x, sol, weight, phi, dphi, ef);

}


void TPZEulerConsLaw2::ContributeAdv(TPZVec<REAL> &x,TPZFMatrix &jacinv,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ek, TPZFMatrix &ef)
{
// contributing the implicit parcell of the residual to the
// rhs.

   // the parcell T1 is always implicit.
   ContributeImplT1(x,sol,dsol,weight,phi,dphi,ek,ef);

   // contributing volume-based quantities
   // diffusive term
   if (fDiff == Implicit_TD)
   {
      // if diffusive term is implicit
      // then the FAD classes must be initialized
      #ifdef _AUTODIFF
         TPZVec<FADREAL> FADsol, FADdsol;
         PrepareFAD(sol, dsol, phi, dphi, FADsol, FADdsol);
	    ContributeImplDiff(x,FADsol,FADdsol, weight, ek, ef);
      #else
         cout << "TPZEulerConsLaw2::Contribute> Implicit diffusive contribution: _AUTODIFF directive not configured -> Using an approximation to the tgMatrix";
         ContributeApproxImplDiff(x,sol,dsol,weight,phi,dphi,ek,ef);
      #endif
   }else
   {
      if (fDiff == ApproxImplicit_TD)
         ContributeApproxImplDiff(x,sol,dsol,weight,phi,dphi,ek,ef);
   }

   // Volume convective term
   if (fConvVol == Implicit_TD)
         ContributeImplConvVol(x,sol,dsol,weight,phi,dphi,ek,ef);

/*
ek(0,0) = 1.;
ek(1,1) = 1.;
ek(2,2) = 1.;
ek(3,3) = 1.;*/
}


void TPZEulerConsLaw2::ContributeInterface(
		TPZVec<REAL> &x,
		TPZVec<REAL> &solL, TPZVec<REAL> &solR,
		TPZFMatrix &dsolL, TPZFMatrix &dsolR,
		REAL weight, TPZVec<REAL> &normal,
		TPZFMatrix &phiL,TPZFMatrix &phiR,
		TPZFMatrix &dphiL,TPZFMatrix &dphiR,
		TPZFMatrix &ek,TPZFMatrix &ef)
{
//return;
   // initial guesses for left and right sol
   // fForcingFunction is null at iterations > 0
   if(fForcingFunction)
   {
      TPZVec<REAL> res;
      int i, nState = NStateVariables();
      fForcingFunction(x, res);
      for(i = 0; i < nState; i++)
         solL[i] = solR[i] = res[i];
   }


   // contributing face-based quantities
   if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
      {
      // if face contribution is implicit,
      // then the FAD classes must be initialized
      #ifdef _AUTODIFF
         TPZVec<FADREAL> FADsolL, FADsolR;
         PrepareInterfaceFAD(solL, solR, phiL, phiR, FADsolL, FADsolR);
         ContributeImplConvFace(x,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef);
      #else
      // forcint explicit contribution and issueing an warning
         cout << "TPZEulerConsLaw2::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
         ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef);
      #endif
      }

   if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
   {
         ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef);
   }

}


void TPZEulerConsLaw2::ContributeBC(TPZVec<REAL> &/*x*/,TPZVec<REAL> &sol,REAL weight,
                                   TPZFMatrix &/*axes*/,TPZFMatrix &phi,TPZFMatrix &ek,
                                   TPZFMatrix &ef,TPZBndCond &bc)
{
  int phr = phi.Rows();
  short in,jn,i,j;
  int nstate = NStateVariables();
  REAL v2[5];//máximo nstate
  for(i=0;i<nstate;i++) v2[i] = bc.Val2()(i,0);

  switch (bc.Type()) {
  case 0 :// Dirichlet condition
    for(in = 0 ; in < phr; in++) {
      for(i = 0 ; i < nstate; i++)
        ef(in*nstate+i,0) += gBigNumber * weight * v2[i] * phi(in,0);
      for (jn = 0 ; jn < phr; jn++) {
        for(i = 0 ; i < nstate; i++)
          ek(in*nstate+i,jn*nstate+i) += gBigNumber * weight * phi(in,0) * phi(jn,0);
      }
    }
    break;
  case 1 :// Neumann condition
    for(in = 0 ; in < phi.Rows(); in++) {
      for(i = 0 ; i < nstate; i++)
        ef(in*nstate+i,0) += v2[i] * phi(in,0) * weight;
    }
    break;
  case 2 :// condiçao mista
    for(in = 0 ; in < phi.Rows(); in++) {
      for(i = 0 ; i < nstate; i++)
        ef(in*nstate+i, 0) += weight * v2[i] * phi(in, 0);
      for (jn = 0 ; jn < phi.Rows(); jn++) {
        for(i = 0 ; i < nstate; i++) for(j = 0 ; j < nstate; j++)
          ek(in*nstate+i,jn*nstate+j) += weight * bc.Val1()(i,j) * phi(in,0) * phi(jn,0);
      }
    }
  }
}


//------------------internal contributions

void TPZEulerConsLaw2::ContributeApproxImplDiff(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   fArtDiff.ContributeApproxImplDiff(fDim, sol, dsol, dphi,
			ek, ef, weight,
			fTimeStep);
}

void TPZEulerConsLaw2::ContributeExplDiff(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
   fArtDiff.ContributeExplDiff(fDim, sol, dsol, dphi,
			ef, weight,
			fTimeStep);
}


#ifdef _AUTODIFF
void TPZEulerConsLaw2::ContributeImplDiff(TPZVec<REAL> &x,
			TPZVec<FADREAL> &sol,TPZVec<FADREAL> &dsol,
			REAL weight,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   fArtDiff.ContributeImplDiff(fDim, sol, dsol,
			ek, ef, weight, fTimeStep);
}
#endif

void TPZEulerConsLaw2::ContributeExplConvFace(TPZVec<REAL> &x,
			TPZVec<REAL> &solL,TPZVec<REAL> &solR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ef)
{
   int nState = NStateVariables();
   TPZVec<REAL > flux(nState,0.);
   Roe_Flux/*Test_Flux*/<REAL>(solL, solR, normal, fGamma, flux);
   int nShapeL = phiL.Rows(),
       nShapeR = phiR.Rows(),
       i_shape, i_state;
   REAL constant = - TimeStep() * weight; // - deltaT * weight

   // Contribution referring to the left element
   for(i_shape = 0; i_shape < nShapeL; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
         ef(i_shape*nState + i_state,0) +=
	    flux[i_state] * phiL(i_shape,0) * constant;

   // Contribution referring to the right element
   // REM: The contributions are negative in comparison
   // to the left elements: opposed normals
   for(i_shape = 0; i_shape < nShapeR; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
         ef((nShapeL + i_shape)*nState + i_state,0) -=
	    flux[i_state] * phiR(i_shape,0) * constant;
}

#ifdef _AUTODIFF

void TPZEulerConsLaw2::ContributeImplConvFace(TPZVec<REAL> &x,
			TPZVec<FADREAL> &solL,TPZVec<FADREAL> &solR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   int nState = NStateVariables();
   TPZVec<FADREAL > flux(nState,0.);
   Roe_Flux/*Test_Flux*/(solL, solR, normal, fGamma, flux);
   int nShapeL = phiL.Rows(),
       nShapeR = phiR.Rows(),
       i_shape, i_state, j,
       nDer = (nShapeL + nShapeR) * nState;
   REAL constant = - TimeStep() * weight; // - deltaT * weight

   // Contribution referring to the left element
   for(i_shape = 0; i_shape < nShapeL; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_shape*nState + i_state;
         ef(index,0) +=
	    flux[i_state].val() * phiL(i_shape,0) * constant;
	 for(j = 0; j < nDer; j++)
	    ek(index, j) += flux[i_state].dx(j) *
	       phiL(i_shape,0) * constant;
      }

   // Contribution referring to the right element
   // REM: The contributions are negative in comparison
   // to the left elements: opposed normals
   for(i_shape = 0; i_shape < nShapeR; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = (nShapeL + i_shape)*nState + i_state;
         ef(index,0) -=
	    flux[i_state].val() * phiR(i_shape,0) * constant;
	 for(j = 0; j < nDer; j++)
	    ek(index, j) -= flux[i_state].dx(j) *
	       phiR(i_shape,0) * constant;
      }
}

#endif

void TPZEulerConsLaw2::ContributeExplConvVol(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,
			REAL weight,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
   TPZVec< TPZVec<REAL> > F(3);
   Flux(sol, F[0], F[1], F[2]);
   REAL constant = -TimeStep() * weight;

   int i_state, i_shape, nShape = phi.Rows(), k;
   int nState = NStateVariables();

   for(i_shape=0; i_shape<nShape; i_shape++)
      for(i_state = 0; i_state<nState; i_state++)
      {
         int index = i_state + i_shape * nState;
         for(k=0; k<fDim; k++)
         ef(index,0) += (F[k])[i_state] * dphi(k, i_shape) * constant;
         // ef(index) += F<scalar>GradPhi
      }

   /*}else{
      cout << "TPZEulerConsLaw2::ContributeConvVol> Unhandled time discretization";
   }*/
}


void TPZEulerConsLaw2::ContributeImplConvVol(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   TPZVec< TPZVec<REAL> > F(3);
   Flux(sol, F[0], F[1], F[2]);
   REAL constant = -TimeStep() * weight;

   TPZVec< TPZDiffMatrix<REAL> > Ai(3);
   JacobFlux(sol, Ai);
   int j_shape, j_state;
   int i_state, i_shape, nShape = phi.Rows(), k;
   int nState = NStateVariables();

   for(i_shape=0; i_shape<nShape; i_shape++)
      for(i_state = 0; i_state<nState; i_state++)
      {
         int index = i_state + i_shape * nState;
         for(k=0; k<fDim; k++)
         {
            // ef(index) += F<scalar>GradPhi
            ef(index,0) += (F[k])[i_state] * dphi(k, i_shape)
                                           * constant;
            for(j_shape = 0; j_shape < nShape; j_shape++)
                for(j_state = 0; j_state < nState; j_state++)
                   ek(index, j_state + j_shape * nState) +=
                                Ai[k](i_state, j_state) *
                                dphi(k, i_shape) *
                                phi(j_shape,0) *
                                constant;
			// dConvVol/dU* =
			// dConvVol/dU * dU/dU*
			// dConvVol/dU = Ai<scalar>GradPhi
			// dU/dU* = phi(j)
         }
      }
}

void TPZEulerConsLaw2::ContributeImplT1(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   int i_shape, ij_state, j_shape;
   int nState = NStateVariables();
   int nShape = phi.Rows();

   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(ij_state = 0; ij_state < nState; ij_state++)
      {
          int index = i_shape * nState + ij_state;
	  // ef += sol*phi(i)
          ef(index, 0) +=
              sol[ij_state] * phi(i_shape,0) *
              weight;
	  // ek += phi(i)*phi(j)
          for(j_shape = 0; j_shape < nShape; j_shape++)
             ek(index, j_shape * nState + ij_state) +=
                phi(i_shape,0) *
                phi(j_shape,0) *
                weight;
      }
}

void TPZEulerConsLaw2::ContributeExplT2(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,
			REAL weight,
			TPZFMatrix &phi,
			TPZFMatrix &ef)
{
   int i_shape, i_state;
   int nState = NStateVariables();
   int nShape = phi.Rows();

   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(i_state = 0; i_state < nState; i_state++)
      {
          int index = i_shape * nState + i_state;
	  // ef += sol*phi(i)
          ef(index, 0) -=
              sol[i_state] * phi(i_shape,0) *
              weight; // the T2 parcell is negative
      }
}


//--------------test purposes

//método para testes de implementa¢ão, verifica¢ão das integrais elementares
void TPZEulerConsLaw2::ContributeTESTE(TPZVec<REAL> &x,TPZFMatrix &jacinv,TPZVec<REAL> &sol,TPZFMatrix &dsol,
				 REAL weight,TPZFMatrix &axes,TPZFMatrix &phi,TPZFMatrix &dphi,
				      TPZFMatrix &ek,TPZFMatrix &ef) {

//   está chamada deve ser feita do Contribute(..)
//   if(0){
//     ContributeTESTE(x,jacinv,sol,dsol,weight,axes,phi,dphi,ek,ef);//PARA TESTES
//     return;
//   }

  static int firsttime = 1;
  static REAL EK[2],EF[3];
  if(firsttime && 0){
    EF[0]=1.; EF[1]=1.; EF[2]=1.;
    cout << "TPZEulerConsLaw2::Contribute EK1+EK2: [0:ambas][1:EK 1][2:EK 2]";
    int par;
    cin >> par;
    if(par==0){EK[0] = 1.0; EK[1] = 1.0;} 
    if(par==1){EK[0] = 1.0; EK[1] = 0.0;}
    if(par==2){EK[0] = 0.0; EK[1] = 1.0;}
    firsttime = 0;
  }
  if(firsttime && 1){
    EK[0] = 1.0; EK[1] = 1.0;
    cout << "TPZEulerConsLaw2::Contribute EF1+EF2+EF3:"
	 << "[0:todas][1:EF 1][2:EF 2][3:EF 3]";
    int par;
    cin >> par;
    if(par==0){EF[0] = 1.0; EF[1] = 1.0; EF[2] = 1.0;} 
    if(par==1){EF[0] = 1.0; EF[1] = 0.0; EF[2] = 0.0;}
    if(par==2){EF[0] = 0.0; EF[1] = 1.0; EF[2] = 0.0;}
    if(par==3){EF[0] = 0.0; EF[1] = 0.0; EF[2] = 1.0;}
    firsttime = 0;
  }

  //int phr = phi.Rows();// phi(in, 0) = phi_in  ,  dphi(i,jn) = dphi_jn/dxi
  int nstate = NStateVariables();//3, 4 ou 5
  if(fForcingFunction) {
    //na 2a iteração deve-se ter fForcingFunction = 0
    TPZManVector<REAL> res(nstate);
    fForcingFunction(x,res);
    for(int i=0;i<nstate;i++) sol[i] = res[i];
  }
  int dim = dphi.Rows();//dx, dy ou dz
  if(Dimension() != dim)
    PZError << "TPZEulerConsLaw2::Contribute dimension error, dimension = " << dim << endl;

//   TPZDiffusionConsLaw diffusion(sol,fGamma,dim,fArtificialDiffusion);
//   TPZVec<REAL> Fx(nstate),Fy(nstate),Fz(nstate);
//   Flux(sol,Fx,Fy,Fz);
//   TPZVec<REAL> gradphi(3,0.),prodpoint(nstate),divF(nstate);
//   diffusion.GradientOfTheFlow(DF1,DF2,DF3);
//   REAL timestep = TimeStep();
//   //REAL delta = diffusion.Delta();
//   REAL delta = diffusion.DeltaOtimo(),Lpl,Hkp,Pkl;
//   TPZFMatrix divF(nstate,1),prodpoint(nstate,nstate);
//   TPZVec<REAL> sum1(nstate,0.),sum2(nstate,0.);

//   for( int in = 0; in < phr; in++ ) {

//     // w * Un
//     for(i=0;i<nstate;i++) sum1[i] = phi(in, 0) * sol[i];

//     // grad(w) . F
//     for(i=0;i<nstate;i++){
//       if(dim>0) sum2[i]  = Fx[i] * dphi(0,in);
//       if(dim>1) sum2[i] += Fy[i] * dphi(1,in);
//       if(dim>2) sum2[i] += Fz[i] * dphi(2,in);
//     }

//     //EF : w * Un + deltaT * (grad(w) . F)
//     for(i=0;i<nstate;i++)
//       ef(in * nstate + i, 0) += weight * (EF[0] * sum1[i] + EF[1] * timestep * sum2[i]);

//     //EK
//     // w * Un+1 + (grad(w) @ T) * div F(Un+1)    
//     for( int jn = 0; jn < phr; jn++ ) {
      
//       // DIFUSÃO + w * Un+1
//       for(k=0;k<nstate;k++){
// 	for(l=0;l<nstate;l++){
	  
// 	  Pkl = 0.0;
// 	  for(p=0;p<nstate;p++){
	    
// 	    if(dim>0) {Hkp  = dphi(0,jn)*Tx(k,p); Lpl  = dphi(0,in)*DF1(p,l);}
// 	    if(dim>1) {Hkp += dphi(1,jn)*Ty(k,p); Lpl += dphi(1,in)*DF2(p,l);}
// 	    if(dim>2) {Hkp += dphi(2,jn)*Tz(k,p); Lpl += dphi(2,in)*DF3(p,l);}
// 	    Pkl += Hkp * Lpl;
	    
// 	  }

// 	  // Dt * delta * (grad(w) o T) * div F (nstatex1)
// 	  ek(nstate  * in + k, nstate  * jn + l) += EK[1] * weight * timestep * delta * Pkl;

// 	  // w * Un+1 (nstatex1)
// 	  if(l == k) ek(nstate * in + k, nstate * jn + l) += EK[0] * weight * phi(in,0) * phi(jn,0);
// 	}
//       }
//     }//jn
//   }//in
}

void TPZEulerConsLaw2::TestOfRoeFlux(REAL &tetainit,REAL &tetamax,REAL &tol,REAL &increment){
/*  //Problema teste choque refletido estacionário de três estados constantes
  //OS VALORES ENCONTRADOS SÃO:
  //61 GRAUS
  //-66.7169 GRAUS (-1.16443 RADIANOS)
  //região R1
  //  REAL r1M = 2.9;//Mach
  REAL r1ro = 1.0;
  REAL r1u = 2.9;
  REAL r1v = 0.0;
  REAL r1w = 0.0;
  REAL r1p = 0.714286;
  REAL r1vel2 = r1u*r1u+r1v*r1v+r1w*r1w;
  //região R2
  //REAL r2M = 2.3781;//Mach
  REAL r2ro = 1.7;
  REAL r2u = 2.61934;
  REAL r2v = -0.50632;
  REAL r2w = 0.0;
  REAL r2p = 1.52819;
  REAL r2vel2 = r2u*r2u+r2v*r2v+r2w*r2w;
  //região R3
  //REAL r3M = 1.94235;//Mach
  REAL r3ro = 2.68728;
  REAL r3u = 2.40140;
  REAL r3v = 0.0;
  REAL r3w = 0.0;
  REAL r3p = 2.93407;
  REAL r3vel2 = r3u*r3u+r3v*r3v+r3w*r3w;
  //procurando a normal à frente estacionaria: aproximadamente 23.28 graus ~ 0.406313 radianos
  //entre as regiões R1 e R2
  //reta (cos ø , sen ø) de ângulo ø
  //(sen ø , -cos ø) normal à reta de ângulo ø apontando da região R1 para a região R2
  REAL teta=0.0,flux_rho,flux_rhou,flux_rhov,flux_rhoE,gama = 1.4;//flux_rhow,
  //REAL increment=0.00001;
  TPZVec<REAL> U1(5,0.),U2(5,0.),U3(5,0.),r1Fx(0),r1Fy(0),r1Fz(0),r2Fx(0),r2Fy(0),r2Fz(0),r3Fx(0),r3Fy(0),r3Fz(0),n(3,0.);
  int i,enter=0;
  U1[0] = r1ro;
  U1[1] = r1ro*r1u;
  U1[2] = r1ro*r1v;
  U1[3] = r1p/(gama-1.0)+r1ro*r1vel2*0.5;
  U1[4] = 0.0;//
  U2[0] = r2ro;
  U2[1] = r2ro*r2u;
  U2[2] = r2ro*r2v;
  U2[3] = r2p/(gama-1.0)+r2ro*r2vel2*0.5;
  U2[4] = 0.0;//
  U3[0] = r3ro;
  U3[1] = r3ro*r3u;
  U3[2] = r3ro*r3v;
  U3[3] = r3p/(gama-1.0)+r3ro*r3vel2*0.5;
  U3[4] = 0.0;
  //REAL tetamax = 2.0*asin(1.0)+0.1;
  REAL soma = 1.0;
  REAL suma = 1.0;
  teta = tetainit;
  while(teta < tetamax){
    teta += increment;
    n[0] = cos(teta);
    n[1] = sin(teta);
    Flux(U1,r1Fx,r1Fy,r1Fz);
    Flux(U2,r2Fx,r2Fy,r2Fz);
    Flux(U3,r3Fx,r3Fy,r3Fz);
    soma = 0.0;
    suma = 0.0;
    for(i=0;i<4;i++){
      soma += (r1Fx[i] - r2Fx[i])*n[0] + (r1Fy[i] - r2Fy[i])*n[1];
    }
    for(i=0;i<4;i++){
      suma += (r2Fx[i] - r3Fx[i])*n[0] + (r2Fy[i] - r3Fy[i])*n[1];
    }
    if(fabs(soma) < tol){
      cout << "TPZEulerConsLaw2::TestOfRoeFlux found angle, angle = " << teta << "\tradians\n";
      //a normal aponta da região R1 para a região R2: U1 é esquerdo e U2 direito
      TPZDiffusionConsLaw::Roe_Flux(U1[0],U1[1],U1[2],U1[3],
				    U2[0],U2[1],U2[2],U2[3],
				    n[0],n[1],gama,
				    flux_rho,flux_rhou,flux_rhov,flux_rhoE);
      cout << "TPZEulerConsLaw2::TestOfRoeFlux flow in the R1 region\n";
      cout << "density : " << r1Fx[0]*n[0]+r1Fy[0]*n[1] << endl
	   << "u*ro    : " << r1Fx[1]*n[0]+r1Fy[1]*n[1] << endl
	   << "v*ro    : " << r1Fx[2]*n[0]+r1Fy[2]*n[1] << endl
	   << "energy  : " << r1Fx[3]*n[0]+r1Fy[3]*n[1] << endl << endl;
      cout << "TPZEulerConsLaw2::TestOfRoeFlux flow in the R2 region\n";
      cout << "density : " << r2Fx[0]*n[0]+r2Fy[0]*n[1] << endl
	   << "u*ro    : " << r2Fx[1]*n[0]+r2Fy[1]*n[1] << endl
	   << "v*ro    : " << r2Fx[2]*n[0]+r2Fy[2]*n[1] << endl
	   << "energy  : " << r2Fx[3]*n[0]+r2Fy[3]*n[1] << endl << endl;
      cout << "TPZEulerConsLaw2::TestOfRoeFlux the calculation of the flow of Roe is\n";
      cout << "density : " << flux_rho << endl
	   << "u*ro    : " << flux_rhou << endl
	   << "v*ro    : " << flux_rhov << endl
	   << "energy  : " << flux_rhoE << endl << endl;
      cout << "\nTPZEulerConsLaw2::TestOfRoeFlux norma da diferenca |F1*n - F2*n| = " << fabs(soma) << "\n\n";
      enter = 1;
    }
    if(fabs(suma) < tol){
      cout << "TPZEulerConsLaw2::TestOfRoeFlux found angle, angle = " << teta << "\tradians\n";
      //a normal aponta da região R1 para a região R2: 
      TPZDiffusionConsLaw::Roe_Flux(U2[0],U2[1],U2[2],U2[3],
				    U3[0],U3[1],U3[2],U3[3],
				    n[0],n[1],gama,
				    flux_rho,flux_rhou,flux_rhov,flux_rhoE);
      cout << "TPZEulerConsLaw2::TestOfRoeFlux flow in the R1 region\n";
      cout << "density : " << r2Fx[0]*n[0]+r2Fy[0]*n[1] << endl
	   << "u*ro    : " << r2Fx[1]*n[0]+r2Fy[1]*n[1] << endl
	   << "v*ro    : " << r2Fx[2]*n[0]+r2Fy[2]*n[1] << endl
	   << "energy  : " << r2Fx[3]*n[0]+r2Fy[3]*n[1] << endl << endl;
      cout << "TPZEulerConsLaw2::TestOfRoeFlux flow in the R2 region\n";
      cout << "density : " << r3Fx[0]*n[0]+r3Fy[0]*n[1] << endl
	   << "u*ro    : " << r3Fx[1]*n[0]+r3Fy[1]*n[1] << endl
	   << "v*ro    : " << r3Fx[2]*n[0]+r3Fy[2]*n[1] << endl
	   << "energy  : " << r3Fx[3]*n[0]+r3Fy[3]*n[1] << endl << endl;
      cout << "TPZEulerConsLaw2::TestOfRoeFlux the calculation of the flow of Roe is\n";
      cout << "density : " << flux_rho << endl
	   << "u*ro    : " << flux_rhou << endl
	   << "v*ro    : " << flux_rhov << endl
	   << "energy  : " << flux_rhoE << endl << endl;
      cout << "\nTPZEulerConsLaw2::TestOfRoeFlux norma da diferenca |F2*n - F3*n| = " << fabs(suma) << "\n\n";
      enter = 1;
      //break;
//       if(fabs(soma) > 1.e-9){
// 	tetamax = teta + 10.0*increment;
// 	tetainit = teta - 10.0*increment;
// 	tol /= 10.0;
// 	increment /=1000.0;
// 	TestOfRoeFlux(tetainit,tetamax,tol,increment);
//       }
    }
  }
  if(enter) cout << "\nTPZEulerConsLaw2::TestOfRoeFlux angle not found, The End\n";
  else cout << "\nTPZEulerConsLaw2::TestOfRoeFlux norma da diferenca |F1*n - F2*n| = " << fabs(soma) << "\n\n";
  cout << "\nTPZEulerConsLaw2::TestOfRoeFlux concluded test\n";
  //OS VALORES ENCONTRADOS SÃO:
  //61 GRAUS 
  //-66.7169 GRAUS (-1.16443 RADIANOS)*/
}

//       if(diffright){// DIFUSÃO NA CARGA : 2
// 	// w * Un+1 
// 	for(i=0;i<nstate;i++)
// 	  ek(nstate * in + i, nstate * jn + i) += weight * phi(in,0) * phi(jn,0);
//       }

//     if(diffright){// DIFUSÃO NA CARGA : 1
//       // grad(w)
//       for(i=0;i<dim;i++) gradphi[i] = dphi(i,in);
//       // grad(w) @ T (point product)
//       diffusion.PointOperator(gradphi,prodpoint);
//       for(i=0;i<dim;i++) gradphi[i] = dsol(i,0);
//       diffusion.Divergence(gradphi,divF);
//       TPZVec<REAL> diff_term(nstate,0.);
//       for(i=0;i<nstate;i++) for(j=0;j<nstate;j++) diff_term[i] += prodpoint(i,j) * divF(j,0);
//       for(i=0;i<nstate;i++)
// 	 ef(in * nstate + i, 0) += weight * timestep * delta * diff_term[i];
//     }

