//$Id: pzeulerconslaw.cpp,v 1.32 2004-06-13 23:33:30 erick Exp $

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
#include <pzsave.h>

#define FASTEST_IMPLICIT

//#define DIVTIMESTEP




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
  if(CFL < 0.0) CFL = OptimalCFL(degree);

  REAL deltaT = CFL*deltax/maxveloc;
  //cout << "TPZCompMesh::Delta Time : " << deltaT << endl;
  TPZConservationLaw2::SetTimeStep(deltaT);
}

int TPZEulerConsLaw2::NStateVariables(int dim) {
  return (2 + dim);//U = (rho, rhou, rhov, rhow, rhoe)
}

int TPZEulerConsLaw2::NStateVariables() {
  return NStateVariables(Dimension());//U = (rho, rhou, rhov, rhow, rhoe)
}

REAL TPZEulerConsLaw2::Pressure(TPZVec<REAL> &U)
{
   REAL press;
   TPZEulerConsLaw2::Pressure(fGamma, fDim, press, U);
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
  if( !strcmp(name,"Mach") )         return 7;//sqrt(u²+v²+w²)/c
  cout << "TPZEulerConsLaw2::VariableIndex not defined\n";
  return TPZMaterial::VariableIndex(name);
}

int TPZEulerConsLaw2::NSolutionVariables(int var){

  if(var == 1 || var == 3 || var == 4 || var == 6 || var == 7) return 1;
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

void TPZEulerConsLaw2::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){

  if(fabs(Sol[0]) < 1.e-10) {
    PZError << "\nTPZEulerConsLaw2::Solution: Density almost null\n"
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
    REAL cspeed;
    REAL us;
    TPZEulerConsLaw2::cSpeed(Sol, fGamma, cspeed);
    TPZEulerConsLaw2::uRes(Sol, us);
    Solout[0] = us / cspeed;
    return;
  } else {
    //cout << "TPZEulerConsLaw2::Solution variable in the base class\n";
    TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
  }
}


void TPZEulerConsLaw2::SetDelta(REAL delta)
{
   fArtDiff.SetDelta(delta);
}

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

template <class T>
void TPZEulerConsLaw2::PrepareFastestInterfaceFAD(
		TPZVec<REAL> &solL,TPZVec<REAL> &solR,
		TPZVec<T> & FADsolL,
		TPZVec<T> & FADsolR)
{
   int nState = solL.NElements();
   int nVars = nState * 2;

   FADsolL.Resize(nState);
   FADsolR.Resize(nState);

   for(int i = 0; i < nState; i++)
   {
      FADsolL[i] = solL[i];
      FADsolL[i].diff(i, nVars);

      FADsolR[i] = solR[i];
      FADsolR[i].diff(i + nState, nVars);
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

void TPZEulerConsLaw2::Contribute(TPZVec<REAL> &x, TPZFMatrix &jacinv,
			      TPZVec<REAL> &sol, TPZFMatrix &dsol, REAL weight,
			      TPZFMatrix &axes, TPZFMatrix &phi,
			      TPZFMatrix &dphi, TPZFMatrix &ef)
{

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
			ef);
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
   if(fResidualType == Residual_RT)
   {
      ContributeExplT2(x,sol,weight,phi,ef);
   }

   // contributing volume-based quantities
   // diffusive term
   if (fDiff == Explicit_TD)
         ContributeExplDiff(x, jacinv, sol,dsol,weight, phi, dphi, ef);

   // Volume Convective term
   if (fConvVol == Explicit_TD)
         ContributeExplConvVol(x, sol, weight, phi, dphi, ef);


}


void TPZEulerConsLaw2::ContributeAdv(TPZVec<REAL> &x,TPZFMatrix &jacinv,
			TPZVec<REAL> &sol, TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ek, TPZFMatrix &ef)
{
// contributing the implicit parcell of the residual to the
// rhs.


   if(fResidualType == Residual_RT)
   {
      // the parcell T1 is always implicit.
      ContributeImplT1(x,sol,dsol,weight, phi,dphi,ek,ef);
   }

   // contributing volume-based quantities
   // diffusive term
   if (fDiff == Implicit_TD)
   {
      // if diffusive term is implicit
      // then the FAD classes must be initialized
      #ifdef _AUTODIFF
         #ifdef FASTEST_IMPLICIT
             ContributeFastestImplDiff(fDim, x, jacinv, sol, dsol,
                                phi, dphi, weight, ek, ef);
         #else
         TPZVec<FADREAL> FADsol, FADdsol;
         PrepareFAD(sol, dsol, phi, dphi, FADsol, FADdsol);
	    ContributeImplDiff(x, jacinv, FADsol,FADdsol, weight, ek, ef);
	 #endif
      #else
         cout << "TPZEulerConsLaw2::Contribute> Implicit diffusive contribution: _AUTODIFF directive not configured -> Using an approximation to the tgMatrix";
         ContributeApproxImplDiff(x, jacinv, sol,dsol,weight,phi,dphi,ek,ef);
      #endif
   }else
   {
         if (fDiff == ApproxImplicit_TD)
            ContributeApproxImplDiff(x, jacinv, sol,dsol,weight,phi,dphi,ek,ef);
   }

   // Volume convective term
   if (fConvVol == Implicit_TD)
      ContributeImplConvVol(x,sol,dsol,weight,phi,dphi,ek,ef);
   /*}else
   {
      // Flux_RT -> contribution only to the residual vector
      if (fDiff == Implicit_TD)
         ContributeExplDiff(x, jacinv, sol,dsol,weight, phi, dphi, ef);
      if (fConvVol == Implicit_TD)
         ContributeExplConvVol(x, sol, weight, phi, dphi, ef);
   }*/
}

void TPZEulerConsLaw2::ContributeAdv(TPZVec<REAL> &x,TPZFMatrix &jacinv,
			TPZVec<REAL> &sol, TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi, TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
// contributing the implicit parcell of the residual to the
// rhs.

   if(fResidualType == Residual_RT)
   {
      ContributeExplT1(x,sol,dsol,weight, phi,dphi,ef);
   }
   // Flux_RT -> contribution only to the residual vector
   if (fDiff == Implicit_TD)
      ContributeExplDiff(x, jacinv, sol,dsol,weight, phi, dphi, ef);
   if (fConvVol == Implicit_TD)
      ContributeExplConvVol(x, sol, weight, phi, dphi, ef);

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

   // contributing face-based quantities
   if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
      {
      // if face contribution is implicit,
      // then the FAD classes must be initialized
      #ifdef _AUTODIFF
         #ifdef FASTEST_IMPLICIT
            ContributeFastestImplConvFace(fDim, x, solL, solR,
                         weight, normal, phiL, phiR, ek, ef);
	 #else
         TPZVec<FADREAL> FADsolL, FADsolR;
         PrepareInterfaceFAD(solL, solR, phiL, phiR, FADsolL, FADsolR);
         ContributeImplConvFace(x,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef);
	 #endif
      #else
      // forcing explicit contribution and issueing an warning
         cout << "TPZEulerConsLaw2::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
         ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef);
      #endif
      }

   if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
   {
      ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef);
   }
}


void TPZEulerConsLaw2::ContributeInterface(
		TPZVec<REAL> &x,
		TPZVec<REAL> &solL, TPZVec<REAL> &solR,
		TPZFMatrix &dsolL, TPZFMatrix &dsolR,
		REAL weight, TPZVec<REAL> &normal,
		TPZFMatrix &phiL,TPZFMatrix &phiR,
		TPZFMatrix &dphiL,TPZFMatrix &dphiR,
		TPZFMatrix &ef)
{

   // contributing face-based quantities
   if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT
                     ||
      fConvFace == Explicit_TD && fContributionTime == Last_CT)
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
          ek(in*nstate+i,jn*nstate+i) -= gBigNumber * weight * phi(in,0) * phi(jn,0);
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
          ek(in*nstate+i,jn*nstate+j) -= weight * bc.Val1()(i,j) * phi(in,0) * phi(jn,0);
      }
    }
  }
}


void TPZEulerConsLaw2::ContributeBCInterface(TPZVec<REAL> &x,
			TPZVec<REAL> &solL, TPZFMatrix &dsolL,
			REAL weight, TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &dphiL,
			TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
   int nstate = NStateVariables();
   TPZVec<REAL> solR(nstate,0.);
   TPZFMatrix dsolR(dsolL.Rows(), dsolL.Cols(),0.);
   TPZFMatrix phiR(0,0), dphiR(0,0);
   int entropyFix;

   // contributing face-based quantities
   if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
      {
      // if face contribution is implicit,
      // then the FAD classes must be initialized
      #ifdef _AUTODIFF
         if(bc.Type() ==5)
         {
	    int entropyFix2;
            TPZManVector<REAL,5 > flux(nstate,0.);
            ComputeGhostState(solL, solR, normal, bc, entropyFix2);
            Roe_Flux<REAL>(solL, solR, normal, fGamma, flux, entropyFix2);
            REAL norflux = flux[1]*normal[1]-flux[2]*normal[0];
            REAL err = fabs(flux[0])+fabs(norflux)+fabs(flux[3]);
            if(err > 1.e-5)
            {
              cout << "fluxo de parede errado 1 err " << err << endl;
            }
         }

         #ifdef FASTEST_IMPLICIT
            ContributeFastestBCInterface(fDim, x, solL, dsolL,
	                      weight, normal, phiL, phiR, ek, ef, bc);
	 #else
         TPZVec<FADREAL> FADsolL, FADsolR;
         PrepareInterfaceFAD(solL, solR, phiL, phiR, FADsolL, FADsolR);
	 ComputeGhostState(FADsolL, FADsolR, normal, bc, entropyFix);
         ContributeImplConvFace(x,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef, entropyFix);
	 #endif
      #else
      // forcint explicit contribution and issueing an warning
         cout << "TPZEulerConsLaw2::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
         ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef);
      #endif
      }

   if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
   {
         ComputeGhostState(solL, solR, normal, bc, entropyFix);
         ContributeExplConvFace(x,solL,solR,weight,normal,phiL,phiR,ef, entropyFix);
   }
} 


void TPZEulerConsLaw2::ContributeBCInterface(TPZVec<REAL> &x,
			TPZVec<REAL> &solL, TPZFMatrix &dsolL,
			REAL weight, TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &dphiL,
			TPZFMatrix &ef,TPZBndCond &bc)
{
   int nstate = NStateVariables();
   TPZVec<REAL> solR(nstate,0.);
   TPZFMatrix dsolR(dsolL.Rows(), dsolL.Cols(),0.);
   TPZFMatrix phiR(0,0), dphiR(0,0);
   int entropyFix;

   if(fConvFace == Implicit_TD && fContributionTime == Advanced_CT
                  ||
      fConvFace == Explicit_TD && fContributionTime == Last_CT)
   {
         ComputeGhostState(solL, solR, normal, bc, entropyFix);

         {// flux tests
            TPZManVector<REAL,3> normal2(2,0.);
            TPZManVector<REAL,5> flux2(nstate,0.);
            normal2[0] = -normal[0];
            normal2[1] = -normal[1];
            TPZManVector<REAL,5 > flux(nstate,0.);
            Roe_Flux<REAL>(solL, solR, normal, fGamma, flux);
            Roe_Flux<REAL>(solR, solL, normal2, fGamma, flux2);
            REAL fluxs = fabs(flux[0]+flux2[0])+fabs(flux[1]+flux2[1])+fabs(flux[2]+flux2[2])+fabs(flux[3]+flux2[3]);
            if(fluxs > 1.e-10)
            {
               cout << "Fluxo nao simetrico fluxs = " << fluxs << endl;
            }

            if(bc.Type() ==5)
            {
               REAL norflux = flux[1]*normal[1]-flux[2]*normal[0];
               REAL err = fabs(flux[0])+fabs(norflux)+fabs(flux[3]);
               if(err > 1.e-5)
               {
                  cout << "fluxo de parede errado 2 err " << err << endl;
                 Roe_Flux<REAL>(solL,solR,normal,fGamma,flux);
               } else
               {
                 Roe_Flux<REAL>(solL,solR,normal,fGamma,flux);
               }
            }
	 } // end of tests

         ContributeExplConvFace(x,solL,solR,weight,normal,
	                        phiL,phiR,ef,entropyFix);
   }
}

#ifdef _AUTODIFF

void TPZEulerConsLaw2::ContributeFastestBCInterface(int dim,
			TPZVec<REAL> &x,
			TPZVec<REAL> &solL, TPZFMatrix &dsolL,
			REAL weight, TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &dphiL,
			TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{

   switch(dim)
   {
      case(1):
      ContributeFastestBCInterface_dim<1>(x, solL, dsolL,
                        weight, normal,
			phiL, dphiL,
			ek, ef, bc);
      break;
      case(2):
      ContributeFastestBCInterface_dim<2>(x, solL, dsolL,
                        weight, normal,
			phiL, dphiL,
			ek, ef, bc);
      break;
      case(3):
      ContributeFastestBCInterface_dim<3>(x, solL, dsolL,
                        weight, normal,
			phiL, dphiL,
			ek, ef, bc);
      break;
      default:
      PZError << "\nTPZEulerConsLaw2::ContributeFastestBCInterface unhandled dimension\n";
      exit(-1);
   }
};


template <int dim>
void TPZEulerConsLaw2::ContributeFastestBCInterface_dim(TPZVec<REAL> &x,
			TPZVec<REAL> &solL, TPZFMatrix &dsolL,
			REAL weight, TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &dphiL,
			TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
#ifdef _TFAD
   typedef TFad<2*(dim+2), REAL> TFADREALInterface;
#endif
#ifdef _FAD
   typedef Fad<REAL> TFADREALInterface;
#endif
#ifdef _TINYFAD
   typedef TinyFad<2*(dim+2), REAL> TFADREALInterface;
#endif

   int entropyFix;

   int nstate = NStateVariables();
   TPZVec<REAL> solR(nstate,0.);
   TPZFMatrix phiR(0,0);

// initial guesses for left and right sol
   // fForcingFunction is null at iterations > 0
   if(fForcingFunction)
   {
      TPZVec<REAL> res;
      fForcingFunction(x, res);
      for(int i = 0; i < nstate; i++)
         solL[i] = solR[i] = res[i];
   }

   TPZVec<TFADREALInterface > FADsolL(nstate),
                                   FADsolR(nstate);

   PrepareFastestInterfaceFAD(solL, solR, FADsolL, FADsolR);

   ComputeGhostState(FADsolL, FADsolR, normal, bc, entropyFix);

   ContributeFastestImplConvFace_T(x, FADsolL, FADsolR,
                                 weight, normal,
				 phiL, phiR,
				 ek, ef,
				 entropyFix);
};

#endif

template <class T>
void TPZEulerConsLaw2:: ComputeGhostState(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> &normal, TPZBndCond &bc, int & entropyFix)
{
  entropyFix = 1;

  int nstate = NStateVariables();
  T vpn=0.;
  T us, un, c;
  REAL Mach, temp;
  int i;
  //Riemann Invariants
  T w1, w2, w5, uninf, usinf,
    cghost, usghost, unghost, p;
  REAL cinf;

  switch (bc.Type()){
  case 3://Dirichlet: nada a fazer a CC é a correta
    for(i=0;i<nstate; i++) solR[i] = bc.Val2().operator()(i,0);
    break;
  case 4://recuperar valor da solu¢ão MEF esquerda: saida livre
    for(i=0;i<nstate; i++) solR[i] = solL[i];
    break;
  case 5://condi¢ão de parede
    for(i=1;i<nstate-1;i++) vpn += solL[i]*T(normal[i-1]);//v.n
    for(i=1;i<nstate-1;i++) solR[i] = solL[i] - T(2.0*normal[i-1])*vpn;
    solR[0] = solL[0];
    solR[nstate-1] = solL[nstate-1];
    entropyFix = 0;
    break;
  case 6://não refletivas (campo distante)
    for(i=0;i<nstate;i++) solR[i] = solL[i];
    break;
  case 7:// INLET (Dirichlet using Mach)
    Mach = bc.Val2().operator()(1,0);
    solR[0] = bc.Val2().operator()(0,0);//solL[0];
    solR[nstate-1] = bc.Val2().operator()(nstate - 1,0);

    temp = Mach * Mach * fGamma * (fGamma - 1);
    us = sqrt(2 * temp * solR[nstate-1] /
              ( solR[0] * (2 + temp)) );

    for(i=1;i<nstate-1;i++) solR[i] = - us * normal[i-1];

/*    solR[nstate-1] = bc.Val2().operator()(nstate - 1,0)/(fGamma -1.) +
                     solR[0] * us * us / 2.;*/

    break;
  case 8:// OUTLET (Dirichlet using Mach)
    Mach = bc.Val2().operator()(1,0);
    if(bc.Val2().operator()(0,0) == 0.)
    {
      solR[0] = solL[0];
    }else
    {
      solR[0] = bc.Val2().operator()(0,0);//solL[0];
    }

    temp = Mach * Mach * fGamma * (fGamma - 1);
    us = sqrt(T(2.) * temp * solR[nstate-1] /
              ( solR[0] * (T(2.) + temp)) );

    for(i=1;i<nstate-1;i++) solR[i] = us * normal[i-1];
    break;
  case 9:// INFLOW/OUTFLOW (dedpending on direction of internal
         // velocity vector.
         // Inputs are in terms of primitive variables
         // rho, Mach, p

    // computing normal velocity and speed norm
    un = 0.;
    us = 0.;
    Mach = 0.;
    for(i = 1; i < nstate-1; i++)
    {
       un += solL[i]/solL[0]*normal[i-1];
       us += solL[i] * solL[i] / solL[0] / solL[0];
       Mach += bc.Val2()(i,0)*bc.Val2()(i,0);
    }
    us = sqrt(us);
    Mach = sqrt(Mach);

    // cinf = sqrt(gamma p / rho)
    cinf = sqrt(fGamma * bc.Val2()(nstate-1,0)/bc.Val2()(0,0));

    //computing the pressure
    // p = (gamma - 1.) * (rhoe - rho*vel^2 / 2)
    p = (fGamma - 1.) * (solL[nstate - 1] - solL[0] * us * us / T(2.));
    //c speed
    c = sqrt(fGamma * p / solL[0]);

    usinf = /*bc.Val2()(1,0)*/ Mach * cinf;
    uninf = un / us * usinf;

    if(un < 0.)// Inflow
    {
       if(/*us>c*/ Mach >= 1.)
       {//supersonic
        // all Riemann invariants retain their imposed values
          solR[0] = bc.Val2()(0,0);
          // rho vel = rho * Mach * c * u_directioni / us
          for(i = 1; i < nstate-1; i++)
             solR[i] = bc.Val2()(0,0) *
                       usinf *
		       solL[i]/solL[0] / us; // versor (direction)
          // rhoe = p / (gamma - 1) + rho * vel*vel/2
          solR[nstate-1] = bc.Val2()(nstate-1,0) / T(fGamma - 1.) +
                         bc.Val2()(0,0) * usinf * usinf / T(2.);
       }
       else
       {//subsonic
        // Invariants w1 and w2 are imposed, w5 computed
          w1 = uninf - T(2.) * cinf/ T(fGamma - 1.);
	  // Modified w2 invariant: w2 = p/rho^(gamma-1)
	  // or w2 = c^2/(gamma * rho^(gamma-1))
	  w2 = cinf * cinf / T(fGamma * pow(bc.Val2()(0,0), fGamma - 1.));
	  // w5 computed based on flow state
	  w5 = un + T(2.) * c / T(fGamma - 1.);

          // computing ghost values
	  cghost = (w5 - w1) * T((fGamma - 1.)/4.);
	  solR[0] = pow(cghost * cghost / (T(fGamma) * w2), 1./(fGamma - 1.));
	  unghost = (w1 + w5) / T(2.);
	  for(i = 1; i < nstate - 1; i++)
	     solR[i] = solR[0]  // rho
	               * unghost / un *  // scale factor
	               solL[i] / solL[0]; // element velocity component
          usghost = us / un * unghost;

	  // rhoe = rho * (c^2 / (gamma(gamma -1)) + vel^2/2)
	  solR[nstate - 1] = solR[0] * (
	                     cghost * cghost /T(fGamma * (fGamma - 1.)) +
			     usghost * usghost / T(2.));
       }
    }else
    { // Outflow
       if(us>c)
       { // supersonic: no BC at all are required
          for(i = 0; i < nstate; i++)
	     solR[i] = solL[i];
       }else
       { // subsonic outlet
         // only the condition w1 referring to the first Riemann invariant
	 // is imposed. As a rule, the imposition of pressure is applied
	 // instead.

         solR[0] = solL[0] * pow(bc.Val2()(nstate-1,0)/p, 1./fGamma);

	 cghost = sqrt(fGamma * bc.Val2()(nstate-1,0)/ solR[0]);

	 unghost = (c - cghost) * T(2./(fGamma - 1.));
         usghost = 0.;
	 // ughost = u + 2.*(c - cghost)/(fGamma - 1.) * normal
	 for(i = 1; i < nstate - 1; i++)
	 {
	    solR[i] = solR[0] * // rho
	              (solL[i] / solL[0] + // element vel
		       unghost * normal[i-1]); // ghost correction
	    usghost += solR[i] * solR[i] / solR[0] / solR[0];
	 }
	 usghost = sqrt(usghost);

         // rhoe = p / (gamma - 1) + rho * vel*vel/2
	 solR[nstate-1] = T(bc.Val2()(nstate-1,0)/(fGamma - 1.)) +
	                solR[0] * usghost * usghost / T(2.);

       }
/*
       if(fabs(val(un)) < .01*val(us))
       {
          if(un < 0.)cout << "\ntangent inlet";
	  if(un > 0.)cout << "\ntangent outlet";
	  if(un == 0.) cout << "\n tangent pure";
       }
*/
    }
  break;


  case 10:// Directional INFLOW
         // velocity vector.
         // Inputs are in terms of primitive variables
         // rho, Machx, Machy, Machz, p

    // computing normal velocity and speed norm
    un = 0.;
    us = 0.;
    Mach = 0.;
    for(i = 1; i < nstate-1; i++)
    {
       un += solL[i]/solL[0]*normal[i-1];
       us += solL[i] * solL[i] / solL[0] / solL[0];
       Mach += bc.Val2()(i,0)*bc.Val2()(i,0);
    }
    us = sqrt(us);
    Mach = sqrt(Mach);

    // cinf = sqrt(gamma p / rho)
    cinf = sqrt(fGamma * bc.Val2()(nstate-1,0)/bc.Val2()(0,0));

    //computing the pressure
    // p = (gamma - 1.) * (rhoe - rho*vel^2 / 2)
    p = (fGamma - 1.) * (solL[nstate - 1] - solL[0] * us * us / T(2.));
    //c speed
    c = sqrt(fGamma * p / solL[0]);

    usinf = /*bc.Val2()(1,0)*/ Mach * cinf;
    uninf = un / us * usinf;

    if(un < 0.)// Inflow
    {
       if(/*us>c*/ Mach >= 1.)
       {//supersonic
        // all Riemann invariants retain their imposed values
          solR[0] = bc.Val2()(0,0);
          // rho vel = rho * Mach * c * u_directioni / us
          for(i = 1; i < nstate-1; i++)
             solR[i] = bc.Val2()(0,0) *
                       usinf *
		       solL[i]/solL[0] / us; // versor (direction)
          // rhoe = p / (gamma - 1) + rho * vel*vel/2
          solR[nstate-1] = bc.Val2()(nstate-1,0) / T(fGamma - 1.) +
                         bc.Val2()(0,0) * usinf * usinf / T(2.);
       }
       else
       {//subsonic
        // Invariants w1 and w2 are imposed, w5 computed
          w1 = uninf - T(2.) * cinf/ T(fGamma - 1.);
	  // Modified w2 invariant: w2 = p/rho^(gamma-1)
	  // or w2 = c^2/(gamma * rho^(gamma-1))
	  w2 = cinf * cinf / T(fGamma * pow(bc.Val2()(0,0), fGamma - 1.));
	  // w5 computed based on flow state
	  w5 = un + T(2.) * c / T(fGamma - 1.);

          // computing ghost values
	  cghost = (w5 - w1) * T((fGamma - 1.)/4.);
	  solR[0] = pow(cghost * cghost / (T(fGamma) * w2), 1./(fGamma - 1.));
	  unghost = (w1 + w5) / T(2.);
	  usghost = us / un * unghost;
	  for(i = 1; i < nstate - 1; i++)
	     solR[i] = solR[0]  // rho
	               * usghost *  // velocity
	               bc.Val2()(i,0) / Mach; // element velocity component


	  // rhoe = rho * (c^2 / (gamma(gamma -1)) + vel^2/2)
	  solR[nstate - 1] = solR[0] * (
	                     cghost * cghost /T(fGamma * (fGamma - 1.)) +
			     usghost * usghost / T(2.));
       }
    }
  break;

  case 11:// SOLID WALL - No input needed

    entropyFix = 0;
    // Invariants w1 and w2 are imposed, w5 computed
    // This gives a set of closed BCs.

    // computing normal velocity and speed norm
    un = 0.;
    us = 0.;
    for(i = 1; i < nstate-1; i++)
    {
       un += solL[i]/solL[0]*normal[i-1];
       us += solL[i] * solL[i] / solL[0] / solL[0];
    }
    us = sqrt(us);
    //computing the pressure
    // p = (gamma - 1.) * (rhoe - rho*vel^2 / 2)
    p = (fGamma - 1.) * (solL[nstate - 1] - solL[0] * us * us / T(2.));
    //c speed
    c = sqrt(fGamma * p / solL[0]);

    // computing the ghost sound speed
    // cghost = c + (gamma - 1) * un / 2
    cghost = c + T((fGamma - 1.)/2.)*un;
    // ghost density
    // rhoghost = (cghost^2*rho^gamma/(gamma * p))^(1/gamma -1)
    solR[0] = pow(cghost * cghost * pow(solL[0],fGamma)/ (p * fGamma), 1./(fGamma - 1.));

    // computing velocity vector
    usghost = 0.;
    for(i = 1; i < nstate-1; i++)
    {
       // vghost = u - (un * n)
       solR[i] = solR[0] * (solL[i]/ solL[0] - un * normal[i-1]);
       usghost += solR[i] * solR[i] / solR[0] / solR[0];
    }
    usghost = sqrt(usghost);

    // rhoe = rho * (c^2 / (gamma(gamma -1)) + vel^2/2)
    solR[nstate - 1] = solR[0] * (
                  cghost * cghost /T(fGamma * (fGamma - 1.)) +
                  usghost * usghost / T(2.));
  break;
  default:
    for(i=0;i<nstate;i++) solR[i] = 0.;
  }
}

//------------------internal contributions

void TPZEulerConsLaw2::ContributeApproxImplDiff(TPZVec<REAL> &x,
			TPZFMatrix &jacinv,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   fArtDiff.ContributeApproxImplDiff(fDim, jacinv, sol, dsol, dphi,
			ek, ef, weight,
			#ifdef DIVTIMESTEP
			   1.
			#else
			   TimeStep()
			#endif
			);
}

void TPZEulerConsLaw2::ContributeExplDiff(TPZVec<REAL> &x,
			TPZFMatrix &jacinv,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
   fArtDiff.ContributeExplDiff(fDim, jacinv, sol, dsol, dphi,
			ef, weight,
			#ifdef DIVTIMESTEP
			   1.
			#else
			   TimeStep()
			#endif
			);
}


#ifdef _AUTODIFF
void TPZEulerConsLaw2::ContributeImplDiff(TPZVec<REAL> &x,
			TPZFMatrix &jacinv,
			TPZVec<FADREAL> &sol,TPZVec<FADREAL> &dsol,
			REAL weight,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   fArtDiff.ContributeImplDiff(fDim, jacinv, sol, dsol,
			ek, ef, weight,
			#ifdef DIVTIMESTEP
			   1.
			#else
			   TimeStep()
			#endif
			);
}


void TPZEulerConsLaw2::ContributeFastestImplDiff(int dim, TPZVec<REAL> &x, TPZFMatrix &jacinv, TPZVec<REAL> &sol, TPZFMatrix &dsol, TPZFMatrix &phi, TPZFMatrix &dphi, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
      fArtDiff.ContributeFastestImplDiff(dim, jacinv, sol, dsol, phi, dphi,
			ek, ef, weight,
			#ifdef DIVTIMESTEP
			   1.
			#else
			   TimeStep()
			#endif
			);
}

#endif

void TPZEulerConsLaw2::ContributeExplConvFace(TPZVec<REAL> &x,
			TPZVec<REAL> &solL,TPZVec<REAL> &solR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ef, int entropyFix)
{
   int nState = NStateVariables();
   TPZVec<REAL > flux(nState,0.);
   Roe_Flux<REAL>(solL, solR, normal, fGamma, flux, entropyFix);
   int nShapeL = phiL.Rows(),
       nShapeR = phiR.Rows(),
       i_shape, i_state;

   #ifdef DIVTIMESTEP
   REAL constant = weight; // weight
   #else
   REAL constant = TimeStep() * weight; // deltaT * weight
   #endif


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
			TPZFMatrix &ek,TPZFMatrix &ef,
			int entropyFix)
{
   int nState = NStateVariables();
   TPZVec<FADREAL > flux(nState,0.);
   Roe_Flux(solL, solR, normal, fGamma, flux, entropyFix);
   
   // Testing whether Roe_Flux<REAL> gives the same result as Roe_Flux<FADREAL>
   
/*   TPZVec<REAL> solL2(nState,0.),solR2(nState,0.),flux2(nState,0.);
   int i;
   for(i=0; i<nState; i++)
   {
      solL2[i] = solL[i].val();
      solR2[i] = solR[i].val();
   }
   Roe_Flux(solL2, solR2, normal, fGamma, flux2,with_entropy_fix);
   REAL diff = fabs(flux[0].val()-flux2[0])+fabs(flux[1].val()-flux2[1])+fabs(flux[2].val()-flux2[2]);
   if(diff != 0.)
   {
      cout << "Roe<FADREAL> is different from Roe<REAL> diff " << diff << endl;
   }*/
   int nShapeL = phiL.Rows(),
       nShapeR = phiR.Rows(),
       i_shape, i_state, j,
       nDer = (nShapeL + nShapeR) * nState;


   #ifdef DIVTIMESTEP
   REAL constant = weight; // weight
   #else
   REAL constant = TimeStep() * weight; // deltaT * weight
   #endif

   // Contribution referring to the left element
   for(i_shape = 0; i_shape < nShapeL; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_shape*nState + i_state;
         ef(index,0) +=
	    flux[i_state].val() * phiL(i_shape,0) * constant;
	 for(j = 0; j < nDer; j++)
	    ek(index, j) -= flux[i_state].dx/*fastAccessDx*/(j) *
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
	    ek(index, j) += flux[i_state].dx/*fastAccessDx*/(j) *
	       phiR(i_shape,0) * constant;
      }
}

void TPZEulerConsLaw2::ContributeFastestImplConvFace(int dim,
			TPZVec<REAL> &x,
			TPZVec<REAL> &solL,TPZVec<REAL> &solR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ek,TPZFMatrix &ef,
			int entropyFix)
{
   switch(dim)
   {
      case(1):
      ContributeFastestImplConvFace_dim<1>(x, solL, solR, weight, normal,
			phiL, phiR, ek, ef, entropyFix);
      break;
      case(2):
      ContributeFastestImplConvFace_dim<2>(x, solL, solR, weight, normal,
			phiL, phiR, ek, ef, entropyFix);
      break;
      case(3):
      ContributeFastestImplConvFace_dim<3>(x, solL, solR, weight, normal,
			phiL, phiR, ek, ef, entropyFix);
      break;
      default:
      PZError << "\nTPZEulerConsLaw2::ContributeFastestImplConvFace unhandled dimension\n";
      exit(-1);
   }
}

template <int dim>
void TPZEulerConsLaw2::ContributeFastestImplConvFace_dim(
			TPZVec<REAL> &x,
			TPZVec<REAL> &solL,TPZVec<REAL> &solR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ek,TPZFMatrix &ef,
			int entropyFix)
{
#ifdef _TFAD
   typedef TFad<2*(dim+2), REAL> TFADREALInterface;
#endif
#ifdef _FAD
   typedef Fad<REAL> TFADREALInterface;
#endif
#ifdef _TINYFAD
   typedef TinyFad<2*(dim+2), REAL> TFADREALInterface;
#endif

   int nstate = NStateVariables(dim);
   TPZVec< TFADREALInterface > FADsolL(nstate),
                                   FADsolR(nstate);
   PrepareFastestInterfaceFAD(solL, solR, FADsolL, FADsolR);

   ContributeFastestImplConvFace_T(x, FADsolL, FADsolR, weight, normal,
			phiL, phiR, ek, ef, entropyFix);
}

template <class T>
void TPZEulerConsLaw2::ContributeFastestImplConvFace_T(TPZVec<REAL> &x,
			TPZVec<T> &FADsolL,TPZVec<T> &FADsolR,
			REAL weight,TPZVec<REAL> &normal,
			TPZFMatrix &phiL,TPZFMatrix &phiR,
			TPZFMatrix &ek,TPZFMatrix &ef,
			int entropyFix)
{
  const int nState = NStateVariables();

   TPZVec<T> FADflux(nState);

   int nShapeL = phiL.Rows(),
       nShapeR = phiR.Rows(),
       i_shape, i_state,  j, k,
       nDerL = nShapeL * nState;


   #ifdef DIVTIMESTEP
   REAL constant = weight; // weight
   #else
   REAL constant = TimeStep() * weight; // deltaT * weight
   #endif


   Roe_Flux(FADsolL, FADsolR, normal, fGamma, FADflux, entropyFix);

   // Contribution referring to the left element
   for(i_shape = 0; i_shape < nShapeL; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = i_shape*nState + i_state;
         ef(index,0) +=
	    FADflux[i_state].val() * phiL(i_shape,0) * constant;
	 for(k = 0; k < nState; k++)
	 {
	    for(j = 0; j < nShapeL; j++)
	       ek(index, j * nState + k) -=
	                       FADflux[i_state].dx/*fastAccessDx*/(k) *
	                       phiL(j) * //df/dUl
	                       phiL(i_shape,0) * //test function
			       constant;
	    for(j = 0; j < nShapeR; j++)
	       ek(index, j*nState + k + nDerL) -=
	                       FADflux[i_state]./*fastAccessDx*/dx(k + nState) *
	                       phiR(j) * //df/dUl
	                       phiL(i_shape,0) * //test function
			       constant;
	  }
      }

   // Contribution referring to the right element
   // REM: The contributions are negative in comparison
   // to the left elements: opposed normals
   for(i_shape = 0; i_shape < nShapeR; i_shape ++)
      for(i_state = 0; i_state < nState; i_state++)
      {
         int index = (nShapeL + i_shape) * nState + i_state;
         ef(index,0) -=
	    FADflux[i_state].val() * phiR(i_shape,0) * constant;
	 for(k = 0; k < nState; k++)
	 {
	    for(j = 0; j < nShapeL; j++)
	       ek(index, j * nState + k) +=
	                       FADflux[i_state]./*fastAccessDx*/dx(k) *
	                       phiL(j) * //df/dUl
	                       phiR(i_shape,0) * //test function
			       constant;
	    for(j = 0; j < nShapeR; j++)
	       ek(index, j * nState + k + nDerL) +=
	                       FADflux[i_state].dx/*fastAccessDx*/(k + nState) *
	                       phiR(j) * //df/dUl
	                       phiR(i_shape,0) * //test function
			       constant;
	 }
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


   #ifdef DIVTIMESTEP
   REAL constant = -weight; // -weight
   #else
   REAL constant = -TimeStep() * weight; // -deltaT * weight
   #endif

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

}


void TPZEulerConsLaw2::ContributeImplConvVol(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   TPZVec< TPZVec<REAL> > F(3);
   Flux(sol, F[0], F[1], F[2]);

   #ifdef DIVTIMESTEP
   REAL constant = -weight; // -weight
   #else
   REAL constant = -TimeStep() * weight; // -deltaT * weight
   #endif

   TPZVec< TPZDiffMatrix<REAL> > Ai(3);
   JacobFlux(fGamma, fDim, sol, Ai);
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
                   ek(index, j_state + j_shape * nState) -=
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


void TPZEulerConsLaw2::ContributeExplT1(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ef)
{
   if(TimeStep()==0.)return;

   int i_shape, ij_state;
   int nState = NStateVariables();
   int nShape = phi.Rows();

   #ifdef DIVTIMESTEP
   REAL constant = weight / TimeStep();
   #else
   REAL constant = weight;
   #endif

   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(ij_state = 0; ij_state < nState; ij_state++)
      {
          int index = i_shape * nState + ij_state;
	  // ef += sol*phi(i)
          ef(index, 0) +=
              sol[ij_state] * phi(i_shape,0) *
              constant;
      }
}

void TPZEulerConsLaw2::ContributeImplT1(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,TPZFMatrix &dsol,
			REAL weight,
			TPZFMatrix &phi,TPZFMatrix &dphi,
			TPZFMatrix &ek,TPZFMatrix &ef)
{
   if(TimeStep()==0.)return;

   int i_shape, ij_state, j_shape;
   int nState = NStateVariables();
   int nShape = phi.Rows();

   #ifdef DIVTIMESTEP
   REAL constant = weight / TimeStep();
   #else
   REAL constant = weight;
   #endif

   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(ij_state = 0; ij_state < nState; ij_state++)
      {
          int index = i_shape * nState + ij_state;
	  // ef += sol*phi(i)
          ef(index, 0) +=
              sol[ij_state] * phi(i_shape,0) *
              constant;
	  // ek += phi(i)*phi(j)
          for(j_shape = 0; j_shape < nShape; j_shape++)
             ek(index, j_shape * nState + ij_state) -=
                phi(i_shape,0) *
                phi(j_shape,0) *
                constant;
      }
}

void TPZEulerConsLaw2::ContributeExplT2(TPZVec<REAL> &x,
			TPZVec<REAL> &sol,
			REAL weight,
			TPZFMatrix &phi,
			TPZFMatrix &ef)
{
   if(TimeStep()==0.)return;

   int i_shape, i_state;
   int nState = NStateVariables();
   int nShape = phi.Rows();

   #ifdef DIVTIMESTEP
   REAL constant = weight / TimeStep();
   #else
   REAL constant = weight;
   #endif

   for(i_shape = 0; i_shape < nShape; i_shape++)
      for(i_state = 0; i_state < nState; i_state++)
      {
          int index = i_shape * nState + i_state;
	  // ef += sol*phi(i)
          ef(index, 0) -=
              sol[i_state] * phi(i_shape,0) *
              constant; // the T2 parcell is negative
      }
}


void TPZEulerConsLaw2::Write(TPZStream &buf, int withclassid)
{
   TPZSaveable::Write(buf, 1);
   TPZConservationLaw2::Write(buf, 0);
   fArtDiff.Write(buf, 0);
   int tmp = static_cast < int > (fDiff);
   buf.Write(& tmp,1);
   tmp = static_cast<int>(fConvVol);
   buf.Write(& tmp,1);
   tmp = static_cast<int>(fConvFace);
   buf.Write(&tmp,1);
}

void TPZEulerConsLaw2::Read(TPZStream &buf, void *context)
{
   TPZSaveable::Read(buf, context);
   TPZConservationLaw2::Read(buf, context);
   fArtDiff.Read(buf, context);
   int diff;
   buf.Read(&diff,1);
   fDiff = static_cast<TPZTimeDiscr>(diff);
   buf.Read(&diff,1);
   fConvVol = static_cast<TPZTimeDiscr>(diff);
   buf.Read(&diff,1);
   fConvFace = static_cast<TPZTimeDiscr>(diff);
}

