#include "pzartdiff.h"
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

REAL TPZArtDiff::fGamma = 1.4;
REAL TPZArtDiff::fDelta = 1.0;
REAL TPZArtDiff::fCFL = 0.0;
//char *TPZArtDiff::fArtificialDiffusion = "LS";


REAL TPZArtDiff::CFL(int degree){

  if(fCFL) return fCFL;
  return (1.0/(2.0*(REAL)degree+1.0));
}

REAL TPZArtDiff::Delta(){
  return fDelta;
}

REAL TPZArtDiff::OptimalDelta(){

  int degree = TPZCompElDisc::gDegree;
  REAL cfl = CFL(degree);
  REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
  return delta;
}

template <class T>
void TPZArtDiff::ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZDiffMatrix<T> > &M, TPZDiffMatrix<T> &Result){

  int dim = M.NElements();
  int size = dphi.NElements();
  if(size<1 || size>3){
    PZError << "TPZArtDiff::PointOperator: error data size";
  }

  Result.Redim(M[0].fRows, M[0].fCols());

  int i;
  for (i=0;i<dim;i++)Result.Add(M[i], dphi[i]);
}

template <class T>
void TPZArtDiff::ComputeTau(int dim, char * ArtDiff,
		 TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau)
{
  Tau.Resize(dim);

  if(!strcmp(ArtDiff,"SUPG")){
    SUPG(Ai, Tau);
    return;
  }
  else if(!strcmp(ArtDiff,"LS")){
    LS(Ai, Tau);
    return;
  }
  else if(!strcmp(ArtDiff,"SUPG")){
    Bornhaus(Ai, Tau);
    return;
  }else
  {
     PZError << "\nTPZArtDiff::Diff: diffusion type '" << ArtDiff << "' not implemented\n";
  }
}

template <class T>
void TPZArtDiff::PrepareDiff(int dim, char * ArtDiff, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau)
{
  JacobFlux(dim, U, Ai);
  ComputeTau(dim, ArtDiff, Ai, Tau);
}

void TPZArtDiff::PrepareFastDiff(int dim, char * ArtDiff, TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv, TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv)
{
  TPZVec<TPZDiffMatrix<REAL> > Ai;
  TPZVec<TPZDiffMatrix<REAL> > Tau;

  JacobFlux(dim, sol, Ai);
  ComputeTau(dim, ArtDiff, Ai, Tau);

  TPZVec<REAL> Div;

  TPZDiffMatrix<REAL> * pdDiv = NULL;
  TPZDiffMatrix<REAL> dDiv;
  if(pTaudDiv) pdDiv = & dDiv;

  //Computing the divergent
  Divergent(dsol, dphi, Ai, Div, pdDiv);


  cout << "\n\ndiv\n" << Div;

  cout << "\n\nddiv\n" << dDiv;

  TauDiv.Resize(dim);
  if(pTaudDiv)pTaudDiv->Resize(dim);

  // computing Tau.Div = {Tx.Div, Ty.Div, Tz.Div}
  // and Tau.dDiv = {Tx.dDiv, Ty.dDiv, Tz.dDiv}, if requested
  int k;
  for(k=0;k<dim;k++)
  {
     Tau[k].Multiply(Div, TauDiv[k]);
     Tau[k].Multiply(dDiv, pTaudDiv->operator[](k));

  }
}

#ifdef _AUTODIFF
void TPZArtDiff::PrepareFastDiff(int dim, char * ArtDiff, TPZVec<FADREAL> &sol,
                 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv)
{
  TPZVec<TPZDiffMatrix<FADREAL> > Ai;
  TPZVec<TPZDiffMatrix<FADREAL> > Tau;

  JacobFlux(dim, sol, Ai);
  ComputeTau(dim, ArtDiff, Ai, Tau);

  //Zeroeing the derivatives... comparison tests
  int i, j, t, l;
  for(t=0;t<3;t++)
     for(i=0;i<5;i++)
        for(j=0;j<5;j++)
	{ Ai[t](i,j).diff(0,30);
          Tau[t](i,j).diff(0,30);
	   for(l=0;l<30;l++)
	   {
              Ai[t](i,j).fastAccessDx(l)=0.;
	      Tau[t](i,j).fastAccessDx(l)=0.;
	   }
	}

  TPZVec<FADREAL> Div;

  //Computing the divergent
  Divergent(dsol, Ai, Div);

  cout << "\n\nFADDiv \n" << Div;
/*  cout << "\n\nAx \n" << Ai[0];
  cout << "\n\nTaux \n" << Tau[0];
*/
  TauDiv.Resize(dim);

  // computing Tau.Div = {Tx.Div, Ty.Div, Tz.Div}
  int k;
  for(k=0;k<dim;k++)
     Tau[k].Multiply(Div, TauDiv[k]);
}

#endif

void TPZArtDiff::Divergent(TPZFMatrix &dsol,
			   TPZFMatrix & dphi,
			   TPZVec<TPZDiffMatrix<REAL> > & Ai,
			   TPZVec<REAL> & Div,
			   TPZDiffMatrix<REAL> * dDiv)
{
   int nstate = Ai[0].Cols();
   int dim = nstate - 2;
   int nshape = dphi.Cols();
   Div.Resize(nstate);
   Div = 0.;

   int i, j, k;

// computing the divergent:
// A.du/dx + B.du/dy + C.du/dz
   for(k=0;k<dim;k++)
      for(i=0;i<nstate; i++)
         for(j=0;j<nstate;j++)
	    {
	       Div[i]+=Ai[k](i,j)*dsol(k,j);
	    }


   if(!dDiv)return;
// computing an approximation to the divergent derivative:

// dDiv/dUj ~= A.d2U/dUidx + B.d2U/dUidy + C.d2U/dUidz

   dDiv->Redim(nstate, nstate * nshape);
   int l;
   REAL buff;
   for(l=0;l<nshape;l++)
      for(j=0;j<nstate;j++)
         for(i=0;i<nstate; i++)
	 {
            buff =0.;
            for(k=0;k<dim;k++)
	    {
               buff+=Ai[k](i,j)*dphi(k,l);
	    }
         dDiv->operator()(i,j+l*nstate)=buff;
         }
}

#ifdef _AUTODIFF
void TPZArtDiff::Divergent(TPZVec<FADREAL> &dsol,
			   TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
			   TPZVec<FADREAL> & Div)
{
   int nstate = Ai[0].Cols();
   int dim = nstate - 2;
   Div.Resize(nstate);

   int i, j, k;

// computing the divergent:
// A.du/dx + B.du/dy + C.du/dz

   for(i=0;i<nstate; i++)
   {
      Div[i]=0.;
      for(j=0;j<nstate;j++)
	 for(k=0;k<dim;k++)
	    Div[i]+=Ai[k](i,j)*dsol[j*dim+k];
   }
}

#endif
/*
void TPZArtDiff::Diff(int dim,
		 TPZVec<TPZDiffMatrix<REAL> > & Ai,
		 TPZVec<TPZDiffMatrix<REAL> > & Tau, TPZVec<REAL> &U,
		 TPZFMatrix & dphi, const int nthphi,
		 TPZVec<REAL> & diff,
		 TPZDiffMatrix<REAL> &approxDiffDer){

  int nstate = dim+2;

  TPZDiffMatrix<REAL> ODotProduct;
  TPZVec<REAL> DPhi(dim);
  int i;
  for(i=0;i<dim;i++)DPhi[i]=dphi(i, nthphi);// aligning derivatives

  ODotOperator(DPhi, Tau, ODotProduct); //gradv [OperatorDot] Tau

  TPZVec<REAL> Divergent(nstate, 0.);
  TPZDiffMatrix<REAL> dDivergent(nstate, nstate);

  for(i=0;i<dim;i++)
     Ai[i].AddDiv(DPhi[i], U, dim, Divergent, dDivergent); //contributing to the divergent

  REAL delta = OptimalDelta();

  ODotProduct.Multiply(Divergent, diff, delta);// computing divergence vector

  ODotProduct.Multiply(dDivergent, approxDiffDer, delta);// computing approximate derivatives
}

#ifdef _AUTODIFF
void TPZArtDiff::Diff(int dim,
		TPZVec<TPZDiffMatrix<FADREAL> > & Ai,
		TPZVec<TPZDiffMatrix<FADREAL> > & Tau,
		TPZVec<FADREAL> &U,
		TPZVec<FADREAL> &dU, int nthphi,
		TPZVec<FADREAL> &diff)
{
  int nstate = dim+2;
  TPZDiffMatrix<FADREAL> ODotProduct;
  TPZVec<REAL> DPhi(dim);
  int i;
  for(i=0;i<dim;i++)DPhi[i]=dU[i].dx(nthphi);

  ODotOperator(DPhi, Tau, ODotProduct);

  TPZVec<FADREAL> Divergent(nstate, 0.);//?is this the initialization order?

  for(i=0;i<dim;i++)Ai[i].AddAlignDiv(dU, i, dim, Divergent);

  ODotProduct.Multiply(Divergent, diff, OptimalDelta());
}

#endif
*/
template <class T>
void TPZArtDiff::SUPG(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
  cout << "TPZDiffusionConsLaw:: SUPG artificial diffusion SUPG not implemented\n";
}

template <class T>
void TPZArtDiff::LS(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){

  Ai[0].Transpose(Tau[0]);
  Ai[1].Transpose(Tau[1]);
  Ai[2].Transpose(Tau[2]);
}

template <class T>
void TPZArtDiff::Bornhaus(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
    cout << "TPZArtDiff::Bornhaus artificial diffusion Bornhaus not implemented\n";
}

template <class T>
void TPZArtDiff::JacobFlux(int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai)
{//OK

  if(U[0] < 1.e-6) {
    PZError << "\nTPZArtDiff::JacobFlux: Density almost null or negative, jacobian evaluation fails\n"
       << "Density = " << U[0] << endl;
       return;
  }

  // dimension must lie between 1 and 3
  if(!(1<=dim && dim<=3)){
    PZError << "\nTPZArtDiff::JacobFlux: Unhandled dimension" << dim;
    return;
  }

  Ai.Resize(dim);
  int i;
  for(i=0;i<dim;i++)Ai[i].Redim(dim+2, dim+2);

  T    u,v,w,e;
  REAL gamma1 = fGamma-1.;
  REAL gamma2 = gamma1/2.;
  REAL gamma3 = fGamma-3;

  if(dim == 3){

    u = U[1]/U[0];
    v = U[2]/U[0];
    w = U[3]/U[0];
    e = U[4]/U[0];

    T    u2 = u*u;
    T    v2 = v*v;
    T    w2 = w*w;
    T    vel = u2+v2+w2;

    Ai[0](0,0) = 0.;
    Ai[0](0,1) = 1.;
    Ai[0](0,2) = 0.;
    Ai[0](0,3) = 0.;
    Ai[0](0,4) = 0.;

    Ai[0](1,0) =  gamma2*vel-u2;
    Ai[0](1,1) = -gamma3*u;
    Ai[0](1,2) = -gamma1*v;
    Ai[0](1,3) = -gamma1*w;
    Ai[0](1,4) =  gamma1;

    Ai[0](2,0) = -u*v;
    Ai[0](2,1) =  v;
    Ai[0](2,2) =  u;
    Ai[0](2,3) =  0.;
    Ai[0](2,4) =  0.;

    Ai[0](3,0) = -u*w;
    Ai[0](3,1) =  w;
    Ai[0](3,2) =  0.;
    Ai[0](3,3) =  u;
    Ai[0](3,4) =  0.;

    Ai[0](4,0) = -fGamma*e*u + gamma1*u*vel;
    Ai[0](4,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
    Ai[0](4,2) = -gamma1*u*v;
    Ai[0](4,3) = -gamma1*u*w;
    Ai[0](4,4) =  fGamma*u;

    Ai[1](0,0) = 0.;
    Ai[1](0,1) = 0.;
    Ai[1](0,2) = 1.;
    Ai[1](0,3) = 0.;
    Ai[1](0,4) = 0.;

    Ai[1](1,0) = -u*v;
    Ai[1](1,1) =  v;
    Ai[1](1,2) =  u;
    Ai[1](1,3) =  0.;
    Ai[1](1,4) =  0.;

    Ai[1](2,0) =  gamma2*vel-v2;
    Ai[1](2,1) = -gamma1*u;
    Ai[1](2,2) = -gamma3*v;
    Ai[1](2,3) = -gamma1*w;
    Ai[1](2,4) =  gamma1;

    Ai[1](3,0) = -v*w;
    Ai[1](3,1) =  0.;
    Ai[1](3,2) =  w;
    Ai[1](3,3) =  v;
    Ai[1](3,4) =  0.;

    Ai[1](4,0) = -fGamma*e*v + gamma1*v*vel;
    Ai[1](4,1) = -gamma1*u*v;
    Ai[1](4,2) =  fGamma*e - gamma1*v2 - gamma2*vel;
    Ai[1](4,3) = -gamma1*v*w;
    Ai[1](4,4) =  fGamma*v;

    Ai[2](0,0) = 0.;
    Ai[2](0,1) = 0.;
    Ai[2](0,2) = 0.;
    Ai[2](0,3) = 1.;
    Ai[2](0,4) = 0.;

    Ai[2](1,0) = -u*w;
    Ai[2](1,1) =  w;
    Ai[2](1,2) =  0.;
    Ai[2](1,3) =  u;
    Ai[2](1,4) =  0.;

    Ai[2](2,0) = -v*w;
    Ai[2](2,1) =  0.;
    Ai[2](2,2) =  w;
    Ai[2](2,3) =  v;
    Ai[2](2,4) =  0.;

    Ai[2](3,0) =  gamma2*vel-w2;
    Ai[2](3,1) = -gamma1*u;
    Ai[2](3,2) = -gamma1*v;
    Ai[2](3,3) = -gamma3*w;
    Ai[2](3,4) =  gamma1;

    Ai[2](4,0) = -fGamma*e*w + gamma1*w*vel;
    Ai[2](4,1) = -gamma1*u*w;
    Ai[2](4,2) = -gamma1*v*w;
    Ai[2](4,3) =  fGamma*e - gamma1*w2 - gamma2*vel;
    Ai[2](4,4) =  fGamma*w;

  } else if(dim == 2){

    u = U[1]/U[0];
    v = U[2]/U[0];
    e = U[3]/U[0];

    T    u2 = u*u;
    T    v2 = v*v;
    T    vel = u2+v2;

    Ai[0](0,0) = 0.;
    Ai[0](0,1) = 1.;
    Ai[0](0,2) = 0.;
    Ai[0](0,3) = 0.;

    Ai[0](1,0) =  gamma2*vel-u2;
    Ai[0](1,1) = -gamma3*u;
    Ai[0](1,2) = -gamma1*v;
    Ai[0](1,3) =  gamma1;

    Ai[0](2,0) = -u*v;
    Ai[0](2,1) =  v;
    Ai[0](2,2) =  u;
    Ai[0](2,3) =  0.;

    Ai[0](3,0) = -fGamma*e*u + gamma1*u*vel;
    Ai[0](3,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
    Ai[0](3,2) = -gamma1*u*v;
    Ai[0](3,3) =  fGamma*u;

    Ai[1](0,0) = 0.;
    Ai[1](0,1) = 0.;
    Ai[1](0,2) = 1.;
    Ai[1](0,3) = 0.;

    Ai[1](1,0) = -u*v;
    Ai[1](1,1) =  v;
    Ai[1](1,2) =  u;
    Ai[1](1,3) =  0.;

    Ai[1](2,0) =  gamma2*vel-v2;
    Ai[1](2,1) = -gamma1*u;
    Ai[1](2,2) = -gamma3*v;
    Ai[1](2,3) =  gamma1;

    Ai[1](3,0) = -fGamma*e*v + gamma1*v*vel;
    Ai[1](3,1) = -gamma1*u*v;
    Ai[1](3,2) =  fGamma*e - gamma1*v2 - gamma2*vel;
    Ai[1](3,3) =  fGamma*v;

  } else if(dim == 1){

    u = U[1]/U[0];
    e = U[2]/U[0];

    T    u2 = u*u;
    T    vel = u2;

    Ai[0](0,0) = 0.;
    Ai[0](0,1) = 1.;
    Ai[0](0,2) = 0.;

    Ai[0](1,0) =  gamma2*vel-u2;
    Ai[0](1,1) = -gamma3*u;
    Ai[0](1,2) =  gamma1;

    Ai[0](2,0) = -fGamma*e*u + gamma1*u*vel;
    Ai[0](2,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
    Ai[0](2,2) =  fGamma*u;
  }
}

//left = **_f    right = **_t
template <class T>
void TPZArtDiff::Roe_Flux(
	       const T & rho_f, const T & rhou_f, const T & rhov_f, const T & rhow_f,
	       const T & rhoE_f, const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhow_t,
	       const T & rhoE_t, const REAL nx, const REAL ny, const REAL nz, const REAL gam,
	       T & flux_rho, T &flux_rhou, T &flux_rhov,
	       T & flux_rhow, T &flux_rhoE){

  T    alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
  T    a1,a2,a3,a4,a5,b1,b2,b3,b4,b5;
  T    ep_t, ep_f, p_t, p_f;
  T    rhouv_t, rhouv_f, rhouw_t, rhouw_f, rhovw_t, rhovw_f;
  T    lambda_f, lambda_t;
  T    delta_rho, delta_rhou, delta_rhov, delta_rhow, delta_rhoE;
  T    hnx, hny, hnz;
  T    tempo11, usc;

  flux_rho = 0;
  flux_rhou = 0;
  flux_rhov = 0;
  flux_rhow = 0;
  flux_rhoE = 0;

  REAL gam1 = gam - 1.0;
  T    irho_f = 1.0/rho_f;
  T    irho_t = 1.0/rho_t;

  //
  //.. Compute the ROE Averages
  //
  //.... some useful quantities
  T    coef1 = sqrt(rho_f);
  T    coef2 = sqrt(rho_t);
  T    somme_coef = coef1 + coef2;
  T    isomme_coef = 1.0/somme_coef;
  T    u_f = rhou_f*irho_f;
  T    v_f = rhov_f*irho_f;
  T    w_f = rhow_f*irho_f;
  T    h_f = (gam * rhoE_f*irho_f) - (.5*gam1) * (u_f * u_f + v_f * v_f + w_f * w_f);
  T    u_t = rhou_t*irho_t;
  T    v_t = rhov_t*irho_t;
  T    w_t = rhow_t*irho_t;
  T    h_t = (gam * rhoE_t*irho_t) - (.5*gam1) * (u_t * u_t + v_t * v_t + w_t * w_t);

  //.... averages
  //REAL rho_ave = coef1 * coef2;
  T    u_ave = (coef1 * u_f + coef2 * u_t) * isomme_coef;
  T    v_ave = (coef1 * v_f + coef2 * v_t) * isomme_coef;
  T    w_ave = (coef1 * w_f + coef2 * w_t) * isomme_coef;
  T    h_ave = (coef1 * h_f + coef2 * h_t) * isomme_coef;
  //
  //.. Compute Speed of sound
  T    scal = u_ave * nx + v_ave * ny + w_ave * nz;
  T    norme = sqrt(nx * nx + ny * ny + nz * nz);
  T    inorme = 1.0/norme;
  T    u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
  T    c_speed = gam1 * (h_ave - 0.5 * u2pv2pw2);
  if(c_speed < 1e-6) c_speed = 1e-6;// <!> zeroes the derivatives?   // avoid division by 0 if critical
  c_speed = sqrt(c_speed);
  T    c_speed2 = c_speed * norme;
  //
  //.. Compute the eigenvalues of the Jacobian matrix
  T    eig_val1 = scal - c_speed2;
  T    eig_val2 = scal;
  T    eig_val3 = scal + c_speed2;
  //
  //.. Compute the ROE flux
  //.... In this part many tests upon the eigenvalues
  //.... are done to simplify calculations
  //.... Here we use the two formes of the ROE flux :
  //.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
  //.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
  //
  if(eig_val2() <= 0.0) {
    p_t = gam1 * (rhoE_t - 0.5 * (rhou_t * rhou_t +
				  rhov_t * rhov_t + rhow_t * rhow_t) * irho_t);
    ep_t = rhoE_t + p_t;
    rhouv_t = rhou_t * v_t;
    rhouw_t = rhou_t * w_t;
    rhovw_t = rhov_t * w_t;
    flux_rho  = rhou_t * nx + rhov_t * ny + rhow_t * nz;
    flux_rhou = (rhou_t * u_t + p_t) * nx + rhouv_t * ny + rhouw_t * nz;
    flux_rhov = rhouv_t * nx + (rhov_t * v_t + p_t) * ny + rhovw_t * nz;
    flux_rhow = rhouw_t * nx + rhovw_t * ny + (rhow_t * w_t + p_t) * nz;
    flux_rhoE = ep_t * (u_t * nx + v_t * ny + w_t * nz);
    //
    //.... A Entropic modification
    //
    p_f = gam1 * (rhoE_f - 0.5 * (rhou_f * rhou_f + rhov_f * rhov_f
				  + rhow_f * rhow_f) * irho_f);
    lambda_f = u_f * nx + v_f * ny + w_f * nz + norme
      * sqrt(gam * p_f * irho_f);
    lambda_t = u_t * nx + v_t * ny + w_t * nz + norme
      * sqrt(gam * p_t * irho_t);
    if ((lambda_f < 0.) && (lambda_t > 0.)) {
      eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
    }
    //
    if (eig_val3 > 0.0) {
      //.. In this case A+ is obtained by multiplying the last
      //.. colomne of T-1 with the last row of T with eig_val3                //Cedric
      delta_rho  = rho_t - rho_f;                                             //right - left
      delta_rhou = rhou_t - rhou_f;                                           //**_t  - **_f
      delta_rhov = rhov_t - rhov_f;
      delta_rhow = rhow_t - rhow_f;
      delta_rhoE = rhoE_t - rhoE_f;
      //
      scal = scal * inorme;
      hnx = nx * inorme;
      hny = ny * inorme;
      hnz = nz * inorme;
      usc = 1.0/c_speed;
      tempo11 = gam1 * usc;
      //.. Last columne of the matrix T-1
      a1 = usc;
      a2 = u_ave * usc + hnx;
      a3 = v_ave * usc + hny;
      a4 = w_ave * usc + hnz;
      a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed + scal;
      //.. Last row of the matrix T * eig_val3
      b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 - scal);
      b2 = 0.5 * (hnx - tempo11 * u_ave);
      b3 = 0.5 * (hny - tempo11 * v_ave);
      b4 = 0.5 * (hnz - tempo11 * w_ave);
      b5 = 0.5 * tempo11;
      //
      alpha1 = b1 * delta_rho;
      alpha2 = b2 * delta_rhou;
      alpha3 = b3 * delta_rhov;
      alpha4 = b4 * delta_rhow;
      alpha5 = b5 * delta_rhoE;
      alpha  = eig_val3 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
      //
      flux_rho  -= a1 * alpha;
      flux_rhou -= a2 * alpha;
      flux_rhov -= a3 * alpha;
      flux_rhow -= a4 * alpha;
      flux_rhoE -= a5 * alpha;
    }
  }
  //
  if(eig_val2 > 0.0) {
    p_f = gam1 * (rhoE_f - 0.5 * (rhou_f * rhou_f +
				  rhov_f * rhov_f + rhow_f * rhow_f) * irho_f);
    ep_f = rhoE_f + p_f;
    rhouv_f = rhou_f * v_f;
    rhouw_f = rhou_f * w_f;
    rhovw_f = rhov_f * w_f;
    flux_rho  = rhou_f * nx + rhov_f * ny + rhow_f * nz;
    flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny + rhouw_f * nz;
    flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny + rhovw_f * nz;
    flux_rhow = rhouw_f * nx + rhovw_f * ny + (rhow_f * w_f + p_f) * nz;
    flux_rhoE = ep_f * (u_f * nx + v_f * ny + w_f * nz);
    //
    // A Entropic modification
    //
    p_t = gam1 * (rhoE_t - 0.5 * (rhou_t * rhou_t +
				  + rhov_t * rhov_t + rhow_t * rhow_t) * irho_t);
    lambda_f = u_f * nx + v_f * ny + w_f * nz - norme
      * sqrt(gam * p_f * irho_f);
    lambda_t   = u_t * nx + v_t * ny + w_t * nz - norme
      * sqrt(gam * p_t * irho_t);
    if ((lambda_f < 0.) && (lambda_t > 0.)) {
      eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
    }
    //
    if (eig_val1 < 0.0) {
      //.. In this case A+ is obtained by multiplying the first
      //.. columne of T-1 with the first row of T with eig_val1
      delta_rho  = rho_t - rho_f;
      delta_rhou = rhou_t - rhou_f;
      delta_rhov = rhov_t - rhov_f;
      delta_rhow = rhow_t - rhow_f;
      delta_rhoE = rhoE_t - rhoE_f;
      //
      scal = scal * inorme;
      hnx = nx * inorme;
      hny = ny * inorme;
      hnz = nz * inorme;
      usc = 1.0/c_speed;
      tempo11 = gam1 * usc;
      //.. First colomne of the matrix T-1
      a1 = usc;
      a2 = u_ave * usc - hnx;
      a3 = v_ave * usc - hny;
      a4 = w_ave * usc - hnz;
      a5 = 0.5 * u2pv2pw2 * usc + 2.5 * c_speed - scal;
      //.. First row of the matrix T * eig_val1
      b1 = 0.5 * (0.5 * tempo11 * u2pv2pw2 + scal);
      b2 = -0.5 * (hnx + tempo11 * u_ave);
      b3 = -0.5 * (hny + tempo11 * v_ave);
      b4 = -0.5 * (hnz + tempo11 * w_ave);
      b5 = 0.5 * tempo11;
      //
      alpha1 = b1 * delta_rho;
      alpha2 = b2 * delta_rhou;
      alpha3 = b3 * delta_rhov;
      alpha4 = b4 * delta_rhow;
      alpha5 = b5 * delta_rhoE;
      alpha  = eig_val1 * (alpha1 + alpha2 + alpha3 + alpha4 + alpha5);
      //
      flux_rho  += a1 * alpha;
      flux_rhou += a2 * alpha;
      flux_rhov += a3 * alpha;
      flux_rhow += a4 * alpha;
      flux_rhoE += a5 * alpha;
    }
  }
}

//left = **_f    right = **_t
template <class T>
void TPZArtDiff::Roe_Flux(const T & rho_f, const T & rhou_f, const T & rhov_f, const T & rhoE_f,
				   const T & rho_t, const T & rhou_t, const T & rhov_t, const T & rhoE_t,
				   const REAL nx, const REAL ny, const REAL gam,
				   T & flux_rho, T & flux_rhou,T & flux_rhov, T & flux_rhoE){

  T    alpha1,alpha2,alpha3,alpha4,a1,a2,a3,a4,b1,b2,b3,b4,alpha;
  T    ep_t, ep_f, p_t, p_f;
  T    rhouv_t, rhouv_f;
  T    lambda_f, lambda_t;
  T    delta_rho, delta_rhou,delta_rhov, delta_rhoE;
  T    hnx, hny;
  T    tempo11, usc;

  flux_rho = 0;
  flux_rhou = 0;
  flux_rhov = 0;	
  flux_rhoE = 0;		
  
  REAL gam1 = gam - 1.0;
  //REAL gam2 = gam * (gam - 1.0);
  //REAL igam = 1.0 / (gam - 1.0);
  
  //
  //.. Compute the ROE Averages
  //
  //.... some useful quantities
  T    coef1 = sqrt(rho_f);
  T    coef2 = sqrt(rho_t);
  T    somme_coef = coef1 + coef2;
  T    u_f = rhou_f/rho_f;
  T    v_f = rhov_f/rho_f;
  T    h_f = (gam * rhoE_f/rho_f) - (gam1 / 2.0) * (u_f * u_f + v_f * v_f);
  T    u_t = rhou_t/rho_t;
  T    v_t = rhov_t/rho_t;
  T    h_t = (gam * rhoE_t/rho_t) - (gam1 / 2.0) * (u_t * u_t + v_t * v_t);
  
  //.... averages
  //REAL rho_ave = coef1 * coef2;
  T    u_ave = (coef1 * u_f + coef2 * u_t) / somme_coef;
  T    v_ave = (coef1 * v_f + coef2 * v_t) / somme_coef;
  T    h_ave = (coef1 * h_f + coef2 * h_t) / somme_coef;  
  //
  //.. Compute Speed of sound
  T    scal = u_ave * nx + v_ave * ny;
  T    norme = sqrt(nx * nx + ny * ny);
  T    u2pv2 = u_ave * u_ave + v_ave * v_ave;
  T    c_speed = gam1 * (h_ave - 0.5 * u2pv2);
  if(c_speed < 1e-6) c_speed = 1e-6;    // avoid division by 0 if critical
  c_speed = sqrt(c_speed);                    
  T    c_speed2 = c_speed * norme;
  //
  //.. Compute the eigenvalues of the Jacobian matrix
  T    eig_val1 = scal - c_speed2;
  T    eig_val2 = scal;
  T    eig_val3 = scal + c_speed2;
  //
  //.. Compute the ROE flux 
  //.... In this part many tests upon the eigenvalues
  //.... are done to simplify calculations                    
  //.... Here we use the two formes of the ROE flux :
  //.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
  //.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
  //
  if(eig_val2 <= 0.0) {
    p_t = gam1 * (rhoE_t - 0.5 * (rhou_t * rhou_t + rhov_t * rhov_t) / rho_t);
    ep_t = rhoE_t + p_t;
    rhouv_t = rhou_t * v_t;
    flux_rho  = rhou_t * nx + rhov_t * ny;
    flux_rhou = (rhou_t * u_t + p_t) * nx + rhouv_t * ny;
    flux_rhov = rhouv_t * nx + (rhov_t * v_t + p_t) * ny;
    flux_rhoE = ep_t * (u_t * nx + v_t * ny);
    //
    //.... A Entropic modification
    //
    p_f = gam1 * (rhoE_f - 0.5 * (rhou_f * rhou_f + rhov_f * rhov_f) / rho_f);
    lambda_f = u_f * nx + v_f * ny + norme * sqrt(gam * p_f / rho_f);
    lambda_t   = u_t * nx + v_t * ny + norme
      * sqrt(gam * p_t / rho_t);
    if ((lambda_f < 0.) && (lambda_t > 0.)) {
      eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
    }
    //           
    if (eig_val3 > 0.0) {
      //.. In this case A+ is obtained by multiplying the last
      //.. colomne of T-1 with the last row of T with eig_val3
      delta_rho  = rho_t - rho_f;
      delta_rhou = rhou_t - rhou_f;
      delta_rhov = rhov_t - rhov_f;
      delta_rhoE = rhoE_t - rhoE_f;
      //
      scal = scal / norme;
      hnx = nx / norme;
      hny = ny / norme;
      usc = 1.0/c_speed;
      tempo11 = gam1 * usc;
      //.. Last columne of the matrix T-1
      a1 = usc;
      a2 = u_ave * usc + hnx;
      a3 = v_ave * usc + hny;
      a4 = 0.5 * u2pv2 * usc + 2.5 * c_speed + scal;
      //.. Last row of the matrix T * eig_val3
      b1 = 0.5 * eig_val3 * (0.5 * tempo11 * u2pv2 - scal);
      b2 = 0.5 * eig_val3 * (hnx - tempo11 * u_ave);
      b3 = 0.5 * eig_val3 * (hny - tempo11 * v_ave);
      b4 = 0.5 * eig_val3 * tempo11;
      //
      alpha1 = a1 * b1 * delta_rho;
      alpha2 = a1 * b2 * delta_rhou;
      alpha3 = a1 * b3 * delta_rhov;
      alpha4 = a1 * b4 * delta_rhoE; 
      alpha = alpha1 + alpha2 + alpha3 + alpha4;
      //
      flux_rho  -= alpha;
      flux_rhou -= a2 * b1 * delta_rho + a2 * b2 * delta_rhou +
	           a2 * b3 * delta_rhov + a2 * b4 * delta_rhoE;
      flux_rhov -= a3 * b1 * delta_rho + a3 * b2 * delta_rhou +
	           a3 * b3 * delta_rhov + a3 * b4 * delta_rhoE;                         
      flux_rhoE -= a4 * b1 * delta_rho + a4 * b2 * delta_rhou +
	           a4 * b3 * delta_rhov + a4 * b4 * delta_rhoE;
    }
  }
  //
  if(eig_val2 > 0.0) {
    p_f = gam1 * (rhoE_f - 0.5 * (rhou_f * rhou_f +
				  + rhov_f * rhov_f) / rho_f);
    ep_f = rhoE_f + p_f;
    rhouv_f = rhou_f * v_f;
    flux_rho  = rhou_f * nx + rhov_f * ny;
    flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny;
    flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny;
    flux_rhoE = ep_f * (u_f * nx + v_f * ny);
    //
    // A Entropic modification
    //
    p_t = gam1 * (rhoE_t - 0.5 * (rhou_t * rhou_t +
				+ rhov_t * rhov_t) / rho_t);        
    lambda_f = u_f * nx + v_f * ny - norme * sqrt(gam * p_f / rho_f);
    lambda_t   = u_t * nx + v_t * ny - norme * sqrt(gam * p_t / rho_t);
    if ((lambda_f < 0.) && (lambda_t > 0.)) {
      eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
    }
    //           
    if (eig_val1 < 0.0) {
      //.. In this case A+ is obtained by multiplying the first
      //.. columne of T-1 with the first row of T with eig_val1
      delta_rho  = rho_t - rho_f;
      delta_rhou = rhou_t - rhou_f;
      delta_rhov = rhov_t - rhov_f;
      delta_rhoE = rhoE_t - rhoE_f;
      //
      scal = scal / norme;
      hnx = nx / norme;
      hny = ny / norme;
      usc = 1.0/c_speed;
      tempo11 = gam1 * usc;
      //.. First colomne of the matrix T-1
      a1 = usc;
      a2 = u_ave * usc - hnx;
      a3 = v_ave * usc - hny;
      a4 = 0.5 * u2pv2 * usc + 2.5 * c_speed - scal;
      //.. First row of the matrix T * eig_val1
      b1 = 0.5 * eig_val1 * (0.5 * tempo11 * u2pv2 + scal);
      b2 = -0.5 * eig_val1 * (hnx + tempo11 * u_ave);
      b3 = -0.5 * eig_val1 * (hny + tempo11 * v_ave);
      b4 = 0.5 * eig_val1 * tempo11;
      //
      alpha1 = a1 * b1 * delta_rho;
      alpha2 = a1 * b2 * delta_rhou;
      alpha3 = a1 * b3 * delta_rhov;
      alpha4 = a1 * b4 * delta_rhoE; 
      alpha = alpha1 + alpha2 + alpha3 + alpha4;
      //
      flux_rho  += alpha;
      flux_rhou += a2 * b1 * delta_rho + a2 * b2 * delta_rhou +
	           a2 * b3 * delta_rhov + a2 * b4 * delta_rhoE;
      flux_rhov += a3 * b1 * delta_rho + a3 * b2 * delta_rhou +
	           a3 * b3 * delta_rhov + a3 * b4 * delta_rhoE;                         
      flux_rhoE += a4 * b1 * delta_rho + a4 * b2 * delta_rhou +
	           a4 * b3 * delta_rhov + a4 * b4 * delta_rhoE;
    }
  }        
}

