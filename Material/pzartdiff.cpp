#include "pzartdiff.h"
#include "TPZDiffusionConsLaw.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

/*REAL TPZArtDiff::fGamma = 1.4;
REAL TPZArtDiff::fDelta = 1.0;
REAL TPZArtDiff::fCFL = 0.0;*/
//char *TPZArtDiff::fArtificialDiffusion = "LS";



TPZArtDiff::TPZArtDiff(TPZArtDiffType type, REAL gamma, REAL CFL, REAL delta):
fArtDiffType(type),
fGamma(gamma),
fDelta(delta),
fCFL(CFL)
{

}

TPZArtDiff::~TPZArtDiff()
{
}

TPZArtDiffType TPZArtDiff::ArtDiffType()
{
   return fArtDiffType;
}

void TPZArtDiff::SetArtDiffType(TPZArtDiffType type)
{
   fArtDiffType = type;
}

TPZString TPZArtDiff::DiffusionName()
{
  TPZString rtstr;
  switch(fArtDiffType)
  {
     case SUPG_AD:
        rtstr = "SUPG";
	return rtstr;
     break;
     case LeastSquares_AD:
        rtstr = "LeastSquares";
	return rtstr;
     break;
     case Bornhaus_AD:
        rtstr = "Bornhaus";
	return rtstr;
     break;
     default:
     PZError << "Unknown artificial diffision term (" << fArtDiffType << ")";
     return rtstr;
  }
}

REAL TPZArtDiff::OptimalDelta()
{
   //int degree = TPZCompElDisc::gDegree;
   REAL cfl = OptimalCFL(/*degree*/);
   REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
   return delta;
}

REAL TPZArtDiff::Delta()
{
   if(fDelta)return fDelta;
   return OptimalDelta();
}

void TPZArtDiff::SetDelta(REAL delta)
{
   fDelta=delta;
}

REAL TPZArtDiff::OptimalCFL(int degree)
{
  if(fCFL) return fCFL;
  return (1.0/(2.0*(REAL)degree+1.0));
}

//-------------------A B C matrices and operations

template <class T>
void TPZArtDiff::JacobFlux(int dim, TPZVec<T> & U,TPZVec<TPZDiffMatrix<T> > &Ai)
{//OK

  Ai.Resize(dim);
  int i;
  for(i=0;i<dim;i++)Ai[i].Redim(dim+2, dim+2);

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
void TPZArtDiff::ODotOperator(TPZVec<REAL> &dphi, TPZVec<TPZVec<T> > &TauDiv, TPZVec<T> &Result){

  int dim = TauDiv.NElements();
  int size = dphi.NElements();
  int neq = TauDiv[0].NElements();
  if(size<1 || size>3){
    PZError << "TPZArtDiff::PointOperator: error data size";
  }

  Result.Resize(neq);
  Result.Fill(0.);

  int i, k;
  for(k=0;k<dim;k++)
     for(i=0;i<neq;i++)Result[i] += TauDiv[k][i] * dphi[k];
}


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

//----------------------Tau tensor

template <class T>
void TPZArtDiff::ComputeTau(int dim, TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau)
{
  Tau.Resize(dim);

  switch(fArtDiffType)
  {
     case SUPG_AD:
        SUPG(Ai, Tau);
     break;
     case LeastSquares_AD:
        LS(Ai, Tau);
     break;
     case Bornhaus_AD:
        Bornhaus(Ai, Tau);
     break;
     default:
     PZError << "Unknown artificial diffision term (" << fArtDiffType << ")";
  }
}


template <class T>
void TPZArtDiff::SUPG(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
  cout << "TPZDiffusionConsLaw:: SUPG artificial diffusion SUPG not implemented\n";
}

template <class T>
void TPZArtDiff::LS(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
  int i, dim = Ai.NElements();
  for(i = 0; i < dim; i++)
  //  Tau[i]=Ai[i];
     Ai[i].Transpose(Tau[i]);
}

template <class T>
void TPZArtDiff::Bornhaus(TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
    cout << "TPZArtDiff::Bornhaus artificial diffusion Bornhaus not implemented\n";
}



//------------------ Diff setup

template <class T>
void TPZArtDiff::PrepareDiff(int dim, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau)
{
  JacobFlux(dim, U, Ai);
  ComputeTau(dim, Ai, Tau);
}

void TPZArtDiff::PrepareFastDiff(int dim, TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv, TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv)
{
  TPZVec<TPZDiffMatrix<REAL> > Ai;
  TPZVec<TPZDiffMatrix<REAL> > Tau;

  JacobFlux(dim, sol, Ai);
  ComputeTau(dim, Ai, Tau);

  TPZVec<REAL> Div;

  TPZDiffMatrix<REAL> * pdDiv = NULL;
  TPZDiffMatrix<REAL> dDiv;
  if(pTaudDiv) pdDiv = & dDiv;

  //Computing the divergent
  Divergent(dsol, dphi, Ai, Div, pdDiv);
/*
  cout << "\n\ndiv\n" << Div;

  cout << "\n\nddiv\n" << dDiv;
*/
  TauDiv.Resize(dim);
  if(pTaudDiv)pTaudDiv->Resize(dim);

  // computing Tau.Div = {Tx.Div, Ty.Div, Tz.Div}
  // and Tau.dDiv = {Tx.dDiv, Ty.dDiv, Tz.dDiv}, if requested
  int k;
  for(k=0;k<dim;k++)
  {
     Tau[k].Multiply(Div, TauDiv[k]);
     if(pTaudDiv)Tau[k].Multiply(dDiv, pTaudDiv->operator[](k));

  }
}

#ifdef _AUTODIFF
void TPZArtDiff::PrepareFastDiff(int dim, TPZVec<FADREAL> &sol,
                 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv)
{
  TPZVec<TPZDiffMatrix<FADREAL> > Ai;
  TPZVec<TPZDiffMatrix<FADREAL> > Tau;

  JacobFlux(dim, sol, Ai);
  ComputeTau(dim, Ai, Tau);

//  #define TEST_PARTIAL_DIFF
  // Uncomment line above to zero the derivatives of
  // Tensors Ai and Tau, so that partial and complete
  // derivatives should be the same.
  #ifdef TEST_PARTIAL_DIFF
  //Zeroeing the derivatives... comparison tests for FAD and Approximate derivative methods
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
  #endif
  TPZVec<FADREAL> Div;

  //Computing the divergent
  Divergent(dsol, Ai, Div);
/*
  cout << "\n\nFADDiv \n" << Div;

  cout << "\n\nA0\n" << Ai[0];
*/
  TauDiv.Resize(dim);

  // computing Tau.Div = {Tx.Div, Ty.Div, Tz.Div}
  int k;
  for(k=0;k<dim;k++)
     Tau[k].Multiply(Div, TauDiv[k]);
}

#endif


//-----------------Contribute

void TPZArtDiff::ContributeApproxImplDiff(int dim, TPZVec<REAL> &sol, TPZFMatrix &dsol,  TPZFMatrix &dphix, TPZFMatrix &ek, TPZFMatrix &ef, REAL weight, REAL timeStep)
{
    REAL delta = Delta();
    REAL constant = /*-*/ weight * delta * timeStep;

    TPZVec<TPZVec<REAL> > TauDiv;
    TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv = NULL;
    TPZVec<TPZDiffMatrix<REAL> > TaudDiv;

    pTaudDiv = & TaudDiv;

    PrepareFastDiff(dim, sol, dsol, dphix, TauDiv, pTaudDiv);

    int i, j, k, l;
    int nshape = dphix.Cols();
    int nstate = dim + 2;
    int neq = nstate*nshape;

    REAL buff;

    // ODotProduct speeded up
    for(l=0;l<nshape;l++)
       for(i=0;i<nstate;i++)
          for(k=0;k<dim;k++)
	     {
	     buff = dphix(k,l) * constant;
	     ef(i+l*nstate,0) += buff * TauDiv[k][i];
             for(j=0;j<neq;j++)
	        ek(i+l*nstate,j) += buff * TaudDiv[k](i,j);
	     }
}

void TPZArtDiff::ContributeExplDiff(int dim, TPZVec<REAL> &sol, TPZFMatrix &dsol,  TPZFMatrix &dphix, TPZFMatrix &ef, REAL weight, REAL timeStep)
{
    REAL delta = Delta();
    REAL constant = /*-*/ weight * delta * timeStep;

    TPZVec<TPZVec<REAL> > TauDiv;

    PrepareFastDiff(dim, sol, dsol, dphix, TauDiv, NULL);

    int i, k, l;
    int nshape = dphix.Cols();
    int nstate = dim + 2;

    // ODotProduct speeded up
    for(l=0;l<nshape;l++)
       for(i=0;i<nstate;i++)
	  for(k=0;k<dim;k++)
	     ef(i+l*nstate,0) += dphix(k,l) * TauDiv[k][i] * constant;
}

#ifdef _AUTODIFF

void TPZArtDiff::ContributeImplDiff(int dim, TPZVec<FADREAL> &sol, TPZVec<FADREAL> &dsol, TPZFMatrix &ek, TPZFMatrix &ef, REAL weight,  REAL timeStep)
{
    REAL delta = Delta();
    REAL constant = /*-*/ delta * weight * timeStep;

    TPZVec<TPZVec<FADREAL> > TauDiv;

    PrepareFastDiff(dim, sol, dsol, TauDiv);

    TPZVec<FADREAL> Diff;
    TPZVec<REAL> gradv(dim);

    int i, j, k, l;
    int nstate = dim + 2;
    int neq = sol[0].size();
    int nshape = neq/nstate;

    for(l=0;l<nshape;l++)
       {
           for(k=0;k<dim;k++)
               gradv[k] = dsol[k].fastAccessDx(/*k+*/l*nstate);// always retrieving this information from the first state variable...
           ODotOperator(gradv, TauDiv, Diff);
	   for(i=0;i<nstate;i++)
	      {
	      ef(i+l*nstate,0) += constant * Diff[i].val();
	      for(j=0;j<neq;j++)
	         ek(i+l*nstate, j) += constant * Diff[i].fastAccessDx(j);
	      }
       }
}

#endif

