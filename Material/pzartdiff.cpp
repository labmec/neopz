#include "pzartdiff.h"
#include "TPZDiffusionConsLaw.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzeulerconslaw.h"

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
void TPZArtDiff::ComputeTau(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > &Ai,
		 TPZVec<TPZDiffMatrix<T> > &Tau)
{
  Tau.Resize(dim);

  switch(fArtDiffType)
  {
     case SUPG_AD:
        SUPG(dim, sol, Ai, Tau);
     break;
     case LeastSquares_AD:
        LS(dim, sol, Ai, Tau);
     break;
     case Bornhaus_AD:
        Bornhaus(dim, sol, Ai, Tau);
     break;
     default:
     PZError << "Unknown artificial diffision term (" << fArtDiffType << ")";
  }
}


template <class T>
void TPZArtDiff::SUPG(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){

   TPZDiffMatrix<T> Rot, RotT,
			 X, Xi,
			 M, Mi,
			 Temp, INVA2B2,
			 LambdaSUPG;
   T us, c;

   TPZEulerConsLaw2::uRes(sol, us);
   TPZEulerConsLaw2::cSpeed(sol, 1.4, c);

   RotMatrix(sol, us, Rot, RotT);
   MMatrix(sol, us, fGamma, M, Mi);
   EigenSystemSUPG(sol, us, c, fGamma, X, Xi, LambdaSUPG);

   RotT.   Multiply(M, Temp);
   Temp.   Multiply(X, INVA2B2);
   INVA2B2.Multiply(LambdaSUPG, Temp);
   Temp.   Multiply(Xi, INVA2B2);
   INVA2B2.Multiply(Mi, Temp);
   Temp.   Multiply(Rot, INVA2B2);

   for(int i = 0; i < Ai.NElements();i++)
   {
      Ai[i].Multiply(INVA2B2, Tau[i]);
     //Ai[i].Multiply(INVA2B2, Temp);
     //Temp.Transpose(Tau[i]);//*/
   }
}

template <class T>
void TPZArtDiff::LS(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
  int i;
  for(i = 0; i < dim; i++)
     Tau[i] = Ai[i];
     //Ai[i].Transpose(Tau[i]);
}

template <class T>
void TPZArtDiff::Bornhaus(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
    cout << "TPZArtDiff::Bornhaus artificial diffusion Bornhaus not implemented\n";
}



//------------------ Diff setup

template <class T>
void TPZArtDiff::PrepareDiff(int dim, TPZVec<T> &U,
		 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau)
{
  TPZEulerConsLaw2::JacobFlux(fGamma, dim, U, Ai);
  ComputeTau(dim, sol, Ai, Tau);
}

void TPZArtDiff::PrepareFastDiff(int dim, TPZVec<REAL> &sol,
                 TPZFMatrix &dsol, TPZFMatrix & dphi,
		 TPZVec<TPZVec<REAL> > & TauDiv, TPZVec<TPZDiffMatrix<REAL> > * pTaudDiv)
{
  TPZVec<TPZDiffMatrix<REAL> > Ai;
  TPZVec<TPZDiffMatrix<REAL> > Tau;

  TPZEulerConsLaw2::JacobFlux(fGamma, dim, sol, Ai);
  ComputeTau(dim, sol, Ai, Tau);

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

  TPZEulerConsLaw2::JacobFlux(fGamma, dim, sol, Ai);
  ComputeTau(dim, sol, Ai, Tau);

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
	        ek(i+l*nstate,j) -= buff * TaudDiv[k](i,j);
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
	         ek(i+l*nstate, j) -= constant * Diff[i].fastAccessDx(j);
	      }
       }
}

template< class T >
void TPZArtDiff::Pressure(REAL gamma, int dim, T& press, TPZVec<T> &U)
{
   TPZEulerConsLaw2::Pressure(gamma, dim, press, U);
}

#endif
