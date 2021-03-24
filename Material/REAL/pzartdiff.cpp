/**
 * @file
 * @brief Contains implementations of the TPZArtDiff methods.
 */

#include "pzartdiff.h" 
#include "TPZDiffusionConsLaw.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"
#include "pzeulerconslaw.h"

#define FASTESTDIFF

TPZArtDiff::TPZArtDiff(TPZArtDiffType type, REAL gamma, REAL CFL, REAL delta):
TPZRegisterClassId(&TPZArtDiff::ClassId),
fArtDiffType(type),
fGamma(gamma),
fDelta(delta),
fCFL(CFL)
{
	
}

TPZArtDiff::TPZArtDiff():
TPZRegisterClassId(&TPZArtDiff::ClassId),
fArtDiffType(None_AD),
fGamma(1.4),
fDelta(0.),
fCFL(1.)
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
			//     break;
		case LeastSquares_AD:
			rtstr = "LeastSquares";
			return rtstr;
			//     break;
		case Bornhaus_AD:
			rtstr = "Bornhaus";
			return rtstr;
			//     break;
		case TrnLeastSquares_AD:
			rtstr = "TrnLeastSquares";
			return rtstr;
			//    break;
		default:
			PZError << "Unknown artificial diffusion term (" << fArtDiffType << ")";
			return rtstr;
	}
}

REAL TPZArtDiff::OptimalDelta()
{
	REAL cfl = OptimalCFL(/*degree*/);
	REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
	return delta;
}

REAL TPZArtDiff::Delta(REAL deltax, TPZVec<STATE> & sol)
{
	REAL dX = 0.;
	
	if(fDelta > 0.)dX = deltax * fDelta;
	if(fDelta < 0.)dX =-fDelta;
	if(fDelta == 0.)dX = OptimalDelta();
	
	STATE us, c, lambdaMax;
	
	switch(fArtDiffType)
	{
		case SUPG_AD:
		case Bornhaus_AD:
			return dX / 2.;
			//     break;
		case LeastSquares_AD:
		case TrnLeastSquares_AD:
			TPZEulerConsLaw::uRes(sol, us);
			TPZEulerConsLaw::cSpeed(sol, fGamma, c);
			lambdaMax = us+c;
			return dX /2. / lambdaMax;
		default:
			PZError << "Unknown artificial diffusion term (" << fArtDiffType << ")";
			return 0.;
	}
}

void TPZArtDiff::SetDelta(REAL delta)
{
	fDelta=delta;
}

REAL TPZArtDiff::OptimalCFL(int degree)
{
	if(fCFL != 0.) return fCFL;
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
	Result.Fill(REAL(0.));
	
	int i, k;
	for(k=0;k<dim;k++)
		for(i=0;i<neq;i++)Result[i] += TauDiv[k][i] * dphi[k];
}


template <class T>
void TPZArtDiff::Divergent(TPZFMatrix<STATE> &dsol,
						   TPZFMatrix<REAL> & phi,
						   TPZFMatrix<REAL> & dphi,
						   TPZVec<TPZDiffMatrix<T> > & Ai,
						   TPZVec<STATE> & Div,
						   TPZDiffMatrix<STATE> * dDiv)
{
	int nstate = Ai[0].Cols();
	int dim = nstate - 2;
	int nshape = dphi.Cols();
	Div.Resize(nstate);
	Div = (STATE(0.));
	
	int i, j, k;
	
	// computing the divergent:
	// A.du/dx + B.du/dy + C.du/dz
	for(k=0;k<dim;k++)
		for(i=0;i<nstate; i++)
			for(j=0;j<nstate;j++)
			{
				Div[i] += Ai[k](i,j).val() * dsol(k,j);
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
					buff+=Ai[k](i,j).val()*dphi(k,l);
				}
				dDiv->operator()(i,j+l*nstate)=buff;
			}
	
	TPZVec<T> ADiv(nstate);
	T temp;
	for( k = 0; k < dim; k++)
	{
		//Ai[k].Multiply(Div, ADiv);
		for(i = 0; i < nstate; i++)
		{
			temp = T(0.);
			for(j = 0; j < nstate; j++)
				temp += Ai[k](i,j) * (T)dsol(k,j);//[j];
			ADiv[i] = temp;
		}
		for(l=0;l<nshape;l++)
			for(j=0;j<nstate;j++)
				for(i=0;i<nstate; i++)
					dDiv->operator()(i,j+l*nstate) +=
					ADiv[i]./*fastAccessDx*/dx(j) * phi(l,0);
	}
}


void TPZArtDiff::Divergent(TPZFMatrix<STATE> &dsol,
						   TPZFMatrix<REAL> & dphi,
						   TPZVec<TPZDiffMatrix<STATE> > & Ai,
						   TPZVec<STATE> & Div,
						   TPZDiffMatrix<STATE> * dDiv)
{
	int nstate = Ai[0].Cols();
	int dim = nstate - 2;
	int nshape = dphi.Cols();
	Div.Resize(nstate);
	Div.Fill(0.);
	
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


//----------------------Tau tensor

template <class T>
void TPZArtDiff::ComputeTau(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > &Ai,
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
			Bornhaus(dim, jacinv, sol, Ai, Tau);
			break;
		case TrnLeastSquares_AD:
			LST(dim, sol, Ai, Tau);
			break;
		default:
			PZError << "Unknown artificial diffision term (" << fArtDiffType << ")";
	}
}


template <class T>
void TPZArtDiff::SUPG(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
	
#ifdef FASTESTDIFF
	
	TPZDiffMatrix<T> RTM, RMi,
	X, Xi,
	Temp, INVA2B2,
	LambdaSUPG;
	T us, c;
	
	TPZEulerConsLaw::uRes(sol, us);
	TPZEulerConsLaw::cSpeed(sol, 1.4, c);
	
	RMMatrix(sol, us, fGamma, RTM, RMi);
	
	EigenSystemSUPG(sol, us, c, fGamma, X, Xi, LambdaSUPG);
	
	
	RTM.     Multiply(X, Temp);
	Temp.   Multiply(LambdaSUPG, INVA2B2);
	INVA2B2.Multiply(Xi, Temp);
	Temp.   Multiply(RMi, INVA2B2);
	
	for(int i = 0; i < Ai.NElements();i++)
	{
		Ai[i].Multiply(INVA2B2, Tau[i]);
	}
#else
	
	TPZDiffMatrix<T> Rot, RotT,
	X, Xi,
	M, Mi,
	Temp, INVA2B2,
	LambdaSUPG;
	T us, c;
	
	TPZEulerConsLaw::uRes(sol, us);
	TPZEulerConsLaw::cSpeed(sol, 1.4, c);
	
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
	}
	
#endif
}

template <class T>
void TPZArtDiff::LS(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
	int i;
	for(i = 0; i < dim; i++)
	{
		Ai[i].Transpose(Tau[i]);
	}
}

template <class T>
void TPZArtDiff::LST(int dim, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
	int i;
	for(i = 0; i < dim; i++)
	{
		Tau[i] = Ai[i];
	}
}

template <class T>
void TPZArtDiff::Bornhaus(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<T> & sol, TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau){
	
#ifdef FASTESTDIFF
	
	int i, j;
	int nstate = Ai[0].Rows();
	
	TPZDiffMatrix<T> RTM, RMi,
	Y, Yi,
	Temp, Temp2,
	BornhausTau(nstate, nstate),
	LambdaBornhaus;
	T us, c;
	TPZVec<REAL> alphas(dim,0.);
	
	TPZEulerConsLaw::uRes(sol, us);
	TPZEulerConsLaw::cSpeed(sol, 1.4, c);
	
	RMMatrix(sol, us, fGamma, RTM, RMi);
	
	BornhausTau.Redim(nstate, nstate);
	
	for( i = 0; i < dim; i++)
	{
		for( j = 0; j < dim; j++)
			alphas[j] = jacinv(i,j);
		ContributeBornhaus(sol, us, c, fGamma, alphas, BornhausTau);
	}
	
	RTM. Multiply(BornhausTau, Temp);
	Temp.Multiply(RMi, BornhausTau);
	
	for( i = 0; i < dim;i++)
	{
		Ai[i].Multiply(BornhausTau, Tau[i]);
	}
	
#else
	
	int i, j;
	int nstate = Ai[0].Rows();
	
	TPZDiffMatrix<T> Rot, RotT,
	Y, Yi,
	M, Mi,
	Temp, Temp2,
	BornhausTau(nstate, nstate),
	LambdaBornhaus;
	T us, c;
	TPZVec<REAL> alphas(dim,0.);
	
	TPZEulerConsLaw::uRes(sol, us);
	TPZEulerConsLaw::cSpeed(sol, 1.4, c);
	
	RotMatrix(sol, us, Rot, RotT);
	MMatrix(sol, us, fGamma, M, Mi);
	
	for( i = 0; i < dim; i++)
	{
		for( j = 0; j < dim; j++)
			alphas[j] = jacinv(i,j);
		
		EigenSystemBornhaus(sol, us, c, fGamma, alphas,
							Y, Yi, LambdaBornhaus);
		
		Y.   Multiply(LambdaBornhaus, Temp);
		Temp.Multiply(Yi, Temp2);
		
		BornhausTau.Add(Temp2);
	}
	
	RotT.   Multiply(M, Temp);
	Temp.   Multiply(BornhausTau, Temp2);
	Temp2.  Multiply(Mi, Temp);
	Temp.   Multiply(Rot, BornhausTau);
	
	BornhausTau.Inverse();
	for( i = 0; i < dim;i++)
	{
		Ai[i].Multiply(BornhausTau, Tau[i]);
	}
#endif
	
}

//------------------ Diff setup
template <class T>
void TPZArtDiff::PrepareDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<T> &U,
							 TPZVec<TPZDiffMatrix<T> > & Ai, TPZVec<TPZDiffMatrix<T> > & Tau)
{
	TPZEulerConsLaw::JacobFlux(fGamma, dim, U, Ai);
	ComputeTau(dim, jacinv, U, Ai, Tau);
}

void TPZArtDiff::PrepareFastDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol,
								 TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> & dphi,
								 TPZVec<TPZVec<STATE> > & TauDiv, TPZVec<TPZDiffMatrix<STATE> > * pTaudDiv)
{
	TPZVec<TPZDiffMatrix<STATE> > Ai;
	TPZVec<TPZDiffMatrix<STATE> > Tau;
	
	TPZEulerConsLaw::JacobFlux(fGamma, dim, sol, Ai);
	ComputeTau(dim, jacinv, sol, Ai, Tau);
	
	TPZVec<STATE> Div;
	
	TPZDiffMatrix<STATE> * pdDiv = NULL;
	TPZDiffMatrix<STATE> dDiv;
	if(pTaudDiv) pdDiv = & dDiv;
	
	//Computing the divergent
	Divergent(dsol, dphi, Ai, Div, pdDiv);
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

void TPZArtDiff::PrepareFastDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<FADREAL> &sol,
								 TPZVec<FADREAL> &dsol, TPZVec<TPZVec<FADREAL> > & TauDiv)
{
	TPZVec<TPZDiffMatrix<FADREAL> > Ai;
	TPZVec<TPZDiffMatrix<FADREAL> > Tau;
	
	TPZEulerConsLaw::JacobFlux(fGamma, dim, sol, Ai);
	ComputeTau(dim, jacinv, sol, Ai, Tau);
	
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
					Ai[t](i,j).dx/*fastAccessDx*/(l)=0.;
					Tau[t](i,j).dx/*fastAccessDx*/(l)=0.;
				}
			}
#endif
	TPZVec<FADREAL> Div;
	
	//Computing the divergent
	Divergent(dsol, Ai, Div);
	TauDiv.Resize(dim);
	
	// computing Tau.Div = {Tx.Div, Ty.Div, Tz.Div}
	int k;
	for(k=0;k<dim;k++)
		Tau[k].Multiply(Div, TauDiv[k]);
}

template <int dim>
void TPZArtDiff::PrepareFastestDiff(TPZFMatrix<REAL> &jacinv,
									TPZVec<STATE> &sol,
									TPZFMatrix<STATE> &dsol,
									TPZFMatrix<REAL> &phi,
									TPZFMatrix<REAL> &dphi,
									TPZVec<TPZVec<STATE> > & TauDiv,
									TPZVec<TPZDiffMatrix<STATE> > & dTauDiv)
{
#ifdef _TFAD
	typedef TFad<dim+2, REAL> TFADREALdim;
#endif
#ifdef _FAD
	typedef Fad<REAL> TFADREALdim;
#endif
#ifdef _TINYFAD
	typedef TinyFad<dim+2, REAL> TFADREALdim;
#endif
	
	const int nstate = sol.NElements();
	const int nshape = phi.Rows();
	
	TPZVec<TFADREALdim >                FADsol(nstate);
	TPZVec<TPZDiffMatrix<TFADREALdim> > FADAi(dim);
	TPZVec<TPZDiffMatrix<REAL> >           Ai(dim);
	TPZVec<TPZDiffMatrix<TFADREALdim> > FADTau(dim);
	TPZVec<TPZDiffMatrix<STATE> >           Tau(dim);
	TPZVec<TPZVec<TFADREALdim> >        FADTauDiv(dim);
	TFADREALdim temp;
	
	int i, j, k, l;
	for(i = 0; i < nstate; i++)
	{
		FADsol[i] = sol[i];
		FADsol[i].diff(i, nstate);
	}
	
	TPZEulerConsLaw::JacobFlux(fGamma, dim, FADsol, FADAi);
	ComputeTau(dim, jacinv, FADsol, FADAi, FADTau);
	
	for( k = 0; k < dim; k++)
	{
		Tau[k].Redim(nstate, nstate);
		Ai [k].Redim(nstate, nstate);
		for(i = 0; i < nstate; i++)
			for( j = 0; j < nstate; j++)
			{
				Tau[k](i,j) = FADTau[k](i,j).val();
				Ai [k](i,j) = FADAi [k](i,j).val();
			}
	}
	
	TPZVec<STATE> Div;
	TPZDiffMatrix<STATE> dDiv;
	
	//Computing the divergent with derivatives
	Divergent(dsol, phi, dphi, FADAi, Div, &dDiv);
	
	TauDiv. Resize(dim);
	dTauDiv.Resize(dim);
	
	//Computing Tau * Div and DTau * Div
	for(k=0;k<dim;k++)
	{
		TauDiv [k].Resize(nstate);
		dTauDiv[k].Redim(nstate, nstate);
		//FADTau[k].Multiply(Div, FADTauDiv[k]);
		FADTauDiv[k].Resize(nstate);
		for(i = 0; i < nstate; i++)
		{
			temp = 0.;
			for(j = 0; j < nstate; j++)
				temp += FADTau[k](i,j) * ((REAL)Div[j]);
			FADTauDiv[k][i] = temp;
		}//
		
		// copying data using REAL
		dTauDiv[k].Redim(nstate, nstate * nshape);
		for(i = 0; i < nstate; i++)
		{
			TauDiv[k][i] = FADTauDiv[k][i].val();
			for(j = 0; j < nstate; j++)
			{
				for(l = 0; l < nshape; l++)
					dTauDiv[k](i,l * nstate + j) =
					FADTauDiv[k][i].dx/*fastAccessDx*/(j) * phi(l,0);
			}
		}
		
	}
	
	TPZDiffMatrix<STATE>  TaudDiv_k; // temporary storage
	
	for(k = 0; k < dim; k++)
	{
		Tau[k].Multiply(dDiv, TaudDiv_k);
		dTauDiv[k].Add(TaudDiv_k, 1.);
	}
	
	//Computing DTauDiv = DTau * Div + Tau * DDiv
}




//-----------------Contribute

void TPZArtDiff::ContributeApproxImplDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,  TPZFMatrix<REAL> &dphix, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, REAL weight, REAL timeStep, REAL deltaX)
{
    REAL delta = Delta(deltaX, sol);
    REAL constant = /*-*/ weight * delta * timeStep;
	
    TPZVec<TPZVec<STATE> > TauDiv;
    TPZVec<TPZDiffMatrix<STATE> > * pTaudDiv = NULL;
    TPZVec<TPZDiffMatrix<STATE> > TaudDiv;
	
    pTaudDiv = & TaudDiv;
	
    PrepareFastDiff(dim, jacinv, sol, dsol, dphix, TauDiv, pTaudDiv);
	
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

void TPZArtDiff::ContributeExplDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol,  TPZFMatrix<REAL> &dphix, TPZFMatrix<STATE> &ef, REAL weight, REAL timeStep, REAL deltaX)
{
    REAL delta = Delta(deltaX, sol);
    REAL constant = /*-*/ weight * delta * timeStep;
	
    TPZVec<TPZVec<STATE> > TauDiv;
	
    PrepareFastDiff(dim, jacinv, sol, dsol, dphix, TauDiv, NULL);
	
    int i, k, l;
    int nshape = dphix.Cols();
    int nstate = dim + 2;
	
    // ODotProduct speeded up
    for(l=0;l<nshape;l++)
		for(i=0;i<nstate;i++)
			for(k=0;k<dim;k++)
				ef(i+l*nstate,0) += dphix(k,l) * TauDiv[k][i] * constant;
}



void TPZArtDiff::ContributeImplDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<FADREAL> &sol, TPZVec<FADREAL> &dsol, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, REAL weight,  REAL timeStep, REAL deltaX)
{
    TPZVec<STATE> solReal(sol.NElements());
    int i;
    for(i = 0; i < sol.NElements(); i++)
		solReal[i] = sol[i].val();
	
    REAL delta = Delta(deltaX, solReal);
    REAL constant = /*-*/ delta * weight * timeStep;
	
    TPZVec<TPZVec<FADREAL> > TauDiv;
	
    PrepareFastDiff(dim, jacinv, sol, dsol, TauDiv);
	
    TPZVec<FADREAL> Diff;
    TPZVec<REAL> gradv(dim);
	
    int j, k, l;
    int nstate = dim + 2;
    int neq = sol[0].size();
    int nshape = neq/nstate;
	
    for(l=0;l<nshape;l++)
	{
		for(k=0;k<dim;k++)
			gradv[k] = dsol[k].dx/*fastAccessDx*/(/*k+*/l*nstate);// always retrieving this information from the first state variable...
		ODotOperator(gradv, TauDiv, Diff);
		for(i=0;i<nstate;i++)
		{
			ef(i+l*nstate,0) += constant * Diff[i].val();
			for(j=0;j<neq;j++)
				ek(i+l*nstate, j) -= constant * Diff[i].dx/*fastAccessDx*/(j);
		}
	}
}

template <int dim>
void TPZArtDiff::ContributeFastestImplDiff_dim(TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, REAL weight, REAL timeStep, REAL deltaX)
{
    REAL delta = Delta(deltaX, sol);
    REAL constant = /*-*/ weight * delta * timeStep;
    REAL buff;
	
    TPZVec<TPZVec<STATE> > TauDiv;
    TPZVec<TPZDiffMatrix<STATE> > dTauDiv;
	
    PrepareFastestDiff<dim>( jacinv, sol, dsol, phi, dphi, TauDiv, dTauDiv);
	
    int i, j, k, l;
    int nshape = dphi.Cols();
    int nstate = dim + 2;
    int neq = nstate * nshape;
	
    // ODotProduct speeded up
	
    for(l=0;l<nshape;l++)
		for(i=0;i<nstate;i++)
			for(k=0;k<dim;k++)
			{
				buff = dphi(k,l) * constant;
				ef(i+l*nstate,0) += buff * TauDiv[k][i];
				for(j=0;j<neq;j++)
					ek(i+l*nstate,j) -= buff * dTauDiv[k](i,j);
			}
}

void TPZArtDiff::ContributeFastestImplDiff(int dim, TPZFMatrix<REAL> &jacinv, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, REAL weight, REAL timeStep, REAL deltaX)
{
	switch(dim)
	{
		case(1):
			ContributeFastestImplDiff_dim<1>(jacinv, sol, dsol,
											 phi, dphi, ek, ef,
											 weight, timeStep, deltaX);
			break;
		case(2):
			ContributeFastestImplDiff_dim<2>(jacinv, sol, dsol,
											 phi, dphi, ek, ef,
											 weight, timeStep, deltaX);
			break;
		case(3):
			ContributeFastestImplDiff_dim<3>(jacinv, sol, dsol,
											 phi, dphi, ek, ef,
											 weight, timeStep, deltaX);
			break;
		default:
			PZError << "\nTPZArtDiff::ContributeFastestImplDiff unhandled dimension\n";
			exit(-1);
	}
}

template< class T >
void TPZArtDiff::Pressure(REAL gamma, int dim, T& press, TPZVec<T> &U)
{
	TPZEulerConsLaw::Pressure(gamma, dim, press, U);
}


void TPZArtDiff::Write(TPZStream &buf, int withclassid) const
{
	int tmp = static_cast<int>(fArtDiffType);
	buf.Write(&tmp,1);
	buf.Write(&fGamma, 1);
	buf.Write(&fDelta, 1);
	buf.Write(&fCFL, 1);
}

void TPZArtDiff::Read(TPZStream &buf, void *context)
{
	int ArtDiffType;
	buf.Read(&ArtDiffType, 1);
	fArtDiffType = static_cast<TPZArtDiffType>(ArtDiffType);
	buf.Read(&fGamma, 1);
	buf.Read(&fDelta, 1);
	buf.Read(&fCFL, 1);
}

int TPZArtDiff::ClassId() const{
    return Hash("TPZArtDiff");
}

