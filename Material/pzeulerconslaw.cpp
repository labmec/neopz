/**
 * \file
 * @brief Contains implementations of the TPZEulerConsLaw methods.
 */
//$Id: pzeulerconslaw.cpp,v 1.54 2011-05-13 20:46:51 phil Exp $

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
#include "pzsave.h"
#include "pzerror.h"

#include <cmath>
using namespace std;

#include "pzlog.h"

#ifdef LOG4CXX

LoggerPtr fluxroe(Logger::getLogger("pz.fluxroe"));
LoggerPtr fluxappr(Logger::getLogger("pz.fluxappr"));

#endif



TPZEulerConsLaw::~TPZEulerConsLaw(){
	
}

TPZEulerConsLaw::TPZEulerConsLaw(int nummat,REAL timeStep,
								   REAL gamma,int dim,
								   TPZArtDiffType artdiff) :
TPZConservationLaw(nummat,timeStep,dim),
fArtDiff(artdiff, gamma),
fDiff(Explicit_TD),
fConvVol(Explicit_TD),
fConvFace(Explicit_TD)
{
	fGamma = gamma;
}

TPZEulerConsLaw::TPZEulerConsLaw() :
TPZConservationLaw(-1, 0, 3),
fArtDiff(LeastSquares_AD, 1.4),
fDiff(Explicit_TD),
fConvVol(Explicit_TD),
fConvFace(Explicit_TD)
{
	fGamma = 1.4;
}

void TPZEulerConsLaw::SetTimeDiscr(TPZTimeDiscr Diff, TPZTimeDiscr ConvVol, TPZTimeDiscr ConvFace)
{
	fDiff = Diff;
	fConvVol = ConvVol;
	fConvFace = ConvFace;
}

REAL TPZEulerConsLaw::OptimalCFL(int degree)
{
	return 1./((2.0*(REAL)degree) + 1.0);
}

REAL TPZEulerConsLaw::SetTimeStep(REAL maxveloc,REAL deltax,int degree)
{
	REAL CFL = fCFL;
	// Notice that fCFL remains 0, so that optimal CFL will
	// be computed unless CFL is redefined.
	if(CFL < 0.0) CFL = OptimalCFL(degree);
	
	REAL deltaT = CFL*deltax/maxveloc;
	//cout << "TPZCompMesh::Delta Time : " << deltaT << endl;
	TPZConservationLaw::SetTimeStep(deltaT);
	
	return deltaT;
}

int TPZEulerConsLaw::NStateVariables(int dim) {
	return (2 + dim);//U = (rho, rhou, rhov, rhow, rhoe)
}

int TPZEulerConsLaw::NStateVariables() {
	return NStateVariables(Dimension());//U = (rho, rhou, rhov, rhow, rhoe)
}

REAL TPZEulerConsLaw::Pressure(TPZVec<REAL> &U)
{
	REAL press;
	TPZEulerConsLaw::Pressure(fGamma, fDim, press, U);
	return press;
}

void TPZEulerConsLaw::Print(std::ostream &out) {
	
	TPZMaterial::Print(out);
	
	TPZConservationLaw::Print(out);
	out << "Artificial Diffusion: " <<
	fArtDiff.DiffusionName().Str() << std::endl;
	out << "Number of State Variables: " << NStateVariables() << std::endl;
	out << "Number of Fluxes: " << NFluxes() << std::endl;
	
	switch(fDiff)
	{
		case(Explicit_TD):
			out << "Explicit Diffusive term\n";
			break;
		case(ApproxImplicit_TD):
			out << "ApproxImplicit Diffusive term\n";
			break;
		case(Implicit_TD):
			out << "Implicit Diffusive term\n";
			break;
		default:
			out << "No Diffusive term\n";
	}
	
	switch(fConvVol)
	{
		case(Explicit_TD):
			out << "Explicit Volume Convective term\n";
			break;
		case(Implicit_TD):
			out << "Implicit Volume Convective term\n";
			break;
		default:
			out << "No Volume Convective term\n";
	}
	
	switch(fConvFace)
	{
		case(Explicit_TD):
			out << "Explicit Face Convective term\n";
			break;
		case(Implicit_TD):
			out << "Implicit Face Convective term\n";
			break;
		case(ApproxImplicit_TD):
			out << "Approximate Implicit Face Convective term\n";
			break;
		default:
			std::cout << "No Face Convective term\n";
	}
	
	
}

int TPZEulerConsLaw::VariableIndex(const std::string &name) {
	if( !strcmp(name.c_str(),"density")  )     return 1;//rho
	if( !strcmp(name.c_str(),"velocity") )     return 2;//(u,v,w)
	if( !strcmp(name.c_str(),"energy")   )     return 3;//E
	if( !strcmp(name.c_str(),"pressure") )     return 4;//p
	if( !strcmp(name.c_str(),"solution") )     return 5;//(ro,u,v,w,E)
	if( !strcmp(name.c_str(),"normvelocity") ) return 6;//sqrt(u+v+w)
	if( !strcmp(name.c_str(),"Mach") )         return 7;//sqrt(u+v+w)/c
	std::cout << "TPZEulerConsLaw::VariableIndex not defined\n";
	return TPZMaterial::VariableIndex(name);
}

int TPZEulerConsLaw::NSolutionVariables(int var){
	
	if(var == 1 || var == 3 || var == 4 || var == 6 || var == 7) return 1;
	if(var == 2) return Dimension();
	if(var == 5) return NStateVariables();
	
	std::cout << "TPZEulerConsLaw::NSolutionVariables not defined\n";
	return 0;
}

REAL TPZEulerConsLaw::DeltaX(REAL detJac)
{
	return REAL(2.) * pow(fabs(detJac),(REAL)(1./fDim));
}

REAL TPZEulerConsLaw::Det(TPZFMatrix<REAL> & Mat)
{
	switch(Mat.Rows())
	{
		case 1:
			return Mat(0,0);
			//	break;
		case 2:
			return Mat(0,0) * Mat(1,1) -
	        Mat(1,0) * Mat(0,1);
			//        break;
		case 3:
			return Mat(0,0) * Mat(1,1) * Mat(2,2) +
			Mat(1,0) * Mat(2,1) * Mat(0,2) +
			Mat(2,0) * Mat(0,1) * Mat(1,2) -
			Mat(0,2) * Mat(1,1) * Mat(2,0) -
			Mat(1,2) * Mat(2,1) * Mat(0,0) -
			Mat(2,2) * Mat(0,1) * Mat(1,0);
			//        break;
		default:
			PZError << "TPZEulerConsLaw::Det error: unhandled matrix size: " <<
			Mat.Rows() << std::endl;
	}
	
	return 0.;
}

int TPZEulerConsLaw::NFluxes()
{
	return Dimension();
}

//-----------------Solutions

void TPZEulerConsLaw::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout){
	
	if(fabs(Sol[0]) < 1.e-10) {
		PZError << "\nTPZEulerConsLaw::Solution: Density almost null\n"
		<< "Density = " << Sol[0] << std::endl;
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
		TPZEulerConsLaw::cSpeed(Sol, fGamma, cspeed);
		TPZEulerConsLaw::uRes(Sol, us);
		Solout[0] = us / cspeed;
		return;
	} else {
		//cout << "TPZEulerConsLaw::Solution variable in the base class\n";
		TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
	}
}


void TPZEulerConsLaw::SetDelta(REAL delta)
{
	fArtDiff.SetDelta(delta);
}

//------------------Differentiable variables setup

#ifdef _AUTODIFF

void TPZEulerConsLaw::PrepareFAD(
								  TPZVec<REAL> & sol, TPZFMatrix<REAL> & dsol,
								  TPZFMatrix<REAL> & phi, TPZFMatrix<REAL> & dphi,
								  TPZVec<FADREAL> & FADsol,
								  TPZVec<FADREAL> & FADdsol)
{
	int nState = NStateVariables();
	int nShape = phi.Rows();
	int i_state, i_shape, k;
	int nDer = nState * nShape;
	
	// initializing the differentiable variables
	FADREAL defaultFAD(nDer, REAL(0.), REAL(0.));
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

void TPZEulerConsLaw::PrepareInterfaceFAD(
										   TPZVec<REAL> &solL,TPZVec<REAL> &solR,
										   TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
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
	FADREAL defaultFAD(nDerL + nDerR, REAL(0.), REAL(0.));
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
void TPZEulerConsLaw::PrepareFastestInterfaceFAD(
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

void TPZEulerConsLaw::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	// initial guesses for sol
	// fForcingFunction is null at iterations > 0
	if(fForcingFunction)
	{
		TPZVec<REAL> res;
		int i, nState = NStateVariables();
		fForcingFunction->Execute(data.x, res);
		for(i = 0; i < nState; i++)
			data.sol[0][i] = res[i];
	}
	
	if(fContributionTime == Last_CT)
	{
		ContributeLast(data.x, data.jacinv,
					   data.sol[0], data.dsol[0],
					   weight,
					   data.phi, data.dphix,
					   ef);
		return;
	}
	
	if(fContributionTime == Advanced_CT)
	{
		ContributeAdv(data.x, data.jacinv,
					  data.sol[0], data.dsol[0],
					  weight,
					  data.phi, data.dphix,
					  ek, ef);
		return;
	}
	
	PZError << "TPZEulerConsLaw::Contribute> Unhandled Contribution Time";
}

void TPZEulerConsLaw::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef)
{
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	// initial guesses for sol
	// fForcingFunction is null at iterations > 0
	if(fForcingFunction)
	{
		TPZVec<REAL> res;
		int i, nState = NStateVariables();
		fForcingFunction->Execute(data.x, res);
		for(i = 0; i < nState; i++)
			data.sol[i] = res[i];
	}
	
	if(fContributionTime == Last_CT)
	{
		ContributeLast(data.x, data.jacinv,
					   data.sol[0], data.dsol[0],
					   weight,
					   data.phi, data.dphix,
					   ef);
		return;
	}
	
	if(fContributionTime == Advanced_CT)
	{
		ContributeAdv(data.x, data.jacinv,
					  data.sol[0], data.dsol[0],
					  weight,
					  data.phi, data.dphix,
					  ef);
		return;
	}
	
	PZError << "TPZEulerConsLaw::Contribute> Unhandled Contribution Time";
}

void TPZEulerConsLaw::ContributeLast(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
									  TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
									  REAL weight,
									  TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
									  TPZFMatrix<REAL> &ef)
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


void TPZEulerConsLaw::ContributeAdv(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
									 TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol,
									 REAL weight,
									 TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
									 TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
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
		std::cout << "TPZEulerConsLaw::Contribute> Implicit diffusive contribution: _AUTODIFF directive not configured -> Using an approximation to the tgMatrix";
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

void TPZEulerConsLaw::ContributeAdv(TPZVec<REAL> &x,TPZFMatrix<REAL> &jacinv,
									 TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol,
									 REAL weight,
									 TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
									 TPZFMatrix<REAL> &ef)
{
	// contributing the implicit parcell of the residual to the
	// rhs.
	
	if(fResidualType == Residual_RT)
	{
		ContributeExplT1(x,sol,dsol,weight, phi,dphi,ef);
	}
	// Flux_RT -> contribution only to the residual vector
	if (fDiff == Implicit_TD || fDiff == ApproxImplicit_TD)
		ContributeExplDiff(x, jacinv, sol,dsol,weight, phi, dphi, ef);
	if (fConvVol == Implicit_TD)
		ContributeExplConvVol(x, sol, weight, phi, dphi, ef);
	
}



// void TPZEulerConsLaw::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
// {
//
//    // contributing face-based quantities
//    if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
//       {
//       // if face contribution is implicit,
//       // then the FAD classes must be initialized
//       #ifdef _AUTODIFF
//          #ifdef FASTEST_IMPLICIT
//             ContributeFastestImplConvFace(fDim, data.x, dataleft.sol, dataright.sol,
//                          weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
// 	 #else
//          TPZVec<FADREAL> FADsolL, FADsolR;
//          PrepareInterfaceFAD(dataleft.sol, dataright.sol, dataleft.phi, dataright.phi, FADsolL, FADsolR);
//          ContributeImplConvFace(data.x,FADsolL,FADsolR, weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
// #ifdef LOG4CXX
//           if(fluxroe->isDebugEnabled()) {
//             std::stringstream sout;
//             ek.Print("computed tangent matrix",sout);
//             ef.Print("computed rhs",sout);
//             LOGPZ_DEBUG(fluxroe,sout.str().c_str());
//           }
// #endif
// 	 #endif
//       #else
//       // forcing explicit contribution and issueing an warning
//          cout << "TPZEulerConsLaw::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
//          ContributeExplConvFace(data.x,dataleft.sol,dataright.sol,weight,data.normal,dataleft.phi,dataright.phi,ef);
//       #endif
//       }
//       else if (fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
//       {
// #ifdef _AUTODIFF
//         REAL facesize = 0.;
//         TPZVec<FADREAL> FADsolL, FADsolR;
//         PrepareInterfaceFAD(dataleft.sol, dataright.sol, dataleft.phi, dataright.phi, FADsolL, FADsolR);
//         ContributeApproxImplConvFace(data.x,facesize,FADsolL,FADsolR, weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
// #endif
// #ifdef LOG4CXX
//         if(fluxroe->isDebugEnabled()){
//           std::stringstream sout;
//           ek.Print("computed tangent matrix",sout);
//           ef.Print("computed rhs",sout);
//           LOGPZ_DEBUG(fluxroe,sout.str().c_str());
//         }
// #endif
//       }
//
//
//    if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
//    {
//       ContributeExplConvFace(data.x,dataleft.sol,dataright.sol,weight,data.normal,dataleft.phi,dataright.phi,ef);
//    }
// }

void TPZEulerConsLaw::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                          REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
#ifdef LOG4CXX
    if(fluxroe->isDebugEnabled()){
		std::stringstream sout;
		sout << "solL " << dataleft.sol << std::endl << "solR " << dataright.sol << std::endl;
		LOGPZ_DEBUG(fluxroe,sout.str().c_str());
		LOGPZ_DEBUG(fluxappr,sout.str().c_str());
    }
#endif
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	// contributing face-based quantities
	if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
	{
		// if face contribution is implicit,
		// then the FAD classes must be initialized
#ifdef _AUTODIFF
#ifdef FASTEST_IMPLICIT
		ContributeFastestImplConvFace(fDim, data.x, dataleft.sol[0], dataright.sol[0],
									  weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
#else
		TPZVec<FADREAL> FADsolL, FADsolR;
		PrepareInterfaceFAD(dataleft.sol[0], dataright.sol[0], dataleft.phi, dataright.phi, FADsolL, FADsolR);
		ContributeImplConvFace(data.x,FADsolL,FADsolR, weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
#ifdef LOG4CXX
        if(fluxroe->isDebugEnabled()){
			std::stringstream sout;
			ek.Print("computed tangent matrix",sout);
			ef.Print("computed rhs",sout);
			LOGPZ_DEBUG(fluxroe,sout.str().c_str());
        }
#endif
#endif
#else
		// forcing explicit contribution and issueing an warning
		std::cout << "TPZEulerConsLaw::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
		//         ContributeApproxImplConvFace(x,data.HSize,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef);
		ContributeExplConvFace(data.x,dataleft.sol[0],dataright.sol[0],weight,data.normal,dataleft.phi,dataright.phi,ef);
#endif
	}
	else if(fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
	{
#ifdef _AUTODIFF
		TPZVec<FADREAL> FADsolL, FADsolR;
		PrepareInterfaceFAD(dataleft.sol[0], dataright.sol[0], dataleft.phi, dataright.phi, FADsolL, FADsolR);
		ContributeApproxImplConvFace(data.x,data.HSize,FADsolL,FADsolR, weight, data.normal, dataleft.phi, dataright.phi, ek, ef);
#endif
#ifdef LOG4CXX
		if(fluxappr->isDebugEnabled()){
			std::stringstream sout;
			ek.Print("computed tangent matrix",sout);
			ef.Print("computed rhs",sout);
			LOGPZ_DEBUG(fluxappr,sout.str().c_str());
		}
#endif
	}
	
	if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
	{
		ContributeExplConvFace(data.x,dataleft.sol[0],dataright.sol[0],weight,data.normal,dataleft.phi,dataright.phi,ef);
	}
}


void TPZEulerConsLaw::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                          REAL weight, TPZFMatrix<REAL> &ef)
{
	
	// contributing face-based quantities
	if ((fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
		||
		(fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
		||
		(fConvFace == Explicit_TD && fContributionTime == Last_CT))
	{
		ContributeExplConvFace(data.x,dataleft.sol[0],dataright.sol[0],weight,data.normal,dataleft.phi,dataright.phi,ef);
	}
	
}

void TPZEulerConsLaw::ContributeBC(TPZMaterialData &data,
                                    REAL weight,
                                    TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
                                    TPZBndCond &bc)
{
	int phr = data.phi.Rows();
	short in,jn,i,j;
	int nstate = NStateVariables();
	REAL v2[5];//m�imo nstate
	for(i=0;i<nstate;i++) v2[i] = bc.Val2()(i,0);
	
	switch (bc.Type()) {
		case 0 :// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i,0) += gBigNumber * weight * v2[i] * data.phi(in,0);
				for (jn = 0 ; jn < phr; jn++) {
					for(i = 0 ; i < nstate; i++)
						ek(in*nstate+i,jn*nstate+i) -= gBigNumber * weight * data.phi(in,0) * data.phi(jn,0);
				}
			}
			break;
		case 1 :// Neumann condition
			for(in = 0 ; in < data.phi.Rows(); in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i,0) += v2[i] * data.phi(in,0) * weight;
			}
			break;
		case 2 :// condi�o mista
			for(in = 0 ; in < data.phi.Rows(); in++) {
				for(i = 0 ; i < nstate; i++)
					ef(in*nstate+i, 0) += weight * v2[i] * data.phi(in, 0);
				for (jn = 0 ; jn < data.phi.Rows(); jn++) {
					for(i = 0 ; i < nstate; i++) for(j = 0 ; j < nstate; j++)
						ek(in*nstate+i,jn*nstate+j) -= weight * bc.Val1()(i,j) * data.phi(in,0) * data.phi(jn,0);
				}
			}
	}
}

void TPZEulerConsLaw::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
											 TPZBndCond &bc)
{
	int nstate = NStateVariables();
	TPZVec<REAL> solR(nstate,0.);
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFMatrix<REAL> dsolR(dataleft.dsol[0].Rows(), dataleft.dsol[0].Cols(),0.);
	TPZFMatrix<REAL> phiR(0,0), dphiR(0,0);
	int entropyFix;
	
	// contributing face-based quantities
	if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
	{
		// if face contribution is implicit,
		// then the FAD classes must be initialized
#ifdef _AUTODIFF
		
		//          if(bc.Type() ==5 && fDim == 2)
		//           {
		//             int entropyFix2;
		//             TPZManVector<REAL,5 > flux(nstate,0.);
		//             ComputeGhostState(solL, solR, normal, bc, entropyFix2);
		//             Roe_Flux<REAL>(solL, solR, normal, fGamma, flux, entropyFix2);
		//             REAL norflux = flux[1]*normal[1]-flux[2]*normal[0];
		//             REAL err = fabs(flux[0])+fabs(norflux)+fabs(flux[3]);
		//             if(err > 1.e-5)
		//             {
		//               cout << "fluxo de parede errado 1 err " << err << endl;
		//             }
		//           }
		
#ifdef FASTEST_IMPLICIT
		ContributeFastestBCInterface(fDim, data.x, dataleft.sol[0], dataleft.dsol[0],
									 weight, data.normal, dataleft.phi, data.dphixl, data.axesleft, ek, ef, bc);
#else
		TPZVec<FADREAL> FADsolL, FADsolR;
		PrepareInterfaceFAD(dataleft.sol[0], solR, dataleft.phi, phiR, FADsolL, FADsolR);
		ComputeGhostState(FADsolL, FADsolR, data.normal, bc, entropyFix);
		ContributeImplConvFace(data.x,FADsolL,FADsolR, weight, data.normal, dataleft.phi, phiR, ek, ef, entropyFix);
#ifdef LOG4CXX
        if(fluxroe->isDebugEnabled()){
			std::stringstream sout;
			ek.Print("computed tangent matrix",sout);
			ef.Print("computed rhs",sout);
			LOGPZ_DEBUG(fluxroe,sout.str().c_str());
        }
#endif
#endif
#else
		// forcint explicit contribution and issueing an warning
		std::cout << "TPZEulerConsLaw::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
		//         ComputeGhostState(FADsolL, FADsolR, normal, bc, entropyFix);
		//         ContributeApproxImplConvFace(x,solL,solR,weight,normal,phiL,phiR,ek,ef,entropyFix);
#endif
	}
	else if (fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
	{
#ifdef _AUTODIFF
		TPZVec<FADREAL> FADsolL, FADsolR;
		PrepareInterfaceFAD(dataleft.sol[0], solR, dataleft.phi, phiR, FADsolL, FADsolR);
		ComputeGhostState(FADsolL, FADsolR, data.normal, bc, entropyFix);
		ContributeApproxImplConvFace(data.x,data.HSize,FADsolL,FADsolR, weight, data.normal, dataleft.phi, phiR, ek, ef, entropyFix);
#ifdef LOG4CXX
		if(fluxappr->isDebugEnabled()){
			std::stringstream sout;
			ek.Print("computed tangent matrix",sout);
			ef.Print("computed rhs",sout);
			LOGPZ_DEBUG(fluxappr,sout.str().c_str());
		}
#endif
#else
		//    ComputeGhostState(solL, solR, normal, bc, entropyFix);
		//    ContributeApproxImplConvFace(x,data.HSize,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef, entropyFix);
#endif
	}
	
	if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
	{
		ComputeGhostState(dataleft.sol[0], solR, data.normal, bc, entropyFix);
		ContributeExplConvFace(data.x,dataleft.sol[0],solR,weight,data.normal,dataleft.phi,phiR,ef, entropyFix);
	}
}

// void TPZEulerConsLaw::ContributeBCInterface(TPZMaterialData &data,
//                                               REAL weight,
//                                               TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
//                                               TPZBndCond &bc)
// {
//    int nstate = NStateVariables();
//    TPZVec<REAL> solR(nstate,0.);
//    TPZFMatrix<REAL> dsolR(dataleft.dsol.Rows(), dataleft.dsol.Cols(),0.);
//    TPZFMatrix<REAL> phiR(0,0), dphiR(0,0);
//    int entropyFix;
//
//    // contributing face-based quantities
//    if (fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
//       {
//       // if face contribution is implicit,
//       // then the FAD classes must be initialized
//       #ifdef _AUTODIFF
// #ifdef DEBUG
// /*         if(bc.Type() ==5 && fDim == 2)
//          {
// 	    int entropyFix2;
//             TPZManVector<REAL,5 > flux(nstate,0.);
//             ComputeGhostState(solL, solR, normal, bc, entropyFix2);
//             Roe_Flux<REAL>(solL, solR, normal, fGamma, flux, entropyFix2);
//             REAL norflux = flux[1]*normal[1]-flux[2]*normal[0];
//             REAL err = fabs(flux[0])+fabs(norflux)+fabs(flux[3]);
//             if(err > 1.e-5)
//             {
//               cout << "fluxo de parede errado 1 err " << err << endl;
//             }
//          }*/
// #endif
//
//          #ifdef FASTEST_IMPLICIT
//             ContributeFastestBCInterface(fDim, x, solL, dsolL,
// 	                      weight, normal, phiL, phiR, ek, ef, bc);
// 	 #else
//          TPZVec<FADREAL> FADsolL, FADsolR;
//          PrepareInterfaceFAD(solL, solR, phiL, phiR, FADsolL, FADsolR);
// 	 ComputeGhostState(FADsolL, FADsolR, normal, bc, entropyFix);
//          ContributeImplConvFace(x,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef, entropyFix);
// #ifdef LOG4CXX
//         if(fluxroe->isDebugEnabled()){
//           std::stringstream sout;
//           ek.Print("computed tangent matrix",sout);
//           ef.Print("computed rhs",sout);
//           LOGPZ_DEBUG(fluxroe,sout.str().c_str());
//         }
// #endif
// #endif
//       #else
//       // forcint explicit contribution and issueing an warning
//          cout << "TPZEulerConsLaw::ContributeInterface> Implicit face convective contribution: _AUTODIFF directive not configured";
//          ContributeExplConvFace(data.x,dataleft.sol,solR,weight,data.normal,dataleft.phi,phiR,ef);
//       #endif
//       }
//       else if (fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
//       {
// #ifdef _AUTODIFF
//         TPZVec<FADREAL> FADsolL, FADsolR;
//         PrepareInterfaceFAD(solL, solR, phiL, phiR, FADsolL, FADsolR);
//         ComputeGhostState(FADsolL, FADsolR, normal, bc, entropyFix);
//         ContributeImplConvFace(x,FADsolL,FADsolR, weight, normal, phiL, phiR, ek, ef, entropyFix);
// #endif
// #ifdef LOG4CXX
//         if(fluxappr->isDebugEnabled()){
//           std::stringstream sout;
//           ek.Print("computed tangent matrix",sout);
//           ef.Print("computed rhs",sout);
//           LOGPZ_DEBUG(fluxappr,sout.str().c_str());
//         }
// #endif
//       }
//
//
//    if(fConvFace == Explicit_TD && fContributionTime == Last_CT)
//    {
//          ComputeGhostState(dataleft.sol, solR, data.normal, bc, entropyFix);
//          ContributeExplConvFace(data.x,dataleft.sol,solR,weight,data.normal,dataleft.phi,phiR,ef, entropyFix);
//    }
// }


void TPZEulerConsLaw::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZFMatrix<REAL> &ef,
											 TPZBndCond &bc)
{
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	int nstate = NStateVariables();
	TPZVec<REAL> solR(nstate,0.);
	TPZFMatrix<REAL> dsolR(dataleft.dsol[0].Rows(), dataleft.dsol[0].Cols(),0.);
	TPZFMatrix<REAL> phiR(0,0), dphiR(0,0);
	int entropyFix;
	
	if((fConvFace == Implicit_TD && fContributionTime == Advanced_CT)
	   ||
	   (fConvFace == ApproxImplicit_TD && fContributionTime == Advanced_CT)
	   ||
	   (fConvFace == Explicit_TD && fContributionTime == Last_CT))
	{
		ComputeGhostState(dataleft.sol[0], solR, data.normal, bc, entropyFix);
		
		if(fDim == 2)
		{// flux tests
            TPZManVector<REAL,3> normal2(2,0.);
            TPZManVector<REAL,5> flux2(nstate,0.);
            normal2[0] = -data.normal[0];
            normal2[1] = -data.normal[1];
            TPZManVector<REAL,5 > flux(nstate,0.);
            Roe_Flux<REAL>(dataleft.sol[0], solR, data.normal, fGamma, flux);
            Roe_Flux<REAL>(solR, dataleft.sol[0], normal2, fGamma, flux2);
            REAL fluxs = fabs(flux[0]+flux2[0])+fabs(flux[1]+flux2[1])+fabs(flux[2]+flux2[2])+fabs(flux[3]+flux2[3]);
            if(fluxs > 1.e-10)
            {
				std::cout << "Fluxo nao simetrico fluxs = " << fluxs << std::endl;
            }
			
            if(bc.Type() ==5)
            {
				REAL norflux = flux[1]*data.normal[1]-flux[2]*data.normal[0];
				REAL err = fabs(flux[0])+fabs(norflux)+fabs(flux[3]);
				if(err > 1.e-5)
				{
					std::cout << "fluxo de parede errado 2 err " << err << std::endl;
					Roe_Flux<REAL>(dataleft.sol[0],solR,data.normal,fGamma,flux);
				} else
				{
					Roe_Flux<REAL>(dataleft.sol[0],solR,data.normal,fGamma,flux);
				}
            }
		} // end of tests
		
		ContributeExplConvFace(data.x,dataleft.sol[0],solR,weight,data.normal,
							   dataleft.phi,phiR,ef,entropyFix);
	}
}

#ifdef _AUTODIFF

void TPZEulerConsLaw::ContributeFastestBCInterface(int dim,
													TPZVec<REAL> &x,
													TPZVec<REAL> &solL, TPZFMatrix<REAL> &dsolL,
													REAL weight, TPZVec<REAL> &normal,
													TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &dphiL, TPZFMatrix<REAL> &axesleft,
													TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
{
	
	switch(dim)
	{
		case(1):
			ContributeFastestBCInterface_dim<1>(x, solL, dsolL,
												weight, normal,
												phiL, dphiL,axesleft,
												ek, ef, bc);
			break;
		case(2):
			ContributeFastestBCInterface_dim<2>(x, solL, dsolL,
												weight, normal,
												phiL, dphiL,axesleft,
												ek, ef, bc);
			break;
		case(3):
			ContributeFastestBCInterface_dim<3>(x, solL, dsolL,
												weight, normal,
												phiL, dphiL,axesleft,
												ek, ef, bc);
			break;
		default:
			PZError << "\nTPZEulerConsLaw::ContributeFastestBCInterface unhandled dimension\n";
			exit(-1);
	}
};


template <int dim>
void TPZEulerConsLaw::ContributeFastestBCInterface_dim(TPZVec<REAL> &x,
														TPZVec<REAL> &solL, TPZFMatrix<REAL> &dsolL,
														REAL weight, TPZVec<REAL> &normal,
														TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &dphiL,TPZFMatrix<REAL> &axesleft,
														TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc)
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
	TPZFMatrix<REAL> phiR(0,0);
	
	// initial guesses for left and right sol
	// fForcingFunction is null at iterations > 0
	if(fForcingFunction)
	{
		TPZVec<REAL> res;
		fForcingFunction->Execute(x, res);
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
void TPZEulerConsLaw:: ComputeGhostState(TPZVec<T> &solL, TPZVec<T> &solR, TPZVec<REAL> &normal, TPZBndCond &bc, int & entropyFix)
{
	entropyFix = 1;
	
	int nstate = NStateVariables();
	T vpn=REAL(0.);
	T us, un, c;
	REAL Mach, temp;
	int i;
	//Riemann Invariants
	T w1, w2, w5, uninf, usinf,
    cghost, usghost, unghost, p;
	REAL cinf;
	
	switch (bc.Type()){
		case 3://Dirichlet: nada a fazer a CC �a correta
			for(i=0;i<nstate; i++) solR[i] = bc.Val2().operator()(i,0);
			break;
		case 4://recuperar valor da solu� MEF esquerda: saida livre
			for(i=0;i<nstate; i++) solR[i] = solL[i];
			break;
		case 5://condi� de parede
			for(i=1;i<nstate-1;i++) vpn += solL[i]*T(normal[i-1]);//v.n
			for(i=1;i<nstate-1;i++) solR[i] = solL[i] - T(2.0*normal[i-1])*vpn;
			solR[0] = solL[0];
			solR[nstate-1] = solL[nstate-1];
			entropyFix = 0;
			break;
		case 6://n� refletivas (campo distante)
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
			
			if(un < REAL(0.))// Inflow
			{
				//if(normal[0]>0.) cout << "\ndirection error\n ";
				
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
					w2 = cinf * cinf / T(fGamma * pow(bc.Val2()(0,0), (fGamma - 1.)));
					// w5 computed based on flow state
					w5 = un + T(2.) * c / T(fGamma - 1.);
					
					// computing ghost values
					cghost = (w5 - w1) * T((fGamma - 1.)/4.);
					solR[0] = pow(cghost * cghost / (T(fGamma) * w2), (1./(fGamma - 1.)));
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
					/*
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
					 usghost * usghost / T(2.));*/
				}
			}else
			{ // Outflow
				//if(normal[0]<0.) cout << "\ndirection error 2\n ";
				if(us>c)
				{ // supersonic: no BC at all are required
					for(i = 0; i < nstate; i++)
						solR[i] = solL[i];
				}else
				{ // subsonic outlet
					// only the condition w1 referring to the first Riemann invariant
					// is imposed. As a rule, the imposition of pressure is applied
					// instead.
					
					solR[0] = solL[0] * pow(bc.Val2()(nstate-1,0)/p, (1./fGamma));
					
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
			
			if(un < REAL(0.))// Inflow
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
					w2 = cinf * cinf / T(fGamma * pow(bc.Val2()(0,0), (fGamma - 1.)));   // Unbelived - only is defined for pow(float, float) and pow(long double, long double) !!! Jorge
					// w5 computed based on flow state
					w5 = un + T(2.) * c / T(fGamma - 1.);
					
					// computing ghost values
					cghost = (w5 - w1) * T((fGamma - 1.)/4.);
					solR[0] = pow(cghost * cghost / (T(fGamma) * w2),(1./(fGamma - 1.)));
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
			solR[0] = pow(cghost * cghost * pow(solL[0],fGamma)/ (p * fGamma), (1./(fGamma - 1.)));
			
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
		case 12: //Symmetry
			for(i=1;i<nstate-1;i++) vpn += solL[i]*T(normal[i-1]);//v.n
			for(i=1;i<nstate-1;i++) solR[i] = solL[i] - T(2.0*normal[i-1])*vpn;
			solR[0] = solL[0];
			solR[nstate-1] = solL[nstate-1];
			entropyFix = 1;
			break;
		default:
			for(i=0;i<nstate;i++) solR[i] = 0.;
	}
}

//------------------internal contributions

void TPZEulerConsLaw::ContributeApproxImplDiff(TPZVec<REAL> &x,
												TPZFMatrix<REAL> &jacinv,
												TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
												REAL weight,
												TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
												TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef)
{
	// computing the determinant of the jacobian
	REAL deltaX = DeltaX(1./Det(jacinv));
	
	fArtDiff.ContributeApproxImplDiff(fDim, jacinv, sol, dsol, dphi,
									  ek, ef, weight,
#ifdef DIVTIMESTEP
									  1.,
#else
									  TimeStep(),
#endif
									  deltaX
									  );
}

void TPZEulerConsLaw::ContributeExplDiff(TPZVec<REAL> &x,
										  TPZFMatrix<REAL> &jacinv,
										  TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
										  REAL weight,
										  TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
										  TPZFMatrix<REAL> &ef)
{
	// computing the determinant of the jacobian
	REAL deltaX = DeltaX(1./Det(jacinv));
	
	fArtDiff.ContributeExplDiff(fDim, jacinv, sol, dsol, dphi,
								ef, weight,
#ifdef DIVTIMESTEP
								1.,
#else
								TimeStep(),
#endif
								deltaX
								);
}


#ifdef _AUTODIFF
void TPZEulerConsLaw::ContributeImplDiff(TPZVec<REAL> &x,
										  TPZFMatrix<REAL> &jacinv,
										  TPZVec<FADREAL> &sol,TPZVec<FADREAL> &dsol,
										  REAL weight,
										  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef)
{
	// computing the determinant of the jacobian
	REAL deltaX = DeltaX(1./Det(jacinv));
	
	fArtDiff.ContributeImplDiff(fDim, jacinv, sol, dsol,
								ek, ef, weight,
#ifdef DIVTIMESTEP
								1.,
#else
								TimeStep(),
#endif
								deltaX
								);
}


void TPZEulerConsLaw::ContributeFastestImplDiff(int dim, TPZVec<REAL> &x, TPZFMatrix<REAL> &jacinv, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	// computing the determinant of the jacobian
	REAL deltaX = DeltaX(1./Det(jacinv));
	
	fArtDiff.ContributeFastestImplDiff(dim, jacinv, sol, dsol, phi, dphi,
									   ek, ef, weight,
#ifdef DIVTIMESTEP
									   1.,
#else
									   TimeStep(),
#endif
									   deltaX
									   );
}

#endif

void TPZEulerConsLaw::ContributeExplConvFace(TPZVec<REAL> &x,
											  TPZVec<REAL> &solL,TPZVec<REAL> &solR,
											  REAL weight,TPZVec<REAL> &normal,
											  TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
											  TPZFMatrix<REAL> &ef, int entropyFix)
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

void TPZEulerConsLaw::ContributeApproxImplConvFace(TPZVec<REAL> &x, REAL faceSize,
													TPZVec<FADREAL> &solL,TPZVec<FADREAL> &solR,
													REAL weight,TPZVec<REAL> &normal,
													TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
													TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix)
{
	int nState = NStateVariables();
	TPZVec<FADREAL > flux(nState,REAL(0.));
	ApproxRoe_Flux(solL, solR, normal, fGamma, flux, entropyFix);
	
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

#endif

void TPZEulerConsLaw::ContributeApproxImplConvFace(TPZVec<REAL> &x, REAL faceSize,
													TPZVec<REAL> &solL,TPZVec<REAL> &solR,
													REAL weight,TPZVec<REAL> &normal,
													TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
													TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef, int entropyFix
													)
{
	int nState = NStateVariables();
	TPZManVector< TPZManVector<REAL,5> ,3> FL(3),FR(3);
	TPZManVector< REAL,5> FN(nState,0.);
	Flux(solL, FL[0], FL[1], FL[2]);
	Flux(solR, FR[0], FR[1], FR[2]);
	TPZFNMatrix<36> DFNL(nState,nState,0.),DFNR(nState,nState,0.);
	
	TPZManVector<REAL,7 > fluxroe(nState,0.);
	Roe_Flux(solL, solR, normal, fGamma, fluxroe,entropyFix);
	TPZManVector< TPZDiffMatrix<REAL> ,3> AL(3),AR(3);
	JacobFlux(fGamma, fDim, solL, AL);
	JacobFlux(fGamma, fDim, solR, AR);
	
	int nShapeL = phiL.Rows(),
	nShapeR = phiR.Rows(),
	i_shape, i_state, i, j_state, j_shape;
	//  nDer = (nShapeL + nShapeR) * nState;
	
	
#ifdef DIVTIMESTEP
	REAL constant = weight; // weight
#else
	REAL constant = TimeStep() * weight; // deltaT * weight
#endif
	
	for(i_state=0; i_state<nState; i_state++)
	{
		FN[i_state]= -(faceSize/constant)*(solR[i_state]-solL[i_state]);
		DFNL(i_state,i_state) = (faceSize/constant);
		DFNR(i_state,i_state) = -(faceSize/constant);
		for(i=0; i<fDim; i++)
		{
			FN[i_state]+=0.5*normal[i]*(FL[i][i_state]+FR[i][i_state]);
			for(j_state=0; j_state<nState; j_state++)
			{
				DFNL(i_state,j_state) +=0.5*normal[i]*(AL[i](i_state,j_state));
				DFNR(i_state,j_state) +=0.5*normal[i]*(AR[i](i_state,j_state));
			}
		}
	}
	// Contribution referring to the left element
	for(i_shape = 0; i_shape < nShapeL; i_shape ++)
		for(i_state = 0; i_state < nState; i_state++)
		{
			int index = i_shape*nState + i_state;
			ef(index,0) +=
			fluxroe[i_state] * phiL(i_shape,0) * constant;
			for(j_shape = 0; j_shape < nShapeL; j_shape++)
			{
				
				for(j_state = 0; j_state < nState; j_state++)
				{
					int jndex = j_shape*nState + j_state;
					ek(index,jndex) -= DFNL(i_state,j_state) *phiL(i_shape,0) * phiL(j_shape,0) * constant;
				}
			}
			for(j_shape = 0; j_shape < nShapeR; j_shape++)
			{
				
				for(j_state = 0; j_state < nState; j_state++)
				{
					int jndex = (nShapeL+j_shape)*nState + j_state;
					ek(index,jndex) -= DFNR(i_state,j_state) *phiL(i_shape,0) * phiR(j_shape,0) * constant;
				}
			}
			/*    ek(index, j) -= flux[i_state].dx(j) *
			 phiL(i_shape,0) * constant;*/
		}
	
	// Contribution referring to the right element
	// REM: The contributions are negative in comparison
	// to the left elements: opposed normals
	for(i_shape = 0; i_shape < nShapeR; i_shape ++)
		for(i_state = 0; i_state < nState; i_state++)
		{
			int index = (nShapeL + i_shape)*nState + i_state;
			ef(index,0) -=
			fluxroe[i_state] * phiR(i_shape,0) * constant;
			for(j_shape = 0; j_shape < nShapeL; j_shape++)
			{
				
				for(j_state = 0; j_state < nState; j_state++)
				{
					int jndex = j_shape*nState + j_state;
					ek(index,jndex) += DFNL(i_state,j_state) *phiR(i_shape,0) * phiL(j_shape,0) * constant;
				}
			}
			for(j_shape = 0; j_shape < nShapeR; j_shape++)
			{
				
				for(j_state = 0; j_state < nState; j_state++)
				{
					int jndex = (nShapeL+j_shape)*nState + j_state;
					ek(index,jndex) += DFNR(i_state,j_state) *phiR(i_shape,0) * phiR(j_shape,0) * constant;
				}
			}
			/*    ek(index, j) += flux[i_state].dx(j) *
			 phiR(i_shape,0) * constant;*/
		}
}

#ifdef _AUTODIFF

void TPZEulerConsLaw::ContributeImplConvFace(TPZVec<REAL> &x,
											  TPZVec<FADREAL> &solL,TPZVec<FADREAL> &solR,
											  REAL weight,TPZVec<REAL> &normal,
											  TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
											  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
											  int entropyFix)
{
	int nState = NStateVariables();
	TPZVec<FADREAL > flux(nState,REAL(0.));
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

void TPZEulerConsLaw::ContributeFastestImplConvFace(int dim,
													 TPZVec<REAL> &x,
													 TPZVec<REAL> &solL,TPZVec<REAL> &solR,
													 REAL weight,TPZVec<REAL> &normal,
													 TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
													 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
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
			PZError << "\nTPZEulerConsLaw::ContributeFastestImplConvFace unhandled dimension\n";
			exit(-1);
	}
}

template <int dim>
void TPZEulerConsLaw::ContributeFastestImplConvFace_dim(
														 TPZVec<REAL> &x,
														 TPZVec<REAL> &solL,TPZVec<REAL> &solR,
														 REAL weight,TPZVec<REAL> &normal,
														 TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
														 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
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
void TPZEulerConsLaw::ContributeFastestImplConvFace_T(TPZVec<REAL> &x,
													   TPZVec<T> &FADsolL,TPZVec<T> &FADsolR,
													   REAL weight,TPZVec<REAL> &normal,
													   TPZFMatrix<REAL> &phiL,TPZFMatrix<REAL> &phiR,
													   TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
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

void TPZEulerConsLaw::ContributeExplConvVol(TPZVec<REAL> &x,
											 TPZVec<REAL> &sol,
											 REAL weight,
											 TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi,
											 TPZFMatrix<REAL> &ef)
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


void TPZEulerConsLaw::ContributeImplConvVol(TPZVec<REAL> &x,
											 TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
											 REAL weight,
											 TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
											 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef)
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


void TPZEulerConsLaw::ContributeExplT1(TPZVec<REAL> &x,
										TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
										REAL weight,
										TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
										TPZFMatrix<REAL> &ef)
{
	if(TimeStep()==0.)
	{
		std::cout << __PRETTY_FUNCTION__ << " Zero timestep bailing out\n";
		return;
	}
	
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

void TPZEulerConsLaw::ContributeImplT1(TPZVec<REAL> &x,
										TPZVec<REAL> &sol,TPZFMatrix<REAL> &dsol,
										REAL weight,
										TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi,
										TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef)
{
	if(TimeStep()==0.)
	{
		std::cout << __PRETTY_FUNCTION__ << " Zero timestep bailing out\n";
		return;
	}
	
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

void TPZEulerConsLaw::ContributeExplT2(TPZVec<REAL> &x,
										TPZVec<REAL> &sol,
										REAL weight,
										TPZFMatrix<REAL> &phi,
										TPZFMatrix<REAL> &ef)
{
	if(TimeStep()==0.)
	{
		std::cout << __PRETTY_FUNCTION__ << " Zero timestep bailing out\n";
		return;
	}
	
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


void TPZEulerConsLaw::Write(TPZStream &buf, int withclassid)
{
	TPZSaveable::Write(buf, 1);
	TPZConservationLaw::Write(buf, 0);
	fArtDiff.Write(buf, 0);
	int tmp = static_cast < int > (fDiff);
	buf.Write(& tmp,1);
	tmp = static_cast<int>(fConvVol);
	buf.Write(& tmp,1);
	tmp = static_cast<int>(fConvFace);
	buf.Write(&tmp,1);
}

void TPZEulerConsLaw::Read(TPZStream &buf, void *context)
{
	TPZSaveable::Read(buf, context);
	TPZConservationLaw::Read(buf, context);
	fArtDiff.Read(buf, context);
	int diff;
	buf.Read(&diff,1);
	fDiff = static_cast<TPZTimeDiscr>(diff);
	buf.Read(&diff,1);
	fConvVol = static_cast<TPZTimeDiscr>(diff);
	buf.Read(&diff,1);
	fConvFace = static_cast<TPZTimeDiscr>(diff);
}

/**
 * This flux encapsulates the two and three dimensional fluxes
 * acquired from the Mouse program
 * This function is called Approx because it evaluates the derivative of the Roe Flux without
 * using automatic differentiation
 *
 * @param solL [in]
 * @param solR [in]
 * @param normal [in]
 * @param gamma [in]
 * @param flux [in]
 */
template <class T>
void TPZEulerConsLaw::ApproxRoe_Flux(TPZVec<T> &solL, TPZVec<T> &solR,
									  TPZVec<REAL> & normal, REAL gamma,
									  TPZVec<T> & flux, int entropyFix)
{
	// Normals outgoing from the BC elements into the
	// mesh elements -> all the normals are opposited to
	// the common convention -> changing the left/right
	// elements and normals.
	int nState = solL.NElements();
	if(nState == 5)
	{
		ApproxRoe_Flux<T>(solL[0], solL[1], solL[2], solL[3], solL[4],
						  solR[0], solR[1], solR[2], solR[3], solR[4],
						  normal[0], normal[1], normal[2],
						  gamma,
						  flux[0], flux[1], flux[2], flux[3], flux[4], entropyFix);
		
	}else if(nState == 4)
	{
		ApproxRoe_Flux<T>(solL[0], solL[1], solL[2], solL[3],
						  solR[0], solR[1], solR[2], solR[3],
						  normal[0], normal[1],
						  gamma,
						  flux[0], flux[1], flux[2], flux[3], entropyFix);
	}else if(nState == 3)
	{
		//using the 2D expression for 1d problem
		T auxL = REAL(0.),
		auxR = REAL(0.),
		fluxaux = REAL(0.);
		auxL = flux[0];
		ApproxRoe_Flux<T>(solL[0], solL[1], auxL, solL[2],
						  solR[0], solR[1], auxR, solR[2],
						  normal[0], 0,
						  gamma,
						  flux[0], flux[1], fluxaux, flux[2], entropyFix);
	}else
	{
		PZError << "No flux on " << nState << " state variables.\n";
	}
}

#ifdef _AUTODIFF
/**
 * Flux of Roe (MOUSE program)
 */
template <>
void TPZEulerConsLaw::ApproxRoe_Flux<FADREAL>(const FADREAL & rho_f,
											   const FADREAL & rhou_f,
											   const FADREAL & rhov_f,
											   const FADREAL & rhow_f,
											   const FADREAL & rhoE_f,
											   const FADREAL & rho_t,
											   const FADREAL & rhou_t,
											   const FADREAL & rhov_t,
											   const FADREAL & rhow_t,
											   const FADREAL & rhoE_t,
											   const REAL nx,
											   const REAL ny,
											   const REAL nz,
											   const REAL gam,
											   FADREAL & flux_rho,
											   FADREAL & flux_rhou,
											   FADREAL & flux_rhov,
											   FADREAL & flux_rhow,
											   FADREAL & flux_rhoE, int entropyFix)
{
	
	typedef FADREAL T;
	//  REAL    alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
	REAL    a1,a2,a3,a4,a5,b1,b2,b3,b4,b5;
	// REAL    ep_t, ep_f, p_t, p_f;
	// REAL    rhouv_t, rhouv_f, rhouw_t, rhouw_f, rhovw_t, rhovw_f;
	// REAL    lambda_f, lambda_t;
	T    delta_rho, delta_rhou, delta_rhov, delta_rhow, delta_rhoE;
	REAL    hnx, hny, hnz;
	REAL    tempo11, usc;
	
	flux_rho = 0;
	flux_rhou = 0;
	flux_rhov = 0;
	flux_rhow = 0;
	flux_rhoE = 0;
	
	REAL gam1 = gam - 1.0;
	REAL    irho_f = REAL(1.0)/rho_f.val();
	REAL    irho_t = REAL(1.0)/rho_t.val();
	
	//
	//.. Compute the ROE Averages
	//
	//.... some useful quantities
	REAL    coef1 = sqrt(rho_f.val());
	REAL    coef2 = sqrt(rho_t.val());
	REAL    somme_coef = coef1 + coef2;
	REAL    isomme_coef = REAL(1.0)/somme_coef;
	REAL    u_f = rhou_f.val()*irho_f;
	REAL    v_f = rhov_f.val()*irho_f;
	REAL    w_f = rhow_f.val()*irho_f;
	REAL    h_f = (gam * rhoE_f.val()*irho_f) - (.5*gam1) * (u_f * u_f + v_f * v_f + w_f * w_f);
	REAL    u_t = rhou_t.val()*irho_t;
	REAL    v_t = rhov_t.val()*irho_t;
	REAL    w_t = rhow_t.val()*irho_t;
	REAL    h_t = (gam * rhoE_t.val()*irho_t) - (.5*gam1) * (u_t * u_t + v_t * v_t + w_t * w_t);
	
	//.... averages
	//REAL rho_ave = coef1 * coef2;
	REAL    u_ave = (coef1 * u_f + coef2 * u_t) * isomme_coef;
	REAL    v_ave = (coef1 * v_f + coef2 * v_t) * isomme_coef;
	REAL    w_ave = (coef1 * w_f + coef2 * w_t) * isomme_coef;
	REAL    h_ave = (coef1 * h_f + coef2 * h_t) * isomme_coef;
	//
	//.. Compute Speed of sound
	REAL    scal = u_ave * nx + v_ave * ny + w_ave * nz;
	REAL    norme = sqrt(nx * nx + ny * ny + nz * nz);
	REAL    inorme = REAL(1.0)/norme;
	REAL    u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
	REAL    c_speed = gam1 * (h_ave - REAL(0.5) * u2pv2pw2);
	if(c_speed < REAL(1e-6)) c_speed = 1e-6;// <!> zeroes the derivatives?   // avoid division by 0 if critical
	c_speed = sqrt(c_speed);
	REAL    c_speed2 = c_speed * norme;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	REAL    eig_val1 = scal - c_speed2;
	REAL    eig_val2 = scal;
	REAL    eig_val3 = scal + c_speed2;
	//
	//.. Compute the ROE flux
	//.... In this part many tests upon the eigenvalues
	//.... are done to simplify calculations
	//.... Here we use the two formes of the ROE flux :
	//.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
	//.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
	//
	if(eig_val2 <= REAL(0.0)) {
		T    irho_t = REAL(1.0)/rho_t;
		T    u_t = rhou_t*irho_t;
		T    v_t = rhov_t*irho_t;
		T    w_t = rhow_t*irho_t;
		T    h_t = (gam * rhoE_t*irho_t) - (.5*gam1) * (u_t * u_t + v_t * v_t + w_t * w_t);
		
		T p_t,ep_t,rhouv_t,rhouw_t,rhovw_t;
		REAL lambda_f,lambda_t;
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t +
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
		REAL p_f = gam1 * (rhoE_f.val() - REAL(0.5) * (rhou_f.val() * rhou_f.val() + rhov_f.val() * rhov_f.val()
													   + rhow_f.val() * rhow_f.val()) * irho_f);
		lambda_f = u_f * nx + v_f * ny + w_f * nz + norme
        * sqrt(gam * p_f * irho_f);
		lambda_t = u_t.val() * nx + v_t.val() * ny + w_t.val() * nz + norme
        * sqrt(gam * p_t.val() * irho_t.val());
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
		}
		//
		if (eig_val3 > REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the last
			//.. colomne of T-1 with the last row of T with eig_val3                //Cedric
			T    alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
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
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. Last columne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc + hnx;
			a3 = v_ave * usc + hny;
			a4 = w_ave * usc + hnz;
			a5 = REAL(0.5) * u2pv2pw2 * usc + REAL(2.5) * c_speed + scal;
			//.. Last row of the matrix T * eig_val3
			b1 = REAL(0.5) * (REAL(0.5) * tempo11 * u2pv2pw2 - scal);
			b2 = REAL(0.5) * (hnx - tempo11 * u_ave);
			b3 = REAL(0.5) * (hny - tempo11 * v_ave);
			b4 = REAL(0.5) * (hnz - tempo11 * w_ave);
			b5 = REAL(0.5) * tempo11;
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
	if(eig_val2 > REAL(0.0)) {
		T p_f,ep_f,rhouv_f,rhovw_f,rhouw_f;
		T    irho_f = REAL(1.0)/rho_f.val();
		
		T    u_f = rhou_f.val()*irho_f;
		T    v_f = rhov_f.val()*irho_f;
		T    w_f = rhow_f.val()*irho_f;
		T    h_f = (gam * rhoE_f.val()*irho_f) - (.5*gam1) * (u_f * u_f + v_f * v_f + w_f * w_f);
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f +
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
		REAL p_t = gam1 * (rhoE_t.val() - REAL(0.5) * (rhou_t.val() * rhou_t.val() +
													   rhov_t.val() * rhov_t.val() + rhow_t.val() * rhow_t.val()) * irho_t);
		REAL lambda_f = u_f.val() * nx + v_f.val() * ny + w_f.val() * nz - norme
        * sqrt(gam * p_f.val() * irho_f.val());
		REAL lambda_t   = u_t * nx + v_t * ny + w_t * nz - norme
        * sqrt(gam * p_t * irho_t);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
		}
		//
		if (eig_val1 < REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the first
			//.. columne of T-1 with the first row of T with eig_val1
			T    alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
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
			usc = REAL(1.0)/c_speed;
			tempo11 = gam1 * usc;
			//.. First colomne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc - hnx;
			a3 = v_ave * usc - hny;
			a4 = w_ave * usc - hnz;
			a5 = REAL(0.5) * u2pv2pw2 * usc + REAL(2.5) * c_speed - scal;
			//.. First row of the matrix T * eig_val1
			b1 = REAL(0.5) * (REAL(0.5) * tempo11 * u2pv2pw2 + scal);
			b2 = -REAL(0.5) * (hnx + tempo11 * u_ave);
			b3 = -REAL(0.5) * (hny + tempo11 * v_ave);
			b4 = -REAL(0.5) * (hnz + tempo11 * w_ave);
			b5 = REAL(0.5) * tempo11;
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
/*
 REAL operator() (const FADREAL &a)
 {
 return a.val();
 }*/

template <>
void TPZEulerConsLaw::ApproxRoe_Flux<FADREAL>(const FADREAL & rho_f,
											   const FADREAL & rhou_f,
											   const FADREAL & rhov_f,
											   const FADREAL & rhoE_f,
											   const FADREAL & rho_t,
											   const FADREAL & rhou_t,
											   const FADREAL & rhov_t,
											   const FADREAL & rhoE_t,
											   const REAL nx,
											   const REAL ny,
											   const REAL gam,
											   FADREAL &flux_rho,
											   FADREAL &flux_rhou,
											   FADREAL &flux_rhov,
											   FADREAL &flux_rhoE, int entropyFix)
{
	typedef FADREAL T;
	typedef REAL locREAL;
	locREAL rho_fv, rhou_fv, rhov_fv, rhoE_fv, rho_tv, rhou_tv, rhov_tv, rhoE_tv;
	
	rho_fv = rho_f.val();
	rhou_fv = rhou_f.val();
	rhov_fv = rhov_f.val();
	rhoE_fv = rhoE_f.val();
	rho_tv = rho_t.val();
	rhou_tv = rhou_t.val();
	rhov_tv = rhov_t.val();
	rhoE_tv = rhoE_t.val();
	
	
	//   rho_fv = rho_f;
	//   rhou_fv = rhou_f;
	//   rhov_fv = rhov_f;
	//   rhoE_fv = rhoE_f;
	//   rho_tv = rho_t;
	//   rhou_tv = rhou_t;
	//   rhov_tv = rhov_t;
	//   rhoE_tv = rhoE_t;
	
	locREAL    c_speed2;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	locREAL    eig_val1;
	locREAL    eig_val2;
	locREAL    eig_val3;
	locREAL    u_fv = rhou_fv/rho_fv;
	locREAL    v_fv = rhov_fv/rho_fv;
	locREAL    u_tv = rhou_tv/rho_tv;
	locREAL    v_tv = rhov_tv/rho_tv;
	REAL gam1 = gam - REAL(1.0);
	locREAL    p_tv = gam1 * (rhoE_tv - REAL(0.5) * (rhou_tv * rhou_tv + rhov_tv * rhov_tv) / rho_tv);
	locREAL    p_fv = gam1 * (rhoE_fv - REAL(0.5) * (rhou_fv * rhou_fv + rhov_fv * rhov_fv) / rho_fv);
	REAL norme = sqrt(nx * nx + ny * ny);
	locREAL scal,c_speed;
	locREAL u_ave,v_ave,u2pv2;
	{
		
		flux_rho = 0;
		flux_rhou = 0;
		flux_rhov = 0;
		flux_rhoE = 0;
		
		//REAL gam2 = gam * (gam - 1.0);
		//REAL igam = 1.0 / (gam - 1.0);
		
		//
		//.. Compute the ROE Averages
		//
		//.... some useful quantities
		locREAL    coef1 = sqrt(rho_fv);
		locREAL    coef2 = sqrt(rho_tv);
		locREAL    somme_coef = coef1 + coef2;
		locREAL    h_f = (gam * rhoE_fv/rho_fv) -  (u_fv * u_fv + v_fv * v_fv) * (gam1 / REAL(2.0));
		locREAL    h_t = (gam * rhoE_tv/rho_tv) -  (u_tv * u_tv + v_tv * v_tv) * (gam1 / REAL(2.0));
		
		//.... averages
		//REAL rho_ave = coef1 * coef2;
		//    cout << "Approx coef1 " << coef1 << "coef2 " << coef2 << "h_f " << h_f << "h_t " << h_t << endl;
		u_ave = (coef1 * u_fv + coef2 * u_tv) / somme_coef;
		v_ave = (coef1 * v_fv + coef2 * v_tv) / somme_coef;
		locREAL    h_ave = (coef1 * h_f + coef2 * h_t) / somme_coef;
		//
		//.. Compute Speed of sound
		scal = u_ave * nx + v_ave * ny;
		u2pv2 = u_ave * u_ave + v_ave * v_ave;
		c_speed = gam1 * (h_ave - REAL(0.5) * u2pv2);
		//    cout << "Approx c_speed " << c_speed << endl << "h_ave " << h_ave << endl << "u2pv2 " << u2pv2 << endl;
		if(c_speed < REAL(1e-6)) c_speed = REAL(1e-6);    // avoid division by 0 if critical
		c_speed = sqrt(c_speed);
		c_speed2 = c_speed * norme;
		//
		//.. Compute the eigenvalues of the Jacobian matrix
		eig_val1 = scal - c_speed2;
		eig_val2 = scal;
		eig_val3 = scal + c_speed2;
		//    cout << "Eigenvalues appox" << eig_val1 << endl << eig_val2 << endl << eig_val3 << endl;
	}
	//
	//.. Compute the ROE flux
	//.... In this part many tests upon the eigenvalues
	//.... are done to simplify calculations
	//.... Here we use the two formes of the ROE flux :
	//.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
	//.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
	//
	if(eig_val2 <= REAL(0.0)) {
		T p_t,ep_t;
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t + rhov_t * rhov_t) / rho_t);
		ep_t = rhoE_t + p_t;
		T    u_t = rhou_t/rho_t;
		T    v_t = rhov_t/rho_t;
		T rhouv_t = rhou_t * v_t;
		flux_rho  = rhou_t * nx + rhov_t * ny;
		flux_rhou = (rhou_t * u_t + p_t) * nx + rhouv_t * ny;
		flux_rhov = rhouv_t * nx + (rhov_t * v_t + p_t) * ny;
		flux_rhoE = ep_t * (u_t * nx + v_t * ny);
		//
		//.... A Entropic modification
		//
		locREAL lambda_f = u_fv * nx + v_fv * ny + norme * sqrt(gam * p_fv / rho_fv);
		locREAL lambda_t   = u_tv * nx + v_tv * ny + norme
        * sqrt(gam * p_tv / rho_tv);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val3 = lambda_t * (eig_val3 - lambda_f) / (lambda_t - lambda_f);
		}
		//
		if (eig_val3 > REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the last
			//.. colomne of T-1 with the last row of T with eig_val3
			T    alpha1,alpha2,alpha3,alpha4,alpha;
			locREAL a1,a2,a3,a4,b1,b2,b3,b4;
			T    delta_rho, delta_rhou,delta_rhov, delta_rhoE;
			delta_rho  = rho_t - rho_f;
			delta_rhou = rhou_t - rhou_f;
			delta_rhov = rhov_t - rhov_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal / norme;
			REAL hnx = nx / norme;
			REAL hny = ny / norme;
			locREAL usc = REAL(1.0)/c_speed;
			locREAL tempo11 = gam1 * usc;
			//.. Last columne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc + hnx;
			a3 = v_ave * usc + hny;
			a4 = REAL(0.5) * u2pv2 * usc + REAL(2.5) * c_speed + scal;
			//.. Last row of the matrix T * eig_val3
			b1 = REAL(0.5) * eig_val3 * (REAL(0.5) * tempo11 * u2pv2 - scal);
			b2 = REAL(0.5) * eig_val3 * (hnx - tempo11 * u_ave);
			b3 = REAL(0.5) * eig_val3 * (hny - tempo11 * v_ave);
			b4 = REAL(0.5) * eig_val3 * tempo11;
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
	if(eig_val2 > REAL(0.0)) {
		T p_f,ep_f;
		T    u_f = rhou_f/rho_f;
		T    v_f = rhov_f/rho_f;
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f +
											rhov_f * rhov_f) / rho_f);
		ep_f = rhoE_f + p_f;
		T rhouv_f = rhou_f * v_f;
		flux_rho  = rhou_f * nx + rhov_f * ny;
		flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny;
		flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny;
		flux_rhoE = ep_f * (u_f * nx + v_f * ny);
		//
		// A Entropic modification
		//
		locREAL lambda_f = u_fv * nx + v_fv * ny - norme * sqrt(gam * p_fv / rho_fv);
		locREAL lambda_t   = u_tv * nx + v_tv * ny - norme * sqrt(gam * p_tv / rho_tv);
		if (entropyFix && (lambda_f < REAL(0.)) && (lambda_t > REAL(0.))) {
			eig_val1 = lambda_f * (lambda_t - eig_val1) / (lambda_t - lambda_f);
		}
		//
		if (eig_val1 < REAL(0.0)) {
			//.. In this case A+ is obtained by multiplying the first
			//.. columne of T-1 with the first row of T with eig_val1
			T    alpha1,alpha2,alpha3,alpha4,alpha;
			locREAL a1,a2,a3,a4,b1,b2,b3,b4;
			T    delta_rho, delta_rhou,delta_rhov, delta_rhoE;
			delta_rho  = rho_t - rho_f;
			delta_rhou = rhou_t - rhou_f;
			delta_rhov = rhov_t - rhov_f;
			delta_rhoE = rhoE_t - rhoE_f;
			//
			scal = scal / norme;
			REAL hnx = nx / norme;
			REAL hny = ny / norme;
			locREAL usc = REAL(1.0)/c_speed;
			locREAL tempo11 = gam1 * usc;
			//.. First colomne of the matrix T-1
			a1 = usc;
			a2 = u_ave * usc - hnx;
			a3 = v_ave * usc - hny;
			a4 = REAL(0.5) * u2pv2 * usc + REAL(2.5) * c_speed - scal;
			//.. First row of the matrix T * eig_val1
			b1 = REAL(0.5) * eig_val1 * (REAL(0.5) * tempo11 * u2pv2 + scal);
			b2 = -REAL(0.5) * eig_val1 * (hnx + tempo11 * u_ave);
			b3 = -REAL(0.5) * eig_val1 * (hny + tempo11 * v_ave);
			b4 = REAL(0.5) * eig_val1 * tempo11;
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
#endif
/**
 Class identificator
 */
int TPZEulerConsLaw::ClassId() const {
    return TPZEULERCONSLAW2ID;
}

template class
TPZRestoreClass< TPZEulerConsLaw, TPZEULERCONSLAW2ID>;

