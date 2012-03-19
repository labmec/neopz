/**
 * \file
 * @brief Contains implementations of the TPZDiffusionConsLaw methods.
 */
#include "TPZDiffusionConsLaw.h"
#include "TPZCompElDisc.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"
using namespace std;

REAL TPZDiffusionConsLaw::fGamma = 1.4;
REAL TPZDiffusionConsLaw::fDelta = -1.0;
REAL TPZDiffusionConsLaw::fCFL = 0.0;
std::string TPZDiffusionConsLaw::fArtificialDiffusion = "LS";

TPZDiffusionConsLaw::~TPZDiffusionConsLaw(){}

TPZDiffusionConsLaw::TPZDiffusionConsLaw(){
	
	fDimension = 0;
	fA.Redim(0,0);
	fB.Redim(0,0);
	fC.Redim(0,0);
}

TPZDiffusionConsLaw::TPZDiffusionConsLaw(TPZVec<REAL> U,REAL gamma,int dim,const std::string &diff) :
fA(0,0), fB(0,0), fC(0,0) {
	fGamma = gamma;
	fDimension = dim;
	fArtificialDiffusion = diff;
	int nstate = 2 + dim;
	fA.Redim(nstate,nstate);
	fB.Redim(nstate,nstate);
	fC.Redim(nstate,nstate);
	JacobFlux(U,fA,fB,fC);
}

void TPZDiffusionConsLaw::GradientOfTheFlow(TPZFMatrix<REAL> &DF1,TPZFMatrix<REAL> &DF2,TPZFMatrix<REAL> &DF3){
	
	DF1 = fA;
	DF2 = fB;
	DF3 = fC;
}

REAL TPZDiffusionConsLaw::CFL(int degree){
	
	if(fCFL != 0.) return fCFL;
	return (1.0/(2.0*(REAL)degree+1.0));
}

REAL TPZDiffusionConsLaw::Delta(){
	return fDelta;
}

REAL TPZDiffusionConsLaw::DeltaOtimo(){
	
	int degree = TPZCompEl::GetgOrder();
	REAL cfl = CFL(degree);
	REAL delta = ( (10./3.)*cfl*cfl - (2./3.)*cfl + 1./10. );
	return delta;
}

void TPZDiffusionConsLaw::Divergence(TPZVec<REAL> &dphi,TPZFMatrix<REAL> &diverg){
	
	fArtificialDiffusion = "LS";
	std::string type = fArtificialDiffusion;
	TPZFMatrix<REAL> diffterm(0,0);
	PointOperator(dphi,diffterm);
	int i,j,nstate=2+fDimension;
	diverg.Redim(nstate,1);
	diverg.Zero();
	// soma das linhas de diffterm da : div(F) = �1/� [�/�] + �2/� [�/�] + �3/� [�/�]
	for(i=0;i<nstate;i++) for(j=0;j<nstate;j++) diverg(i,0) += diffterm(i,j);// 5x1 (3D)
	fArtificialDiffusion = type;
}

void TPZDiffusionConsLaw::PointOperator(TPZVec<REAL> &dphi,TPZFMatrix<REAL> &diff_term){
	
	diff_term.Zero();
	int size = dphi.NElements();
	int nstate = 2+fDimension;
	if(size != 3){
		//esta forma permite trabalhar com 1, 2 e 3 dimens�s
		cout << "TPZDiffusionConsLaw::PointOperator error data size";
	}
	TPZFMatrix<REAL> Tx(nstate,nstate),Ty(nstate,nstate),Tz(nstate,nstate);
	TPZFMatrix<REAL> Trx(nstate,nstate),Try(nstate,nstate),Trz(nstate,nstate);
	diff_term.Redim(nstate,nstate);
	Tau(Tx,Ty,Tz);
	Tx.Transpose(&Trx);
	Ty.Transpose(&Try);
	Tz.Transpose(&Trz);
	int i,j;
	// �/� *  (T1)'  + �/� * (T2)'  + �/� * (T3)'
	for(i=0;i<nstate;i++){
		for(j=0;j<nstate;j++){
			diff_term(i,j) = dphi[0] * Trx(i,j) + dphi[1] * Try(i,j) + dphi[2] * Trz(i,j);
		}
	}
}

void TPZDiffusionConsLaw::Tau(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz){
	
	if(!strcmp(fArtificialDiffusion.c_str(),"SUPG")){
		SUPG(Tx,Ty,Tz);
		return;
	}
	if(!strcmp(fArtificialDiffusion.c_str(),"LS")){
		LS(Tx,Ty,Tz);
		return;
	}
	if(!strcmp(fArtificialDiffusion.c_str(),"BORNHAUS")){
		Bornhaus(Tx,Ty,Tz);
		return;
	}
	cout << "\nTPZDiffusionConsLaw::Tau case not implemented\n";
}

void TPZDiffusionConsLaw::SUPG(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz){
	cout << "TPZDiffusionConsLaw:: SUPG artificial diffusion SUPG not implemented\n";
}

void TPZDiffusionConsLaw::LS(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz){
	
	fA.Transpose(&Tx);
	fB.Transpose(&Ty);
	fC.Transpose(&Tz);
}

void TPZDiffusionConsLaw::Bornhaus(TPZFMatrix<REAL> &Tx,TPZFMatrix<REAL> &Ty,TPZFMatrix<REAL> &Tz){
    cout << "TPZDiffusionConsLaw::Bornhaus artificial diffusion Bornhaus not implemented\n";
}

void TPZDiffusionConsLaw::JacobFlux(TPZVec<REAL> U,TPZFMatrix<REAL> &A,TPZFMatrix<REAL> &B,TPZFMatrix<REAL> &C) {//OK
	
	if(fabs(U[0]) < 1.e-6) {
		cout << "\nTPZDiffusionConsLaw::JacobFlux: Densidade quase nula, o jacobiano falha\n";
		cout << "Densidade = " << U[0] << endl;
		U[0] = 1.e-6;
		cout << "Nova densidade = " << U[0] << endl;
		//cout << "Programa abortado\n\n";
		//exit(-1);
	}
	
	int cap = A.Rows();
	if(cap < 3 || cap > 5){
		cout << "\n\nTPZDiffusionConsLaw::JacobFlux case not trated\n\n";
		A.Redim(0,0);
		B.Redim(0,0);
		C.Redim(0,0);
		exit(-1);
	}
	
	REAL u,v,w,e;
	REAL gamma1 = fGamma-1.;
	REAL gamma2 = gamma1/2.;
	REAL gamma3 = fGamma-3;
	
	if(cap == 5){
		
		A.Redim(5,5);
		B.Redim(5,5);
		C.Redim(5,5);
		
		u = U[1]/U[0];
		v = U[2]/U[0];
		w = U[3]/U[0];
		e = U[4]/U[0];
		
		REAL u2 = u*u;
		REAL v2 = v*v;
		REAL w2 = w*w;
		REAL vel = u2+v2+w2;
		
		A(0,0) = 0.;
		A(0,1) = 1.;
		A(0,2) = 0.;
		A(0,3) = 0.;
		A(0,4) = 0.;
		
		A(1,0) =  gamma2*vel-u2;
		A(1,1) = -gamma3*u;
		A(1,2) = -gamma1*v;
		A(1,3) = -gamma1*w;
		A(1,4) =  gamma1;
		
		A(2,0) = -u*v;
		A(2,1) =  v;
		A(2,2) =  u;
		A(2,3) =  0.;
		A(2,4) =  0.;
		
		A(3,0) = -u*w;
		A(3,1) =  w;
		A(3,2) =  0.;
		A(3,3) =  u;
		A(3,4) =  0.;
		
		A(4,0) = -fGamma*e*u + gamma1*u*vel;
		A(4,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
		A(4,2) = -gamma1*u*v;
		A(4,3) = -gamma1*u*w;
		A(4,4) =  fGamma*u;
		
		B(0,0) = 0.;
		B(0,1) = 0.;
		B(0,2) = 1.;
		B(0,3) = 0.;
		B(0,4) = 0.;
		
		B(1,0) = -u*v;
		B(1,1) =  v;
		B(1,2) =  u;
		B(1,3) =  0.;
		B(1,4) =  0.;
		
		B(2,0) =  gamma2*vel-v2;
		B(2,1) = -gamma1*u;
		B(2,2) = -gamma3*v;
		B(2,3) = -gamma1*w;
		B(2,4) =  gamma1;
		
		B(3,0) = -v*w;
		B(3,1) =  0.;
		B(3,2) =  w;
		B(3,3) =  v;
		B(3,4) =  0.;
		
		B(4,0) = -fGamma*e*v + gamma1*v*vel;
		B(4,1) = -gamma1*u*v;
		B(4,2) =  fGamma*e - gamma1*v2 - gamma2*vel;
		B(4,3) = -gamma1*v*w;
		B(4,4) =  fGamma*v;
		
		C(0,0) = 0.;
		C(0,1) = 0.;
		C(0,2) = 0.;
		C(0,3) = 1.;
		C(0,4) = 0.;
		
		C(1,0) = -u*w;
		C(1,1) =  w;
		C(1,2) =  0.;
		C(1,3) =  u;
		C(1,4) =  0.;
		
		C(2,0) = -v*w;
		C(2,1) =  0.;
		C(2,2) =  w;
		C(2,3) =  v;
		C(2,4) =  0.;
		
		C(3,0) =  gamma2*vel-w2;
		C(3,1) = -gamma1*u;
		C(3,2) = -gamma1*v;
		C(3,3) = -gamma3*w;
		C(3,4) =  gamma1;
		
		C(4,0) = -fGamma*e*w + gamma1*w*vel;
		C(4,1) = -gamma1*u*w;
		C(4,2) = -gamma1*v*w;
		C(4,3) =  fGamma*e - gamma1*w2 - gamma2*vel;
		C(4,4) =  fGamma*w;
		
	} else if(cap == 4){
		
		A.Redim(4,4);
		B.Redim(4,4);
		C.Redim(4,4);
		
		u = U[1]/U[0];
		v = U[2]/U[0];
		e = U[3]/U[0];
		
		REAL u2 = u*u;
		REAL v2 = v*v;
		REAL vel = u2+v2;
		
		A(0,0) = 0.;
		A(0,1) = 1.;
		A(0,2) = 0.;
		A(0,3) = 0.;
		
		A(1,0) =  gamma2*vel-u2;
		A(1,1) = -gamma3*u;
		A(1,2) = -gamma1*v;
		A(1,3) =  gamma1;
		
		A(2,0) = -u*v;
		A(2,1) =  v;
		A(2,2) =  u;
		A(2,3) =  0.;
		
		A(3,0) = -fGamma*e*u + gamma1*u*vel;
		A(3,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
		A(3,2) = -gamma1*u*v;
		A(3,3) =  fGamma*u;
		
		B(0,0) = 0.;
		B(0,1) = 0.;
		B(0,2) = 1.;
		B(0,3) = 0.;
		
		B(1,0) = -u*v;
		B(1,1) =  v;
		B(1,2) =  u;
		B(1,3) =  0.;
		
		B(2,0) =  gamma2*vel-v2;
		B(2,1) = -gamma1*u;
		B(2,2) = -gamma3*v;
		B(2,3) =  gamma1;
		
		B(3,0) = -fGamma*e*v + gamma1*v*vel;
		B(3,1) = -gamma1*u*v;
		B(3,2) =  fGamma*e - gamma1*v2 - gamma2*vel;
		B(3,3) =  fGamma*v;
		
	} else if(cap == 3){
		
		A.Redim(3,3);
		B.Redim(3,3);
		C.Redim(3,3);
		
		u = U[1]/U[0];
		e = U[2]/U[0];
		
		REAL u2 = u*u;
		REAL vel = u2;
		
		A(0,0) = 0.;
		A(0,1) = 1.;
		A(0,2) = 0.;
		
		A(1,0) =  gamma2*vel-u2;
		A(1,1) = -gamma3*u;
		A(1,2) =  gamma1;
		
		A(2,0) = -fGamma*e*u + gamma1*u*vel;
		A(2,1) =  fGamma*e - gamma1*u2 - gamma2*vel;
		A(2,2) =  fGamma*u;
	}
}

//left = **_f    right = **_t
void TPZDiffusionConsLaw::Roe_Flux(
								   REAL rho_f, REAL rhou_f, REAL rhov_f, REAL rhow_f,
								   REAL rhoE_f,REAL rho_t, REAL rhou_t, REAL rhov_t, REAL rhow_t,
								   REAL rhoE_t, REAL nx, REAL ny, REAL nz, REAL gam,
								   REAL &flux_rho, REAL &flux_rhou, REAL &flux_rhov,
								   REAL &flux_rhow, REAL &flux_rhoE){
	
	REAL alpha1,alpha2,alpha3,alpha4,alpha5,alpha;
	REAL a1,a2,a3,a4,a5,b1,b2,b3,b4,b5;
	REAL ep_t, ep_f, p_t, p_f;
	REAL rhouv_t, rhouv_f, rhouw_t, rhouw_f, rhovw_t, rhovw_f;
	REAL lambda_f, lambda_t;
	REAL delta_rho, delta_rhou, delta_rhov, delta_rhow, delta_rhoE;
	REAL hnx, hny, hnz;
	REAL tempo11, usc;
	
	flux_rho = 0;
	flux_rhou = 0;
	flux_rhov = 0;
	flux_rhow = 0;
	flux_rhoE = 0;
	
	REAL gam1 = gam - 1.0;
	REAL irho_f = 1.0/rho_f;
	REAL irho_t = 1.0/rho_t;
	
	//
	//.. Compute the ROE Averages
	//
	//.... some useful quantities
	REAL coef1 = sqrt(rho_f);
	REAL coef2 = sqrt(rho_t);
	REAL somme_coef = coef1 + coef2;
	REAL isomme_coef = 1.0/somme_coef;
	REAL u_f = rhou_f*irho_f;
	REAL v_f = rhov_f*irho_f;
	REAL w_f = rhow_f*irho_f;
	REAL h_f = (gam * rhoE_f*irho_f) - (.5*gam1) * (u_f * u_f + v_f * v_f + w_f * w_f);
	REAL u_t = rhou_t*irho_t;
	REAL v_t = rhov_t*irho_t;
	REAL w_t = rhow_t*irho_t;
	REAL h_t = (gam * rhoE_t*irho_t) - (.5*gam1) * (u_t * u_t + v_t * v_t + w_t * w_t);
	
	//.... averages
	//REAL rho_ave = coef1 * coef2;
	REAL u_ave = (coef1 * u_f + coef2 * u_t) * isomme_coef;
	REAL v_ave = (coef1 * v_f + coef2 * v_t) * isomme_coef;
	REAL w_ave = (coef1 * w_f + coef2 * w_t) * isomme_coef;
	REAL h_ave = (coef1 * h_f + coef2 * h_t) * isomme_coef;
	//
	//.. Compute Speed of sound
	REAL scal = u_ave * nx + v_ave * ny + w_ave * nz;
	REAL norme = sqrt(nx * nx + ny * ny + nz * nz);
	REAL inorme = 1.0/norme;
	REAL u2pv2pw2 = u_ave * u_ave + v_ave * v_ave + w_ave * w_ave;
	REAL c_speed = gam1 * (h_ave - 0.5 * u2pv2pw2);
	if(c_speed < 1e-6) c_speed = 1e-6;    // avoid division by 0 if critical
	c_speed = sqrt(c_speed);
	REAL c_speed2 = c_speed * norme;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	REAL eig_val1 = scal - c_speed2;
	REAL eig_val2 = scal;
	REAL eig_val3 = scal + c_speed2;
	//
	//.. Compute the ROE flux
	//.... In this part many tests upon the eigenvalues
	//.... are done to simplify calculations
	//.... Here we use the two formes of the ROE flux :
	//.... phi(Wl,Wr) = F(Wl) + A-(Wroe)(Wr - Wl)
	//.... phi(Wl,Wr) = F(Wr) - A+(Wroe)(Wr - Wl)
	//
	if(eig_val2 <= 0.0) {
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
									  rhov_t * rhov_t + rhow_t * rhow_t) * irho_t);
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
void TPZDiffusionConsLaw::Roe_Flux(REAL rho_f, REAL rhou_f, REAL rhov_f, REAL rhoE_f,
								   REAL rho_t, REAL rhou_t, REAL rhov_t, REAL rhoE_t,
								   REAL nx, REAL ny, REAL gam,
								   REAL &flux_rho, REAL &flux_rhou,REAL &flux_rhov, REAL &flux_rhoE){
	
	REAL alpha1,alpha2,alpha3,alpha4,a1,a2,a3,a4,b1,b2,b3,b4,alpha;
	REAL ep_t, ep_f, p_t, p_f;
	REAL rhouv_t, rhouv_f;
	REAL lambda_f, lambda_t;
	REAL delta_rho, delta_rhou,delta_rhov, delta_rhoE;
	REAL hnx, hny;
	REAL tempo11, usc;
	
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
	REAL coef1 = sqrt(rho_f);
	REAL coef2 = sqrt(rho_t);
	REAL somme_coef = coef1 + coef2;
	REAL u_f = rhou_f/rho_f;
	REAL v_f = rhov_f/rho_f;
	REAL h_f = (gam * rhoE_f/rho_f) - (gam1 / 2.0) * (u_f * u_f + v_f * v_f);
	REAL u_t = rhou_t/rho_t;
	REAL v_t = rhov_t/rho_t;
	REAL h_t = (gam * rhoE_t/rho_t) - (gam1 / 2.0) * (u_t * u_t + v_t * v_t);
	
	//.... averages
	//REAL rho_ave = coef1 * coef2;
	REAL u_ave = (coef1 * u_f + coef2 * u_t) / somme_coef;
	REAL v_ave = (coef1 * v_f + coef2 * v_t) / somme_coef;
	REAL h_ave = (coef1 * h_f + coef2 * h_t) / somme_coef;
	//
	//.. Compute Speed of sound
	REAL scal = u_ave * nx + v_ave * ny;
	REAL norme = sqrt(nx * nx + ny * ny);
	REAL u2pv2 = u_ave * u_ave + v_ave * v_ave;
	REAL c_speed = gam1 * (h_ave - 0.5 * u2pv2);
	if(c_speed < 1e-6) c_speed = 1e-6;    // avoid division by 0 if critical
	c_speed = sqrt(c_speed);
	REAL c_speed2 = c_speed * norme;
	//
	//.. Compute the eigenvalues of the Jacobian matrix
	REAL eig_val1 = scal - c_speed2;
	REAL eig_val2 = scal;
	REAL eig_val3 = scal + c_speed2;
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
		p_f = gam1 * (rhoE_f - REAL(0.5) * (rhou_f * rhou_f +
											rhov_f * rhov_f) / rho_f);
		ep_f = rhoE_f + p_f;
		rhouv_f = rhou_f * v_f;
		flux_rho  = rhou_f * nx + rhov_f * ny;
		flux_rhou = (rhou_f * u_f + p_f) * nx + rhouv_f * ny;
		flux_rhov = rhouv_f * nx + (rhov_f * v_f + p_f) * ny;
		flux_rhoE = ep_f * (u_f * nx + v_f * ny);
		//
		// A Entropic modification
		//
		p_t = gam1 * (rhoE_t - REAL(0.5) * (rhou_t * rhou_t +
											rhov_t * rhov_t) / rho_t);
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
