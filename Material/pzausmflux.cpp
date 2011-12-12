/**
 * \file
 * @brief Contains implementations of the TPZAUSMFlux methods.
 */
//$Id: pzausmflux.cpp,v 1.3 2009-08-12 21:05:31 fortiago Exp $

#include "pzausmflux.h"

TPZAUSMFlux::TPZAUSMFlux(REAL gamma){
	this->fGamma = gamma;
	this->fAlpha = 3./16.; 
	this->fBeta = 1./8.;
}//method

TPZAUSMFlux::TPZAUSMFlux(const TPZAUSMFlux &cp){
	this->fGamma = cp.fGamma;
	this->fAlpha = cp.fAlpha;
	this->fBeta = cp.fBeta;
}

void TPZAUSMFlux::ComputeFlux(TPZVec<REAL> &solL, TPZVec<REAL> &solR, TPZVec<REAL> &normal, TPZVec<REAL> & F){
	REAL LeftSoundSpeed,RightSoundSpeed,
	LeftSpeed,RightSpeed,
	LeftNormalSpeed, RightNormalSpeed,
	LeftEnthalpy,RightEnthalpy,
	LeftPress, RightPress;
	
	this->ComputeInitialData(solL,normal,LeftSoundSpeed,LeftSpeed,LeftNormalSpeed,LeftEnthalpy,LeftPress);
	this->ComputeInitialData(solR,normal,RightSoundSpeed,RightSpeed,RightNormalSpeed,RightEnthalpy,RightPress);
	
	const REAL NumericalSoundSpeed = NumSoundSpeed(LeftSoundSpeed,RightSoundSpeed);
	const REAL LeftNumMach = LeftNormalSpeed/NumericalSoundSpeed;
	const REAL RightNumMach = RightNormalSpeed/NumericalSoundSpeed;
	
	const REAL FacePressure = this->FacePressure(LeftPress,RightPress,LeftNumMach,RightNumMach);
	const REAL FaceMach = this->FaceMachNumber(LeftNumMach,RightNumMach);
	const REAL MassFlux = this->MassFlux(NumericalSoundSpeed,solL[0], solR[0], FaceMach);
	
	//sol={rho,rhou,rhov,rhow,rhoe}
	const REAL uL = solL[1]/solL[0];
	const REAL vL = solL[2]/solL[0];
	const REAL wL = solL[3]/solL[0];
	
	const REAL uR = solR[1]/solR[0];
	const REAL vR = solR[2]/solR[0];
	const REAL wR = solR[3]/solR[0];
	
	F[0] = 0.5*MassFlux* (1.+1.) -0.5*fabs(MassFlux)* (1.-1.) + 0.;
	F[1] = 0.5*MassFlux* (uL+uR) -0.5*fabs(MassFlux)* (uR-uL) + FacePressure*normal[0];
	F[2] = 0.5*MassFlux* (vL+vR) -0.5*fabs(MassFlux)* (vR-vL) + FacePressure*normal[1];
	F[3] = 0.5*MassFlux* (wL+wR) -0.5*fabs(MassFlux)* (wR-wL) + FacePressure*normal[2];
	F[4] = 0.5*MassFlux* (LeftEnthalpy+RightEnthalpy) -0.5*fabs(MassFlux)* (RightEnthalpy-LeftEnthalpy) + 0.;
	
}//void

void TPZAUSMFlux::ComputeInitialData(TPZVec<REAL>&sol,TPZVec<REAL> &normal, REAL&soundSpeed, 
                                     REAL &Speed, REAL &NormalSpeed, REAL &Enthalpy, REAL &press){
	press = this->Pressure(sol);
	soundSpeed = this->SoundSpeed(sol,press);
	Speed = this->Speed(sol,normal,NormalSpeed);
	Enthalpy = this->Enthalpy(soundSpeed,Speed);
}//void

REAL TPZAUSMFlux::SoundSpeed(TPZVec<REAL> &sol,REAL press){
	const REAL rho = sol[0];
	if(rho < 1e-10){
		PZError << "TPZEulerEquation(::cSpeed Too small or negative density\n";
		DebugStop();
	}
	
//	const REAL temp = this->fGamma * press;	
//	if(temp < 1e-10){ // too low or negative
//		PZError << "TPZEulerEquation::cSpeed Too low or negative numerator\n";
//	}
	
	const REAL c = sqrt(this->fGamma * press/ sol[0]);
	return c;
}//method

REAL TPZAUSMFlux::Pressure(TPZVec<REAL> &sol){
	if(fabs(sol[0]) < 1.e-6){
		PZError << "TPZEulerEquation::Pressure - Negative or too small density"
		<< sol[0] << std::endl;
		DebugStop();
	}
	
	REAL press = 0.0;
	
	//sol = (rho , rho u , rho v , rho w , rho e)
	const REAL rho_velocity2 = ( sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3] )/sol[0];
	press = ((this->fGamma-1.)*( sol[4] - 0.5 * rho_velocity2 ));
	
	if(press < 0){
		REAL temp = (this->fGamma-1.)*sol[4];
		PZError << "TPZEulerEquation::Pressure Negative pressure: " << press << " (gama-1)*E = " << temp << std::endl;
		DebugStop();
	}
	return press;
}//method

REAL TPZAUSMFlux::Speed(TPZVec<REAL> &sol, TPZVec<REAL> &normal, REAL &NormalSpeed){
	const REAL rho = sol[0];
	const REAL rhoU = sol[1];
	const REAL rhoV = sol[2];
	const REAL rhoW = sol[3];
	const REAL u = rhoU/rho;
	const REAL v = rhoV/rho;
	const REAL w = rhoW/rho;
	const REAL speed = sqrt(u*u+v*v+w*w);
	NormalSpeed = u*normal[0]+v*normal[1]+w*normal[2];
	return speed;
}//method

REAL TPZAUSMFlux::Enthalpy(REAL soundSpeed, REAL speed){
	const REAL enthalpy = soundSpeed*soundSpeed/(this->fGamma-1.)+0.5*speed*speed;
	return enthalpy;
}//method

REAL TPZAUSMFlux::FacePressure(REAL pL, REAL pR, REAL Ml, REAL Mr){
	REAL Pplus;
	if(Ml >= 1.) Pplus = 1.;
	else if(Ml <= -1.) Pplus = 0.;
	else{
		Pplus = 0.25*(Ml+1.)*(Ml+1.)*(2.-Ml)+this->fAlpha*Ml*(Ml*Ml-1.)*(Ml*Ml-1.);
	}
	
	REAL Pminus;
	if(Mr >= 1.) Pminus = 0.;
	else if(Mr <= -1.) Pminus = 1.;
	else{
		Pminus = 0.25*(Mr-1.)*(Mr-1.)*(2.+Mr)-this->fAlpha*Mr*(Mr*Mr-1.)*(Mr*Mr-1.);
	}
	
	const REAL pFace = pL*Pplus+pR*Pminus;
	return pFace;
	
}//method

REAL TPZAUSMFlux::FaceMachNumber(REAL Ml, REAL Mr){
	REAL Mplus;
	if(Ml >= 1.) Mplus = Ml;
	else if(Ml <= -1.) Mplus = 0.;
	else{
		Mplus = 0.25*(Ml+1.)*(Ml+1.)+this->fBeta*(Ml*Ml-1.)*(Ml*Ml-1.);
	}
	
	REAL Mminus;
	if(Mr >= 1.) Mminus = 0.;
	else if(Mr <= -1.) Mminus = Mr;
	else{
		Mminus = -0.25*(Mr-1.)*(Mr-1.)-this->fBeta*(Mr*Mr-1.)*(Mr*Mr-1.);
	}
	
	const REAL faceMach = Mplus+Mminus;
	return faceMach;
	
}//method

REAL TPZAUSMFlux::NumSoundSpeed(REAL LeftSoundSpeed,REAL RightSoundSpeed){
	return sqrt(LeftSoundSpeed*RightSoundSpeed);
}//method

REAL TPZAUSMFlux::MassFlux(REAL NumericalSoundSpeed, REAL rhoL, REAL rhoR, REAL FaceMach){
	REAL min = 0., max = 0.;
	if(FaceMach > 0.) max = FaceMach;
	if(FaceMach < 0.) min = FaceMach;
	
	const REAL m = NumericalSoundSpeed*(rhoL*max+rhoR*min);
	return m;
}//method

/*  {
 const REAL rho = sol[0];
 const REAL rhoU = sol[1];
 const REAL rhoV = sol[2];
 const REAL rhoW = sol[3];
 const REAL rhoE = sol[4];
 const REAL u = rhoU/rho;
 const REAL v = rhoV/rho;
 const REAL w = rhoW/rho;
 const REAL p = this->Pressure(sol);
 F.Resize(5,3);
 F(0,0) = rhoU;        F(0,1) = rhoV;        F(0,2) = rhoW;
 F(1,0) = rhoU*u+p;    F(1,1) = rhoU*v;      F(1,2) = rhoU*w;
 F(2,0) = rhoV*u;      F(2,1) = rhoV*v+p;    F(2,2) = rhoV*w;
 F(3,0) = rhoW*u;      F(3,1) = rhoW*v;      F(3,2) = rhoW*w+p;
 F(4,0) = (rhoE+p)*u;  F(4,1) = (rhoE+p)*v;  F(4,2) = (rhoE+p)*w;
 }//void
 */
