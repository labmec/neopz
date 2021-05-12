/**
 * \file
 * @brief Contains implementations of the TPZAUSMFlux methods.
 */

#include "pzausmflux.h"

TPZAUSMFlux::TPZAUSMFlux(STATE gamma) {
	this->fGamma = gamma;
	this->fAlpha = 3./16.; 
	this->fBeta = 1./8.;
}

TPZAUSMFlux::TPZAUSMFlux(const TPZAUSMFlux &cp) {
	this->fGamma = cp.fGamma;
	this->fAlpha = cp.fAlpha;
	this->fBeta = cp.fBeta;
}

void TPZAUSMFlux::ComputeFlux(TPZVec<STATE> &solL, TPZVec<STATE> &solR, TPZVec<REAL> &normal, TPZVec<STATE> & F) {
	STATE LeftSoundSpeed,RightSoundSpeed,
	LeftSpeed,RightSpeed,
	LeftNormalSpeed, RightNormalSpeed,
	LeftEnthalpy,RightEnthalpy,
	LeftPress, RightPress;
	
	this->ComputeInitialData(solL,normal,LeftSoundSpeed,LeftSpeed,LeftNormalSpeed,LeftEnthalpy,LeftPress);
	this->ComputeInitialData(solR,normal,RightSoundSpeed,RightSpeed,RightNormalSpeed,RightEnthalpy,RightPress);
	
	const STATE NumericalSoundSpeed = NumSoundSpeed(LeftSoundSpeed,RightSoundSpeed);
	const STATE LeftNumMach = LeftNormalSpeed/NumericalSoundSpeed;
	const STATE RightNumMach = RightNormalSpeed/NumericalSoundSpeed;
	
	const STATE FacePressure = this->FacePressure(LeftPress,RightPress,LeftNumMach,RightNumMach);
	const STATE FaceMach = this->FaceMachNumber(LeftNumMach,RightNumMach);
	const STATE MassFlux = this->MassFlux(NumericalSoundSpeed,solL[0], solR[0], FaceMach);
	
	const STATE uL = solL[1]/solL[0];
	const STATE vL = solL[2]/solL[0];
	const STATE wL = solL[3]/solL[0];
	
	const STATE uR = solR[1]/solR[0];
	const STATE vR = solR[2]/solR[0];
	const STATE wR = solR[3]/solR[0];
	
	F[0] = ((STATE)0.5)*(MassFlux* ((STATE)(2.)));
	F[1] = ((STATE)0.5)*MassFlux* (uL+uR) -((STATE)0.5)*fabs(MassFlux)* (uR-uL) + FacePressure*normal[0];
	F[2] = 0.5*MassFlux* (vL+vR) -0.5*fabs(MassFlux)* (vR-vL) + FacePressure*((STATE)normal[1]);
	F[3] = 0.5*MassFlux* (wL+wR) -0.5*fabs(MassFlux)* (wR-wL) + FacePressure*normal[2];
	F[4] = 0.5*MassFlux* (LeftEnthalpy+RightEnthalpy) -0.5*fabs(MassFlux)* (RightEnthalpy-LeftEnthalpy) + 0.;
	
}//void

void TPZAUSMFlux::ComputeInitialData(TPZVec<STATE>&sol,TPZVec<REAL> &normal, STATE&soundSpeed, 
                                     STATE &Speed, STATE &NormalSpeed, STATE &Enthalpy, STATE &press){
	press = this->Pressure(sol);
	soundSpeed = this->SoundSpeed(sol,press);
	Speed = this->Speed(sol,normal,NormalSpeed);
	Enthalpy = this->Enthalpy(soundSpeed,Speed);
}//void

STATE TPZAUSMFlux::SoundSpeed(TPZVec<STATE> &sol,STATE press){
	const STATE rho = sol[0];
	if(rho < 1e-10){
		PZError << "TPZEulerEquation(::cSpeed Too small or negative density\n";
		DebugStop();
	}
	
	const STATE c = sqrt(this->fGamma * press/ sol[0]);
	return c;
}

STATE TPZAUSMFlux::Pressure(TPZVec<STATE> &sol){
	if(fabs(sol[0]) < 1.e-6){
		PZError << "TPZEulerEquation::Pressure - Negative or too small density"
		<< sol[0] << std::endl;
		DebugStop();
	}
	
	STATE press = 0.0;
	
	const STATE rho_velocity2 = ( sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3] )/sol[0];
	press = ((this->fGamma-1.)*( sol[4] - 0.5 * rho_velocity2 ));
	
	if(press < 0){
		STATE temp = (this->fGamma-1.)*sol[4];
		PZError << "TPZEulerEquation::Pressure Negative pressure: " << press << " (gama-1)*E = " << temp << std::endl;
		DebugStop();
	}
	return press;
}//method

STATE TPZAUSMFlux::Speed(TPZVec<STATE> &sol, TPZVec<REAL> &normal, STATE &NormalSpeed){
	const STATE rho = sol[0];
	const STATE rhoU = sol[1];
	const STATE rhoV = sol[2];
	const STATE rhoW = sol[3];
	const STATE u = rhoU/rho;
	const STATE v = rhoV/rho;
	const STATE w = rhoW/rho;
	const STATE speed = sqrt(u*u+v*v+w*w);
	NormalSpeed = u*normal[0]+v*normal[1]+w*normal[2];
	return speed;
}//method

STATE TPZAUSMFlux::Enthalpy(STATE soundSpeed, STATE speed){
	const STATE enthalpy = soundSpeed*soundSpeed/(this->fGamma-1.)+0.5*speed*speed;
	return enthalpy;
}//method

STATE TPZAUSMFlux::FacePressure(STATE pL, STATE pR, STATE Ml, STATE Mr){
	STATE Pplus;
	if(Ml >= 1.) Pplus = 1.;
	else if(Ml <= -1.) Pplus = 0.;
	else{
		Pplus = 0.25*(Ml+1.)*(Ml+1.)*(2.-Ml)+this->fAlpha*Ml*(Ml*Ml-1.)*(Ml*Ml-1.);
	}
	
	STATE Pminus;
	if(Mr >= 1.) Pminus = 0.;
	else if(Mr <= -1.) Pminus = 1.;
	else{
		Pminus = 0.25*(Mr-1.)*(Mr-1.)*(2.+Mr)-this->fAlpha*Mr*(Mr*Mr-1.)*(Mr*Mr-1.);
	}
	
	const STATE pFace = pL*Pplus+pR*Pminus;
	return pFace;
	
}//method

STATE TPZAUSMFlux::FaceMachNumber(STATE Ml, STATE Mr){
	STATE Mplus;
	if(Ml >= 1.) Mplus = Ml;
	else if(Ml <= -1.) Mplus = 0.;
	else{
		Mplus = 0.25*(Ml+1.)*(Ml+1.)+this->fBeta*(Ml*Ml-1.)*(Ml*Ml-1.);
	}
	
	STATE Mminus;
	if(Mr >= 1.) Mminus = 0.;
	else if(Mr <= -1.) Mminus = Mr;
	else{
		Mminus = -0.25*(Mr-1.)*(Mr-1.)-this->fBeta*(Mr*Mr-1.)*(Mr*Mr-1.);
	}
	
	const STATE faceMach = Mplus+Mminus;
	return faceMach;
	
}//method

STATE TPZAUSMFlux::NumSoundSpeed(STATE LeftSoundSpeed,STATE RightSoundSpeed){
	return sqrt(LeftSoundSpeed*RightSoundSpeed);
}//method

STATE TPZAUSMFlux::MassFlux(STATE NumericalSoundSpeed, STATE rhoL, STATE rhoR, STATE FaceMach){
	STATE min = 0., max = 0.;
	if(FaceMach > 0.) max = FaceMach;
	if(FaceMach < 0.) min = FaceMach;
	
	const STATE m = NumericalSoundSpeed*(rhoL*max+rhoR*min);
	return m;
}//method
