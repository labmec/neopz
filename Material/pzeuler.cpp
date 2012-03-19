/**
 * \file
 * @brief Contains implementations of the TPZEulerEquation methods.
 */
//$Id: pzeuler.cpp,v 1.8 2011-02-11 17:27:11 phil Exp $

#include "pzeuler.h"

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

TPZEulerEquation::CALCType TPZEulerEquation::gType = EFlux;

REAL TPZEulerEquation::gGamma = 1.4;

#ifdef LinearConvection
void TPZEulerEquation::SetLinearConvection(TPZCompMesh * cmesh, TPZVec<REAL> &Celerity){
	std::map<int, TPZAutoPointer<TPZMaterial> >::iterator w;
	for(w = cmesh->MaterialVec().begin(); w != cmesh->MaterialVec().end(); w++){
		TPZMaterial * mat = w->second.operator->();
		TPZEulerEquation * matcast = dynamic_cast< TPZEulerEquation * >(mat);
		if(matcast){
			cout << "old celerity = ";
			for(int i = 0; i < matcast->fCelerity.NElements(); i++)  cout << matcast->fCelerity[i] << "  ";
			cout << "\n";
			matcast->fCelerity = Celerity;
			cout << "new celerity = ";
			for(int i = 0; i < matcast->fCelerity.NElements(); i++)  cout << matcast->fCelerity[i] << "  ";
			cout << "\n";
		}
	}
}
#endif

void TPZEulerEquation::FromPrimitiveToConservative(TPZVec<REAL> &sol,REAL gamma){
	
#ifdef LinearConvection
	return;
#endif
	
	double keepP = sol[4];
	//sol = {rho, u, v, w, p}
	double rhoE = 0.5*sol[0]*(sol[1]*sol[1]+sol[2]*sol[2]+sol[3]*sol[3])+sol[4]/(gamma-1.);
	sol[1] = sol[1]*sol[0];
	sol[2] = sol[2]*sol[0];
	sol[3] = sol[3]*sol[0];
	sol[4] = rhoE;
	double p = TPZEulerEquation::Pressure(sol,gamma);
	if(fabs(p-keepP) > 1e-4){
		std::cout << "\np = " << p << "  keepP = " << keepP << "\n";
	}
}//void

void TPZEulerEquation::FromConservativeToPrimitive(TPZVec<REAL> &sol,REAL gamma){
	
#ifdef LinearConvection
	return;
#endif
	// The material must to have not null density
	if(IsZero(sol[0])) {
		PZError << "\nTPZEulerEquation::FromConservativeToPrimitive: Density almost null\n" << "Density = " << sol[0] << std::endl;
		DebugStop();
	}
	
	double p = TPZEulerEquation::Pressure(sol,gamma);
	sol[1] = sol[1]/sol[0];
	sol[2] = sol[2]/sol[0];
	sol[3] = sol[3]/sol[0];
	sol[4] = p;
}//void

TPZEulerEquation::~TPZEulerEquation(){
	
}

TPZEulerEquation::TPZEulerEquation(int nummat, REAL gamma) : 
TPZDiscontinuousGalerkin(nummat),fAUSMFlux(gamma),fGradientFlux(){
	gGamma = gamma;
}

TPZEulerEquation::TPZEulerEquation():TPZDiscontinuousGalerkin(),fAUSMFlux(-1.),fGradientFlux(){
	
}

TPZEulerEquation::TPZEulerEquation(const TPZEulerEquation &cp) : 
TPZDiscontinuousGalerkin(cp),fAUSMFlux(cp.fAUSMFlux),fGradientFlux(cp.fGradientFlux){
	
}

TPZAutoPointer<TPZMaterial> TPZEulerEquation::NewMaterial(){
	return new TPZEulerEquation(*this);
}

int TPZEulerEquation::NStateVariables() {
	return 5; //U = (rho, rhou, rhov, rhow, rhoe)
}

int TPZEulerEquation::Dimension(){
	return 3;
}

void TPZEulerEquation::Print(std::ostream &out) {
	TPZDiscontinuousGalerkin::Print(out);
	out << "gGamma = " << gGamma << "\n";
}

int TPZEulerEquation::VariableIndex(const std::string &name) {
	if( !strcmp(name.c_str(),"density")  )     return 1;//rho
	if( !strcmp(name.c_str(),"velocity") )     return 2;//(u,v,w)
	if( !strcmp(name.c_str(),"energy")   )     return 3;//rhoE
	if( !strcmp(name.c_str(),"pressure") )     return 4;//p
	if( !strcmp(name.c_str(),"solution") )     return 5;//(ro,u,v,w,E)
	if( !strcmp(name.c_str(),"normvelocity") ) return 6;//sqrt(u+v+w)
	if( !strcmp(name.c_str(),"Mach") )         return 7;//sqrt(u+v+w)/c
	std::cout << "TPZEulerEquation::VariableIndex not defined\n";
	return TPZMaterial::VariableIndex(name);
}

int TPZEulerEquation::NSolutionVariables(int var){
	
	if(var == 1 || var == 3 || var == 4 || var == 6 || var == 7) return 1;
	if(var == 2) return Dimension();
	if(var == 5) return NStateVariables();
	
	std::cout << "TPZEulerEquation::NSolutionVariables not defined\n";
	return 0;
}

void TPZEulerEquation::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout){
	
#ifndef LinearConvection
	if(IsZero(Sol[0])){
		PZError << "\nTPZEulerEquation::Solution: Density almost null\n" << "Density = " << Sol[0] << std::endl;
		DebugStop();
	}
#endif
	
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
		Solout[0] = Sol[pos];//energy = rhoE
		return;
	} else if(var == 4) {
		Solout.Resize(1);
		Solout[0] = Pressure(Sol,gGamma);//pressure
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
		const REAL cspeed = this->cSpeed(Sol);
		const REAL us = this->uRes(Sol);
		Solout[0] = us / cspeed;
		return;
	} else {
		//cout << "TPZEulerEquation::Solution variable in the base class\n";
		TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
	}
}

void TPZEulerEquation::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
	//nothing to be done here
	std::cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}

void TPZEulerEquation::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ef){
	//nothing to be done here
	std::cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}

void TPZEulerEquation::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
	this->ContributeInterface(data,dataleft,dataright,weight,ef);
	std::cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}


void TPZEulerEquation::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ef){
#ifdef LinearConvection
	if(gType == EFlux){
		fGradientFlux.ApplyLimiter(data);
		TPZManVector<REAL,15> Flux(5,0.);
		double dot = 0.;
		for(int i = 0; i < 3; i++) dot += data.normal[i]*fCelerity[i];
		if(dot > 0.){
			Flux[0] = dot*dataleft.sol[0];
		}
		else{
			Flux[0] = dot*dataright.sol[0];
		}
		
		for(int i = 0; i < 5; i++) ef(i,0)   += +1. * weight*Flux[i];
		for(int i = 0; i < 5; i++) ef(i+20,0) += -1. * weight*Flux[i];
	}
	if(gType == EGradient){
		TPZManVector<REAL,15> Flux(15);
		fGradientFlux.ComputeFlux(dataleft.sol,dataright.sol,data.normal,Flux);
		
		for(int i = 0; i < 15; i++) ef(i+5,0)    += +1. * weight*Flux[i];
		for(int i = 0; i < 15; i++) ef(i+20+5,0) += -1. * weight*Flux[i];
	}
	return;
#endif
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	if(gType == EFlux){
		fGradientFlux.ApplyLimiter(data,dataleft,dataright);
		TPZEulerEquation::FromPrimitiveToConservative(dataleft.sol[0],gGamma);
		TPZEulerEquation::FromPrimitiveToConservative(dataright.sol[0],gGamma);
		TPZManVector<REAL,15> Flux(5);
		fAUSMFlux.ComputeFlux(dataleft.sol[0],dataright.sol[0],data.normal,Flux);
		
		for(int i = 0; i < 5; i++) ef(i,0)   += +1. * weight*Flux[i];
		for(int i = 0; i < 5; i++) ef(i+20,0) += -1. * weight*Flux[i];
	}
	if(gType == EGradient){
		TPZManVector<REAL,15> Flux(15);
		fGradientFlux.ComputeFlux(dataleft.sol[0],dataright.sol[0],data.normal,Flux);
		
		for(int i = 0; i < 15; i++) ef(i+5,0)    += +1. * weight*Flux[i];
		for(int i = 0; i < 15; i++) ef(i+20+5,0) += -1. * weight*Flux[i];
	}
	
}//void

void TPZEulerEquation::ContributeBC(TPZMaterialData &data,
                                    REAL weight,
                                    TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef,
                                    TPZBndCond &bc){
	std::cout << __PRETTY_FUNCTION__ << " - this method should not be called for Finite Volume Method\n";
	DebugStop();
}

void TPZEulerEquation::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,
											 TPZBndCond &bc){
	this->ContributeBCInterface(data,dataleft,weight,ef,bc);
	std::cout << "\nWarning at " << __PRETTY_FUNCTION__ << " - this method should not be called";
}//void

void TPZEulerEquation::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZFMatrix<REAL> &ef,
											 TPZBndCond &bc){
    
    int numbersol = dataleft.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }

#ifdef LinearConvection
	if(gType == EFlux){
		if (bc.Type() == EFreeSlip){
			TPZFNMatrix<100> fakeef(2*ef.Rows(),ef.Cols(),0);
			dataright.sol = dataleft.sol;
			this->ContributeInterface(data,weight,fakeef);
			for(int i = 0; i < ef.Rows(); i++) ef(i,0) += fakeef(i,0);
		}//if FreeSlip
	}
	if(gType == EGradient){
		if (bc.Type() == EFreeSlip){
			TPZManVector<REAL,15> Flux(5);
			fGradientFlux.ComputeFlux(dataleft.sol[0],dataleft.sol[0],data.normal,Flux);
			for(int i = 0; i < 15; i++) ef(i+5,0)    += +1. * weight*Flux[i];
		}//if FreeSlip
	}
	return;
#endif
	if(gType == EFlux){
		TPZEulerEquation::FromPrimitiveToConservative(dataleft.sol[0],gGamma);
		if (bc.Type() == EFreeSlip){
			TPZManVector<REAL,15> Flux(5);
// #warning One needs to invert the velocity component
			fAUSMFlux.ComputeFlux(dataleft.sol[0],dataleft.sol[0],data.normal,Flux);
			for(int i = 0; i < 5; i++) ef(i,0)   += +1. * weight*Flux[i];
		}//if FreeSlip
	}
	if(gType == EGradient){
		if (bc.Type() == EFreeSlip){
			TPZManVector<REAL,15> Flux(5);
			fGradientFlux.ComputeFlux(dataleft.sol[0],dataleft.sol[0],data.normal,Flux);
			for(int i = 0; i < 15; i++) ef(i+5,0)    += +1. * weight*Flux[i];
		}//if FreeSlip
	}
}

REAL TPZEulerEquation::Pressure(TPZVec<REAL> &U, double gamma){
	
	if(U[0] < REAL(1e-10)){
		PZError << "TPZEulerEquation::Pressure - Negative or too small density "
		<< U[0] << std::endl;
		DebugStop();
	}
	
	REAL press = 0.0;
	
	//U = (U0,U1,U2,U3,U4) = (ro , ro u , ro v , ro w , ro e)
	REAL rho_velocity = ( U[1]*U[1] + U[2]*U[2] + U[3]*U[3] )/U[0];
	press = ((gamma-1.)*( U[4] - 0.5 * rho_velocity ));
	
	if(press < 0){
		REAL temp = (gamma-1.)*U[4];
		PZError << "TPZEulerEquation::Pressure Negative pressure: " << press << " (gama-1)*E = " << temp << std::endl;
		DebugStop();
	}
	return press;
	
}//method


REAL TPZEulerEquation::cSpeed(TPZVec<REAL> & sol){
	
	if(sol[0] < REAL(1e-10)){
		PZError << "TPZEulerEquation(::cSpeed Too small or negative density " << sol[0] << std::endl;
		DebugStop();
	}
	
	const REAL press = this->Pressure(sol,gGamma);
	const REAL temp = gGamma * press;
	
	if(temp < REAL(1e-10)) // too low or negative
	{
		PZError << "TPZEulerEquation::cSpeed Too low or negative numerator\n";
	}
	const REAL c = sqrt(gGamma * press/ sol[0]);
	return c;
	
}//method

REAL TPZEulerEquation::uRes(TPZVec<REAL> & sol){
	const REAL temp = sol[1]*sol[1] + sol[2]*sol[2] + sol[3]*sol[3];
	if(temp < REAL(1e-40)){
		PZError << "TPZEulerEquation::uRes Zero Velocity\n";
		DebugStop();
	}
	const REAL us = sqrt(temp)/sol[0];
	return us;
}

void TPZEulerEquation::ComputeEulerFlux(TPZVec<REAL> &sol, TPZFMatrix<REAL> & F){
	const double rho = sol[0];
	const double rhoU = sol[1];
	const double rhoV = sol[2];
	const double rhoW = sol[3];
	const double rhoE = sol[4];
	const double u = rhoU/rho;
	const double v = rhoV/rho;
	const double w = rhoW/rho;
	const double p = this->Pressure(sol,gGamma);
	F.Resize(5,3);
	F(0,0) = rhoU;        F(0,1) = rhoV;        F(0,2) = rhoW;
	F(1,0) = rhoU*u+p;    F(1,1) = rhoU*v;      F(1,2) = rhoU*w;
	F(2,0) = rhoV*u;      F(2,1) = rhoV*v+p;    F(2,2) = rhoV*w;
	F(3,0) = rhoW*u;      F(3,1) = rhoW*v;      F(3,2) = rhoW*w+p;
	F(4,0) = (rhoE+p)*u;  F(4,1) = (rhoE+p)*v;  F(4,2) = (rhoE+p)*w;
}//void
