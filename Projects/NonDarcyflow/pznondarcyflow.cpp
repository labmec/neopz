/*
 *  pznondarcyflow.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "pznondarcyflow.h"
#include "pzbndcond.h"
#include "pzaxestools.h"

TPZNonDarcyFlow::TPZNonDarcyFlow() : TPZMaterial()
{
	
}

TPZNonDarcyFlow::TPZNonDarcyFlow(int matid) : TPZMaterial(matid)
{
	
}


TPZNonDarcyFlow::TPZNonDarcyFlow(const TPZNonDarcyFlow &mat) : TPZMaterial(mat)
{
	
}

TPZNonDarcyFlow::~TPZNonDarcyFlow()
{
	
}

void TPZNonDarcyFlow::FillDataRequirements(TPZMaterialData &data)
{
	data.SetAllRequirements(false);
	data.fNeedsSol = true;
	data.fNeedsNormal = true;
}

void TPZNonDarcyFlow::FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
{
	data.SetAllRequirements(false);
	data.fNeedsSol = true;
	data.fNeedsNormal = true;
}

void TPZNonDarcyFlow::Print(std::ostream &out) {
	out << "\tBase class print:\n";
	TPZMaterial::Print(out);
	out << "name of material : " << this->Name() << "\n";
}

int TPZNonDarcyFlow::VariableIndex(const std::string &name) {
	if (!strcmp("Pressure", name.c_str())) return 0;
	if (!strcmp("Velocity", name.c_str())) return 1; // in SI m2/s
	return TPZMaterial::VariableIndex(name);
}

int TPZNonDarcyFlow::NSolutionVariables(int var) {
	switch(var) {
		case 0:
			return 1;
		case 1:
			return 3;
		default:
			return TPZMaterial::NSolutionVariables(var);
	}
}

void TPZNonDarcyFlow::Solution(TPZMaterialData &data, int var,
								TPZVec<REAL> &Solout) {
	
	TPZFMatrix<STATE> &dsol = data.dsol[0];
	
	REAL p = data.sol[0][0];
	REAL mu,rho;
	this->ReturnVisc(p,mu);
	this->ReturnRho(p,rho);	
	REAL kh = globData.fKh;
	REAL kv = globData.fKv;
	
	switch(var) {
		case 0:
		{
			Solout[0] = p;	
		}
			break;
		case 1:
		{
			Solout[0] = - (kh * rho) / mu * dsol(0,0) * 2. * M_PI * data.x[0]; // flux in 1r;  // in SI m2/s
			Solout[1] = - (kv * rho) / mu * dsol(1,0) * 2. * M_PI * data.x[0]; // flux in 1z;  // in SI m2/s
		}
			break;
		default:
			TPZMaterial::Solution(data,var,Solout);
	}
}

void TPZNonDarcyFlow::Contribute(TPZMaterialData &data, REAL weight,
								  TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	
	TPZFMatrix<STATE> &phi = data.phi;
	TPZFMatrix<STATE> &dphia = data.dphix;
	TPZManVector<REAL,3> &x = data.x;
	const int nphi = phi.Rows();
	const REAL r = x[0];
	const REAL pir2 = 2.*M_PI*r;
	
	TPZFMatrix<STATE> &axes = data.axes;
	TPZFMatrix<STATE> dphi;
	TPZAxesTools<STATE>::Axes2XYZ(dphia,dphi,axes);
	
	
	if (globData.fState == ELastState){
		REAL dpdr = data.dsol[0](0,0);
		REAL rho, por;
		
		REAL p = data.sol[0][0];
		this->ReturnPor(p,por);
		this->ReturnRho(p,rho);
		REAL dt = globData.fDeltaT;
		
		for (int i = 0 ; i < nphi ; i++){
			REAL res = - rho / dt * phi(i,0) * por * pir2;
			ef(i,0) += - 1. * res * weight;
		}
	}
	else if (globData.fState == ECurrentState){
		

		VarFad pf(data.sol[0][0],0);
		VarFad dpdrf(data.dsol[0](0,0),1);
		VarFad dpdzf(data.dsol[0](1,0),2);
		VarFad rhof, porf, muf;
		REAL kh,kv;
		
		this->ReturnPor(pf,porf);
		this->ReturnRho(pf,rhof);
		this->ReturnVisc(pf,muf);	
		kh = globData.fKh;
		kv = globData.fKv;
		REAL dt = globData.fDeltaT;
		
		//REAL g = globData.Gravity();
		for (int i = 0 ; i < nphi ; i++){
			VarFad resfad = (kh/muf * rhof * dpdrf * dphi(0,i) * pir2) + (kv/muf * rhof * (dpdzf /*+ rhof*g*/) * dphi(1,i) * pir2);
			resfad += rhof / dt * phi(i,0) * porf * pir2;
			REAL res = -1.*weight*resfad.val();
			ef(i,0) += res;
			
			const REAL dRdp = resfad.dx(0);
			const REAL dRdpdr = resfad.dx(1);
			const REAL dRdpdz = resfad.dx(2);
			
			for (int j = 0 ; j < nphi ; j++){
				REAL fadStiff = dRdp * phi(j,0) + dRdpdr * dphi(0,j)+ dRdpdz * dphi(1,j);
				ek(i,j) += fadStiff * weight;
			}
		}
	}
	else {
		DebugStop(); // why no state?????
	}

}

void TPZNonDarcyFlow::ContributeBC(TPZMaterialData &data, REAL weight,
									TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) {
	
	TPZFMatrix<STATE> &phi = data.phi;
	const int nphi = phi.Rows();
	
	if(bc.Val2().Rows() != 1 || bc.Val2().Cols() != 1){
		PZError << "Val2 must be of size (1,1)!\n";
		DebugStop();
	}
	
	if (globData.fState == ELastState) return;
	
	REAL v2 = bc.Val2()(0,0);
	if (bc.HasForcingFunction()){
		TPZManVector <REAL,1> Pbc(1,0.);
		bc.ForcingFunction()->Execute(data.x,Pbc); // here, data.x is useless
		v2 = Pbc[0];
	}
	
	const REAL p = data.sol[0][0];
	const REAL DdirValue = p - v2;
	
	switch(bc.Type()){
		case 0: //Dirichlet
			for (int i = 0 ; i < nphi ; i++){
				ef(i,0) += -1.*phi(i,0)*gBigNumber*DdirValue*weight;
				for (int j = 0 ; j < nphi ; j++){
					ek(i,j) += phi(i,0)*phi(j,0)*gBigNumber*weight;
				}
			}
			break;
		case 1: // Neumann - vazao massica
			for (int i = 0 ; i < nphi ; i++){
				ef(i,0) += -1.*phi(i,0)*weight*v2;
			}
			break;
		default:
			PZError << __PRETTY_FUNCTION__ << "bc type not implemented\n";
			DebugStop();
	}
	
}

int TPZNonDarcyFlow::ClassId() const{
	DebugStop();
	return -6378;
}

// -------------------------------------------------------------------------------------------

void TPZNonDarcyFlow::Write(TPZStream &buf, int withclassid) const{
	PZError << "Method Not Implemented!!\n";
	DebugStop();
	
}

// -------------------------------------------------------------------------------------------

void TPZNonDarcyFlow::Read(TPZStream &buf, void *context) {
	PZError << "Method Not Implemented!!\n";
	DebugStop();
	
}

// --
template<class T>
void TPZNonDarcyFlow::ReturnVisc(T &p, T &visc)
{
	visc = globData.fViscRef*exp(globData.fCVisc*(p - globData.fPref)); 
}

template<class T>
void TPZNonDarcyFlow::ReturnRho(T &p, T &rho)
{
	rho = globData.fRhoRef*exp(globData.fCRho*(p - globData.fPref)); 
}

template<class T>
void TPZNonDarcyFlow::ReturnPor(T &p, T &por)
{
	por = globData.fPorRef*exp(globData.fCPor*(p - globData.fPref)); 
}