/*
 *  ElasticMatInterface2D.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/23/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "PoroElasticMatInterface2D.h"
#include "pzdiscgal.h"
#include "pzelasmat.h"
#include "pzporoelastic2d.h"
#include "pzlog.h"
#include "pzaxestools.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.poroelastic2d"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
#endif

PoroElasticMatInterface2D::PoroElasticMatInterface2D() : TPZPoroElastic2d(){
	fknu = 1000000.0;
	fktu = 1000000.0;
	fknp = 1000000.0;
	fktp = 1000000.0;
	fcontribute = false;
	this->SetDimension(1);
}

PoroElasticMatInterface2D::PoroElasticMatInterface2D(int mat,int dim, bool DoContribute) : TPZPoroElastic2d(mat,dim){
	fknu = 1000000.0;
	fktu = 1000000.0;
	fknp = 1000000.0;
	fktp = 1000000.0;
	fcontribute = DoContribute;
	this->SetDimension(1);	
}

PoroElasticMatInterface2D::~PoroElasticMatInterface2D()
	{
	}

// Contribute Interface implementations

void PoroElasticMatInterface2D::SetPenalty(REAL knu, REAL ktu, REAL knp, REAL ktp)
{
	fknu = knu;
	fktu = ktu;
	fknp = knp;
	fktp = ktp;	
}

void PoroElasticMatInterface2D::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, 
													REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef)
{
	if (this->fcontribute) 
	{
		//	Definition of penalty constants note: this constans for nolinear analysis are funtion of normal and tangencial forces.
		REAL knu = this->fknu;
		REAL ktu = this->fktu;
		int sidn = ek.Rows();
		
		TPZFMatrix<REAL> &phiL = dataleftvec[0].phi;
		TPZFMatrix<REAL> &phiR = datarightvec[0].phi;
		TPZFMatrix<REAL>	&phipL	=	dataleftvec[1].phi;
		TPZFMatrix<REAL>	&dphipL	=	dataleftvec[1].dphix;
		TPZFMatrix<REAL>	&phipR	=	datarightvec[1].phi;
		TPZFMatrix<REAL>	&dphipR	=	datarightvec[1].dphix;	
		TPZFMatrix<REAL>	du(2,2);
		int phrpL = phipL.Rows();	
		int phrpR = phipR.Rows();	
		
		int &LeftPOrder=dataleftvec[0].p;
		int &RightPOrder=datarightvec[0].p;
		
		REAL &faceSize=data.HSize;
		
		
		int nrowl = phiL.Rows();
		int nrowr = phiR.Rows();
		int il,jl,ir,jr,id;	
		
		
		// For elastic part
		for (int in = 0; in < 2*(nrowl + nrowr) + ( phrpL + phrpR) ; in++) {
			ef(in,0) = 0.0;
		}
		REAL t[2] = {-data.normal[1],data.normal[0]};
		REAL n[2] = {data.normal[0],data.normal[1]};
		TPZFNMatrix<4,REAL> nx(2,2),tx(2,2);
		for (int i=0; i<2; i++) {
			for (int j=0; j<2 ; j++) {
				nx(i,j) = n[i]*n[j];
				tx(i,j) = t[i]*t[j];
			}
		}
		
		// Left Left contribution
		for(il=0; il < nrowl; il++) {
			for(jl=0; jl < nrowl; jl++) {
				ek(2*il,2*jl) += weight * (+ knu *phiL(il)*phiL(jl)*nx(0,0) + ktu*phiL(il)*phiL(jl)*tx(0,0));
				ek(2*il+1,2*jl+1) += weight * (+ knu *phiL(il)*phiL(jl)*nx(1,1) + ktu*phiL(il)*phiL(jl)*tx(1,1));
				ek(2*il+1,2*jl) += weight * (+ knu *phiL(il)*phiL(jl)*nx(1,0) + ktu*phiL(il)*phiL(jl)*tx(1,0));
				ek(2*il,2*jl+1) += weight * (+ knu *phiL(il)*phiL(jl)*nx(0,1) + ktu*phiL(il)*phiL(jl)*tx(0,1));
			}
		}
		// Right Left contribution	
		for(ir=0; ir < nrowr; ir++) {
			for(jl=0; jl < nrowl; jl++) {
				ek(2*ir+2*(nrowl)+phrpL,2*jl) += weight * (- knu*phiR(ir)*phiL(jl)*nx(0,0) - ktu*phiR(ir)*phiL(jl)*tx(0,0));
				ek(2*ir+1+2*(nrowl)+phrpL,2*jl+1) += weight * (- knu*phiR(ir)*phiL(jl)*nx(1,1) - ktu*phiR(ir)*phiL(jl)*tx(1,1));
				ek(2*ir+1+2*(nrowl)+phrpL,2*jl) += weight * (- knu*phiR(ir)*phiL(jl)*nx(1,0) - ktu*phiR(ir)*phiL(jl)*tx(1,0));
				ek(2*ir+2*(nrowl)+phrpL,2*jl+1) += weight * (- knu*phiR(ir)*phiL(jl)*nx(0,1) - ktu*phiR(ir)*phiL(jl)*tx(0,1));			
			}
		}
		
		// Left Right contribution		
		for(il=0; il < nrowl; il++) {
			for(jr=0; jr < nrowr; jr++) {
				ek(2*il,2*jr+2*(nrowl)+phrpL) += weight * (- knu*phiR(jr)*phiL(il)*nx(0,0) - ktu*phiR(jr)*phiL(il)*tx(0,0));
				ek(2*il+1,2*jr+1+2*(nrowl)+phrpL) += weight * (- knu*phiR(jr)*phiL(il)*nx(1,1) - ktu*phiR(jr)*phiL(il)*tx(1,1));
				ek(2*il+1,2*jr+2*(nrowl)+phrpL) += weight * (- knu*phiR(jr)*phiL(il)*nx(1,0) - ktu*phiR(jr)*phiL(il)*tx(1,0));
				ek(2*il,2*jr+1+2*(nrowl)+phrpL) += weight * (- knu*phiR(jr)*phiL(il)*nx(0,1) - ktu*phiR(jr)*phiL(il)*tx(0,1));				
			}
		}
		
		// Right Right contribution		
		for(ir=0; ir < nrowr; ir++) {
			for(jr=0; jr < nrowr; jr++) {
				ek(2*ir+2*(nrowl)+phrpL,2*jr+2*(nrowl)+phrpL) += weight * (+ knu *phiR(ir)*phiR(jr)*nx(0,0) + ktu*phiR(ir)*phiR(jr)*tx(0,0));
				ek(2*ir+1+2*(nrowl)+phrpL,2*jr+1+2*(nrowl)+phrpL) += weight * (+ knu *phiR(ir)*phiR(jr)*nx(1,1) + ktu*phiR(ir)*phiR(jr)*tx(1,1));
				ek(2*ir+1+2*(nrowl)+phrpL,2*jr+2*(nrowl)+phrpL) += weight * (+ knu *phiR(ir)*phiR(jr)*nx(1,0) + ktu*phiR(ir)*phiR(jr)*tx(1,0));
				ek(2*ir+2*(nrowl)+phrpL,2*jr+1+2*(nrowl)+phrpL) += weight * (+ knu *phiR(ir)*phiR(jr)*nx(0,1) + ktu*phiR(ir)*phiR(jr)*tx(0,1));	
			}
		}	
	}

}


void PoroElasticMatInterface2D::Solution(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleftvec, TPZVec<TPZMaterialData> &datarightvec, int var, TPZVec<REAL> &Solout)
{
	
	// Loading left and right solution
	
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZManVector<REAL,3> SolULeft, SolPLeft;
	TPZFNMatrix <6> DSolULeft, DSolPLeft;
	TPZFNMatrix <9> axesULeft, axesPLeft;
	
	TPZManVector<REAL,3> SolURight, SolPRight;	
	TPZFNMatrix <6> DSolURight, DSolPRight;
	TPZFNMatrix <9> axesURight, axesPRight;

	REAL t[2] = {-data.normal[1],data.normal[0]};
	REAL n[2] = {data.normal[0],data.normal[1]};	
	
	int StateVar = dataleftvec.size();
	
	for(int istate = 0 ; istate < StateVar ; istate++)
	{
		if (dataleftvec[istate].sol.size() != 1 && datarightvec[istate].sol.size() ) 
		{	
			std::cout << "Data no initialized on dataleftvec or datarightvec " << std::endl;
			DebugStop();
		}
	}

	SolULeft = dataleftvec[0].sol[0];
	DSolULeft = dataleftvec[0].dsol[0];
	axesULeft = dataleftvec[0].axes;
	SolPLeft = dataleftvec[0].sol[0];
	DSolPLeft = dataleftvec[0].dsol[0];
	axesPLeft = dataleftvec[0].axes;
	
	SolURight = datarightvec[0].sol[0];
	DSolURight = datarightvec[0].dsol[0];
	axesURight = datarightvec[0].axes;
	SolPRight = datarightvec[0].sol[0];
	DSolPRight = datarightvec[0].dsol[0];
	axesPRight = datarightvec[0].axes;
		
	REAL epsx;
	REAL epsy;
	REAL epsxy;
	REAL SigX;
	REAL SigY;
	REAL SigZ;
	REAL Tau, DSolxy[2][2];
	REAL DivU;
	REAL TMass;

	
	DSolxy[0][0] = DSolULeft(0,0)*axesULeft(0,0)+DSolULeft(1,0)*axesULeft(1,0); // dUx/dx
	DSolxy[1][0] = DSolULeft(0,0)*axesULeft(0,1)+DSolULeft(1,0)*axesULeft(1,1); // dUx/dy
	
	DSolxy[0][1] = DSolULeft(0,1)*axesULeft(0,0)+DSolULeft(1,1)*axesULeft(1,0); // dUy/dx
	DSolxy[1][1] = DSolULeft(0,1)*axesULeft(0,1)+DSolULeft(1,1)*axesULeft(1,1); // dUy/dy
	
	DivU = DSolxy[0][0]+DSolxy[1][1]+0.0;	
	
	epsx = DSolxy[0][0];// du/dx
	epsy = DSolxy[1][1];// dv/dy
	epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);

//	SigX = -((flambda + 2*fmu)*(epsx) + (flambda)*epsy - falpha*SolP[0]);
//	SigY = -((flambda + 2*fmu)*(epsy) + (flambda)*epsx - falpha*SolP[0]);	
//	Tau = 2.0*fmu*epsxy;		
//	Solout[0] = (SigX/SolP[0])
	
	std::cout << "Pressure Left : " << SolPLeft[0] << std::endl;
	std::cout << "Pressure Right : " << SolPRight[0] << std::endl;	

	
	
}


void PoroElasticMatInterface2D::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();
}

void PoroElasticMatInterface2D::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();	
}

void PoroElasticMatInterface2D::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &left, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	PZError << "\nFATAL ERROR - Method not implemented: " << __PRETTY_FUNCTION__ << "\n";
	DebugStop();	
}