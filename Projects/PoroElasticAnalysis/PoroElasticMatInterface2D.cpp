/*
 *  ElasticMatInterface2D.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/23/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "PoroElasticMatInterface2D.h"

#include "pzelasmat.h"
#include "pzporoelastic2d.h"
#include "pzlog.h"
#include "pzaxestools.h"
#include "pzcompel.h"
#include "TPZMultiphysicsInterfaceEl.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.poroelastic2d");
#endif

#ifdef PZ_LOG
static TPZLogger logdata("pz.material.poroelastic.data");
#endif

PoroElasticMatInterface2D::PoroElasticMatInterface2D() : TPZPoroElastic2d(){
	fknu = 1000000.0;
	fktu = 1000000.0;
	fknp = 1000000.0;
	fktp = 1000000.0;
	fcontribute = false;
	this->SetDimension(1);
}

PoroElasticMatInterface2D::PoroElasticMatInterface2D(int mat,int dim, bool DoContribute, REAL fmu) : TPZPoroElastic2d(mat,dim){
	fknu = 1000000.0;
	fktu = 1000000.0;
	fknp = 1000000.0;
	fktp = 1000000.0;
	fcontribute = DoContribute;
	fFrictionmu = fmu;
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

void PoroElasticMatInterface2D::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	if (this->fcontribute) 
	{
		//	Definition of penalty constants note: this constans for nolinear analysis are funtion of normal and tangencial forces.
		REAL knu = this->fknu;
		REAL ktu = this->fktu;
//		int sidn = ek.Rows();
		
		TPZFMatrix<REAL> &phiL = dataleftvec[0].phi;
		TPZFMatrix<REAL> &phiR = datarightvec[0].phi;
		TPZFMatrix<REAL>	&phipL	=	dataleftvec[1].phi;
//		TPZFMatrix<REAL>	&dphipL	=	dataleftvec[1].dphix;
		TPZFMatrix<REAL>	&phipR	=	datarightvec[1].phi;
//		TPZFMatrix<REAL>	&dphipR	=	datarightvec[1].dphix;
		TPZFMatrix<REAL>	du(2,2);
		int phrpL = phipL.Rows();	
		int phrpR = phipR.Rows();	
		
//		int &LeftPOrder=dataleftvec[0].p;
//		int &RightPOrder=datarightvec[0].p;
		
//		REAL &faceSize=data.HSize;
		
		
		int nrowl = phiL.Rows();
		int nrowr = phiR.Rows();
		int il,jl,ir,jr;
		
		
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

/** Returns the variable index associated with the name */
int PoroElasticMatInterface2D::VariableIndex(const std::string &name)
{
	//	Elasticity Variables
	if(!strcmp("NormalRight",name.c_str()))				return	1;
	if(!strcmp("NormalLeft",name.c_str()))				return	2;
	if(!strcmp("TangentialRight",name.c_str()))				return	3;
	if(!strcmp("TangentialLeft",name.c_str()))				return	4;	
	if(!strcmp("AverageNormal",name.c_str()))						return	5;
	if(!strcmp("AverageTangential",name.c_str()))						return	6;	
	if(!strcmp("FRFactor",name.c_str()))						return	7;
	if(!strcmp("DrivingStress",name.c_str()))						return	8;
	
	return TPZMaterial::VariableIndex(name);
}

int PoroElasticMatInterface2D::NSolutionVariables(int var){
	if(var == 1)	return 3;
	if(var == 2)	return 3;
	if(var == 3)	return 3;
	if(var == 4)	return 3;
	if(var == 5)	return 1;
	if(var == 6)	return 1;
	if(var == 7)	return 1;
	if(var == 8)	return 1;	
	return TPZMaterial::NSolutionVariables(var);
}

void PoroElasticMatInterface2D::Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right)
{
	
	Solout.Resize( this->NSolutionVariables(var));

	int NumberOfStateVar = dataleftvec.size();	
	
    TPZVec<TPZMaterialData> dataleft(NumberOfStateVar),dataright(NumberOfStateVar);
	for(int istate = 0 ; istate < NumberOfStateVar ; istate++)
	{
		if (dataleftvec[istate].sol.size() != 1 && datarightvec[istate].sol.size() ) 
		{	
			std::cout << "Data not initialized on dataleftvec or datarightvec " << std::endl;
			DebugStop();
		}
        dataleft[istate] = dataleftvec[istate];
        dataright[istate] = datarightvec[istate];
	}
	
	//	Using the Solution method
    TPZPoroElastic2d * LeftPoroelastic = dynamic_cast<TPZPoroElastic2d *> (Left->Material());
    TPZPoroElastic2d * RightPoroelastic = dynamic_cast<TPZPoroElastic2d *>(Right->Material());
	
	
	TPZVec<STATE> LeftSigmaX,LeftSigmaY;
	TPZVec<STATE> RightSigmaX,RightSigmaY;
	TPZVec<STATE> LeftTau,RightTau;	
	
	LeftPoroelastic->Solution(dataleft,3,LeftSigmaX);
	RightPoroelastic->Solution(dataright,3,RightSigmaX);
	LeftPoroelastic->Solution(dataleft,4,LeftSigmaY);
	RightPoroelastic->Solution(dataright,4,RightSigmaY);
	LeftPoroelastic->Solution(dataleft,5,LeftTau);
	RightPoroelastic->Solution(dataright,5,RightTau);		

	REAL n[2] = {-data.normal[1],data.normal[0]};
	REAL t[2] = {-data.normal[0],-data.normal[1]};
		
	REAL NormalOnRight	= 0.0;
	REAL NormalOnLeft	= 0.0;
	REAL TangentialRight	= 0.0;
	REAL TangentialLeft	= 0.0;	
	REAL AverageNormal	= 0.0;
	REAL AverageTangential	= 0.0;	
//	REAL StressMagnitude = 0.0;
	REAL Stresstrx = 0.0;
	REAL Stresstry = 0.0;
	REAL Stresstlx = 0.0;
	REAL Stresstly = 0.0;
	REAL FailureStress = 0.0;
	REAL DrivingStress = 0.0;	
	
	// Left Stress vector on plane with normal n and tangent t
	Stresstlx = LeftSigmaX[0]*n[0] + LeftTau[0]*n[1];
	Stresstly = LeftTau[0]*n[0] + LeftSigmaY[0]*n[1];	
	REAL StressVectorL[2] = {Stresstlx,Stresstly};
	
	// Right Stress vector on plane with normal n and tangent t	
	Stresstrx = RightSigmaX[0]*n[0] + RightTau[0]*n[1];
	Stresstry = RightTau[0]*n[0] + RightSigmaY[0]*n[1];	
	REAL StressVectorR[2] = {Stresstrx,Stresstry};	
	
	// Calculation of normal of left part

	NormalOnLeft = StressVectorL[0]*n[0]+StressVectorL[1]*n[1];
	TangentialLeft = StressVectorL[0]*t[0]+StressVectorL[1]*t[1]; 	
	
	// Calculation of normal of Right part	
	
	NormalOnRight = StressVectorR[0]*n[0]+StressVectorR[1]*n[1];
	TangentialRight = StressVectorR[0]*t[0]+StressVectorR[1]*t[1]; 		
	

	
//	StressMagnitude = (StressVectorMagnitudeL+StressVectorMagnitudeR)/2.0;
	AverageNormal = (NormalOnLeft+NormalOnRight)/2.0;
	AverageTangential = (TangentialLeft+TangentialRight)/2.0;
		
	DrivingStress = AverageTangential - this->fFrictionmu*AverageNormal;
	//	Normal Right
	if(var == 1){
		Solout[0] = NormalOnRight*n[0];		
		Solout[1] = NormalOnRight*n[1];
		Solout[2] = 0.0;		
		return;
	}
	
	//	Normal Left
	if(var == 2){
		Solout[0] = -1.0*NormalOnLeft*n[0];		
		Solout[1] = -1.0*NormalOnLeft*n[1];
		Solout[2] = 0.0;		
		return;
	}
	
	//	Tangential Right
	if(var == 3){
		Solout[0] = TangentialRight*t[0];		
		Solout[1] = TangentialRight*t[1];
		Solout[2] = 0.0;
		return;
	}
	
	//	Tangential Left
	if(var == 4){
		Solout[0] = -1.0*TangentialLeft*t[0];		
		Solout[1] = -1.0*TangentialLeft*t[1];
		Solout[2] = 0.0;
		return;
	}	
	
	//	Normal Average
	if(var == 5){
		Solout[0] = AverageNormal;
		return;
	}
	
	//	Tangential Average
	if(var == 6){
		Solout[0] = AverageTangential;
		return;
	}
	
	//	Fault Reactivation Factor
	if(var == 7){
		Solout[0] = FailureStress;
		return;
	}	
	
	//	DrivingStress Average
	if(var == 8){
		Solout[0] = DrivingStress;
		return;
	}	
	
//	std::stringstream filetemp;
//	filetemp << "Output/Fault" << 1 << ".txt";
//	std::string FaultFile = filetemp.str();	
//	std::ofstream out("Output/Fault.txt");
//	TPZFMatrix<REAL> Data(10,5,0.0);
//	
//	Data.Print("FaultDataSolution", out, EMathematicaInput);
	
	
	
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

int PoroElasticMatInterface2D::ClassId() const {
    return Hash("PoroElasticMatInterface2D") ^ TPZPoroElastic2d::ClassId() << 1;
}
