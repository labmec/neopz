/*
 *  pzviscoelastic.cpp
 *  pos_processamento
 *
 *  Created by Pamela Diaz on 7/16/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */


#include "pzviscoelastic.h"

TPZViscoelastic::TPZViscoelastic(TPZMatWithMem<TPZFMatrix<REAL>, TPZElasticity3D> &matwithmem,int id,REAL lambdaE,REAL muE, REAL lambdaV, REAL muV, REAL alphaT):TPZMatWithMem<TPZFMatrix<REAL>, TPZElasticity3D>(matwithmem), flambdaE(lambdaE),fmuE(muE),flambdaV(lambdaV),fmuV(muV),falphaT(alphaT)
{

	REAL lambda = flambdaE-(falphaT*flambdaV)/(1+falphaT);
	REAL mu = fmuE -(falphaT*fmuV)/(1+falphaT);
	
	std::cout<<"lambda"<<std::endl;
	std::cout<<lambda<<std::endl;
	std::cout<<"mu"<<std::endl;
	std::cout<<mu<<std::endl;
	
	fE=mu*(3*lambda+2*mu)/(lambda+mu);
	fPoisson= lambda/(2*(lambda+mu));
	
	std::cout<<"elasticity modulus"<<std::endl;
	std::cout<<fE<<std::endl;
	std::cout<<"poisson"<<std::endl;
	std::cout<<fPoisson<<std::endl;
	TPZMatWithMem<TPZFMatrix<REAL>,TPZElasticity3D>::Print(std::cout);
	
	
}

void TPZViscoelastic::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef)
{
	
	if(fUpdateMem== 0){
		
		TPZMatWithMem<TPZFMatrix<REAL>,TPZElasticity3D>::Contribute(data,weight,ek,ef);
		TPZFMatrix<REAL> dphi = data.dphix;
		TPZFMatrix<REAL> phi = data.phi;
		
		const int phr = phi.Rows();
		//this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
		//{dvdy*dudx, dvdy*dudy, dvdy*dudz},
		//{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
		TPZFNMatrix<9> Deriv(3,3);
		
		int index = data.intGlobPtIndex;
		//int index = 0;    
		TPZFNMatrix<6>  qsi;//
	    qsi = *MemItem(index);
	    int nstate = NStateVariables();
	    int in;
		REAL val;
		//qsi.Print("Valores de qsi ",std::cout);
		
		
		for(in = 0; in < phr; in++) 
		{ 
			//in: test function index	
			val = 0.;
			val -= qsi(_XX_,0) * dphi(0,in); // |
			val -= qsi(_XY_,0) * dphi(1,in); // fk
			val -= qsi(_XZ_,0) * dphi(2,in); // |
			val = val/(1+falphaT);
			ef(in*nstate+0,0) += weight * val;
			
			//Second equation: fb and fk
			val = 0.;
			val -= qsi(_XY_,0) * dphi(0,in); // |
			val -= qsi(_YY_,0) * dphi(1,in); // fk
			val -= qsi(_YZ_,0) * dphi(2,in); // |
			val = val/(1+falphaT);
			ef(in*nstate+1,0) += weight * val;
			
			//third equation: fb and fk
			val = 0.;
			val -= qsi(_XZ_,0) * dphi(0,in); // |
			val -= qsi(_YZ_,0) * dphi(1,in); // fk
			val -= qsi(_ZZ_,0) * dphi(2,in); // |
			val = val/(1+falphaT);
			ef(in*nstate+2,0) += weight * val;
		}
	}
	else{
		
		TPZFMatrix<REAL> dphi = data.dphix;
		TPZFMatrix<REAL> phi = data.phi;
		
		//this matrix will store {{dvdx*dudx, dvdx*dudy, dvdx*dudz},
		//{dvdy*dudx, dvdy*dudy, dvdy*dudz},
		//{dvdz*dudx, dvdz*dudy, dvdz*dudz}}
		TPZFNMatrix<9> Deriv(3,3);
		
		int index = data.intGlobPtIndex;
		TPZFNMatrix<6>  qsi;
		TPZFNMatrix<6>  Strain(6,1);
		TPZFNMatrix<6>  qsin1(6,1);
		
		
	    qsi = *MemItem(index);
		
		TPZFNMatrix<9> DSolXYZ(3,3,0.);
		DSolXYZ = data.dsol[0];
		//data.axes.Multiply(data.dsol,DSolXYZ,1/*transpose*/);
		
		Strain.Redim(6,1);
		Strain(_XX_,0) = DSolXYZ(0,0);
		Strain(_YY_,0) = DSolXYZ(1,1);
		Strain(_ZZ_,0) = DSolXYZ(2,2);
		Strain(_XY_,0) = 0.5 * ( DSolXYZ(1,0) + DSolXYZ(0,1) );
		Strain(_XZ_,0) = 0.5 * ( DSolXYZ(2,0) + DSolXYZ(0,2) );
		Strain(_YZ_,0) = 0.5 * ( DSolXYZ(2,1) + DSolXYZ(1,2) );
		
		qsi = *MemItem(index);
		
		
		REAL tr;
		tr = Strain(_XX_,0)+Strain(_YY_,0)+Strain(_ZZ_,0);
		
		//REAL lambdaE,REAL muE, REAL lambdaV, REAL muV, REAL alphaT
		
		qsin1(_XX_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_XX_,0)*fmuV) + qsi(_XX_,0))/(1 + falphaT);
		qsin1(_YY_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_YY_,0)*fmuV) + qsi(_YY_,0))/(1 + falphaT);
		qsin1(_ZZ_,0) = (falphaT*(-(tr)*flambdaV - 2*Strain(_ZZ_,0)*fmuV) + qsi(_ZZ_,0))/(1 + falphaT);
		qsin1(_XY_,0) = (-2*falphaT*Strain(_XY_,0)*fmuV + qsi(_XY_,0))/(1 + falphaT);
		qsin1(_XZ_,0) = (-2*falphaT*Strain(_XZ_,0)*fmuV + qsi(_XZ_,0))/(1 + falphaT);
		qsin1(_YZ_,0) = (-2*falphaT*Strain(_YZ_,0)*fmuV + qsi(_YZ_,0))/(1 + falphaT);
		
		//qsin1.Print("qsin1",std::cout);
		MemItem(index) = qsin1;
		
		
	}
	
}

