/*
 *  LaplacianosDesacoplados.cpp
 *  PZ
 *
 *  Created by Agnaldo on 10/18/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include "poissondesacoplados.h"
#include "pzbndcond.h"

using namespace std;

TwoUncoupledPoisson::TwoUncoupledPoisson():TPZDiscontinuousGalerkin(), fXf1(0.), fXf2(0.),fDim(1){
	fK1 = 1.;
	fK2 = 1.;
}

TwoUncoupledPoisson::~TwoUncoupledPoisson(){
}

int TwoUncoupledPoisson::NStateVariables() {
	return 2;
}

void TwoUncoupledPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){

	
	int nref =  datavec.size();
	if (nref != 2 ) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
	TPZFMatrix  &phiu =  datavec[0].phi;
	TPZFMatrix &dphiu = datavec[0].dphix;
	int phru = phiu.Rows();

	TPZFMatrix  &phip =  datavec[1].phi;
	TPZFMatrix &dphip = datavec[1].dphix;
	int phrp = phip.Rows();
	
	//Equacao de Poisson
	//primeiro equacao
	int kd, in, jn;
	for(in = 0; in < phru; in++ ) {		
		ef(in, 0) += weight*fXf1*phiu(in,0); 
								 
		for(jn = 0; jn < phru; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in,jn) += weight *fK1*dphiu(kd,in)*dphiu(kd,jn); 
			}
		}
	}
	
	//segundo equacao
	for(in = 0; in < phrp; in++) {
		ef(in+phru, 0) += weight*fXf2*phip(in,0); 
		
		for(jn = 0; jn < phrp; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in+phru, jn+phru) += weight *fK2*dphip(kd,in)*dphip(kd,jn); 
			}
		}
	}

}

void TwoUncoupledPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix &ek,
									   TPZFMatrix &ef,TPZBndCond &bc) {
	
	
	int nref =  datavec.size();
	if (nref != 2 ) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	if (bc.Type() > 1 ) {
		cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
		DebugStop();
	}
	
	TPZFMatrix  &phiu = datavec[0].phi;
	TPZFMatrix  &phip = datavec[1].phi;
		
	int phru = phiu.Rows();
	int phrp = phip.Rows();
	short in,jn;
	REAL v2[1];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			//primeira equacao
			for(in = 0 ; in < phru; in++) {
				ef(in,0) += gBigNumber * v2[0]*phiu(in,0)*weight;
				for (jn = 0 ; jn < phru; jn++) {
					ek(in,jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
				}
			}
			
			//segunda equacao
			for(in = 0 ; in < phrp; in++) {
				ef(in+phru,0) += gBigNumber * v2[1]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+phru,jn+phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			break;
			
		case 1 :			// Neumann condition
			//primeira equacao
			for(in = 0 ; in < phru; in++) {
				ef(in,0) += v2[0] * phiu(in,0) * weight;
			}
			
			//seguna equacao
			for(in = 0 ; in < phrp; in++) {
				ef(in+phru,0) += v2[1] * phip(in,0) * weight;
			}
			break;
	}
		
}
