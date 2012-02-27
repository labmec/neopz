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
#include "pzaxestools.h"

//#include "pzdiscgal.h"
//#include "pzmaterialdata.h"
//#include "pzmaterialid.h"

using namespace std;

TwoUncoupledPoisson::TwoUncoupledPoisson():TPZDiscontinuousGalerkin(), fXf1(0.), fXf2(0.),fDim(1){
	fK1 = 1.;
	fK2 = 1.;
}

TwoUncoupledPoisson::TwoUncoupledPoisson(int matid, int dim):TPZDiscontinuousGalerkin(matid), fXf1(0.), fXf2(0.),fDim(dim){
	fK1 = 1.;
	fK2 = 1.;
}

TwoUncoupledPoisson::~TwoUncoupledPoisson(){
}

int TwoUncoupledPoisson::NStateVariables() {
	return 1;
}

void TwoUncoupledPoisson::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Laplace operator multiplier fK  da primeira equacao  "<< fK1 << endl;
	out << "Laplace operator multiplier fK  da segunda equacao  "<< fK2<< endl;
	out << "Forcing vector fXf da primeira equacao  " << fXf1 << endl;
	out << "Forcing vector fXf da seguna equacao  " << fXf2 << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
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
				ek(in,jn) += -weight *fK1*dphiu(kd,in)*dphiu(kd,jn); 
			}
		}
	}
	
	//segundo equacao
	for(in = 0; in < phrp; in++) {
		ef(in+phru, 0) += weight*fXf2*phip(in,0); 
		
		for(jn = 0; jn < phrp; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in+phru, jn+phru) += -weight *fK2*dphip(kd,in)*dphip(kd,jn); 
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
	REAL v2[2];
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

//void TwoUncoupledPoisson::ContributeInterface(TPZVec<TPZMaterialData> &datavec, REAL weight,
//                                          TPZFMatrix &ek,TPZFMatrix &ef){
//	
//	TPZFMatrix &dphiLdAxesu = datavec[0].dphixl;
//	TPZFMatrix &dphiRdAxesu = datavec[0].dphixr;
//	TPZFMatrix &phiLu = datavec[0].phil;
//	TPZFMatrix &phiRu = datavec[0].phir;
//	TPZManVector<REAL,3> &normalu = datavec[0].normal;
//	
//			
//	TPZFNMatrix<660> dphiLu, dphiRu;
//	TPZAxesTools::Axes2XYZ(dphiLdAxesu, dphiLu, datavec[0].axesleft);
//	TPZAxesTools::Axes2XYZ(dphiRdAxesu, dphiRu, datavec[0].axesright);
//	
//	int &LeftPOrderu=datavec[0].leftp;
//	int &RightPOrderu=datavec[0].rightp;
//	
//	REAL &faceSizeu=datavec[0].HSize;
//	
//	
//	int nrowlu = phiLu.Rows();
//	int nrowru = phiRu.Rows();
//	int il,jl,ir,jr,id;
//	
//		
//	//diffusion term
//	REAL leftK, rightK;
//	leftK  = fK1;
//	rightK = fK1;
//	
//	// 1) phi_I_left, phi_J_left
//	for(il=0; il<nrowlu; il++) {
//		REAL dphiLinormal = 0.;
//		for(id=0; id<fDim; id++) {
//			dphiLinormal += dphiLu(id,il)*normalu[id];
//		}
//		for(jl=0; jl<nrowlu; jl++) {
//			REAL dphiLjnormal = 0.;
//			for(id=0; id<fDim; id++) {
//				dphiLjnormal += dphiLu(id,jl)*normalu[id];
//			}
//			ek(il,jl) += weight * leftK*(this->fSymmetry * 0.5*dphiLinormal*phiL(jl,0)-0.5*dphiLjnormal*phiL(il,0));
//		}
//	}
//	
//	// 2) phi_I_right, phi_J_right
//	for(ir=0; ir<nrowru; ir++) {
//		REAL dphiRinormal = 0.;
//		for(id=0; id<fDim; id++) {
//			dphiRinormal += dphiRu(id,ir)*normalu[id];
//		}
//		for(jr=0; jr<nrowr; jr++) {
//			REAL dphiRjnormal = 0.;
//			for(id=0; id<fDim; id++) {
//				dphiRjnormal += dphiRu(id,jr)*normalu[id];
//			}
//			ek(ir+nrowl,jr+nrowl) += weight * rightK * (
//														this->fSymmetry * (-0.5 * dphiRinormal * phiR(jr) ) + 0.5 * dphiRjnormal * phiR(ir)
//														);
//		}
//	}
//	
//	// 3) phi_I_left, phi_J_right
//	for(il=0; il<nrowl; il++) {
//		REAL dphiLinormal = 0.;
//		for(id=0; id<fDim; id++) {
//			dphiLinormal += dphiLu(id,il)*normalu[id];
//		}
//		for(jr=0; jr<nrowr; jr++) {
//			REAL dphiRjnormal = 0.;
//			for(id=0; id<fDim; id++) {
//				dphiRjnormal += dphiRu(id,jr)*normalu[id];
//			}
//			ek(il,jr+nrowl) += weight * (
//										 this->fSymmetry * (-0.5 * dphiLinormal * leftK * phiR(jr) ) - 0.5 * dphiRjnormal * rightK * phiL(il)
//										 );
//		}
//	}
//	
//	// 4) phi_I_right, phi_J_left
//	for(ir=0; ir<nrowr; ir++) {
//		REAL dphiRinormal = 0.;
//		for(id=0; id<fDim; id++) {
//			dphiRinormal += dphiRu(id,ir)*normalu[id];
//		}
//		for(jl=0; jl<nrowl; jl++) {
//			REAL dphiLjnormal = 0.;
//			for(id=0; id<fDim; id++) {
//				dphiLjnormal += dphiLu(id,jl)*normalu[id];
//			}
//			ek(ir+nrowl,jl) += weight * (
//										 this->fSymmetry * 0.5 * dphiRinormal * rightK * phiL(jl) + 0.5 * dphiLjnormal * leftK * phiR(ir)
//										 );
//		}
//	}
//	
//	
//	
//	if (this->IsSymetric()){
//		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
//	}
//}

/** Returns the variable index associated with the name */
int TwoUncoupledPoisson::VariableIndex(const std::string &name){
	if(!strcmp("SolutionU",name.c_str()))        return  1;
	if(!strcmp("SolutionP",name.c_str()))        return  2;
	if(!strcmp("DerivateU",name.c_str()))        return  3;
	if(!strcmp("DerivateP",name.c_str()))        return  4;
	
	return TPZMaterial::VariableIndex(name);
}

int TwoUncoupledPoisson::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 1;
	if((var == 3) || (var == 4)) return fDim;
	return TPZMaterial::NSolutionVariables(var);
}

void TwoUncoupledPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolU, SolP;
	TPZFMatrix DSolU, DSolP;
	TPZFMatrix axesU, axesP;
	
	SolU=datavec[0].sol[0];
	DSolU=datavec[0].dsol[0];
	axesU=datavec[0].axes;
	SolP=datavec[1].sol[0];
	DSolP=datavec[1].dsol[0];
	axesP=datavec[1].axes;
	
	if(var == 1){
		Solout[0] = SolU[0];//function (state variable u)
		return;
	}
	
	if(var == 2){
		Solout[0] = SolP[0];//function (state variable p)
		return;
	}
	
	if(var == 3) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9> dsoldx;
			TPZAxesTools::Axes2XYZ(DSolU, dsoldx, axesU);
			Solout[id] = dsoldx(id,0);//derivate of u
		}
		return;
	}//var == 3
	
	if(var == 4) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9> dsoldx;
			TPZAxesTools::Axes2XYZ(DSolP, dsoldx, axesP);
			Solout[id] = dsoldx(id,0);//derivate of p
		}
		return;
	}//var == 4
}


void TwoUncoupledPoisson::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
	DebugStop();
}

void TwoUncoupledPoisson::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
	DebugStop();
}
//int IntegrationRuleOrder(TPZVec<int> elPMaxOrder) const
//{
//	 
//}

