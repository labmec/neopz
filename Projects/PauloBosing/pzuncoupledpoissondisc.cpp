//
//  pzuncoupledpoissondisc.cpp
//  PZ
//
//  Created by Agnaldo Farias on 2/19/13.
//
//

#include "pzuncoupledpoissondisc.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"


using namespace std;

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc():TPZMaterial(), fXf1(0.), fXf2(0.),fDim(1){
	fK1 = 1.;
	fK2 = 1.;
}

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc(int matid, int dim):TPZMaterial(matid), fXf1(0.), fXf2(0.),fDim(dim){
	fK1 = 1.;
	fK2 = 1.;
}

TPZMatUncoupledPoissonDisc::~TPZMatUncoupledPoissonDisc(){
}

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc(const TPZMatUncoupledPoissonDisc &copy):TPZMaterial(copy){
    
    this->operator=(copy);
}

TPZMatUncoupledPoissonDisc & TPZMatUncoupledPoissonDisc::operator=(const TPZMatUncoupledPoissonDisc &copy){
    
    TPZMaterial::operator = (copy);
	fXf1  = copy.fXf1;
    fXf2  = copy.fXf2;
	fDim = copy.fDim;
	fK1  = copy.fK1;
    fK2  = copy.fK2;
    
	fSymmetry1 = copy.fSymmetry1;
    fSymmetry2 = copy.fSymmetry2;
    
	fPenaltyConstant1 = copy.fPenaltyConstant1;
    fPenaltyConstant2 = copy.fPenaltyConstant2;
    
	return *this;
}


void TPZMatUncoupledPoissonDisc::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Laplace operator multiplier fK  da primeira equacao  "<< fK1 << endl;
	out << "Laplace operator multiplier fK  da segunda equacao  "<< fK2<< endl;
	out << "Forcing vector fXf da primeira equacao  " << fXf1 << endl;
	out << "Forcing vector fXf da seguna equacao  " << fXf2 << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMatUncoupledPoissonDisc::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
	
	int nref =  datavec.size();
	if (nref != 2 ) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
    
	TPZFMatrix<REAL>  &phiu1 =  datavec[0].phi;
	TPZFMatrix<REAL> &dphiu1 = datavec[0].dphix;
	int phru1 = phiu1.Rows();
    
	TPZFMatrix<REAL>  &phiu2 =  datavec[1].phi;
	TPZFMatrix<REAL> &dphiu2 = datavec[1].dphix;
	int phru2 = phiu2.Rows();
	
    //----- nas duas equacoes ------
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		fXf1 = res[0];
        fXf2 = res[0];
	}
    
	//Equacao de Poisson
	// ------ primeira equacao ------
    
	int kd, in, jn;
	for(in = 0; in < phru1; in++ ) {
		ef(in, 0) += -1.*weight*fXf1*phiu1(in,0);
        
		for(jn = 0; jn < phru1; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in,jn) += weight*fK1*dphiu1(kd,in)*dphiu1(kd,jn);
			}
		}
	}
	
//	//----- so na segunda equacao ------
//    if(fForcingFunction) {
//		TPZManVector<STATE> res(1);
//		fForcingFunction->Execute(datavec[1].x,res);
//		fXf2 = res[0];
//	}
    
	for(in = 0; in < phru2; in++) {
		ef(in+phru1, 0) += -1.*weight*fXf2*phiu2(in,0);
		
		for(jn = 0; jn < phru2; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in+phru1, jn+phru1) += weight*fK2*dphiu2(kd,in)*dphiu2(kd,jn);
			}
		}
	}
    
}

void TPZMatUncoupledPoissonDisc::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,
                                              TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
	
	int nref =  datavec.size();
	if (nref != 2) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
	TPZFMatrix<REAL>  &phiu1 = datavec[0].phi;
	TPZFMatrix<REAL>  &phiu2 = datavec[1].phi;
    
	int phru1 = phiu1.Rows();
	int phru2 = phiu2.Rows();
	short in,jn;
	REAL v2[2];
	v2[0] = bc.Val2()(0,0); //condicao de contorno da primeira equacao
	v2[1] = bc.Val2()(1,0); //condicao de contorno da segunda equaca
	
	switch (bc.Type()) {
		case 0 : // Dirichlet condition
			
            //primeira equacao
			for(in = 0 ; in < phru1; in++) {
				ef(in,0) += gBigNumber * v2[0]*phiu1(in,0)*weight;
				for (jn = 0 ; jn < phru1; jn++) {
					ek(in,jn) += gBigNumber*phiu1(in,0)*phiu1(jn,0)*weight;
				}
			}
			
			//segunda equacao
			for(in = 0 ; in < phru2; in++) {
				ef(in+phru1,0) += gBigNumber * v2[1]*phiu2(in,0)*weight;
				for (jn = 0 ; jn < phru2; jn++) {
					ek(in+phru1,jn+phru1) += gBigNumber*phiu2(in,0)*phiu2(jn,0)*weight;
				}
			}
			break;
			
		case 1 : // Neumann condition
			//primeira equacao
			for(in = 0 ; in < phru1; in++) {
				ef(in,0) += v2[0]*phiu1(in,0)*weight;
			}
			
			//seguna equacao
			for(in = 0 ; in < phru2; in++) {
				ef(in+phru1,0) += v2[1]*phiu2(in,0)*weight;
			}
			break;
            
        case 10 : // Dirichlet-Neumann condition
            
			//primeira equacao
			for(in = 0 ; in < phru1; in++) {
				ef(in,0) += gBigNumber * v2[0]*phiu1(in,0)*weight;
				for (jn = 0 ; jn < phru1; jn++) {
					ek(in,jn) += gBigNumber*phiu1(in,0)*phiu1(jn,0)*weight;
				}
			}
			
			//segunda equacao
			for(in = 0 ; in < phru2; in++) {
				ef(in+phru1,0) += v2[1]*phiu2(in,0)*weight;
			}
			break;
            
        case 11 : // Neumann-Dirichlet condition
            
			//primeira equacao
			for(in = 0 ; in < phru1; in++) {
				ef(in,0) += v2[0]*phiu1(in,0)*weight;
			}
			
			//segunda equacao
			for(in = 0 ; in < phru2; in++) {
				ef(in+phru1,0) += gBigNumber*v2[1]*phiu2(in,0)*weight;
				for (jn = 0 ; jn < phru2; jn++) {
					ek(in+phru1,jn+phru1) += gBigNumber*phiu2(in,0)*phiu2(jn,0)*weight;
				}
			}
			break;
            
            
	}
    
}

void TPZMatUncoupledPoissonDisc::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

	int il,jl,ir,jr,id;
    
//========== Primeira equacao ==========    
{
    //dados da primeira variavel (u)
	TPZFMatrix<REAL> &dphiLdAxes_u1 = dataleft[0].dphix;
	TPZFMatrix<REAL> &dphiRdAxes_u1 = dataright[0].dphix;
	TPZFMatrix<REAL> &phiL_u1 = dataleft[0].phi;
	TPZFMatrix<REAL> &phiR_u1 = dataright[0].phi;
	TPZManVector<REAL,3> &normal_u1 = data.normal;//No PZ a normal eh do elemento de menor indice para o de maior indice
    
	TPZFNMatrix<660> dphiL_u1, dphiR_u1;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes_u1, dphiL_u1, dataleft[0].axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes_u1, dphiR_u1, dataright[0].axes);
    
	int &LeftPOrder_u1=dataleft[0].p;
	int &RightPOrder_u1=dataright[0].p;
    
	REAL &faceSize_u1=data.HSize;
    
	int nrowl_u1 = phiL_u1.Rows();
	int nrowr_u1 = phiR_u1.Rows();
    int first_right = nrowl_u1 + dataleft[1].phi.Rows();
    
    //diffusion term
    REAL leftK1, rightK1;
	leftK1  = fK1;
	rightK1 = fK1;
    
    //Contribuicao na matriz de rigidez
    //No PZ a normal eh do elemento de menor indice para o de maior indice
    //por isso tem-se sinais contrarios nas equacoes  (comparadas com a teoria) da contribuicao da matriz de rigidez

	// 1) phi_I_left, phi_J_left 
	for(il=0; il<nrowl_u1; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u1(id,il)*normal_u1[id];
		}
		for(jl=0; jl<nrowl_u1; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u1(id,jl)*normal_u1[id];
			}
			ek(il,jl) += weight*leftK1*(this->fSymmetry1*0.5*dphiLinormal*phiL_u1(jl,0)-0.5*dphiLjnormal*phiL_u1(il,0));
		}
	}

	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr_u1; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u1(id,ir)*normal_u1[id];
		}
		for(jr=0; jr<nrowr_u1; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u1(id,jr)*normal_u1[id];
			}
			ek(ir+first_right,jr+first_right) += weight*rightK1*(this->fSymmetry1*(-0.5*dphiRinormal*phiR_u1(jr)) + 0.5*dphiRjnormal*phiR_u1(ir));
		}
	}

	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl_u1; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u1(id,il)*normal_u1[id];
		}
		for(jr=0; jr<nrowr_u1; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u1(id,jr)*normal_u1[id];
			}
			ek(il,jr+first_right) += weight*leftK1*(this->fSymmetry1*(-0.5*dphiLinormal*phiR_u1(jr)) - 0.5*dphiRjnormal*phiL_u1(il));
		}
	}


    // 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr_u1; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u1(id,ir)*normal_u1[id];
		}
		for(jl=0; jl<nrowl_u1; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u1(id,jl)*normal_u1[id];
			}
			ek(ir+first_right,jl) += weight*rightK1*(this->fSymmetry1*0.5*dphiRinormal*phiL_u1(jl) + 0.5*dphiLjnormal*phiR_u1(ir));
		}
	}
    
    if(this->fPenaltyConstant1 !=0 ){
        
        REAL penalty = fPenaltyConstant1*(0.5*(abs(leftK1)*LeftPOrder_u1*LeftPOrder_u1 + abs(rightK1)*RightPOrder_u1*RightPOrder_u1))/faceSize_u1;
        
        // 1) left i / left j
		for(il=0; il<nrowl_u1; il++) {
			for(jl=0; jl<nrowl_u1; jl++) {
				ek(il,jl) += weight*penalty*phiL_u1(il,0)*phiL_u1(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr_u1; ir++) {
			for(jr=0; jr<nrowr_u1; jr++) {
				ek(ir+first_right,jr+first_right) += weight*penalty*phiR_u1(ir,0)*phiR_u1(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl_u1; il++) {
			for(jr=0; jr<nrowr_u1; jr++) {
				ek(il,jr+first_right) += -1.0*weight*penalty*phiR_u1(jr,0)*phiL_u1(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr_u1; ir++) {
			for(jl=0; jl<nrowl_u1; jl++) {
				ek(ir+first_right,jl) += -1.0*weight*penalty*phiL_u1(jl,0)*phiR_u1(ir,0);
			}
		}

    }
}
    
    
// =========== Segunda equacao =========
{
    //dados da segunda variavel (p)
    TPZFMatrix<REAL> &dphiLdAxes_u2 = dataleft[1].dphix;
	TPZFMatrix<REAL> &dphiRdAxes_u2 = dataright[1].dphix;
	TPZFMatrix<REAL> &phiL_u2 = dataleft[1].phi;
	TPZFMatrix<REAL> &phiR_u2 = dataright[1].phi;
	TPZManVector<REAL,3> &normal_u2 = data.normal;
    
	TPZFNMatrix<660> dphiL_u2, dphiR_u2;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes_u2, dphiL_u2, dataleft[1].axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes_u2, dphiR_u2, dataright[1].axes);
    
	int &LeftPOrder_u2=dataleft[1].p;
	int &RightPOrder_u2=dataright[1].p;
    
	REAL &faceSize_u2=data.HSize;
    
	int nrowl_u2 = phiL_u2.Rows();
	int nrowr_u2 = phiR_u2.Rows();
    
    int first_left = dataleft[0].phi.Rows();
    int first_right = dataleft[0].phi.Rows()+dataleft[1].phi.Rows()+dataright[0].phi.Rows();
    
    //diffusion term
	REAL leftK2, rightK2;
    leftK2  = fK2;
	rightK2 = fK2;
    
    // 1) phi_I_left, phi_J_left
	for(il=0; il<nrowl_u2; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u2(id,il)*normal_u2[id];
		}
		for(jl=0; jl<nrowl_u2; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u2(id,jl)*normal_u2[id];
			}
			ek(first_left + il,first_left + jl) += weight*leftK2*(this->fSymmetry2*0.5*dphiLinormal*phiL_u2(jl,0)-0.5*dphiLjnormal*phiL_u2(il,0));
		}
	}
    
	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr_u2; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u2(id,ir)*normal_u2[id];
		}
		for(jr=0; jr<nrowr_u2; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u2(id,jr)*normal_u2[id];
			}
			ek(first_right + ir, first_right + jr) += weight*rightK2*(this->fSymmetry2*(-0.5*dphiRinormal*phiR_u2(jr)) + 0.5*dphiRjnormal*phiR_u2(ir));
		}
	}
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl_u2; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u2(id,il)*normal_u2[id];
		}
		for(jr=0; jr<nrowr_u2; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u2(id,jr)*normal_u2[id];
			}
			ek(first_left + il, first_right+ jr) += weight*leftK2*(this->fSymmetry2*(-0.5*dphiLinormal*phiR_u2(jr)) - 0.5*dphiRjnormal*phiL_u2(il));
		}
	}
    
    
    // 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr_u2; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u2(id,ir)*normal_u2[id];
		}
		for(jl=0; jl<nrowl_u2; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u2(id,jl)*normal_u2[id];
			}
			ek(first_right + ir, first_left + jl) += weight*rightK2*(this->fSymmetry2*0.5*dphiRinormal*phiL_u2(jl) + 0.5*dphiLjnormal*phiR_u2(ir));
		}
	}
    
    if(this->fPenaltyConstant2 !=0 ){
        
        REAL penalty = fPenaltyConstant2*(0.5*(abs(leftK2)*LeftPOrder_u2*LeftPOrder_u2 + abs(rightK2)*RightPOrder_u2*RightPOrder_u2))/faceSize_u2;
        
        // 1) left i / left j
		for(il=0; il<nrowl_u2; il++) {
			for(jl=0; jl<nrowl_u2; jl++) {
				ek(first_left+il, first_left+jl) += weight*penalty*phiL_u2(il,0)*phiL_u2(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr_u2; ir++) {
			for(jr=0; jr<nrowr_u2; jr++) {
				ek(first_right+ir, first_right + jr) += weight*penalty*phiR_u2(ir,0)*phiR_u2(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl_u2; il++) {
			for(jr=0; jr<nrowr_u2; jr++) {
				ek(first_left+il, first_right+jr) += -1.0*weight*penalty*phiR_u2(jr,0)*phiL_u2(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr_u2; ir++) {
			for(jl=0; jl<nrowl_u2; jl++) {
				ek(first_right+ir, first_left + jl) += -1.0*weight*penalty*phiL_u2(jl,0)*phiR_u2(ir,0);
			}
		}
        
    }
}

	if (this->IsSymetricOne() && this->IsSymetricTwo()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}

void TPZMatUncoupledPoissonDisc::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    //dados da primeira equacao (variavel u) 
    TPZFMatrix<REAL> &dphiL_u1 = dataleft[0].dphix;
	TPZFMatrix<REAL> &phiL_u1 = dataleft[0].phi;
	TPZManVector<REAL,3> &normal_u1 = data.normal;
	int POrder_u1= dataleft[0].p;
	REAL faceSize_u1=data.HSize;
    
    int nrowl_u1;
	nrowl_u1 = phiL_u1.Rows();
    
    //dados da segunda equacao (variavel p)
    TPZFMatrix<REAL> &dphiL_u2 = dataleft[1].dphix;
	TPZFMatrix<REAL> &phiL_u2 = dataleft[1].phi;
	TPZManVector<REAL,3> &normal_u2 = data.normal;
	int POrder_u2= dataleft[1].p;
	REAL faceSize_u2=data.HSize;
    int nrowl_u2;
	nrowl_u2 = phiL_u2.Rows();
    
    
	int il,jl,id;
    
	switch(bc.Type()) {
		case 0: // Dirichlet nas duas equacoes
			
            //primeira equacao
			for(il=0; il<nrowl_u1; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL_u1(id,il)*normal_u1[id];
				}
				ef(il,0) += weight*fK1*dphiLinormal*bc.Val2()(0,0)*this->fSymmetry1;
				for(jl=0; jl<nrowl_u1; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL_u1(id,jl)*normal_u1[id];
					}
					ek(il,jl) += weight*fK1*(this->fSymmetry1*dphiLinormal*phiL_u1(jl,0) - dphiLjnormal*phiL_u1(il,0));
				}
			}
            
            //segunda equacao
			for(il=0; il<nrowl_u2; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL_u2(id,il)*normal_u2[id];
				}
				ef(il+nrowl_u1,0) += weight*fK2*dphiLinormal*bc.Val2()(1,0)*this->fSymmetry2;
				for(jl=0; jl<nrowl_u2; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL_u2(id,jl)*normal_u2[id];
					}
					ek(il+nrowl_u1,jl+nrowl_u1) += weight*fK2*(this->fSymmetry2*dphiLinormal*phiL_u2(jl,0) - dphiLjnormal*phiL_u2(il,0));
				}
			}
            
            break;
			
		case 1: // Neumann nas duas equacoes
            
            //primeira equacao
			for(il=0; il<nrowl_u1; il++) {
				ef(il,0) += weight*phiL_u1(il,0)*bc.Val2()(0,0);
			}
            
            //segunda equacao
            for(il=0; il<nrowl_u2; il++) {
				ef(il+nrowl_u1,0) += weight*phiL_u2(il,0)*bc.Val2()(1,0);
			}
			break;
			
        case 11: // Neumann nas primeira e Dirichlet na segunda equacao
            
            //primeira equacao
			for(il=0; il<nrowl_u1; il++) {
				ef(il,0) += weight*phiL_u1(il,0)*bc.Val2()(0,0);
			}
            
            //segunda equacao
			for(il=0; il<nrowl_u2; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL_u2(id,il)*normal_u2[id];
				}
				ef(il+nrowl_u1,0) += weight*fK2*dphiLinormal*bc.Val2()(1,0)*this->fSymmetry2;
				for(jl=0; jl<nrowl_u2; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL_u2(id,jl)*normal_u2[id];
					}
					ek(il+nrowl_u1,jl+nrowl_u1) += weight*fK2*(this->fSymmetry2*dphiLinormal*phiL_u2(jl,0) - dphiLjnormal*phiL_u2(il,0));
				}
			}
			break;

		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}
    
    if (this->IsSymetricOne() && this->IsSymetricTwo()){
		if (!ek.VerifySymmetry()) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
	
// ------- Penalty: primeira equacao --------
	if (this->fPenaltyConstant1 != 0.){
        
        nrowl_u1 = phiL_u1.Rows();
        const REAL penalty = fPenaltyConstant1*abs(fK1)*POrder_u1*POrder_u1/faceSize_u1; //Cp^2/h
        
        switch(bc.Type()) {
            case 0: // Dirichlet
                for(il=0; il<nrowl_u1; il++) {
                    ef(il,0) += weight*penalty*phiL_u1(il,0)*bc.Val2()(0,0);
                    for(jl=0; jl<nrowl_u1; jl++) {
                        ek(il,jl) += weight*penalty*phiL_u1(il,0)*phiL_u1(jl,0);
                    }
                }
                
                break;
                
            default:
                PZError << "TPZMatUncoupledPoissonDisc::Wrong boundary condition type\n";
                break;
        }
    }
    
// --------- Penalty: segunda equacao ----------
	if (this->fPenaltyConstant2 != 0.){
        
        nrowl_u1 = phiL_u1.Rows();
        nrowl_u2 = phiL_u2.Rows();
        const REAL penalty = fPenaltyConstant2*abs(fK2)*POrder_u2*POrder_u2/faceSize_u2; //Cp^2/h
        
        switch(bc.Type()) {
            case 0: // Dirichlet
                for(il=0; il<nrowl_u2; il++) {
                    ef(il+nrowl_u1,0) += weight*penalty*phiL_u2(il,0)*bc.Val2()(1,0);
                    for(jl=0; jl<nrowl_u2; jl++) {
                        ek(il+nrowl_u1,jl+nrowl_u1) += weight*penalty*phiL_u2(il,0)*phiL_u2(jl,0);
                    }
                }
                
                break;
                
            default:
                PZError << "TPZMatUncoupledPoissonDisc::Wrong boundary condition type\n";
                break;
        }
    }

}

/** Returns the variable index associated with the name */
int TPZMatUncoupledPoissonDisc::VariableIndex(const std::string &name){
	if(!strcmp("Solution_u1",name.c_str()))        return  1;
	if(!strcmp("Solution_u2",name.c_str()))        return  2;
	if(!strcmp("Derivate_u1",name.c_str()))        return  3;
	if(!strcmp("Derivate_u2",name.c_str()))        return  4;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMatUncoupledPoissonDisc::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 1;
	if((var == 3) || (var == 4)) return fDim;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatUncoupledPoissonDisc::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> Sol_u1, Sol_u2;
	TPZFMatrix<STATE> DSol_u1, DSol_u2;
	TPZFMatrix<REAL> axes_u1, axes_u2;
	
	Sol_u1=datavec[0].sol[0];
	DSol_u1=datavec[0].dsol[0];
	axes_u1=datavec[0].axes;
    
	Sol_u2=datavec[1].sol[0];
	DSol_u2=datavec[1].dsol[0];
	axes_u2=datavec[1].axes;
	
	if(var == 1){
		Solout[0] = Sol_u1[0];//function (state variable u)
		return;
	}
	
	if(var == 2){
		Solout[0] = Sol_u2[0];//function (state variable p)
		return;
	}
	
	if(var == 3) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9,STATE> dsoldx;
			TPZAxesTools<STATE>::Axes2XYZ(DSol_u1, dsoldx, axes_u1);
			Solout[id] = dsoldx(id,0);//derivate of u
		}
		return;
	}//var == 3
	
	if(var == 4) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9,STATE> dsoldx;
			TPZAxesTools<STATE>::Axes2XYZ(DSol_u2, dsoldx, axes_u2);
			Solout[id] = dsoldx(id,0);//derivate of p
		}
		return;
	}//var == 4
}


void TPZMatUncoupledPoissonDisc::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	DebugStop();
}

void TPZMatUncoupledPoissonDisc::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	DebugStop();
}
