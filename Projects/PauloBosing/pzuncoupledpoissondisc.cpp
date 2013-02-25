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


using namespace std;

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc():TPZDiscontinuousGalerkin(), fXf1(0.), fXf2(0.),fDim(1){
	fK1 = 1.;
	fK2 = 1.;
}

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc(int matid, int dim):TPZDiscontinuousGalerkin(matid), fXf1(0.), fXf2(0.),fDim(dim){
	fK1 = 1.;
	fK2 = 1.;
}

TPZMatUncoupledPoissonDisc::~TPZMatUncoupledPoissonDisc(){
}

TPZMatUncoupledPoissonDisc::TPZMatUncoupledPoissonDisc(const TPZMatUncoupledPoissonDisc &copy):TPZDiscontinuousGalerkin(copy){
    
    this->operator=(copy);
}

TPZMatUncoupledPoissonDisc & TPZMatUncoupledPoissonDisc::operator=(const TPZMatUncoupledPoissonDisc &copy){
    
    TPZDiscontinuousGalerkin::operator = (copy);
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

int TPZMatUncoupledPoissonDisc::NStateVariables() {
	return 1;
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

void TPZMatUncoupledPoissonDisc::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
	
	int nref =  datavec.size();
	if (nref != 2 ) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
    
	TPZFMatrix<REAL>  &phiu =  datavec[0].phi;
	TPZFMatrix<REAL> &dphiu = datavec[0].dphix;
	int phru = phiu.Rows();
    
	TPZFMatrix<REAL>  &phip =  datavec[1].phi;
	TPZFMatrix<REAL> &dphip = datavec[1].dphix;
	int phrp = phip.Rows();
	
	//Equacao de Poisson
	// ------ primeira equacao ------
    
	int kd, in, jn;
	for(in = 0; in < phru; in++ ) {
		ef(in, 0) += weight*fXf1*phiu(in,0);
        
		for(jn = 0; jn < phru; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in,jn) += weight*fK1*dphiu(kd,in)*dphiu(kd,jn);
			}
		}
	}
	
	//----- segunda equacao ------
    
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		fXf2 = res[0];
	}
    
	for(in = 0; in < phrp; in++) {
		ef(in+phru, 0) += weight*fXf2*phip(in,0);
		
		for(jn = 0; jn < phrp; jn++ ) {
			for(kd=0; kd<fDim; kd++) {
				ek(in+phru, jn+phru) += weight*fK2*dphip(kd,in)*dphip(kd,jn);
			}
		}
	}
    
}

void TPZMatUncoupledPoissonDisc::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<REAL> &ek,
                                              TPZFMatrix<REAL> &ef,TPZBndCond &bc) {
	
	
	int nref =  datavec.size();
	if (nref != 2) {
		cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
	TPZFMatrix<REAL>  &phiu = datavec[0].phi;
	TPZFMatrix<REAL>  &phip = datavec[1].phi;
    
	int phru = phiu.Rows();
	int phrp = phip.Rows();
	short in,jn;
	REAL v2[2];
	v2[0] = bc.Val2()(0,0); //condicao de contorno da primeira equacao
	v2[1] = bc.Val2()(1,0); //condicao de contorno da segunda equaca
	
	switch (bc.Type()) {
		case 0 : // Dirichlet condition
			
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
			
		case 1 : // Neumann condition
			//primeira equacao
			for(in = 0 ; in < phru; in++) {
				ef(in,0) += v2[0]*phiu(in,0)*weight;
			}
			
			//seguna equacao
			for(in = 0 ; in < phrp; in++) {
				ef(in+phru,0) += v2[1]*phip(in,0)*weight;
			}
			break;
            
        case 10 : // Dirichlet-Neumann condition
            
			//primeira equacao
			for(in = 0 ; in < phru; in++) {
				ef(in,0) += gBigNumber * v2[0]*phiu(in,0)*weight;
				for (jn = 0 ; jn < phru; jn++) {
					ek(in,jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
				}
			}
			
			//segunda equacao
			for(in = 0 ; in < phrp; in++) {
				ef(in+phru,0) += v2[1]*phip(in,0)*weight;
			}
			break;
            
        case 11 : // Neumann-Dirichlet condition
            
			//primeira equacao
			for(in = 0 ; in < phru; in++) {
				ef(in,0) += v2[0]*phiu(in,0)*weight;
			}
			
			//segunda equacao
			for(in = 0 ; in < phrp; in++) {
				ef(in+phru,0) += gBigNumber*v2[1]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+phru,jn+phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			break;
            
            
	}
    
}

void TPZMatUncoupledPoissonDisc::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

	int il,jl,ir,jr,id;
    int dim_eku;
    
//========== Primeira equacao ==========    
{
    //dados da primeira variavel (u)
	TPZFMatrix<REAL> &dphiLdAxes_u = dataleft[0].dphix;
	TPZFMatrix<REAL> &dphiRdAxes_u = dataright[0].dphix;
	TPZFMatrix<REAL> &phiL_u = dataleft[0].phi;
	TPZFMatrix<REAL> &phiR_u = dataright[0].phi;
	TPZManVector<REAL,3> &normal_u = data.normal;//No PZ a normal eh do elemento de menor indice para o de maior indice
    
	TPZFNMatrix<660> dphiL_u, dphiR_u;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes_u, dphiL_u, dataleft[0].axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes_u, dphiR_u, dataright[0].axes);
    
	int &LeftPOrder_u=dataleft[0].p;
	int &RightPOrder_u=dataright[0].p;
    
	REAL &faceSize_u=data.HSize;
    
	int nrowl_u = phiL_u.Rows();
	int nrowr_u = phiR_u.Rows();
    dim_eku = nrowr_u +nrowl_u;
    
    //diffusion term
    REAL leftK1, rightK1;
	leftK1  = fK1;
	rightK1 = fK1;
    
    //symmetry parameter
    REAL theta1;
    theta1 = this->fSymmetry1;
    
    //Contribuicao na matriz de rigidez
    //No PZ a normal eh do elemento de menor indice para o de maior indice
    //por isso tem-se sinais contrarios nas equacoes  (comparadas com a teoria) da contribuicao da matriz de rigidez

	// 1) phi_I_left, phi_J_left 
	for(il=0; il<nrowl_u; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u(id,il)*normal_u[id];
		}
		for(jl=0; jl<nrowl_u; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u(id,jl)*normal_u[id];
			}
			ek(il,jl) += weight*leftK1*(theta1*0.5*dphiLinormal*phiL_u(jl,0)-0.5*dphiLjnormal*phiL_u(il,0));
		}
	}

	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr_u; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u(id,ir)*normal_u[id];
		}
		for(jr=0; jr<nrowr_u; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u(id,jr)*normal_u[id];
			}
			ek(ir+nrowl_u,jr+nrowl_u) += weight*rightK1*(theta1*(-0.5*dphiRinormal*phiR_u(jr)) + 0.5*dphiRjnormal*phiR_u(ir));
		}
	}

	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl_u; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_u(id,il)*normal_u[id];
		}
		for(jr=0; jr<nrowr_u; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_u(id,jr)*normal_u[id];
			}
			ek(il,jr+nrowl_u) += weight*leftK1*(theta1*(-0.5*dphiLinormal*phiR_u(jr)) - 0.5*dphiRjnormal*phiL_u(il));
		}
	}


    // 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr_u; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_u(id,ir)*normal_u[id];
		}
		for(jl=0; jl<nrowl_u; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_u(id,jl)*normal_u[id];
			}
			ek(ir+nrowl_u,jl) += weight*rightK1*(theta1*0.5*dphiRinormal*phiL_u(jl) + 0.5*dphiLjnormal*phiR_u(ir));
		}
	}
    
    if(this->fPenaltyConstant1 !=0 ){
        
        REAL penalty = fPenaltyConstant1*(0.5*(abs(leftK1)*LeftPOrder_u*LeftPOrder_u + abs(rightK1)*RightPOrder_u*RightPOrder_u))/faceSize_u;
        
        // 1) left i / left j
		for(il=0; il<nrowl_u; il++) {
			for(jl=0; jl<nrowl_u; jl++) {
				ek(il,jl) += weight*penalty*phiL_u(il,0)*phiL_u(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr_u; ir++) {
			for(jr=0; jr<nrowr_u; jr++) {
				ek(ir+nrowl_u,jr+nrowl_u) += weight*penalty*phiR_u(ir,0)*phiR_u(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl_u; il++) {
			for(jr=0; jr<nrowr_u; jr++) {
				ek(il,jr+nrowl_u) += -1.0*weight*penalty*phiR_u(jr,0)*phiL_u(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr_u; ir++) {
			for(jl=0; jl<nrowl_u; jl++) {
				ek(ir+nrowl_u,jl) += -1.0*weight*penalty*phiL_u(jl,0)*phiR_u(ir,0);
			}
		}

    }
}
    
// =========== Segunda equacao =========
{
    //dados da segunda variavel (p)
    TPZFMatrix<REAL> &dphiLdAxes_p = dataleft[1].dphix;
	TPZFMatrix<REAL> &dphiRdAxes_p = dataright[1].dphix;
	TPZFMatrix<REAL> &phiL_p = dataleft[1].phi;
	TPZFMatrix<REAL> &phiR_p = dataright[1].phi;
	TPZManVector<REAL,3> &normal_p = data.normal;
    
	TPZFNMatrix<660> dphiL_p, dphiR_p;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes_p, dphiL_p, dataleft[1].axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes_p, dphiR_p, dataright[1].axes);
    
	int &LeftPOrder_p=dataleft[1].p;
	int &RightPOrder_p=dataright[1].p;
    
	REAL &faceSize_p=data.HSize;
    
	int nrowl_p = phiL_p.Rows();
	int nrowr_p = phiR_p.Rows();
    
    //diffusion term
	REAL leftK2, rightK2;
    leftK2  = fK2;
	rightK2 = fK2;
    
    //symmetry parameter
    REAL theta2;
    theta2 = this->fSymmetry2;
    
    
    // 1) phi_I_left, phi_J_left
	for(il=0; il<nrowl_p; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_p(id,il)*normal_p[id];
		}
		for(jl=0; jl<nrowl_p; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_p(id,jl)*normal_p[id];
			}
			ek(dim_eku + il,jl + dim_eku) += weight*leftK2*(theta2*0.5*dphiLinormal*phiL_p(jl,0)-0.5*dphiLjnormal*phiL_p(il,0));
		}
	}
    
	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr_p; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_p(id,ir)*normal_p[id];
		}
		for(jr=0; jr<nrowr_p; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_p(id,jr)*normal_p[id];
			}
			ek(dim_eku + ir+nrowl_p, dim_eku + jr+nrowl_p) += weight*rightK2*(theta2*(-0.5*dphiRinormal*phiR_p(jr)) + 0.5*dphiRjnormal*phiR_p(ir));
		}
	}
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl_p; il++) {
		REAL dphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiLinormal += dphiL_p(id,il)*normal_p[id];
		}
		for(jr=0; jr<nrowr_p; jr++) {
			REAL dphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiRjnormal += dphiR_p(id,jr)*normal_p[id];
			}
			ek(dim_eku + il, dim_eku+ jr+nrowl_p) += weight*leftK2*(theta2*(-0.5*dphiLinormal*phiR_p(jr)) - 0.5*dphiRjnormal*phiL_p(il));
		}
	}
    
    
    // 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr_p; ir++) {
		REAL dphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			dphiRinormal += dphiR_p(id,ir)*normal_p[id];
		}
		for(jl=0; jl<nrowl_p; jl++) {
			REAL dphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				dphiLjnormal += dphiL_p(id,jl)*normal_p[id];
			}
			ek(dim_eku + ir + nrowl_p, dim_eku+ jl) += weight*rightK2*(theta2*0.5*dphiRinormal*phiL_p(jl) + 0.5*dphiLjnormal*phiR_p(ir));
		}
	}
    
    if(this->fPenaltyConstant2 !=0 ){
        
        REAL penalty = fPenaltyConstant2*(0.5*(abs(leftK2)*LeftPOrder_p*LeftPOrder_p + abs(rightK2)*RightPOrder_p*RightPOrder_p))/faceSize_p;
        
        // 1) left i / left j
		for(il=0; il<nrowl_p; il++) {
			for(jl=0; jl<nrowl_p; jl++) {
				ek(dim_eku+il, dim_eku+jl) += weight*penalty*phiL_p(il,0)*phiL_p(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr_p; ir++) {
			for(jr=0; jr<nrowr_p; jr++) {
				ek(dim_eku+ir+nrowl_p, dim_eku + jr+nrowl_p) += weight*penalty*phiR_p(ir,0)*phiR_p(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl_p; il++) {
			for(jr=0; jr<nrowr_p; jr++) {
				ek(dim_eku+il, dim_eku+jr+nrowl_p) += -1.0*weight*penalty*phiR_p(jr,0)*phiL_p(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr_p; ir++) {
			for(jl=0; jl<nrowl_p; jl++) {
				ek(dim_eku+ir+nrowl_p, dim_eku + jl) += -1.0*weight*penalty*phiL_p(jl,0)*phiR_p(ir,0);
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
    TPZFMatrix<REAL> &dphiL_u = dataleft[0].dphix;
	TPZFMatrix<REAL> &phiL_u = dataleft[0].phi;
	TPZManVector<REAL,3> &normal_u = data.normal;
	int POrder_u= dataleft[0].p;
	REAL faceSize_u=data.HSize;
    
    int nrowl_u;
	nrowl_u = phiL_u.Rows();
    
    //dados da segunda equacao (variavel p)
    TPZFMatrix<REAL> &dphiL_p = dataleft[1].dphix;
	TPZFMatrix<REAL> &phiL_p = dataleft[1].phi;
	TPZManVector<REAL,3> &normal_p = data.normal;
	int POrder_p= dataleft[1].p;
	REAL faceSize_p=data.HSize;
    int nrowl_p;
	nrowl_p = phiL_p.Rows();
    
    
	int il,jl,id;
    
	switch(bc.Type()) {
		case 0: // Dirichlet nas duas equacoes
			
            //primeira equacao
			for(il=0; il<nrowl_u; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL_u(id,il)*normal_u[id];
				}
				ef(il,0) += weight*fK1*dphiLinormal*bc.Val2()(0,0)*this->fSymmetry1;
				for(jl=0; jl<nrowl_u; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL_u(id,jl)*normal_u[id];
					}
					ek(il,jl) += weight*fK1*(this->fSymmetry1*dphiLinormal*phiL_u(jl,0) - dphiLjnormal*phiL_u(il,0));
				}
			}
            
            //segunda equacao
			for(il=0; il<nrowl_p; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL_p(id,il)*normal_p[id];
				}
				ef(il+nrowl_u,0) += weight*fK2*dphiLinormal*bc.Val2()(1,0)*this->fSymmetry2;
				for(jl=0; jl<nrowl_p; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL_p(id,jl)*normal_p[id];
					}
					ek(il+nrowl_u,jl+nrowl_u) += weight*fK2*(this->fSymmetry2*dphiLinormal*phiL_p(jl,0) - dphiLjnormal*phiL_p(il,0));
				}
			}
            
            break;
			
		case 1: // Neumann nas duas equacoes
            
            //primeira equacao
			for(il=0; il<nrowl_u; il++) {
				ef(il,0) += weight*phiL_u(il,0)*bc.Val2()(0,0);
			}
            
            //segunda equacao
            for(il=0; il<nrowl_p; il++) {
				ef(il+nrowl_u,0) += weight*phiL_p(il,0)*bc.Val2()(1,0);
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
        
        nrowl_u = phiL_u.Rows();
        const REAL penalty = fPenaltyConstant1*abs(fK1)*POrder_u*POrder_u/faceSize_u; //Cp^2/h
        
        switch(bc.Type()) {
            case 0: // Dirichlet
                for(il=0; il<nrowl_u; il++) {
                    ef(il,0) += weight*penalty*phiL_u(il,0)*bc.Val2()(0,0);
                    for(jl=0; jl<nrowl_u; jl++) {
                        ek(il,jl) += weight*penalty*phiL_u(il,0)*phiL_u(jl,0);
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
        
        nrowl_u = phiL_u.Rows();
        nrowl_p = phiL_p.Rows();
        const REAL penalty = fPenaltyConstant2*abs(fK2)*POrder_p*POrder_p/faceSize_p; //Cp^2/h
        
        switch(bc.Type()) {
            case 0: // Dirichlet
                for(il=0; il<nrowl_p; il++) {
                    ef(il+nrowl_u,0) += weight*penalty*phiL_p(il,0)*bc.Val2()(1,0);
                    for(jl=0; jl<nrowl_p; jl++) {
                        ek(il+nrowl_u,jl+nrowl_u) += weight*penalty*phiL_p(il,0)*phiL_p(jl,0);
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
	if(!strcmp("SolutionU",name.c_str()))        return  1;
	if(!strcmp("SolutionP",name.c_str()))        return  2;
	if(!strcmp("DerivateU",name.c_str()))        return  3;
	if(!strcmp("DerivateP",name.c_str()))        return  4;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMatUncoupledPoissonDisc::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 1;
	if((var == 3) || (var == 4)) return fDim;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatUncoupledPoissonDisc::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolU, SolP;
	TPZFMatrix<REAL> DSolU, DSolP;
	TPZFMatrix<REAL> axesU, axesP;
	
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
			TPZAxesTools<REAL>::Axes2XYZ(DSolU, dsoldx, axesU);
			Solout[id] = dsoldx(id,0);//derivate of u
		}
		return;
	}//var == 3
	
	if(var == 4) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9> dsoldx;
			TPZAxesTools<REAL>::Axes2XYZ(DSolP, dsoldx, axesP);
			Solout[id] = dsoldx(id,0);//derivate of p
		}
		return;
	}//var == 4
}


void TPZMatUncoupledPoissonDisc::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
	DebugStop();
}

void TPZMatUncoupledPoissonDisc::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
	DebugStop();
}
