//
//  pzporoelasticmf2d.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/20/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>
#include <string>

#include "pzporoelasticmf2d.h"
#include "pzelasmat.h" 
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static PZLogger logdata("pz.material.poroelastic2d.data");
#endif

TPZPoroElasticMF2d::EState TPZPoroElasticMF2d::gState = ECurrentState;

TPZPoroElasticMF2d::TPZPoroElasticMF2d():TPZMaterial(){
    fE = 0.;
    fnu = 0.;
    flambda = 0;
    fmu = 0.;
    falpha = 0.;
    fSe = 0.;
    fk = 0.;
    fvisc = 0.;
    fCf = 0.;
    fSaux = 0.;
    
    fmatId = 0.;
	fDim = 0;
    fLref = 0.;
    fkovervisc=0.;
    
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
    fpref = 0.;
    
    fSf = 0.;
    
	fPlaneStress = 1;
    fReturnSolutionDimension = false;
}

TPZPoroElasticMF2d::TPZPoroElasticMF2d(int matid, int dim):TPZMaterial(matid){
    
    fE = 0.;
    fnu = 0.;
    flambda = 0;
    fmu = 0.;
    falpha = 0.;
    fSe = 0.;
    fk = 0.;
    fvisc = 0.;
    fkovervisc=0.;
    fCf = 0.;
    fSaux = 0.;
    
    fmatId = matid;
	fDim = dim;
    fLref = 0.;
    
    fmatIdSourceTerm = 0;
    
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
    fpref = 0.;
    
    fSf = 0.;
    
	fPlaneStress = 1;
    fReturnSolutionDimension = false;
}

TPZPoroElasticMF2d::~TPZPoroElasticMF2d(){
}

TPZPoroElasticMF2d::TPZPoroElasticMF2d(const TPZPoroElasticMF2d &copy){
    this->operator=(copy);
}

TPZPoroElasticMF2d &TPZPoroElasticMF2d::operator=(const TPZPoroElasticMF2d &copy){
    
    ff = copy.ff;
    fnu = copy.fnu;
    falpha = copy.falpha;
    fk = copy.fk;
    fvisc = copy.fvisc;
    fPlaneStress = copy.fPlaneStress;
    fE = copy.fE;
    fDim = copy.fDim;
    fPlaneStress = copy.fPlaneStress;
    fmatId = copy.fmatId;
    fSe = copy.fSe;
    flambda = copy.flambda;
    fmu = copy.fmu;
    fpref = copy.fpref;
    fLref = copy.fLref;
    fkovervisc = copy.fkovervisc;
    fReturnSolutionDimension = copy.fReturnSolutionDimension;
    fCf = copy.fCf;
    fSaux = copy.fSaux;
    
    fSf=copy.fSf;
    fmatIdSourceTerm=copy.fmatIdSourceTerm;
    
    return *this;
}

void TPZPoroElasticMF2d::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	
    
	int nref =  datavec.size();
	if (nref != 3 ) {
		std::cout << " Erro.!! datavec tem que ser de tamanho 3 \n";
		DebugStop();
	}
	
	// Setting the phis
	TPZFMatrix<REAL>  &phiu =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiQ =  datavec[1].phi;
    TPZFMatrix<REAL>  &phip =  datavec[2].phi;
    
    TPZFMatrix<REAL> &dphiu = datavec[0].dphix;
	TPZFMatrix<REAL> &dphiQ = datavec[1].dphix;
    
    int phru, phrq, phrp;
    phru = phiu.Rows();
    phrq = datavec[1].fVecShapeIndex.NElements();
    phrp = phip.Rows();
	
    TPZFMatrix<STATE> du(2,2);
    TPZFMatrix<REAL> &axes = datavec[0].axes;
	//current state (n+1)
	if(gState == ECurrentState)
	{
		REAL fEover1MinNu2 = fE/(1.-fnu*fnu);  ///4G(lamb+G)/(lamb+2G)
		REAL fEover21PlusNu = fE/(2.*(1.+fnu));//2fE/(2.*(1+fnu)); ///2G=2mi
		
		/*
		 * Plain strain materials values
		 * 2G=2mi=nu2*F, lamb=fnu*F, lamb+2G=nu1*F
		 */
		REAL nu1 = 1.-fnu;
		REAL nu2 = (1.-2.*fnu)/2.;//2(1-2*fnu);
		REAL F = fE/((1.+fnu)*(1.-2.*fnu));
		
		//Calculate the matrix contribution for elastic problem. Matrix A
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			ef(2*in, 0) += weight*ff[0]*phiu(in, 0);
			ef(2*in+1, 0) += weight*ff[1]*phiu(in, 0);
			
			for(int jn = 0; jn < phru; jn++)
			{
				du(0,1) = dphiu(0,jn)*axes(0,0)+dphiu(1,jn)*axes(1,0);
				du(1,1) = dphiu(0,jn)*axes(0,1)+dphiu(1,jn)*axes(1,1);
				
				if (fPlaneStress != 1){
					/* Plain Strain State */
					ek(2*in,2*jn) += weight*(nu1*du(0,0)*du(0,1) + nu2*du(1,0)*du(1,1))*F;
					
					ek(2*in,2*jn+1) += weight*(fnu*du(0,0)*du(1,1) + nu2*du(1,0)*du(0,1))*F;
					
					ek(2*in+1,2*jn) += weight*(fnu*du(1,0)*du(0,1) + nu2*du(0,0)*du(1,1))*F;
					
					ek(2*in+1,2*jn+1) += weight*(nu1*du(1,0)*du(1,1) + nu2*du(0,0)*du(0,1))*F;
				}
				else{
					/* Plain stress state */
					ek(2*in,2*jn) += weight*(fEover1MinNu2*du(0,0)*du(0,1) + fEover21PlusNu*du(1,0)*du(1,1));
					
					ek(2*in,2*jn+1) += weight*(fEover1MinNu2*fnu*du(0,0)*du(1,1) + fEover21PlusNu*du(1,0)*du(0,1));
					
					ek(2*in+1,2*jn) += weight*(fEover1MinNu2*fnu*du(1,0)*du(0,1) + fEover21PlusNu*du(0,0)*du(1,1));
					
					ek(2*in+1,2*jn+1) += weight*(fEover1MinNu2*du(1,0)*du(1,1) + fEover21PlusNu*du(0,0)*du(0,1));
				}
			}
		}
        
        // Calculate the matrix contribution for flux. Marix C
        REAL ratiomuk =fvisc/fk;
        
        for(int iq=0; iq<phrq; iq++)
        {
            ef(2*phru+iq, 0) += 0.;//falta colocar termo da gravidade (densidadefluido*phi.campoGravidade)
            
            int ivecind = datavec[1].fVecShapeIndex[iq].first;
            int ishapeind = datavec[1].fVecShapeIndex[iq].second;
            for (int jq=0; jq<phrq; jq++)
            {
                int jvecind = datavec[1].fVecShapeIndex[jq].first;
                int jshapeind = datavec[1].fVecShapeIndex[jq].second;
                REAL prod = datavec[1].fDeformedDirections(0,ivecind)*datavec[1].fDeformedDirections(0,jvecind)+
                datavec[1].fDeformedDirections(1,ivecind)*datavec[1].fDeformedDirections(1,jvecind)+
                datavec[1].fDeformedDirections(2,ivecind)*datavec[1].fDeformedDirections(2,jvecind);//dot product between u and v
                ek(2*phru+iq,2*phru+jq) += fTimeStep*ratiomuk*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
            }
        }
        
        // Calculate the matrix contribution for pressure. Matrix (-1)*D
        if(fForcingFunction) {
            TPZManVector<STATE> res(1);
            fForcingFunction->Execute(datavec[2].x,res);
            fSf = res[0];
        }

		for(int in = 0; in < phrp; in++)
		{
			ef(2*phru+phrq+in,0) += (-1.)*fTimeStep*weight*fSf*phip(in,0);//source term
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(in+2*phru+phrq, jn+2*phru+phrq) += (-1.)*fSe*weight*phip(in,0)*phip(jn,0);
			}
		}
        
		// Coupling terms between displacement and pressure. Matrix B
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
                //Matrix B
				ek(2*in,2*phru+phrq+jn) += (-1.)*falpha*weight*phip(jn,0)*du(0,0);
				ek(2*in+1,2*phru+phrq+jn) += (-1.)*falpha*weight*phip(jn,0)*du(1,0);
                
                //matrix BˆT
                ek(2*phru+phrq+jn,2*in) += (-1.)*falpha*weight*phip(jn,0)*du(0,0);
				ek(2*phru+phrq+jn,2*in+1) += (-1.)*falpha*weight*phip(jn,0)*du(1,0);
			}
		}
        
        // Coupling terms between flux and pressure. Matrix E
        for(int i=0; i<phrq; i++)
        {
            int ivecind = datavec[1].fVecShapeIndex[i].first;
            int ishapeind = datavec[1].fVecShapeIndex[i].second;
            
            TPZFNMatrix<3> ivec(3,1);
            ivec(0,0) = datavec[1].fDeformedDirections(0,ivecind);
            ivec(1,0) = datavec[1].fDeformedDirections(1,ivecind);
            ivec(2,0) = datavec[1].fDeformedDirections(2,ivecind);
            TPZFNMatrix<3> axesvec(3,1);
            datavec[1].axes.Multiply(ivec,axesvec);
            
            REAL divwq = 0.;
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            for (int j=0; j<phrp; j++) {
                REAL fact = (-1.)*fTimeStep*weight*phip(j,0)*divwq;
                
                // Matrix E
                ek(2*phru+i,2*phru+phrq+j) += fact;
                
                // Matrix E^T
                ek(2*phru+phrq+j,2*phru+i) += fact;
            }
        }
		
	}//end if
	
    
	//Last state (n)
	if(gState == ELastState)
	{
        //Coupling terms between displacement and pressure. Matrix BˆT
        for(int in = 0; in < phru; in++ ){
            du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
            
            for(int jn = 0; jn < phrp; jn++){
                ek(2*phru+phrq+jn,2*in) += (-1.)*falpha*weight*phip(jn,0)*du(0,0);
				ek(2*phru+phrq+jn,2*in+1) += (-1.)*falpha*weight*phip(jn,0)*du(1,0);
            }
        }
        
        // Calculate the matrix contribution for pressure. Matrix (-1)*D
        for(int in = 0; in < phrp; in++){
            for(int jn = 0; jn < phrp; jn++) {
				ek(in+2*phru+phrq, jn+2*phru+phrq) += (-1.)*fSe*weight*phip(in,0)*phip(jn,0);
            }
        }
    }
	
#ifdef LOG4CXX
	if(logdata.isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		ef.Print("ef = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
}

void TPZPoroElasticMF2d::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    REAL v2x, v2y;
    
    if(bc.HasForcingFunction())
    {
		TPZManVector<STATE> res(3);
		bc.ForcingFunction()->Execute(datavec[0].x,res);
		v2x = res[0];
        v2y = res[1];
	}else
    {
        v2x = bc.Val2()(0,0);
        v2y = bc.Val2()(1,0);
    }
    for(int in = 0 ; in < phiu.Rows(); in++){
        ef(2*in,0) += gBigNumber*v2x*phiu(in,0)*weight; // x displacement forced v2 displacement
        ef(2*in+1,0) += gBigNumber*v2y*phiu(in,0)*weight; // y displacement forced v2 displacement
        
        for(int jn = 0 ; jn < phiu.Rows(); jn++){
            ek(2*in,2*jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
            ek(2*in+1,2*jn+1) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyDirichlet_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    int phrq = datavec[1].fVecShapeIndex.NElements();
    REAL v2=0.;
    //---
    
    if(bc.HasForcingFunction()){
		TPZManVector<STATE> res(3);
		bc.ForcingFunction()->Execute(datavec[1].x,res);
		v2 = res[2];
	}
    else{
        
       v2 = bc.Val2()(2,0);
    }
    
    for(int iq=0; iq<phrq; iq++)
    {
        //the contribution of the Dirichlet boundary condition appears in the flow equation
        ef(2*phru+iq,0) += (-1.)*fTimeStep*v2*phiQ(iq,0)*weight;
    }
}

void TPZPoroElasticMF2d::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    REAL v2x, v2y;
    v2x = bc.Val2()(0,0);
    v2y = bc.Val2()(1,0);
    for(int in = 0 ; in <phru; in++){
        ef(2*in,0) += v2x*phiu(in,0)*weight;   // traction in x
        ef(2*in+1,0) += v2y*phiu(in,0)*weight; // traction in y
    }
}

void TPZPoroElasticMF2d::ApplyNeumann_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    int phrq = datavec[1].fVecShapeIndex.NElements();
    
    REAL v2;
    if(bc.HasForcingFunction()){
		TPZManVector<STATE> res(4);
		bc.ForcingFunction()->Execute(datavec[1].x,res);
		v2 = res[3];
	}
    else{
        
        v2 = bc.Val2()(2,0);
    }
    
    for(int iq=0; iq<phrq; iq++)
    {
        ef(2*phru+iq,0)+= gBigNumber*fTimeStep*v2*phiQ(iq,0)*weight;
        for (int jq=0; jq<phrq; jq++) {
            
            ek(2*phru+iq,2*phru+jq)+= gBigNumber*fTimeStep*phiQ(iq,0)*phiQ(jq,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    for(int in = 0; in < phru; in++){
        ef(2*in,0) += bc.Val2()(0,0)*phiu(in,0)*weight;   // Neumann , Sigmaij
        ef(2*in+1,0) += bc.Val2()(1,0)*phiu(in,0)*weight; // Neumann
        
        for(int jn = 0 ; jn < phru; jn++){
            ek(2*in,2*jn) += bc.Val1()(0,0)*phiu(in,0)*phiu(jn,0)*weight;
            ek(2*in+1,2*jn) += bc.Val1()(1,0)*phiu(in,0)*phiu(jn,0)*weight;
            ek(2*in+1,2*jn+1) += bc.Val1()(1,1)*phiu(in,0)*phiu(jn,0)*weight;
            ek(2*in,2*jn+1) += bc.Val1()(0,1)*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyMixed_QP(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    int phrq = datavec[1].fVecShapeIndex.NElements();
    for(int iq = 0; iq < phrq; iq++) {
        
        ef(2*phru+iq,0) += bc.Val2()(2,0)*phiQ(iq,0)*weight;
        for (int jq = 0; jq < phrq; jq++) {
            
            ek(2*phru+iq,2*phru+jq) += weight*bc.Val1()(2,0)*phiQ(iq,0)*phiQ(jq,0);
        }
    }
}

void TPZPoroElasticMF2d::ApplyDirichletFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    REAL v2x = bc.Val2()(0,0);
    for(int in = 0 ; in < phiu.Rows(); in++){
        ef(2*in,0) += gBigNumber*v2x*phiu(in,0)*weight; // x displacement forced v2 displacement
        for(int jn = 0 ; jn < phiu.Rows(); jn++){
            ek(2*in,2*jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyDirichletFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    REAL v2y;
    v2y = bc.Val2()(1,0);
    for(int in = 0 ; in < phiu.Rows(); in++){
        ef(2*in+1,0) += gBigNumber*v2y*phiu(in,0)*weight; // y displacement forced v2 displacement
        
        for(int jn = 0 ; jn < phiu.Rows(); jn++){
            ek(2*in+1,2*jn+1) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}


void TPZPoroElasticMF2d::ApplyNeumannFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    for(int in = 0 ; in <phru; in++){
        ef(2*in+1,0) += bc.Val2()(1,0)*phiu(in,0)*weight; // traction in y
    }
}


void TPZPoroElasticMF2d::ApplyNeumannFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    REAL v2x;
    v2x = bc.Val2()(0,0);
    for(int in = 0 ; in <phru; in++){
        ef(2*in,0) += v2x*phiu(in,0)*weight;   // traction in x
    }
}

//void TPZPoroElasticMF2d::ApplySourceTerm_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
//    
//    TPZFMatrix<REAL>  &phip =  datavec[2].phi;
//    int phru = datavec[0].phi.Rows();
//    int phrq = datavec[1].fVecShapeIndex.NElements();
//    int phrp =  phip.Rows();
//    
//    for(int in = 0; in < phrp; in++)
//    {
//        ef(2*phru+phrq+in,0) += weight*fSf*phip(in,0);
//    }
//}

void TPZPoroElasticMF2d::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,
                                      TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	//The Last state (n) not include boundary conditions
	if(gState == ELastState){
		return;
	}
    
	int nref =  datavec.size();
	if (nref != 3) {
		std::cout << " Error.! datavec has to be of size 3 \n";
		DebugStop();
	}
	if (bc.Val2().Rows() != 3) {
		std::cout << " Error.! This material must be the value of the boundary condition for ux, uy, and p\n";
		std::cout << " Therefore, we need to pass an matrix Val2(3,1)\n";
		DebugStop();
	}
	
	if (bc.Val1().Rows() != 3) {
		std::cout << " Error.! This material must be the value of the boundary condition for ux, uy, and p\n";
		DebugStop();
	}
	
    switch (bc.Type())
    {
        case 0: // Dirichlet condition for two equations (elastic and mixed)
            ApplyDirichlet_U(datavec, weight, ek, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 11: // Neumann condition for two equations (elastic and mixed)
            ApplyNeumann_U(datavec, weight, ef, bc);
            ApplyNeumann_QP(datavec, weight, ek, ef, bc);
            break;
            
        case 22:
            ApplyMixed_U(datavec, weight, ek, ef, bc);
            ApplyMixed_QP(datavec, weight, ek, ef, bc);
            break;
            
        case 1: //Dirichlet condition for elastic (0) and  Neumann condition (1) for mixed problem
            ApplyDirichlet_U(datavec, weight, ek, ef, bc);
            ApplyNeumann_QP(datavec, weight, ek, ef, bc);
            break;
            
        case 10: // Neumann condition for elastic and Dirichlet condition for mixed problem
            ApplyNeumann_U(datavec, weight, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 100: // Dirichlet for two equations, but with free boundary in y for elastic equation
            ApplyDirichletFreeY_U(datavec, weight, ek, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 200: //Neumann condition free in x for elastic and Dirichlet condition for mixed problem
            ApplyNeumannFreeX_U(datavec, weight, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 300: //Mixed condition free in y for elastic and Neumann condition for mixed problem
            ApplyDirichletFreeY_U(datavec, weight, ek, ef, bc);
            ApplyNeumann_QP(datavec, weight, ek, ef, bc);
            break;
            
        case 20: //Mixed condition for elastic and Dirichlet condition for mixed problem
            ApplyMixed_U(datavec, weight, ek, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 21: //Mixed condition for elastic and Neumann condition for mixed problem
            ApplyMixed_U(datavec, weight, ek, ef, bc);
            ApplyNeumann_QP(datavec, weight, ek, ef, bc);
            break;
            
        case 400: //Dirichlet free in y e Neumann free in x for elastic and Dirichlet condition for mixed problem
            ApplyDirichletFreeY_U(datavec, weight, ek, ef, bc);
            ApplyNeumannFreeX_U(datavec, weight, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;
            
        case 500: //Dirichlet free in x e Neumann free in y for elastic and Dirichlet condition for mixed problem
            ApplyDirichletFreeX_U(datavec, weight, ek, ef, bc);
            ApplyNeumannFreeY_U(datavec, weight, ef, bc);
            ApplyDirichlet_QP(datavec, weight, ef, bc);
            break;

            
    }
	
}

void TPZPoroElasticMF2d::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\t E   = " << fE   << std::endl;
	out << "\t nu   = " << fnu   << std::endl;
	out << "\t Forcing function F   = " << ff[0] << ' ' << ff[1]   << std::endl;
	out << "Biot constant  = " << falpha << std::endl;
	out << "Permeabilidade da rocha fK "<< fk << std::endl;
	out << "Viscosidade do fluido "<< fvisc << std::endl;
	out << "2D problem " << fPlaneStress << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

/** Returns the variable index associated with the name */
int TPZPoroElasticMF2d::VariableIndex(const std::string &name){
	//variables of solid
	if(!strcmp("Displacement",name.c_str()))        return  1;
    if(!strcmp("DisplacementX",name.c_str()))  return 2;
	if(!strcmp("DisplacementY",name.c_str()))  return 3;
    
	if(!strcmp("SigmaX",name.c_str()))        return  4;
	if(!strcmp("SigmaY",name.c_str()))        return  5;
	if(!strcmp("TauXY",name.c_str()))        return  6;
    if(!strcmp("PressureStress",name.c_str()))        return  7;
	
	//variables of fluid
	if(!strcmp("PorePressure",name.c_str()))        return  8;
	if(!strcmp("Fluxo",name.c_str()))        return  9;//using Hdiv
    if(!strcmp("FluxoX",name.c_str()))        return  17;//using Hdiv
    if(!strcmp("FluxoY",name.c_str()))        return  18;//using Hdiv
    if(!strcmp("MinusKMuGradP",name.c_str()))     return  10;
	
	//Exact soluion
	if(!strcmp("ExactPressure",name.c_str()))  return 11;
    if(!strcmp("ExactDisplacementX",name.c_str()))  return 12;
    if(!strcmp("ExactDisplacementY",name.c_str()))  return 13;
    if(!strcmp("ExactSigmaX",name.c_str()))  return 14;
    if(!strcmp("ExactSigmaY",name.c_str()))  return 15;
    if(!strcmp("ExactFluxo",name.c_str()))  return 16;
    if(!strcmp("ExactDisplacement",name.c_str()))  return 19;
    
	return TPZMaterial::VariableIndex(name);
}

int TPZPoroElasticMF2d::NSolutionVariables(int var){
	if(var == 1) return fDim;
    if(var ==  2) return 1;
	if(var ==  3) return 1;
    if(var ==  4) return 1;
	if(var ==  5) return 1;
	if(var ==  6) return 1;
	if(var ==  7) return 1;
	
	if(var ==  8) return 1;
	if(var ==  9) return fDim;
	if(var ==  10) return fDim;
    
	if(var == 11) return 1;
    if(var == 12) return 1;
    if(var == 13) return 1;
    if(var == 14) return 1;
    if(var == 15) return 1;
    if(var == 16) return fDim;
    if(var == 17) return 1;
    if(var == 18) return 1;
    if(var == 19) return fDim;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZPoroElasticMF2d::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> SolU, SolP;
	TPZFMatrix<STATE> DSolU, DSolP;
	TPZFMatrix<REAL> axesU, axesP;
	
	TPZVec<REAL> ptx(3);
	TPZVec<STATE> solExata(5);
	TPZFMatrix<STATE> flux(2,1);
	
	SolU=datavec[0].sol[0];
	DSolU=datavec[0].dsol[0];
	axesU=datavec[0].axes;
    
	SolP=datavec[2].sol[0];
	DSolP=datavec[2].dsol[0];
    axesP=datavec[2].axes;
    
    //variables of solid
	//function (state variable u)
	if(var == 1){
        if(fReturnSolutionDimension ==true){
            Solout[0] = SolU[0]*fpref*fLref*fSaux;
            Solout[1] = SolU[1]*fpref*fLref*fSaux;
            
        }
        else{
            Solout[0] = SolU[0];
            Solout[1] = SolU[1];
        }
		//Solout[2] = 0.;
		return;
	}//var1
	
	//function (state variable ux)
	if(var == 2){
		if(fReturnSolutionDimension ==true)
            Solout[0] = SolU[0]*fpref*fLref*fSaux;
        else
            Solout[0] = SolU[0];
        
		return;
	}//var2
	
	//function (state variable uy)
	if(var == 3){
        if(fReturnSolutionDimension ==true) Solout[0] = SolU[1]*fpref*fLref*fSaux;
		else Solout[0] = SolU[1];
		return;
	}//var3
	
	
	REAL epsx;
	REAL epsy;
	REAL epsxy;
	REAL SigX;
	REAL SigY;
	REAL Tau, DSolxy[2][2];
    
	// dudx - dudy
	DSolxy[0][0] = DSolU(0,0)*axesU(0,0)+DSolU(1,0)*axesU(1,0);
	DSolxy[1][0] = DSolU(0,0)*axesU(0,1)+DSolU(1,0)*axesU(1,1);
	
	// dvdx - dvdy
	DSolxy[0][1] = DSolU(0,1)*axesU(0,0)+DSolU(1,1)*axesU(1,0);
	DSolxy[1][1] = DSolU(0,1)*axesU(0,1)+DSolU(1,1)*axesU(1,1);
	
	epsx = DSolxy[0][0];// du/dx
	epsy = DSolxy[1][1];// dv/dy
	epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
	REAL Gmodule = fE/(1-fnu*fnu);
	
	if (this->fPlaneStress){
		SigX = Gmodule*(epsx+fnu*epsy);
		SigY = Gmodule*(fnu*epsx+epsy);
	}
	else{
		SigX = fE/((1.-2.*fnu)*(1.+fnu))*((1.-fnu)*epsx+fnu*epsy);
		SigY = fE/((1.-2.*fnu)*(1.+fnu))*(fnu*epsx+(1.-fnu)*epsy);
	}
    
	if(var == 4) {
		if(fReturnSolutionDimension ==true) Solout[0] = SigX*fpref;
		else Solout[0] = SigX;
		return;
	}//var4
	
	if(var == 5) {
        if(fReturnSolutionDimension ==true) Solout[0] = SigY*fpref;
		else Solout[0] = SigY;
		return;
	}//var5
	
	if(var == 6) {
		Tau = fE*epsxy/(1.+fnu);
        if(fReturnSolutionDimension==true) Solout[0] = Tau*fpref;
		else Solout[0] = Tau;
		return;
	}//var6
    
    if(var == 7) {
        if(fReturnSolutionDimension==true) Solout[0] = (SigX+SigY)*fpref;
		else Solout[0] = SigX+SigY;
		return;
	}//var7
    
    
    //-----------------
    
    //variables of fluid
	if(var == 8) {
        if(fReturnSolutionDimension ==true) Solout[0] = SolP[0]*fpref;
		else Solout[0] = SolP[0];
		return;
	}//var8
	
	if (var == 9){
        if(fReturnSolutionDimension ==true){
            Solout[0] = fpref*(fkovervisc)*datavec[1].sol[0][0];
            Solout[1] = fpref*(fkovervisc)*datavec[1].sol[0][1];
        }
        else{
            Solout[0] = datavec[1].sol[0][0];
            Solout[1] = datavec[1].sol[0][1];
        }
        
		return;
	}//var9
    
    if (var == 17){
        if(fReturnSolutionDimension ==true){
            Solout[0] = fpref*(fkovervisc)*datavec[1].sol[0][0];
        }
        else{
            Solout[0] = datavec[1].sol[0][0];
        }
        
		return;
	}//var17

    if (var == 18){
        if(fReturnSolutionDimension ==true){
            Solout[0] = fpref*(fkovervisc)*datavec[1].sol[0][1];
        }
        else{
            Solout[0] = datavec[1].sol[0][1];
        }
        
		return;
	}//var18

    
    if (var == 10){
		int id;
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSolP, dsoldx, axesP);
		for(id=0 ; id<fDim; id++) {
            if(fReturnSolutionDimension==true) Solout[id] = -1.*dsoldx(id,0)*fpref*(fkovervisc);
			else Solout[id] = -1.*dsoldx(id,0);
		}
		return;
	}//var10
    
    
    //Exact soluion
	if(var == 11){
		fExactSol->Execute(datavec[2].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var11
    
    if(var == 12){
		fExactSol->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[1];
		return;
	}//var12
    
    if(var == 13){
		fExactSol->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[2];
		return;
	}//var13
    
    if(var == 19){
		fExactSol->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[1];
        Solout[1] = solExata[2];
		return;
	}//var19

    
    if(var == 14){
		fExactSol->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[3];
		return;
	}//var14
    
    if(var == 15){
		fExactSol->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[4];
		return;
	}//var15
    
    if(var == 16){
		fExactSol->Execute(datavec[1].x, solExata,flux);
		Solout[0] = flux(0,0);
        Solout[1] = flux(1,0);
		return;
	}//var16
}


void TPZPoroElasticMF2d::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                             REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
	DebugStop();
}

void TPZPoroElasticMF2d::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                               REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	DebugStop();
}

void TPZPoroElasticMF2d::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(false);
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = false;
	}
	
}


