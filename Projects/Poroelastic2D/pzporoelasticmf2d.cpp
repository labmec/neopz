/*
 *  pzporoelasticmf2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 03/15/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>

#include "pzporoelasticmf2d.h"
#include "pzelasmat.h" 
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic2d.data"));
#endif

TPZPoroElasticMF2d::EState TPZPoroElasticMF2d::gState = ECurrentState;

TPZPoroElasticMF2d::TPZPoroElasticMF2d():TPZDiscontinuousGalerkin(), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.), fPlaneStress(0) {
	fE = 0.;
	fDim = 2;
	fmatId = 0;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1.;
	
}

TPZPoroElasticMF2d::TPZPoroElasticMF2d(int matid, int dim):TPZDiscontinuousGalerkin(matid), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.),fPlaneStress(0) {
	fE = 0.;
	fDim = dim;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1;
	fmatId = matid;
}

TPZPoroElasticMF2d::~TPZPoroElasticMF2d(){
}


int TPZPoroElasticMF2d::NStateVariables() {
	return 1;
}


void TPZPoroElasticMF2d::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef){
	

	int nref =  datavec.size();
	if (nref != 3 ) {
		std::cout << " Erro.!! datavec tem que ser de tamanho 3 \n";
		DebugStop();
	}
	
	// Setting the phis
	TPZFMatrix<>  &phiu =  datavec[0].phi;
    TPZFMatrix<>  &phiQ =  datavec[1].phi;
    TPZFMatrix<>  &phip =  datavec[2].phi;
    
    TPZFMatrix<> &dphiu = datavec[0].dphix;
	TPZFMatrix<> &dphiQ = datavec[1].dphix;
	TPZFMatrix<> &dphip = datavec[2].dphix;

    int phru, phrq, phrp;
    phru = phiu.Rows();
    phrq = phiQ.Rows();
    phrp = phip.Rows();
    
    if(phrq!=datavec[1].fVecShapeIndex.NElements()){
        PZError << "\n inconsistent input Flux data : \n";
		return;
    }
    
	
    TPZFMatrix<> du(2,2);
    TPZFMatrix<> &axes = datavec[0].axes;
	//current state (n+1)
	if(gState == ECurrentState)
	{	   
		REAL fEover1MinNu2 = fE/(1-fnu*fnu);  ///4G(lamb+G)/(lamb+2G)
		REAL fEover21PlusNu = 2.*fE/(2.*(1+fnu));/*fE/(2.*(1+fnu));*/ ///2G=2mi
		
		/*
		 * Plain strain materials values
		 * 2G=2mi=nu2*F, lamb=fnu*F, lamb+2G=nu1*F 
		 */
		REAL nu1 = 1 - fnu;
		REAL nu2 = (1-2*fnu);//(1-2*fnu)/2;
		REAL F = fE/((1+fnu)*(1-2*fnu));
		
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
            ef(2*phru+iq, 0) += 0.;
            
            int ivecind = datavec[1].fVecShapeIndex[iq].first;
            int ishapeind = datavec[1].fVecShapeIndex[iq].second;
            for (int jq=0; jq<phrq; jq++) 
            {
                int jvecind = datavec[1].fVecShapeIndex[jq].first;
                int jshapeind = datavec[1].fVecShapeIndex[jq].second;
                REAL prod = datavec[1].fNormalVec(0,ivecind)*datavec[1].fNormalVec(0,jvecind)+
                datavec[1].fNormalVec(1,ivecind)*datavec[1].fNormalVec(1,jvecind)+
                datavec[1].fNormalVec(2,ivecind)*datavec[1].fNormalVec(2,jvecind);//dot product between u and v 
                ek(2*phru+iq,2*phru+jq) += fTimeStep*ratiomuk*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
            }
        }
        
        // Calculate the matrix contribution for pressure. Matrix -D
		for(int in = 0; in < phrp; in++)
		{
			ef(2*phru+phrq+in,0) += 0.; 
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
            ivec(0,0) = datavec[1].fNormalVec(0,ivecind);
            ivec(1,0) = datavec[1].fNormalVec(1,ivecind);
            ivec(2,0) = datavec[1].fNormalVec(2,ivecind);
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
        for(int in = 0; in < phru; in++ ){
            du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
            
            for(int jn = 0; jn < phrp; jn++){
                ek(2*phru+phrq+jn,2*in) += (-1.)*falpha*weight*phip(jn,0)*du(0,0);		
				ek(2*phru+phrq+jn,2*in+1) += (-1.)*falpha*weight*phip(jn,0)*du(1,0);
            }
        }
        for(int in = 0; in < phrp; in++){
            for(int jn = 0; jn < phrp; jn++) {
				ek(in+2*phru+phrq, jn+2*phru+phrq) += (-1.)*fSe*weight*phip(in,0)*phip(jn,0); 
            }
        }
    }
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		ef.Print("ef = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
}

void TPZPoroElasticMF2d::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    for(int in = 0 ; in < phiu.Rows(); in++){
        ef(2*in,0) += gBigNumber*bc.Val2()(0,0)*phiu(in,0)*weight; // x displacement forced v2 displacement      
        ef(2*in+1,0) += gBigNumber*bc.Val2()(1,0)*phiu(in,0)*weight; // y displacement forced v2 displacement 
        
        for(int jn = 0 ; jn < phiu.Rows(); jn++){
            ek(2*in,2*jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
            ek(2*in+1,2*jn+1) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyDirichlet_PQ(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    for(int iq=0; iq<phiQ.Rows(); iq++)
    {
        //the contribution of the Dirichlet boundary condition appears in the flow equation
        ef(2*phru+iq,0) += -fTimeStep*bc.Val2()(2,0)*phiQ(iq,0)*weight;
    }
}

void TPZPoroElasticMF2d::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    for(int in = 0 ; in <phru; in++){  
        ef(2*in,0) += bc.Val2()(0,0)*phiu(in,0)*weight;   // traction in x 
        ef(2*in+1,0) += bc.Val2()(1,0)*phiu(in,0)*weight; // traction in y
    }     
}

void TPZPoroElasticMF2d::ApplyNeumann_PQ(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    for(int iq=0; iq<phiQ.Rows(); iq++)
    {
        ef(2*phru+iq,0)+= gBigNumber*bc.Val2()(2,0)*phiQ(iq,0)*weight;
        for (int jq=0; jq<phiQ.Rows(); jq++) {
            
            ek(2*phru+iq,2*phru+jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight; 
        }
    }   
}

void TPZPoroElasticMF2d::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
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

void TPZPoroElasticMF2d::ApplyMixed_PQ(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    int phru = datavec[0].phi.Rows();
    TPZFMatrix<> &phiQ = datavec[1].phi;
    for(int iq = 0; iq < phiQ.Rows(); iq++) {
        
        ef(2*phru+iq,0) += bc.Val2()(2,0)*phiQ(iq,0)*weight;
        for (int jq = 0; jq < phiQ.Rows(); jq++) {
           
            ek(2*phru+iq,2*phru+jq) += weight*bc.Val1()(2,0)*phiQ(iq,0)*phiQ(jq,0);
        }
    }
}

void TPZPoroElasticMF2d::ApplyDirichletFreeY_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    for(int in = 0 ; in < phiu.Rows(); in++){
        ef(2*in,0) += gBigNumber*bc.Val2()(0,0)*phiu(in,0)*weight; // x displacement forced v2 displacement      
        for(int jn = 0 ; jn < phiu.Rows(); jn++){
            ek(2*in,2*jn) += gBigNumber*phiu(in,0)*phiu(jn,0)*weight;
        }
    }
}

void TPZPoroElasticMF2d::ApplyNeumannFreeX_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc)
{
    TPZFMatrix<>  &phiu = datavec[0].phi;
    int phru = phiu.Rows();
    for(int in = 0 ; in <phru; in++){  
        ef(2*in+1,0) += bc.Val2()(1,0)*phiu(in,0)*weight; // traction in y
    }     
}

void TPZPoroElasticMF2d::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<> &ek,
                                      TPZFMatrix<> &ef,TPZBndCond &bc) 
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
            ApplyDirichlet_PQ(datavec, weight, ef, bc);
            break;
            
        case 11: // Neumann condition for two equations (elastic and mixed)
            ApplyNeumann_U(datavec, weight, ef, bc);
            ApplyNeumann_PQ(datavec, weight, ek, ef, bc);
            break;
            
        case 22:
            ApplyMixed_U(datavec, weight, ek, ef, bc);
            ApplyMixed_PQ(datavec, weight, ek, ef, bc);
            break;
        
        case 1: //Dirichlet condition for elastic (0) and  Neumann condition (1) for mixed problem
            ApplyDirichlet_U(datavec, weight, ek, ef, bc);
            ApplyNeumann_PQ(datavec, weight, ek, ef, bc);
            break;
            
        case 10: // Neumann condition for elastic and Dirichlet condition for mixed problem
            ApplyNeumann_U(datavec, weight, ef, bc);
            ApplyDirichlet_PQ(datavec, weight, ef, bc);
            break;
            
        case 100: // Dirichlet for two equations, but with free boundary in y for elastic equation
            ApplyDirichletFreeY_U(datavec, weight, ek, ef, bc);
            ApplyDirichlet_PQ(datavec, weight, ef, bc);
            break;
            
        case 200: //Neumann condition free in x for elastic and Dirichlet condition for mixed problem
            ApplyNeumannFreeX_U(datavec, weight, ef, bc);
            ApplyDirichlet_PQ(datavec, weight, ef, bc);
            break;
            
        case 300: //Dirichlet condition free in y for elastic and Neumann condition for mixed problem
            ApplyDirichletFreeY_U(datavec, weight, ek, ef, bc);
            ApplyNeumann_PQ(datavec, weight, ek, ef, bc);
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
	//variaveis da elasticidade
	if(!strcmp("Displacement",name.c_str()))        return  1;
	if(!strcmp("Pressure",name.c_str()))        return  2;
	if(!strcmp("SigmaX",name.c_str()))        return  3;
	if(!strcmp("SigmaY",name.c_str()))        return  4;
	if(!strcmp("TauXY",name.c_str()))        return  5;
	if(!strcmp("DisplacementX",name.c_str()))  return 8;
	if(!strcmp("DisplacementY",name.c_str()))  return 9;
	
	//variaveis da pressao
	if(!strcmp("SolutionP",name.c_str()))        return  6;
	if(!strcmp("MinusKGradP",name.c_str()))        return  7;
	
	//solucao exata problema teste 1D
	if(!strcmp("PressaoExata",name.c_str()))  return 10;
	if(!strcmp("DeslocamentoYExata",name.c_str()))  return 11;
	if(!strcmp("SigmaYExata",name.c_str()))  return 12;
		
	return TPZMaterial::VariableIndex(name);
}

int TPZPoroElasticMF2d::NSolutionVariables(int var){
	if(var == 1) return 3;
	if(var == 2) return 1;
	if(var == 3) return 1;
	if(var == 4) return 1;
	if(var == 5) return 1;
	if(var == 6) return 1;
	if(var == 7) return fDim;
	if(var == 8) return 1;
	if(var == 9) return 1;
	if(var == 10) return 1;
	if(var == 11) return 1;
	if(var == 12) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZPoroElasticMF2d::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolU, SolP;
	TPZFMatrix<> DSolU, DSolP;
	TPZFMatrix<> axesU, axesP;
	
	TPZVec<REAL> ptx(3), solExata(3);
	TPZFMatrix<> flux(3,1);
    
    if (datavec[0].sol.size() != 1) {
        DebugStop();
    }
	
	SolU=datavec[0].sol[0];
	DSolU=datavec[0].dsol[0];
	axesU=datavec[0].axes;
	SolP=datavec[1].sol[0];
	DSolP=datavec[1].dsol[0];
	axesP=datavec[1].axes;
	
	//function (state variable u)
	if(var == 1){
		Solout[0] = SolU[0];
		Solout[1] = SolU[1];
		Solout[2] = 0.;
		return;
	}//var1
	
	//function (state variable ux)
	if(var == 8){
		Solout[0] = SolU[0];
		return;
	}//var8
	
	//function (state variable uy)
	if(var == 9){
		Solout[0] = SolU[1];
		return;
	}//var9

	//solucao pressao exata
	if(var == 10){
		fForcingFunctionExact->Execute(datavec[1].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var10
	
	//solucao deslocamento y exata
	if(var == 11){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[1];
		return;
	}//var11
	
	//solucao sigmaY exata
	if(var == 12){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[2];
		return;
	}//var12
	
	
	//-----------------
	if(var == 6) {
		Solout[0] = SolP[0];//derivate of P
		return;
	}//var6
	
	if (var == 7){ //MinusKGradU
		int id;
		//REAL val = 0.;
		TPZFNMatrix<9> dsoldx;
		TPZAxesTools<REAL>::Axes2XYZ(DSolP, dsoldx, axesP);
		for(id=0 ; id<fDim; id++) {
			Solout[id] = -1.*(fk/fvisc)*dsoldx(id,0);
		}
		return;
	}//var7
			
	//-----------------
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
			
	if(var == 2) {
		Solout[0] = SigX+SigY;
		return;
	}//var2
	
	if(var ==3) {
		Solout[0] = SigX;
		return;
	}//var3
	
	if(var == 4) {
		Solout[0] = SigY;
		return;
	}//var4
	
	if(var == 5) {
		Tau = fE*epsxy/(1.+fnu);
		Solout[0] = Tau;
		return;
	}//var5

}


void TPZPoroElasticMF2d::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                           REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef){
	DebugStop();
}

void TPZPoroElasticMF2d::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                             REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
	DebugStop();
}



