/*
 *  pzporoelastic2d.cpp
 *  PZ
 *
 *  Created by Agnaldo on 11/28/11.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */

#include <iostream>
#include <string>

#include "pzporoelastic2d.h"
#include "pzelasmat.h" 
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.poroelastic.data"));
#endif

TPZPoroElastic2d::EState TPZPoroElastic2d::gState = ECurrentState;

TPZPoroElastic2d::TPZPoroElastic2d():TPZDiscontinuousGalerkin(), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.), fPlaneStress(0) {
	fE = 0.;
	fDim = 2;
	fmatId = 0;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1.;
	
}

TPZPoroElastic2d::TPZPoroElastic2d(int matid, int dim):TPZDiscontinuousGalerkin(matid), ff(0), fnu(0.), falpha(0.), fk(0.), fvisc(0.),fPlaneStress(0) {
	fE = 0.;
	fDim = dim;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1;
	fmatId = matid;
}

TPZPoroElastic2d::~TPZPoroElastic2d(){
}


int TPZPoroElastic2d::NStateVariables() {
	return 1;
}


void TPZPoroElastic2d::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
	

	int nref =  datavec.size();
	if (nref != 2 ) {
		std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	
	//--------- Matrix size is for each element is  phiu.Rows() plus  phip.Rows() ---------------
	//Setting the size of first block of first problem. Elastic problem
	TPZFMatrix  &phiu =  datavec[0].phi;
	TPZFMatrix &dphiu = datavec[0].dphix;
	TPZFMatrix &axes=datavec[0].axes;
	int phcu, phru, dphcu, dphru;
	phru = phiu.Rows();
	phcu = phiu.Cols();
	dphcu = dphiu.Cols();
	dphru = dphiu.Rows();
	
	if(phcu != 1 || dphru != 2 || phru != dphcu) 
	{
		PZError << "\n inconsistent input Elasticity data : \n" <<
		"phi.Cols() = " << phiu.Cols() << " dphi.Cols() = " << dphiu.Cols() <<
		" phi.Rows = " << phiu.Rows() << " dphi.Rows = " << dphiu.Rows() <<"\n";
		return;
	}
	
	// Setting the size of second block of second problem. transport problem 
	TPZFMatrix  &phip =  datavec[1].phi;
	TPZFMatrix &dphip = datavec[1].dphix;
	int phrp = phip.Rows();
	
	TPZFMatrix du(2,2);
	
	int efr, efc, ekr, ekc;  
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	
	if(ekr != (2*phru + phrp) || ekc != (2*phru + phrp) || efr != (2*phru + phrp) || efc != 1)
	{
		PZError << "\n inconsistent input data : \n" << "\nek.Rows() = " << ek.Rows() <<
		" ek.Cols() = " << ek.Cols() << "\nef.Rows() = " << ef.Rows() << " ef.Cols() = " << ef.Cols() << "\n";
		return;
	}
	
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
		
		//Elastic equation: Calculate the matrix contribution for elastic problem 
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
		
		
		// Coupling terms right upper side of global matrix 
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*in,2*phru+jn) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*in+1,2*phru+jn) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));							
			}
		}
		
		
		// Coupling terms left lower side of global matrix 
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*phru+jn,2*in) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*phru+jn,2*in+1) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));							
			}
		}
		
		
		//Equacao de Poisson: pressao 
		// Calculate the matrix contribution for transport problem from
		const REAL DeltaT = fTimeStep;
		for(int in = 0; in < phrp; in++)
		{
			ef(in+2*phru, 0) += 0.; 
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(in+2*phru, jn+2*phru) += (-1.)*weight*fSe*phip(in,0)*phip(jn,0); 
				for(int kd=0; kd<fDim; kd++) 
				{
					ek(in+2*phru, jn+2*phru) += (-1.)*weight *(fk/fvisc)*DeltaT*dphip(kd,in)*dphip(kd,jn);
				}
			}
		}
	}//end if
	
	//Last state (n)
	if(gState == ELastState)
	{				
		for(int in = 0; in < phru; in++ )
		{
			du(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
			du(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
			
			for(int jn = 0; jn < phrp; jn++)
			{
				ek(2*phru+jn,2*in) += (-1.)*falpha*weight*(phip(jn,0)*du(0,0));		
				ek(2*phru+jn,2*in+1) += (-1.)*falpha*weight*(phip(jn,0)*du(1,0));							
			}
		}
		
		for(int in = 0; in < phrp; in++){
			for(int jn = 0; jn < phrp; jn++) {
				ek(in+2*phru, jn+2*phru) += (-1.)*weight*fSe*phip(in,0)*phip(jn,0); 
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

void TPZPoroElastic2d::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix &ek,
									   TPZFMatrix &ef,TPZBndCond &bc) 
{
	//The Last state (n) not include boundary conditions
	if(gState == ELastState){
		return;
	}
			
	int nref =  datavec.size();
	if (nref != 2) {
		std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	if (bc.Val2().Rows() != 3) {
		std::cout << " Erro.!! Neste material precisa-se do valor da condicao de contorno para ux, uy e p.\n";
		std::cout << "portanto precisa-se passar uma matrix Val2(3,1).\n";
		DebugStop();
	}
	
	if (bc.Val1().Rows() != 3) {
		std::cout << " Erro.!! Neste material precisa-se do valor da condicao de contorno para ux, uy e p.\n";
		DebugStop();
	}
	
	TPZFMatrix  &phiu = datavec[0].phi;
	TPZFMatrix  &phip = datavec[1].phi;
	
	int phru = phiu.Rows();
	int phrp = phip.Rows();
	short in,jn;
	REAL v2[3];
	v2[0] = bc.Val2()(0,0);//referente a elasticidade em x
	v2[1] = bc.Val2()(1,0);//referente a elasticidade em y
	v2[2] = bc.Val2()(2,0);//referente a pressao
	
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition for two equations
			//Equacao da elasticidade
			for(in = 0 ; in < phru; in++) {
				ef(2*in,0) += BIGNUMBER*v2[0]*phiu(in,0)*weight;  /// x displacement  forced v2 displacement      
				ef(2*in+1,0) += BIGNUMBER*v2[1]*	phiu(in,0)*weight;   /// y displacement  forced v2 displacement 
				
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;
					ek(2*in+1,2*jn+1) += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;
				}
			}
						
			//segunda equacao
			//Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += gBigNumber * v2[2]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			break;
			
		case 11 :/// Neumann condition for two equations
		{
			//Equacao da elasticidade
			for(in = 0 ; in <phru; in++) {           // componentes da tracao normal ao contorno
				ef(2*in,0) += v2[0]*phiu(in,0)*weight;   // tracao em x  (ou pressao)
				ef(2*in+1,0) += v2[1]*phiu(in,0)*weight; // tracao em y (ou pressao) , nula se nao h
			}      // ou deslocamento nulo  v2 = 0
			
			
			//Equacao de Poisson: pressao 
			const REAL DeltT = fTimeStep;
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += v2[2]*DeltT*phip(in,0) * weight;
			}
			break;
		}	
		case 22 : /// Mixed condition for two equations
		{
			//Equacao da elasticidade			
			for(in = 0 ; in < phru; in++) {
				ef(2*in, 0) += v2[0] * phiu(in, 0) * weight;   // Neumann , Sigmaij
				ef(2*in+1, 0) += v2[1] * phiu(in, 0) * weight; // Neumann
				
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += bc.Val1()(0,0)*phiu(in,0)*
					phiu(jn,0)*weight;         // peso de contorno => integral de contorno
					ek(2*in+1,2*jn) += bc.Val1()(1,0)*phiu(in,0)*
					phiu(jn,0)*weight;
					ek(2*in+1,2*jn+1) += bc.Val1()(1,1)*phiu(in,0)*
					phiu(jn,0)*weight;
					ek(2*in,2*jn+1) += bc.Val1()(0,1)*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			} 
			
			///Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru, 0) += v2[2] * phip(in, 0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += bc.Val1()(2,0)*phip(in,0)*
					phip(jn,0)*weight;     // peso de contorno => integral de contorno
				}
			}
			
			break;
		}
			
		case 10: //Neumann condition for elastic and Dirichlet condition for pressure
		{			
		
			//Neumann para Equacao da elasticidade
			for(in = 0 ; in <phru; in++) {           // componentes da tracao normal ao contorno
				ef(2*in,0) += v2[0]*phiu(in,0)*weight;   // tracao em x  (ou pressao)
				ef(2*in+1,0) += v2[1]*phiu(in,0)*weight; // tracao em y (ou pressao) , nula se nao h
			}      // ou deslocamento nulo  v2 = 0
			
			//Dirichlet para Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += gBigNumber * v2[2]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			
			break;
		}
			
		case 1: //Dirichlet condition for elastic (0) and  Neumann condition (1) for pressure
		{			
			
			//Dirichlet para Equacao da elasticidade
			for(in = 0 ; in < phru; in++) {
				ef(2*in,0) += BIGNUMBER*v2[0]*phiu(in,0)*weight;    /// x displacement forced v2 displacement
				ef(2*in+1,0) += BIGNUMBER*v2[1]*phiu(in,0)*weight;   /// y displacement  forced v2 displacement      
				
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;
					ek(2*in+1,2*jn+1) += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight;
				}
			}
			
			//Neumann para Equacao de Poisson: pressao 
			const REAL DeltT = fTimeStep;
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += v2[2]*DeltT*phip(in,0) * weight;
			}
						
			break;
		}	
		
		///Condicoes de fronteiras livres
		case 100:///Dirichlet para as duas equacoes mas com fronteira livre em y na equacao da elasticiade
			
			//Equacao da elasticidade
			for(in = 0 ; in < phru; in++) {
				ef(2*in,0) += BIGNUMBER*v2[0] *   // x displacement
				phiu(in,0)*weight;        // forced v2 displacement
				
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += BIGNUMBER*phiu(in,0)*
					phiu(jn,0)*weight;
				}
			}
			
			//Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += gBigNumber * v2[2]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			break;
					
		case 200: //Neumann condition free in x for elastic and Dirichlet condition for pressure
		{			
			//Neumann para Equacao da elasticidade
			for(in = 0 ; in <phru; in++) {           // componentes da tracao normal ao contorno
				ef(2*in+1,0) += v2[1]*phiu(in,0)*weight; // tracao em y (ou pressao) , nula se nao h
			}      // ou deslocamento nulo  v2 = 0
			
			//Dirichlet para Equacao da pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += gBigNumber * v2[2]*phip(in,0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
				}
			}
			
			break;
		}
		case 300: //Dirichlet condition free in y for elastic and Neumann condition for pressure
		{			
			///Dirichlet para Equacao da elasticidade
			for(in = 0 ; in < phru; in++) {
				ef(2*in,0) += BIGNUMBER*v2[0] *phiu(in,0)*weight;     /// x displacement forced v2 displacement
				
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += BIGNUMBER*phiu(in,0)*phiu(jn,0)*weight; /// x displacement
				}
			}
			
			///Neumann para Equacao da pressao 
			const REAL DeltT = fTimeStep;
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += v2[2]*DeltT*phip(in,0) * weight;
			}
				
			break;
		}

	}
	
}

void TPZPoroElastic2d::Print(std::ostream &out) {
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
int TPZPoroElastic2d::VariableIndex(const std::string &name){
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

int TPZPoroElastic2d::NSolutionVariables(int var){
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

void TPZPoroElastic2d::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolU, SolP;
	TPZFMatrix DSolU, DSolP;
	TPZFMatrix axesU, axesP;
	
	TPZVec<REAL> ptx(3), solExata(3);
	TPZFMatrix flux(3,1);
    
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
		fForcingFunctionExact(datavec[1].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var10
	
	//solucao deslocamento y exata
	if(var == 11){
		fForcingFunctionExact(datavec[0].x, solExata,flux);
		Solout[0] = solExata[1];
		return;
	}//var11
	
	//solucao sigmaY exata
	if(var == 12){
		fForcingFunctionExact(datavec[0].x, solExata,flux);
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
		TPZAxesTools::Axes2XYZ(DSolP, dsoldx, axesP);
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


void TPZPoroElastic2d::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, 
                                           REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
	DebugStop();
}

void TPZPoroElastic2d::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                             REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
	DebugStop();
}



