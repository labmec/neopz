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
	
	// Matrix size is for each element is  phiu.Rows() plus  phip.Rows()
	
	// Setting the size of first block of first problem. Elastic problem
	TPZFMatrix  &phiu =  datavec[0].phi;
//	TPZFMatrix &dphiu = datavec[0].dphix;
	int phru = phiu.Rows();
	
	// Setting the size of second block of second problem. transport problem 
	TPZFMatrix  &phip =  datavec[1].phi;
	TPZFMatrix &dphip = datavec[1].dphix;
	int phrp = phip.Rows();
	
	//Equacao da elasticidade
	// Calculate the matrix contribution for elastic problem
	TPZAutoPointer <TPZElasticityMaterial> matelast = new TPZElasticityMaterial(fmatId, fE, fnu, ff[0], ff[1], fPlaneStress);
	TPZFNMatrix<10000> ekelastic(2*phru,2*phru),efelastic(2*phru,1);
	matelast->Contribute(datavec[0], weight, ekelastic, efelastic);
	// It adds Source matrix on current matrix from position (sRow, sCol). 
	ek.AddSub(0,0, ekelastic);
	ef.AddSub(0, 0, efelastic);		
	
	//Equacao de Poisson: pressao 
	// Calculate the matrix contribution for transport problem from
	for(int in = 0; in < phrp; in++) {
		ef(in+2*phru, 0) += 0.; 
		
		for(int jn = 0; jn < phrp; jn++) {
			for(int kd=0; kd<fDim; kd++) {
				ek(in+2*phru, jn+2*phru) += weight *dphip(kd,in)*dphip(kd,jn); 
			}
		}
	}
	
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		ef.Print("ef = ",sout,EMathematicaInput);
		ekelastic.Print("ekelastic = ",sout,EMathematicaInput);
		efelastic.Print("efelastic = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
			
}

void TPZPoroElastic2d::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix &ek,
									   TPZFMatrix &ef,TPZBndCond &bc) 
{
		
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
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	v2[2] = bc.Val2()(2,0);
		
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition for two equations
			//Equacao da elasticidade
			for(in = 0 ; in < phru; in++) {
				ef(2*in,0) += gBigNumber * v2[0] *   // x displacement
				phiu(in,0) * weight;        // forced v2 displacement
				ef(2*in+1,0) += gBigNumber * v2[1] * // x displacement
				phiu(in,0) * weight;        // forced v2 displacement
				for (jn = 0 ; jn < phru; jn++) {
					ek(2*in,2*jn) += gBigNumber * phiu(in,0) *
					phiu(jn,0) * weight;
					ek(2*in+1,2*jn+1) += gBigNumber * phiu(in,0) *
					phiu(jn,0) * weight;
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
			
		case 11 :		// Neumann condition for two equations
			//Equacao da elasticidade
			for(in = 0 ; in <phru; in++) {           // componentes da tracao normal ao contorno
				ef(2*in,0) += v2[0] * phiu(in,0)*weight;   // tracao em x  (ou pressao)
				ef(2*in+1,0) += v2[1] * phiu(in,0)*weight; // tracao em y (ou pressao) , nula se nao h
			}      // ou deslocamento nulo  v2 = 0
			
			
			//Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru,0) += v2[2]*phip(in,0) * weight;
			}
			break;
			
			case 22 : // // Mixed condition for two equations
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
			
			//Equacao de Poisson: pressao 
			for(in = 0 ; in < phrp; in++) {
				ef(in+2*phru, 0) += v2[2] * phip(in, 0)*weight;
				for (jn = 0 ; jn < phrp; jn++) {
					ek(in+2*phru,jn+2*phru) += bc.Val1()(2,0)*phip(in,0)*
					phip(jn,0)*weight;     // peso de contorno => integral de contorno
				}
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
	if(!strcmp("SolutionU",name.c_str()))        return  1;
	if(!strcmp("SolutionP",name.c_str()))        return  2;
	if(!strcmp("DerivateU",name.c_str()))        return  3;
	if(!strcmp("DerivateP",name.c_str()))        return  4;
	
	return TPZMaterial::VariableIndex(name);
}


void TPZPoroElastic2d::ContributeInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef){
	DebugStop();
}

void TPZPoroElastic2d::ContributeBCInterface(TPZMaterialData &data, REAL weight, TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc){
	DebugStop();
}



