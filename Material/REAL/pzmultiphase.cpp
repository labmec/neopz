//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzlog.h"
#include "pzmultiphase.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.pzmultiphase.data"));
#endif

TPZMultiphase::TPZMultiphase(): TPZDiscontinuousGalerkin(){
	ff = 0.;
}

TPZMultiphase::TPZMultiphase(int matid, int dim): TPZDiscontinuousGalerkin(matid){
	ff = 0.;
  
}

TPZMultiphase::~TPZMultiphase(){
}

int TPZMultiphase::NStateVariables() {
	return 1;
}

void TPZMultiphase::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
#ifdef DEBUG
	int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
#endif
    
   
}

void TPZMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
    
#ifdef DEBUG
    int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	if (bc.Type() > 1 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
		DebugStop();
	}
#endif
	
    
}

/** Returns the variable index associated with the name */
int TPZMultiphase::VariableIndex(const std::string &name){
	if(!strcmp("Flux",name.c_str()))        return  1;
	if(!strcmp("Pressure",name.c_str()))    return  2;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMultiphase::NSolutionVariables(int var){
	if(var == 1) return 2;
	if(var == 2) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMultiphase::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolP, SolQ;
    
	SolQ = datavec[0].sol[0];
	SolP = datavec[1].sol[0];

}


void TPZMultiphase::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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




