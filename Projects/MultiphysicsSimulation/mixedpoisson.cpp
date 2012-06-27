//
//  mixedpoisson.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzlog.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
#endif

TPZMixedPoisson::TPZMixedPoisson(): TPZMatPoisson3d(), fDim(1){
	fk = 1.;
	ff = 0.;
}

TPZMixedPoisson::TPZMixedPoisson(int matid, int dim): TPZMatPoisson3d(matid,dim), fDim(dim){
	fk = 1.;
	ff = 0.;
  
}

TPZMixedPoisson::~TPZMixedPoisson(){
}

int TPZMixedPoisson::NStateVariables() {
	return 1;
}

void TPZMixedPoisson::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Coeficient which multiplies the gradient operator "<< fk << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMixedPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
#ifdef DEBUG
	int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
#endif
    
    if(fForcingFunction) {            
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);   
		ff = res[0];
	}
    
    // Setting the phis
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phip =  datavec[1].phi;
	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
	
	//Calculate the matrix contribution for flux. Matrix A
    REAL ratiok = 1./fk;
    for(int iq=0; iq<phrq; iq++)
    {
        ef(iq, 0) += 0.;
        
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        for (int jq=0; jq<phrq; jq++) 
        {
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            REAL prod = datavec[0].fNormalVec(0,ivecind)*datavec[0].fNormalVec(0,jvecind)+
            datavec[0].fNormalVec(1,ivecind)*datavec[0].fNormalVec(1,jvecind)+
            datavec[0].fNormalVec(2,ivecind)*datavec[0].fNormalVec(2,jvecind);//dot product between u and v 
            ek(iq,jq) += ratiok*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
        }
    }
    
	
	// Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fNormalVec(0,ivecind);
        ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
        ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        TPZFNMatrix<3> axesvec(3,1);
        datavec[0].axes.Multiply(ivec,axesvec);
        
        REAL divwq = 0.;
        for(int iloc=0; iloc<fDim; iloc++)
        {
            divwq += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
        }
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;            
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
        }
    }
    
    for(int ip=0; ip<phrp; ip++){
         ef(phrq+ip,0) += (-1.)*weight*ff*phip(ip,0);
    }

//#ifdef LOG4CXX
//    if(logdata->isDebugEnabled())
//	{
//        std::stringstream sout;
//        sout<<"\n\n Matriz ek e vetor fk \n ";
//        ek.Print("ekmph = ",sout,EMathematicaInput);
//        ef.Print("efmph = ",sout,EMathematicaInput);
//        LOGPZ_DEBUG(logdata,sout.str());
//	}
//#endif
    
}

void TPZMixedPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
    
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
	
	TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
	int phrq = datavec[0].fVecShapeIndex.NElements();

	REAL v2[2];
	v2[0] = bc.Val2()(0,0);
	
	switch (bc.Type()) {
		case 0 :		// Dirichlet condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2[0]*phiQ(iq,0)*weight;
            }
            break;
			
		case 1 :			// Neumann condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2[0]*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight; 
                }
            }  
			break;
        
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
				ef(iq,0) += v2[0]*phiQ(iq,0)*weight;
				for (int jq = 0; jq < phrq; jq++) {
					ek(iq,jq) += weight*bc.Val1()(0,0)*phiQ(iq,0)*phiQ(jq,0);
				}
			}
            
            break;
	}
    
}

/** Returns the variable index associated with the name */
int TPZMixedPoisson::VariableIndex(const std::string &name){
	if(!strcmp("Flux",name.c_str()))        return  1;
	if(!strcmp("Pressure",name.c_str()))        return  2;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMixedPoisson::NSolutionVariables(int var){
	if(var == 1) return fDim;
	if(var == 2) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMixedPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<REAL> SolP, SolQ;
    
   // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 1){ //function (state variable Q)
		Solout[0] = datavec[0].sol[0][0];
        Solout[1] = datavec[0].sol[0][1];
		return;
	}
    
    if(var == 2){
		Solout[0] = SolP[0];//function (state variable p)
		return;
	}
    
    
	
	
	
}


void TPZMixedPoisson::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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




