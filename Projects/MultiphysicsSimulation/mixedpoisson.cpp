//
//  mixedpoisson.cpp
//  PZ
//
//  Created by Agnaldo Farias on 5/28/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzlog.h"
#include "mixedpoisson.h"


#include <iostream>


TPZMixedPoisson::TPZMixedPoisson():TPZDiscontinuousGalerkin(){
	fk = 1.;
	ff = 0.;
    fDim=2;
}

TPZMixedPoisson::TPZMixedPoisson(int matid):TPZDiscontinuousGalerkin(matid){
	fk = 1.;
	ff = 0.;
    fDim=2;
}

TPZMixedPoisson::~TPZMixedPoisson(){};

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
    phrq = phiQ.Rows();
    phrp = phip.Rows();

    
	
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
            REAL prod = datavec[0].fNormalVec(0,ivecind)*datavec[1].fNormalVec(0,jvecind)+
            datavec[1].fNormalVec(1,ivecind)*datavec[0].fNormalVec(1,jvecind)+
            datavec[1].fNormalVec(2,ivecind)*datavec[0].fNormalVec(2,jvecind);//dot product between u and v 
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
        datavec[1].axes.Multiply(ivec,axesvec);
        
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
    
}
