//
//  pzelasthybrid.cpp
//  PZ
//
//  Created by Joao on 18/09/12.
//
//
/**
 * @file
 * @brief Contains implementations of the TPZElasticityMaterial methods.
 */

#include "pzelasthybrid.h"

#include "pzelasmat.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logdata("pz.material.elasticity.data");
#endif

#include <fstream>
using namespace std;

TPZElasticityHybridMaterial::TPZElasticityHybridMaterial() : TPZRegisterClassId(&TPZElasticityHybridMaterial::ClassId),
TPZElasticityMaterial() {
}

TPZElasticityHybridMaterial::TPZElasticityHybridMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress) : 
TPZRegisterClassId(&TPZElasticityHybridMaterial::ClassId),
TPZElasticityMaterial(num, E, nu, fx, fy, plainstress) {
	
}

TPZElasticityHybridMaterial::~TPZElasticityHybridMaterial() {
}

int TPZElasticityHybridMaterial::NStateVariables() const {
	return 2;
}

void TPZElasticityHybridMaterial::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    std::cout<<" Passou pelo Contribute"<<std::endl;
    TPZMaterialData::MShapeFunctionType shapetype = data.fShapeType;
    if(shapetype==data.EVecShape){
        ContributeVecShape(data,weight,ek, ef);
        return;
    }
    
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &axes=data.axes;
	
	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if(phc != 1 || dphr != 2 || phr != dphc ||
	   ekr != phr*2 || ekc != phr*2 ||
	   efr != phr*2 ){
		PZError << "\nTPZElasticityMaterial.contr, inconsistent input data : \n" <<
		"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
		"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
		return;
		//		PZError.show();
	}
	if(fForcingFunction) {            // phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(data.x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	TPZFMatrix<STATE> du(2,2);
	/*
	 * Plain strain materials values
	 */
	REAL nu1 = 1 - fnu;//(1-nu)
	REAL nu2 = (1-2*fnu)/2;
	REAL F = fE/((1+fnu)*(1-2*fnu));
	
	for( int in = 0; in < phr; in++ ) {
		du(0,0) = dphi(0,in)*axes(0,0)+dphi(1,in)*axes(1,0);
		du(1,0) = dphi(0,in)*axes(0,1)+dphi(1,in)*axes(1,1);
		
        for (int col = 0; col < efc; col++)
        {
            ef(2*in, col) += weight * (ff[0] * phi(in, 0)- du(0,0)*fPreStressXX - du(1,0)*fPreStressXY);  // dire�o x
            ef(2*in+1, col) += weight * (ff[1] * phi(in, 0)- du(0,0)*fPreStressYY - du(1,0)*fPreStressXY);// dire�o y <<<----
        }
		for( int jn = 0; jn < phr; jn++ ) {
			du(0,1) = dphi(0,jn)*axes(0,0)+dphi(1,jn)*axes(1,0);
			du(1,1) = dphi(0,jn)*axes(0,1)+dphi(1,jn)*axes(1,1);
			
			
			if (fPlaneStress != 1){
				/* Plain Strain State */
				ek(2*in,2*jn) += weight * (
										   nu1 * du(0,0)*du(0,1)+ nu2 * du(1,0)*du(1,1)
										   ) * F;
				
				ek(2*in,2*jn+1) += weight * (
											 fnu*du(0,0)*du(1,1)+ nu2*du(1,0)*du(0,1)
											 ) * F;
				
				ek(2*in+1,2*jn) += weight * (
											 fnu*du(1,0)*du(0,1)+ nu2*du(0,0)*du(1,1)
											 ) * F;
				
				ek(2*in+1,2*jn+1) += weight * (
											   nu1*du(1,0)*du(1,1)+ nu2*du(0,0)*du(0,1)
											   ) * F;
			}
			else{
				/* Plain stress state */
				ek(2*in,2*jn) += weight * (
										   fEover1MinNu2 * du(0,0)*du(0,1)+ fEover21PlusNu * du(1,0)*du(1,1)
										   );
				
				ek(2*in,2*jn+1) += weight * (
											 fEover1MinNu2*fnu*du(0,0)*du(1,1)+ fEover21PlusNu*du(1,0)*du(0,1)
											 );
				
				ek(2*in+1,2*jn) += weight * (
											 fEover1MinNu2*fnu*du(1,0)*du(0,1)+ fEover21PlusNu*du(0,0)*du(1,1)
											 );
				
				ek(2*in+1,2*jn+1) += weight * (
											   fEover1MinNu2*du(1,0)*du(1,1)+ fEover21PlusNu*du(0,0)*du(0,1)
											   );
			}
		}
	}
	
    //#ifdef PZ_LOG
    //	if(logdata.isDebugEnabled())
    //	{
    //		std::stringstream sout;
    //		ek.Print("ek_elastmat = ",sout,EMathematicaInput);
    //		ef.Print("ef_elastmat = ",sout,EMathematicaInput);
    //		LOGPZ_DEBUG(logdata,sout.str())
    //	}
    //#endif
	
}

void TPZElasticityHybridMaterial::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    std::cout<< "Passou ContributeInterface"<< std::endl;
}

void TPZElasticityHybridMaterial::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    std::cout<< "Passou ContributeBCInterface"<< std::endl;
}


TPZElasticityHybridMaterial::TPZElasticityHybridMaterial(const TPZElasticityHybridMaterial &copy) :
TPZRegisterClassId(&TPZElasticityHybridMaterial::ClassId),
TPZElasticityMaterial(copy)
{
}


int TPZElasticityHybridMaterial::ClassId() const{
    return Hash("TPZElasticityHybridMaterial") ^ TPZElasticityMaterial::ClassId() << 1;
}



