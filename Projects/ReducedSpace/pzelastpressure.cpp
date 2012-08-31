//
//  pzelastpressure.cpp
//  PZ
//
//  Created by Agnaldo Farias on 8/27/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>
#include <string>

#include "pzelastpressure.h"
#include "pzelasmat.h" 
#include "pzbndcond.h"
#include "pzaxestools.h"

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.elastpressure"));
#endif


TPZElastPressure::TPZElastPressure():TPZDiscontinuousGalerkin(), ff(0), fnu(0.), fk(0.),fXf(0.), fPlaneStress(0) {
	fE = 0.;
	fDim = 2;
	fmatId = 0;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1.;
	
}

TPZElastPressure::TPZElastPressure(int matid, int dim):TPZDiscontinuousGalerkin(matid), ff(0), fnu(0.), fk(0.), fXf(0.),fPlaneStress(0) {
	fE = 0.;
	fDim = dim;
	ff.resize(2);
	ff[0]=0.;
	ff[1]=0.;
	fPlaneStress = 1;
	fmatId = matid;
}

TPZElastPressure::~TPZElastPressure(){
}


int TPZElastPressure::NStateVariables() {
	return 1;
}

void TPZElastPressure::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\t E   = " << fE   << std::endl;
	out << "\t nu   = " << fnu   << std::endl;
	out << "\t Forcing function F   = " << ff[0] << ' ' << ff[1]   << std::endl; 
	out << "Permeabilidade da rocha fK "<< fk << std::endl;
	out << "2D problem " << fPlaneStress << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}


void TPZElastPressure::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
    int nref =  datavec.size();
	if (nref != 2 ) {
		std::cout << " Error.!! the size of datavec is equal to two\n";
        std::cout << " datavec[0]->elasticity and datavec[1]->pressure\n";
		DebugStop();
	}
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!=0){
        
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
    
    //Calculate the matrix contribution for elastic problem.
    TPZFMatrix<REAL> &dphiu = datavec[0].dphix;
	TPZFMatrix<REAL> &phiu = datavec[0].phi;
	TPZFMatrix<REAL> &axes=datavec[0].axes;
	
	int phcu,efcu;
	phcu = phiu.Cols();
	efcu = ef.Cols();
	
	if(fForcingFunction) {// phi(in, 0) :  node in associated forcing function
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(datavec[0].x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	TPZFNMatrix<4,STATE> dphix_i(2,1),dphiy_i(2,1), dphix_j(2,1), dphiy_j(2,1);
	/*
	 * Plain strain materials values
	 */
    REAL fEover1MinNu2 = fE/(1-fnu*fnu);  ///4G(lamb+G)/(lamb+2G)
    REAL fEover21PlusNu = 2.*fE/(2.*(1+fnu));/*fE/(2.*(1+fnu));*/ ///2G=2mi
	REAL nu1 = 1 - fnu;//(1-nu)
	REAL nu2 = (1-2*fnu)/2;
	REAL F = fE/((1+fnu)*(1-2*fnu));
    
	for(int in = 0; in < phcu; in++) {
		dphix_i(0,0) = dphiu(0,in)*axes(0,0)+dphiu(1,in)*axes(1,0);
		dphix_i(1,0) = dphiu(0,in)*axes(0,1)+dphiu(1,in)*axes(1,1);
		dphiy_i(0,0) = dphiu(2,in)*axes(0,0)+dphiu(3,in)*axes(1,0);
		dphiy_i(1,0) = dphiu(2,in)*axes(0,1)+dphiu(3,in)*axes(1,1);
		
        for (int col = 0; col < efcu; col++) 
        {
            ef(in,col) += weight*(ff[0]*phiu(0, in) + ff[1]*phiu(1, in));
        }		
		for( int jn = 0; jn < phcu; jn++ ) {
            
            dphix_j(0,0) = dphiu(0,jn)*axes(0,0)+dphiu(1,jn)*axes(1,0);
            dphix_j(1,0) = dphiu(0,jn)*axes(0,1)+dphiu(1,jn)*axes(1,1);
            dphiy_j(0,0) = dphiu(2,jn)*axes(0,0)+dphiu(3,jn)*axes(1,0);
            dphiy_j(1,0) = dphiu(2,jn)*axes(0,1)+dphiu(3,jn)*axes(1,1);
			
			
			if (fPlaneStress != 1){/* Plain Strain State */
                ek(in,jn) += weight*(nu1*dphix_i(0,0)*dphix_j(0,0) + nu2*dphix_i(1,0)*dphix_j(1,0) +
                                     
                                     fnu*dphix_i(0,0)*dphiy_j(1,0) + nu2*dphix_i(1,0)*dphiy_j(0,0) +
                 
                                     fnu*dphiy_i(1,0)*dphix_j(0,0) + nu2*dphiy_i(0,0)*dphix_j(1,0) +
                 
                                     nu1*dphiy_i(1,0)*dphiy_j(1,0) + nu2*dphiy_i(0,0)*dphiy_j(0,0))*F;
			}
			else{/* Plain stress state */
                ek(in,jn) += weight*(fEover1MinNu2*dphix_i(0,0)*dphix_j(0,0) + fEover21PlusNu*dphix_i(1,0)*dphix_j(1,0) +
                                     
                                     fEover1MinNu2*dphix_i(0,0)*dphiy_j(1,0) + fEover21PlusNu*dphix_i(1,0)*dphiy_j(0,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dphix_j(0,0) + fEover21PlusNu*dphiy_i(0,0)*dphix_j(1,0) +
                                     
                                     fEover1MinNu2*dphiy_i(1,0)*dphiy_j(1,0) + fEover21PlusNu*dphiy_i(0,0)*dphiy_j(0,0));
            }
		}
	}
    
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek_reduced = ",sout,EMathematicaInput);
		ef.Print("ef_reduced = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
    
     // Calculate the matrix contribution for pressure
    ContributePressure(datavec, weight, ek, ef);
}

void TPZElastPressure::ContributePressure(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef){
    
   if(!datavec[1].phi) return;
    TPZFMatrix<REAL>  &phip =  datavec[1].phi;
    TPZFMatrix<REAL>  &dphip =  datavec[1].dphix;
    int phrp = phip.Rows();
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
    int phcu = phiu.Cols();
    
//    int r =1;// fk.Rows();
//	int c =1;// fk.Cols();
//	TPZFMatrix<STATE> submat(r,c);
//    submat.Redim(r, c);
//    
//	for(int in=0 ; in < phrp; ++in){
//        
//        REAL tmpef = phip(in,0)*weight;
//		TPZFMatrix<REAL> tmpTPZFMatrix1 = fXf*tmpef;
//		ef.AddSub(in*r + phcu,0,tmpTPZFMatrix1);
//        //ef(in*r,0) = fXf[0]*tmpef;
//		
//		for(int jn=0 ; jn<phrp; ++jn){
//            
//			REAL temp = dphip(0,in)*dphip(0,jn)*weight;
//			submat = fk*temp;
//			ek.AddSub(in*r + phcu,jn*c + phcu,submat);
//		}
//	}
    
    for(int in = 0; in<phrp; in++){
        
        ef(in+phcu,0)+=weight*(fXf*phip(in,0));
        
        for(int jn=0; jn<phrp; jn++){
          
            ek(in+phcu, jn+phcu)+=fk*dphip(0,in)*dphip(0,jn)*weight;
        }
    }
}

void TPZElastPressure::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	const REAL big  = TPZMaterial::gBigNumber;
	int phc = phiu.Cols();
    
    for(int in = 0 ; in < phc; in++) {
        for (int il = 0; il <fNumLoadCases; il++) 
        {
            TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
            
            ef(in,il) += big*(v2(0,il)*phiu(0,in) + v2(1,il)*phiu(1,in))*weight;
        }
        for (int jn = 0 ; jn < phc; jn++) {
            
            ek(in,jn) += big*(phiu(0,in)*phiu(0,jn) + phiu(1,in)*phiu(1,jn))*weight;
        }
    }
}

void TPZElastPressure::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ef,TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	int phc = phiu.Cols();
    ef.Print();
    for (int in = 0; in < phc; in++) 
    {
        for (int il = 0; il <fNumLoadCases; il++) 
        {
            TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
            ef(in,il)+= weight*(v2(0,il)*phiu(0,in) + v2(1,il)*phiu(1,in));
        }
    }
    
    ef.Print();
    
}

void TPZElastPressure::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	int phc = phiu.Cols();
    
    for(int in = 0 ; in < phc; in++) 
    {
        for (int il = 0; il <fNumLoadCases; il++) 
        {
            TPZFNMatrix<3,STATE> v2 = bc.Val2(il);
            ef(in,il)+= weight*(v2(0,il)*phiu(0,in) + v2(1,il)*phiu(1,in));
        }
        
        for (int jn = 0; jn <phc; jn++) {
            
            ek(in,jn) += bc.Val1()(0,0)*phiu(0,in)*phiu(0,jn)*weight 
            
            + bc.Val1()(1,0)*phiu(1,in)*phiu(0,jn)*weight
            
            + bc.Val1()(0,1)*phiu(0,in)*phiu(1,jn)*weight
            
            + bc.Val1()(1,1)*phiu(1,in)*phiu(1,jn)*weight;
        }
    }
}

void TPZElastPressure::ApplyDirichlet_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	int c_u = phiu.Cols();
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    int phrp = phip.Rows();
    
//	int r =1;// fk.Rows();
//    int numnod = (ek.Rows()-c_u)/r;
//    
//    for(int in=0 ; in<numnod ; ++in){
//        for(int idf = 0;idf<r;idf++) {
//            (ef)(in*r+idf + c_u,0) += gBigNumber*phip(in,0)*bc.Val2()(2+idf,0)*weight;
//        }
//        for(int jn=0 ; jn<numnod ; ++jn) {
//            for(int idf = 0;idf<r;idf++) {
//                ek(in*r+idf + c_u, jn*r + idf + c_u) += gBigNumber*phip(in,0)*phip(jn,0)*weight;
//            }
//        }
//    }
    
    for(int in = 0; in<phrp; in++){
        ef(in+c_u,0)+=gBigNumber*bc.Val2()(2,0)*phip(in,0)*weight; 
        
        for(int jn = 0; jn<phrp; jn++){
            ek(in+c_u,jn+c_u)+=gBigNumber*phip(in,0)*phip(jn,0)*weight;
        }
    }
}

void TPZElastPressure::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek, TPZFMatrix<> &ef, TPZBndCond &bc){
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	int c_u = phiu.Cols();
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    int  phrp = phip.Rows();
    
//	int r = 1;//fk.Rows();
//    int numnod = (ek.Rows()-c_u)/r;
//    
//    for(int in=0 ; in<numnod ; ++in){
//        for(int idf = 0;idf<r;idf++) {
//            (ef)(in*r+idf+c_u,0) += phip(in,0)*bc.Val2()(2+idf,0)*weight;
//        }
//    }
    
    for(int in=0; in<phrp; in++){
        ef(in+c_u,0) += bc.Val2()(2,0)*phip(in,0)*weight;
    }
}

void TPZElastPressure::ApplyMixed_P(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc){
    
    
    TPZFMatrix<REAL> &phiu = datavec[0].phi;
	int c_u = phiu.Cols();
    TPZFMatrix<REAL> &phip = datavec[1].phi;
    int phrp = phip.Rows();
//    
//	int r =1;// fk.Rows();
//    int numnod = (ek.Rows()-c_u)/r;
//    int in, idf, jn, jdf;
//    for(in=0 ; in<numnod ; ++in){
//        for(idf = 0;idf<r;idf++) {
//            (ef)(in*r+idf,0) += phip(in,0)*bc.Val2()(2+idf,0)*weight;
//        }
//        for(jn=0 ; jn<numnod ; ++jn) {
//            for(idf = 0;idf<r;idf++) {
//                for(jdf = 0;jdf<r;jdf++) {
//                    ek(in*r+idf,jn*r+jdf) += bc.Val1()(2+idf,jdf)*phip(in,0)*phip(jn,0)*weight;
//                }
//            }
//        }
//    }
    
    for(int in=0; in<phrp; in++){
        ef(in+c_u,0)+=bc.Val2()(2,0)*phip(in,0)*weight;
        for(int jn = 0; jn<phrp; jn++){
            ek(in+c_u,jn+c_u) += bc.Val1()(2,0)*phip(in,0)*phip(jn,0)*weight;
        }
    }
}
                                                                                                                                           
void TPZElastPressure::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc){
    
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!=0){
        
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
	bc.Val1().Print();
    bc.Val2().Print();
	switch (bc.Type()) {
		case 0 :			// Dirichlet condition in  both elasticity and pressure equations
			
            ApplyDirichlet_U(datavec, weight, ek,ef,bc);
            ApplyDirichlet_P(datavec, weight, ek,ef,bc);
            break;
			
		case 1 :			// Neumann condition in  both elasticity and pressure equations
            ApplyNeumann_U(datavec, weight,ef,bc);
            ContributePressure(datavec, weight, ek, ef);
            break;
			
		case 2 :            // Mixed condition in  both elasticity and pressure equations
            ApplyMixed_U(datavec, weight, ek,ef,bc);
            ApplyMixed_P(datavec, weight, ek,ef,bc);
            break;
            
        case 10:        // Neumann condition only on the elasticity equation
        {
            // Calculate the matrix contribution for pressure
            ApplyNeumann_U(datavec, weight,ef,bc);
            break;
        }
            
        case 20:        // Mixed condition only on the elasticity equation
             ApplyMixed_U(datavec, weight, ek,ef,bc);
            
        case 21:        // Neumann condition only on the pressure equation pressure
            ApplyNeumann_P(datavec, weight, ek,ef,bc);
            break;
            
        case 22:        // Dirichlet condition only on the pressure equation pressure
            ApplyDirichlet_P(datavec, weight, ek,ef,bc);
            break;
            
    }
}

int TPZElastPressure::VariableIndex(const std::string &name){
    if(!strcmp("Pressure",name.c_str()))        return  1;
    if(!strcmp("MinusKGradP",name.c_str()))     return  2;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZElastPressure::NSolutionVariables(int var){
    if(var == 1) return 1;
    if(var ==  2) return 1;
    return TPZMaterial::NSolutionVariables(var);
}


void TPZElastPressure::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout){
   
    Solout.Resize(this->NSolutionVariables(var));
    
    TPZVec<REAL> SolP;
	TPZFMatrix<> DSolP;
	TPZFMatrix<> axesP;
    SolP=datavec[1].sol[0];
	DSolP=datavec[1].dsol[0];
    axesP=datavec[1].axes;
    
    int ssol = datavec[1].sol[0].size();
    if(ssol==0) return;
    if(var == 1) {
		Solout[0] = SolP[0];
		return;
	}//var1
    
    if (var == 2){
		int id;
		TPZFNMatrix<9,REAL> dsoldx;
		TPZAxesTools<REAL>::Axes2XYZ(DSolP, dsoldx, axesP);
		for(id=0 ; id<1; id++) {
			Solout[id] = -1.*(fk)*dsoldx(id,0);
		}
		return;
	}//var2

    
}

void TPZElastPressure::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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

