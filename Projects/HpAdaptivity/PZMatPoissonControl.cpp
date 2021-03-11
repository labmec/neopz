/**
 * @file
 * @brief Contains implementations of the TPZMatPoisson3d methods.
 */

#include "pzpoisson3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "PZMatPoissonControl.h"
#include <cmath>
#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.tpzmatpoissonControl.HpAdapitivity"));
#endif


//
//  PZMatPoissonD3.cpp
//  PZ
//
//  Created by Douglas Castro on 5/23/14.
//
//

#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.TPZMatPoissonControl.data"));
#endif

using namespace std;

TPZMatPoissonControl::TPZMatPoissonControl():TPZMaterial(){
	
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = 2;
    
    /** Material id not initialized */
    fMatId = -1;
    
    /** Coeficiente que multiplica o gradiente */
    fK = 1.;
    
    falpha = 1.;
    
    
}

TPZMatPoissonControl::TPZMatPoissonControl(int matid, int dim):TPZMaterial(matid){
	
    if(dim<0 || dim >3){
        DebugStop();
    }
    
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = dim;
    
    /** Material id no initialized */
    fMatId = matid;
    
    /** Coeficiente que multiplica o gradiente */
    fK = 1.;
    
    falpha = 1.;
}

TPZMatPoissonControl::~TPZMatPoissonControl(){
}

TPZMatPoissonControl::TPZMatPoissonControl(const TPZMatPoissonControl &copy):TPZMaterial(copy){
    
    this->operator=(copy);
}

TPZMatPoissonControl & TPZMatPoissonControl::operator=(const TPZMatPoissonControl &copy){
    
    TPZMaterial::operator = (copy);
    this->fF = copy.fF; //fF
    this->fDim = copy.fDim;
    this->fMatId = copy.fMatId;
    this->fK = copy.fK;
    this->falpha = copy.falpha;
    
	return *this;
}

void TPZMatPoissonControl::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    out << "Dimesion of problem " << fDim << endl;
	out << "Material ID  "<< fMatId << endl;
    out << "Forcing function  "<< fF << endl;
    out << "Grad Coeficient  "<< fK << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}


// Contribute methods
void TPZMatPoissonControl::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	
    DebugStop();
}


void TPZMatPoissonControl::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
    /** monta a matriz
     |Ch   0    Ah^T |    |f1|
     |0 alphaCh -Ch^T| =  |0 |
     |Ah   -Ch     0 |    |f2|
     
     **/
    
#ifdef PZDEBUG
	int nref =  datavec.size();
	if (nref != 3) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
#endif
    
    REAL fXfLoc = fF;
    REAL fXfadLoc = fFad;
    
    if(fForcingFunction) {
		TPZManVector<STATE> res(2);
		fForcingFunction->Execute(datavec[0].x,res);
		fXfLoc = res[0];
        fXfadLoc = res[1];
	}
    
    // Setting the phis
    TPZFMatrix<REAL>  &phiY =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiU =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiP =  datavec[2].phi;
    
    
	TPZFMatrix<REAL> &dphiY = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[2].dphix;
    
    TPZFNMatrix<200,REAL> dphip(3,datavec[1].dphix.Cols(),0.0);
    
    
    //Numeros de linhas: rows
    int rY, rU, rP;
    rY = phiY.Rows();
    rU = phiU.Rows();
    rP = phiP.Rows();
    
    int in, jn, kd;
    
    for(in = 0; in < rY; in++ ){
        
        ef(in, 0) += weight*fXfadLoc*phiY(in,0);
        
        //BOLCO: Ch
        for(jn = 0; jn < rY; jn++ ){
            ek(in,jn) = weight*phiY(in,0)*phiY(jn,0);
        }
        
        
        for(jn = 0; jn < rP; jn++)
        {
            for(kd=0; kd<fDim; kd++) {
                
                //BLOCO: Ah^T
                ek(in,jn+rY+rU) += weight*dphiY(kd,in)*dphiP(kd,jn);
                
                //BLOCO: Ah
                ek(in+rY+rU,jn) += weight*dphiY(kd,in)*dphiP(kd,jn);
            }
        }
    }
    
    for(in = 0; in < rU; in++ ){
        
        //BOLCO: alpha*Ch
        for(jn = 0; jn < rU; jn++){
            
            ek(in+rY,jn+rY) = falpha*weight*phiU(in,0)*phiU(jn,0);
        }
        
        for(jn = 0; jn < rP; jn++){
            
            
            ef(in+rY+rU, 0) += weight*fXfLoc*phiY(in,0);
            
            
            //BLOCO: -Ch^T
            ek(in+rY,jn+rY+rU) = weight*phiU(in,0)*phiP(jn,0);
            
            //BLOCO: -Ch
            ek(in+rY+rU,jn+rY) = weight*phiU(in,0)*phiP(jn,0);
        }
    }
    
}

void TPZMatPoissonControl::ContributeBC(TPZVec<TPZMaterialData> &datavec,REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    int nref =  datavec.size();
    if (nref != 3) {
        cout << " Erro.!! datavec tem que ser de tamanho 3 \n";
        DebugStop();
    }
    
    TPZFMatrix<REAL>  &phiY = datavec[0].phi;
    TPZFMatrix<REAL>  &phiU = datavec[1].phi;
    TPZFMatrix<REAL>  &phiP = datavec[2].phi;
    
    
    int rY, rU, rP;
    rY = phiY.Rows();
    rU = phiU.Rows();
    rP = phiP.Rows();
    
    short in,jn;
    REAL v2[2];
    v2[0] = bc.Val2()(0,0); //condicao de contorno da primeira equacao
    v2[1] = bc.Val2()(1,0); //condicao de contorno da segunda equaca
    
    switch (bc.Type()) {
        case 0 : // Dirichlet condition
            
            //primeira equacao
            for(in = 0 ; in < rY; in++) {
                
              //ef(in+rY+rU,0) += gBigNumber * v2[0]*phiP(in,0)*weight;
                ef(in,0) += gBigNumber * v2[0]*phiY(in,0)*weight;
                
                
                for (jn = 0 ; jn < rY; jn++) {
                   // ek(in+rY+rU,jn) += gBigNumber*phiP(in,0)*phiY(jn,0)*weight;
                    ek(in,jn) += gBigNumber*phiY(in,0)*phiY(jn,0)*weight;
                }
                
            }
            break;
            
    }
    
}


//Contribute interface methods
void TPZMatPoissonControl::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatPoissonControl::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    DebugStop();
}


void TPZMatPoissonControl::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}


void TPZMatPoissonControl::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}

void TPZMatPoissonControl::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
	}
}

/** Returns the variable index associated with the name */
int TPZMatPoissonControl::VariableIndex(const std::string &name){
	if(!strcmp("State",name.c_str()))           return 1;
	if(!strcmp("Control",name.c_str()))       return 2;
    if(!strcmp("LagrangeMult",name.c_str()))      return 3;
    if(!strcmp("ExactState",name.c_str()))      return 4;
    return TPZMaterial::VariableIndex(name);
}

int TPZMatPoissonControl::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 1;
    if(var == 3) return 1;
    if(var == 4) return 1;
    
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatPoissonControl::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> SolY, SolU, SolP;
    
    SolY = datavec[0].sol[0];
    SolU = datavec[1].sol[0];
    SolP = datavec[2].sol[0];
    
    if(var == 1){ //function (state variable Y)
        Solout[0] = SolY[0];
        
        return;
	}
    
    if(var == 2){
		Solout[0] = SolU[0];//function (state variable U)
		return;
	}
    
    if(var == 3){
        Solout[0] = SolP[0];//function (state variable P)
        return;
    }
    
    
    TPZVec<REAL> ptx(3);
	TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(fDim,1);
    
    //Exact soluion
	if(var == 4){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}
}


void TPZMatPoissonControl::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
		datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
}
