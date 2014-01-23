/**
 * @file
 * @brief Contains the methods of the TPZMixedPoisson class (multiphysics environment)
 * @author Agnaldo Farias
 * @date 2012/05/28
 */

#include "pzlog.h"
#include "mixedpoisson.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.mixedpoisson.data"));
#endif

TPZMixedPoisson::TPZMixedPoisson(): TPZMatPoisson3d(), fDim(1) {
	fk = 1.;
    fvisc = 1.;
	ff = 0.;
    fIsStabilized = false;
    fdelta1 = 0.;
    fdelta2 = 0.;
    fUseHdois = false;
    fh2 = 1.;
}

TPZMixedPoisson::TPZMixedPoisson(int matid, int dim): TPZMatPoisson3d(matid,dim), fDim(dim) {
	fk = 1.;
    fvisc = 1.;
	ff = 0.;
    fIsStabilized = false;
    fdelta1 = 0.;
    fdelta2 = 0.;
    fUseHdois = false;
    fh2 = 1.;
}

TPZMixedPoisson::~TPZMixedPoisson() {
}

TPZMixedPoisson::TPZMixedPoisson(const TPZMixedPoisson &cp){
    fk = cp.fk;
    fvisc = cp.fvisc;
    ff = cp.ff;
    fIsStabilized = cp.fIsStabilized;
    fdelta1 = cp.fdelta1;
    fdelta2 = cp.fdelta2;
    fUseHdois = cp.fUseHdois;
    fh2 = cp.fh2;
}

TPZMixedPoisson & TPZMixedPoisson::operator=(const TPZMixedPoisson &copy){
    fk = copy.fk;
    fvisc = copy.fvisc;
    ff = copy.ff;
    fIsStabilized = copy.fIsStabilized;
    fdelta1 = copy.fdelta1;
    fdelta2 = copy.fdelta2;
    fUseHdois = copy.fUseHdois;
    fh2 = copy.fh2;
    return *this;
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

void TPZMixedPoisson::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
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
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    
    REAL &faceSize = datavec[0].HSize;
    if(fUseHdois==true){
        fh2 = faceSize*faceSize;
    }else fh2 = 1.;
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
	
	//Calculate the matrix contribution for flux. Matrix A
    REAL InvK = fvisc/fk;
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fNormalVec(0,ivecind);
        ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
        ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        
        //Inserindo termo de estabilizacao no termo de fonte
        REAL divqi = 0.;
        if(fIsStabilized==true)
        {
            //calculando div(qi)
            TPZFNMatrix<3> axesvec(3,1);
            datavec[0].axes.Multiply(ivec,axesvec);
            for(int iloc=0; iloc<fDim; iloc++)
            {
                divqi += axesvec(iloc,0)*dphiQ(iloc,ishapeind);
            }
            ef(iq, 0) += weight*(fdelta2*divqi*ff);
        }
        
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3> jvec(3,1);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            jvec(0,0) = datavec[0].fNormalVec(0,jvecind);
            jvec(1,0) = datavec[0].fNormalVec(1,jvecind);
            jvec(2,0) = datavec[0].fNormalVec(2,jvecind);
            
            //dot product between u and v
            REAL prod = ivec(0,0)*jvec(0,0) + ivec(1,0)*jvec(1,0) + ivec(2,0)*jvec(2,0);          
            ek(iq,jq) += InvK*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
            
            
            //Inserindo termos de estabilizacao na matriz do fluxo
            if(fIsStabilized==true)
            {
                //termos de delta1
                ek(iq,jq) += (-1.)*weight*fh2*fdelta1*InvK*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
                
                
                //termos de delta2
                REAL divqj = 0.;
                TPZFNMatrix<3> axesvec(3,1);
                datavec[0].axes.Multiply(jvec,axesvec);
                //calculando div(qj)
                for(int jloc=0; jloc<fDim; jloc++)
                {
                    divqj += axesvec(jloc,0)*dphiQ(jloc,jshapeind);
                }
                ek(iq,jq) += weight*fdelta2*divqi*divqj;
            }

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
            
            
            //Inserindo termo de estabilizacao: delta1
            if(fIsStabilized==true)
            {
                REAL dotVGradP = 0.;
                for(int k =0; k<fDim; k++)
                {
                    dotVGradP += ivec(k,0)*phiQ(ishapeind,0)*dphiP(k,jp);
                }
                
                REAL integration = (-1.)*weight*fh2*fdelta1*dotVGradP;
                
                // Estabilizacao delta1 na Matrix B
                ek(iq, phrq+jp) += integration;
                
                // Estabilizacao delta1 na Matrix BˆT
                ek(phrq+jp,iq) += integration;
            }
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
         ef(phrq+ip,0) += (-1.)*weight*ff*phip(ip,0);
    }
    
    //Contribution for estabilization delta1 for gradPu*gradPv. Matrix D
    if(fIsStabilized==true)
    {
        REAL mK = fk/fvisc;
        for(int ip=0; ip<phrp; ip++)
        {
             for(int jp=0; jp<phrp; jp++)
             {
                 for(int k =0; k<fDim; k++)
                 {
                     ek(phrq+ip, phrq+jp) += (-1.)*weight*fh2*fdelta1*mK*dphiP(k,ip)*dphiP(k,jp);
                 }
                 
             }
        }
    }
//
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

void TPZMixedPoisson::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef DEBUG
    int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
		DebugStop();
	}
	if (bc.Type() > 2 ) {
        std::cout << " Erro.!! Neste material utiliza-se apenas condicoes de Neumann e Dirichlet\n";
		DebugStop();
	}
#endif
	
	TPZFMatrix<REAL>  &phiQ = datavec[0].phi;
	int phrq = phiQ.Rows();

	REAL v2;
    if(bc.HasForcingFunction())
    {
		TPZManVector<STATE> res(3);
		bc.ForcingFunction()->Execute(datavec[0].x,res);
		v2 = res[0];
	}else
    {
        v2 = bc.Val2()(0,0);
    }

	switch (bc.Type()) {
		case 0 :		// Dirichlet condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                //the contribution of the Dirichlet boundary condition appears in the flow equation
                ef(iq,0) += (-1.)*v2*phiQ(iq,0)*weight;
            }
            break;
			
		case 1 :			// Neumann condition
			//primeira equacao
			for(int iq=0; iq<phrq; iq++)
            {
                ef(iq,0)+= gBigNumber*v2*phiQ(iq,0)*weight;
                for (int jq=0; jq<phrq; jq++) {
                    
                    ek(iq,jq)+= gBigNumber*phiQ(iq,0)*phiQ(jq,0)*weight; 
                }
            }  
			break;
        
        case 2 :			// mixed condition
            for(int iq = 0; iq < phrq; iq++) {
                
				ef(iq,0) += v2*phiQ(iq,0)*weight;
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
	if(!strcmp("Pressure",name.c_str()))    return  2;
    if(!strcmp("GradFluxX",name.c_str()))   return  3;
    if(!strcmp("GradFluxY",name.c_str()))   return  4;
    if(!strcmp("DivFlux",name.c_str()))   return  5;
    
    if(!strcmp("ExactPressure",name.c_str()))  return 6;
    if(!strcmp("ExactFlux",name.c_str()))  return 7;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMixedPoisson::NSolutionVariables(int var){
	if(var == 1) return fDim;
	if(var == 2) return 1;
    if(var == 3) return 3;
    if(var == 4) return 3;
    if(var == 5) return 1;
    if(var == 6) return 1;
    if(var == 7) return fDim;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMixedPoisson::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	
	TPZVec<STATE> SolP, SolQ;
    
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
    
    if(var==3){
        Solout[0]=datavec[0].dsol[0](0,0);
        Solout[1]=datavec[0].dsol[0](1,0);
        Solout[2]=datavec[0].dsol[0](2,0);
        return;
    }

    if(var==4){
        Solout[0]=datavec[0].dsol[0](0,1);
        Solout[1]=datavec[0].dsol[0](1,1);
        Solout[2]=datavec[0].dsol[0](2,1);
        return;
    }
    
    if(var==5){
        Solout[0]=datavec[0].dsol[0](0,0)+datavec[0].dsol[0](1,1);
        return;
    }
    
    TPZVec<REAL> ptx(3);
	TPZVec<STATE> solExata(1);
	TPZFMatrix<STATE> flux(2,1);
    
    //Exact soluion
	if(var == 6){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = solExata[0];
		return;
	}//var6
    
    if(var == 7){
		fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
		Solout[0] = flux(0,0);
        Solout[1] = flux(1,0);
		return;
	}//var7

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




