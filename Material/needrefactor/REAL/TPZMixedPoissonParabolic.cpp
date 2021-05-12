/**
 * @file
 * @brief Contains the methods of the TPZMixedPoissonParabolic class (multiphysics environment)
 * @author Philippe Devloo
 * @date 2017/09/04
 */

#include "TPZMixedPoissonParabolic.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef PZ_LOG
static TPZLogger logdata("pz.mixedpoisson.data");
static TPZLogger logerror("pz.mixedpoisson.error");
#endif

TPZMixedPoissonParabolic::TPZMixedPoissonParabolic(): TPZMixedPoisson(), fRho(1.), fDeltaT(1.)
{
}

TPZMixedPoissonParabolic::TPZMixedPoissonParabolic(int matid, int dim): TPZMixedPoisson(matid,dim), fRho(1.), fDeltaT(1.)
{
}

TPZMixedPoissonParabolic::~TPZMixedPoissonParabolic() {
}

TPZMixedPoissonParabolic::TPZMixedPoissonParabolic(const TPZMixedPoissonParabolic &cp) : TPZMixedPoisson(cp), fRho(cp.fRho), fDeltaT(cp.fDeltaT)
{
}

TPZMixedPoissonParabolic & TPZMixedPoissonParabolic::operator=(const TPZMixedPoissonParabolic &copy){
    TPZMixedPoisson::operator=(copy);
    fRho = copy.fRho;
    fDeltaT = copy.fDeltaT;
    return *this;
} 


void TPZMixedPoissonParabolic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    TPZMixedPoisson::Print(out);
    out << "rho = " << fRho << std::endl;
    out << "delta t = " << fDeltaT << fDeltaT << std::endl;
}

void TPZMixedPoissonParabolic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    
#ifdef PZDEBUG
	int nref =  datavec.size();
	if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
		DebugStop();
	}
    if(fIsStabilized || fUseHdois)
    {
        DebugStop();
    }
#endif
    
    STATE force = ff;
    if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(datavec[1].x,res);
		force = res[0];
	}
    
    STATE pressureN = datavec[1].sol[0][0];
    
    TPZFNMatrix<9,STATE> PermTensor;
    this->GetPermeability(PermTensor);
    TPZFNMatrix<9,STATE> InvPermTensor;
    this->GetInvPermeability(InvPermTensor);
    //int rtens = 2*fDim;
    if(fPermeabilityFunction){
        PermTensor.Zero();
        InvPermTensor.Zero();
        TPZFNMatrix<18,STATE> resultMat(2*fDim,fDim,0.);
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        
        for(int id=0; id<fDim; id++){
            for(int jd=0; jd<fDim; jd++){
                
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+fDim,jd);
            }
        }
    }
    
    // Setting the phis
    TPZFMatrix<REAL> &phiQ = datavec[0].phi;
    TPZFMatrix<REAL> &phip = datavec[1].phi;
	TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
	
#ifdef PZDEBUG
    if(phrp+phrq != ek.Rows())
    {
        DebugStop();
    }
#endif
	//Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        //ef(iq, 0) += 0.;
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
            //ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
            //ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        }
        
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<3; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            //jvecZ.Print("mat1 = ");
            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
        }
    }
    
	
	// Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<3; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }
        TPZFNMatrix<3,REAL> axesvec(3,1,0.);
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
    
    for (int ip=0; ip<phrp; ip++) {
        for (int jp=0; jp < phrp; jp++) {
            ek(phrq+ip,phrq+jp) += (-1.)*weight*phip(ip,0)*phip(jp,0)*fRho/fDeltaT;
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*fRho*pressureN/fDeltaT*phip(ip,0);
    }
    
    
}

void TPZMixedPoissonParabolic::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) {


#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    if(fIsStabilized || fUseHdois)
    {
        DebugStop();
    }
#endif
    TPZFMatrix<REAL> &phip = datavec[1].phi;

    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();

    STATE pressureN = datavec[1].sol[0][0];

    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*fRho*pressureN/fDeltaT*phip(ip,0);
    }
    

}
void TPZMixedPoissonParabolic::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
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
        TPZFNMatrix<9,STATE> gradu(Dimension(),1);
		bc.ForcingFunction()->Execute(datavec[0].x,res,gradu);
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


void TPZMixedPoissonParabolic::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(false);
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = false;
        datavec[i].fNeedsHSize = false;
	}
    datavec[1].fNeedsSol = true;
}


