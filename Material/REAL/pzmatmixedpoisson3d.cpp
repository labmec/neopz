//
//  pzmatmixedpoisson3d.cpp
//  PZ
//
//  Created by Agnaldo Farias on 03/11/14.
//
//

#include "pzmatmixedpoisson3d.h"
#include "pzlog.h"
#include "pzfmatrix.h"
#include "pzmaterialdata.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzlog.h"

#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.TPZMatMixedPoisson3D.data"));
#endif

using namespace std;

TPZMatMixedPoisson3D::TPZMatMixedPoisson3D():TPZMaterial(){
    
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = 2;
    
    /** Material id not initialized */
    fMatId = -1;
    
    /** Coeficiente que multiplica o gradiente */
    fK = 1.;
    
    fvisc = 1.;
    fInvK.Resize(1, 1);
    fTensorK.Resize(1, 1);
    fTensorK.Identity();
    fInvK.Identity();
    fPermeabilityFunction = NULL;
    fmatLagr = new TPZLagrangeMultiplier();
    
}

TPZMatMixedPoisson3D::TPZMatMixedPoisson3D(int matid, int dim):TPZMaterial(matid){
    
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
    
    fvisc = 1.;
    
    fInvK.Redim(dim, dim);
    fTensorK.Resize(dim, dim);
    fInvK.Identity();
    fTensorK.Identity();
    fPermeabilityFunction = NULL;
    fmatLagr =  new TPZLagrangeMultiplier();
}

TPZMatMixedPoisson3D::~TPZMatMixedPoisson3D(){
}

TPZMatMixedPoisson3D::TPZMatMixedPoisson3D(const TPZMatMixedPoisson3D &copy):TPZMaterial(copy){
    
    this->operator=(copy);
}

TPZMatMixedPoisson3D & TPZMatMixedPoisson3D::operator=(const TPZMatMixedPoisson3D &copy){
    
    TPZMaterial::operator = (copy);
    this->fF = copy.fF; //fF
    this->fDim = copy.fDim;
    this->fMatId = copy.fMatId;
    this->fK = copy.fK;
    this->fmatLagr = copy.fmatLagr;
    
    return *this;
}

int TPZMatMixedPoisson3D::NStateVariables() {
    return 1;//(1+fDim);
}

void TPZMatMixedPoisson3D::Print(std::ostream &out) {
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
// esse metodo esta ok
void TPZMatMixedPoisson3D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
#ifdef DEBUG
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    REAL force = fF;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    
    // Setting the phis
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phip =  datavec[1].phi;
    TPZFMatrix<REAL> &dphipLoc = datavec[1].dphix;
    
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    //Matrix B: the contribution of skeletal elements is done by TPZLagrangeMultiplier material
    if(phrq==0) {
        int nr = phiQ.Rows();
        if(nr<1) DebugStop();
        REAL multiplier = -1.;
        fmatLagr->SetMultiplier(multiplier);
        fmatLagr->Contribute(datavec, weight, ek, ef);
        return;
    }
    
    //Calculate the matrix contribution for flux. Matrix A
    // A matriz de rigidez é tal que A{ij}=\int_\Omega K^{-1} \varphi_j\cdot\varphi_i d\Omega
    // K, futuramente sera uma matriz ou funcao, deve-se ter cuidado com essa parte da inversao de K
    
    TPZFNMatrix<3,REAL> PermTensor = fTensorK;
    TPZFNMatrix<3,REAL> InvPermTensor = fInvK;
    
    if(fPermeabilityFunction){
        PermTensor.Redim(fDim,fDim);
        InvPermTensor.Redim(fDim,fDim);
        TPZFNMatrix<3,REAL> resultMat;
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(datavec[1].x,res,resultMat);
        
        for(int id=0; id<fDim; id++){
            for(int jd=0; jd<fDim; jd++){
                
                PermTensor(id,jd) = resultMat(id,jd);
                InvPermTensor(id,jd) = resultMat(id+fDim,jd);
            }
        }
    }
    
    
    //REAL InvK = 1./fK;
    for (int iq = 0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fNormalVec(0,ivecind);
        ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
        ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        
        TPZFNMatrix<3,REAL> jvecZ(fDim,1,0.);
        
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3> jvec(3,1);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            jvec(0,0) = datavec[0].fNormalVec(0,jvecind);
            jvec(1,0) = datavec[0].fNormalVec(1,jvecind);
            jvec(2,0) = datavec[0].fNormalVec(2,jvecind);
            
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<fDim; id++){
                for(int jd=0; jd<fDim; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            REAL prod = 0.;
            for(int id=0; id < fDim;id++) prod += ivec(id,0)*jvecZ(id,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod;
            
        }
    }
    
    // Coupling terms between flux and pressure inside the element. Matrix B
    // A matriz de rigidez é tal que B{ij}=\int_\Omega \nabla\phi_j\cdot\varphi_i d\Omega
    
    TPZFNMatrix<200,REAL> dphip(3,datavec[1].dphix.Cols(),0.0);
    for (int ip = 0; ip<dphip.Cols(); ip++) {
        for (int d = 0; d<dphipLoc.Rows(); d++) {
            for (int j=0; j< 3; j++) {
                dphip(j,ip)+=datavec[1].axes(d,j)*dphipLoc(d,ip);
            }
        }
    }
    
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3> ivec(3,1);
        ivec(0,0) = datavec[0].fNormalVec(0,ivecind);
        ivec(1,0) = datavec[0].fNormalVec(1,ivecind);
        ivec(2,0) = datavec[0].fNormalVec(2,ivecind);
        
        for (int jp=0; jp<phrp; jp++)
        {
            //dot product between  varphi and grad phi
            REAL prod = 0.;
            for(int iloc=0; iloc<3; iloc++)
            {
                prod += (ivec(iloc,0)*phiQ(ishapeind,0))*dphip(iloc,jp);
            }
            
            REAL fact = weight*prod;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
            
        }
    }
    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
}


void TPZMatMixedPoisson3D::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef DEBUG
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 2 \n";
        DebugStop();
    }
//    if (!bc.Type() ) {
//        std::cout << " Erro.!!\n";
//        DebugStop();
//    }
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
            
//        case 3 :          //Contribution of skeletal elements.
//        {
//            REAL multiplier = -1.;
//            fmatLagr->SetMultiplier(multiplier);
//            fmatLagr->Contribute(datavec, weight, ek, ef);
//            break;
//        }
        default:
        {
            std::cout << "Boundary condition not found!";
            DebugStop();
            break;
        }


    }
    
}

void TPZMatMixedPoisson3D::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = false;
        datavec[i].fNeedsNormal = false;
    }
}

/** Returns the variable index associated with the name */
int TPZMatMixedPoisson3D::VariableIndex(const std::string &name){
    if(!strcmp("Flux",name.c_str()))           return 1;
    if(!strcmp("Pressure",name.c_str()))       return 2;
    if(!strcmp("GradFluxX",name.c_str()))      return 3;
    if(!strcmp("GradFluxY",name.c_str()))      return 4;
    if(!strcmp("GradFluxZ",name.c_str()))      return 5;
    if(!strcmp("ExactPressure",name.c_str()))  return 6;
    if(!strcmp("ExactFlux",name.c_str()))      return 7;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZMatMixedPoisson3D::NSolutionVariables(int var){
    if(var == 1) return fDim;
    if(var == 2) return 1;
    if(var == 3) return 3;
    if(var == 4) return 3;
    if(var == 5) return 3;
    if(var == 6) return 1;
    if(var == 7) return fDim;
    return TPZMaterial::NSolutionVariables(var);
}


void TPZMatMixedPoisson3D::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZVec<STATE> SolP, SolQ;
    
    // SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 1){ //function (state variable Q)
        for (int ip = 0; ip<Dimension(); ip++)
        {
            Solout[ip] = datavec[0].sol[0][ip];
        }
        return;
    }
    
    if(var == 2){
        Solout[0] = SolP[0];//function (state variable p)
        return;
    }
    
    if(var==3){
        for (int ip = 0; ip<Dimension(); ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,0);
        }
        return;
    }
    
    if(var==4){
        for (int ip = 0; ip<Dimension(); ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,1);
        }
        return;
    }
    
    if(var==5){
        for (int ip = 0; ip<Dimension(); ip++)
        {
            Solout[ip] = datavec[0].dsol[0](ip,2);
        }
        return;
    }
    
    TPZVec<REAL> ptx(3);
    TPZVec<STATE> solExata(1);
    TPZFMatrix<STATE> flux(fDim,1);
    
    //Exact soluion
    if(var == 6){
        fForcingFunctionExact->Execute(datavec[1].x, solExata,flux);
        Solout[0] = solExata[0];
        return;
    }//var6
    
    if(var == 7){
        fForcingFunctionExact->Execute(datavec[0].x, solExata,flux);
        for (int ip = 0; ip<Dimension(); ip++)
        {
            Solout[ip] = flux(ip,0);
        }
        return;
    }//var7
    
}


void TPZMatMixedPoisson3D::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = false;
        datavec[i].fNeedsNeighborSol = false;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = false;
    }
}
