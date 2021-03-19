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

#ifdef PZ_LOG
static TPZLogger logdata("pz.TPZMatMixedPoisson3D.data");
#endif

using namespace std;

TPZMatMixedPoisson3D::
TPZMatMixedPoisson3D():
TPZRegisterClassId(&TPZMatMixedPoisson3D::ClassId),
TPZMaterial(){
    
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = 2;
    
    /** Material id not initialized */
    fMatId = -1;
    
    falpha = 0.;
    fvisc = 1.;
    fInvK.Resize(1, 1);
    fTensorK.Resize(1, 1);
    fTensorK.Identity();
    fInvK.Identity();
    fPermeabilityFunction = NULL;
    fReactionTermFunction = NULL;
    fmatLagr = new TPZLagrangeMultiplier();
    fSecondIntegration = false;
    fReactionTerm = false;
    
}

TPZMatMixedPoisson3D::TPZMatMixedPoisson3D(int matid, int dim):
TPZRegisterClassId(&TPZMatMixedPoisson3D::ClassId),
TPZMaterial(matid){
    
    
    /** Valor da funcao de carga */
    fF = 0.; //fF
    
    /** Dimensao do dominio */
    fDim = dim;
    
    /** Material id no initialized */
    fMatId = matid;
    
    fvisc = 1.;
    falpha = 0.;
    fInvK.Redim(dim, dim);
    fTensorK.Resize(dim, dim);
    fInvK.Identity();
    fTensorK.Identity();
    fPermeabilityFunction = NULL;
    fReactionTermFunction = NULL;
    fmatLagr =  new TPZLagrangeMultiplier();
    fSecondIntegration = false;
    fReactionTerm = false;
}

TPZMatMixedPoisson3D::~TPZMatMixedPoisson3D(){
}

TPZMatMixedPoisson3D::TPZMatMixedPoisson3D(const TPZMatMixedPoisson3D &copy):TPZRegisterClassId(&TPZMatMixedPoisson3D::ClassId),
TPZMaterial(copy){
    
    this->operator=(copy);
    this->fF = copy.fF;
    this->falpha = copy.falpha;
    this->fDim = copy.fDim;
    this->fMatId = copy.fMatId;
    this->fmatLagr = copy.fmatLagr;
    this->fTensorK = copy.fTensorK;
    this->fInvK = copy.fInvK;
    this->fSecondIntegration = copy.fSecondIntegration;
    this->fReactionTerm = copy.fReactionTerm;
}

TPZMatMixedPoisson3D & TPZMatMixedPoisson3D::operator=(const TPZMatMixedPoisson3D &copy){
    
    TPZMaterial::operator = (copy);
    this->fF = copy.fF;
    this->falpha = copy.falpha;
    this->fDim = copy.fDim;
    this->fMatId = copy.fMatId;
    this->fmatLagr = copy.fmatLagr;
    this->fTensorK = copy.fTensorK;
    this->fInvK = copy.fInvK;
    this->fSecondIntegration = copy.fSecondIntegration;
    this->fReactionTerm = copy.fReactionTerm;
    
    return *this;
}

int TPZMatMixedPoisson3D::NStateVariables() const {
    return 1;
}

void TPZMatMixedPoisson3D::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Dimesion of problem " << fDim << endl;
    out << "Material ID  "<< fMatId << endl;
    out << "Forcing function  "<< fF << endl;
    out << "Base Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}


// Contribute methods
// esse metodo esta ok
//This method use normalized piola contravariant mapping for nonlinear mappings. With second integration by parts
void TPZMatMixedPoisson3D::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    int phrq = datavec[0].fVecShapeIndex.NElements();
    
    if(!fSecondIntegration){
        
        if(/*HDivPiola==0 ||*/ phrq == 0){
            std::cout << "\n Error.\n";
            std::cout << "If you want to use the variational formulation with second integration by parts, you need to set fSecondIntegration==true\n";
            DebugStop();
        }
        
        ContributeWithoutSecondIntegration(datavec,weight,ek,ef);
        return;
    }
    
#ifdef PZDEBUG
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
    
    TPZFNMatrix<660> GradofPhiL2;
    TPZAxesTools<REAL>::Axes2XYZ(datavec[1].dphix, GradofPhiL2, datavec[1].axes);
    
    
    int phrp = phip.Rows();
    
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
//        PermTensor.Redim(fDim,fDim);
//        InvPermTensor.Redim(fDim,fDim);
        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(datavec[0].x,res,resultMat);
        
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
        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3> jvec(3,1);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            jvec(0,0) = datavec[0].fDeformedDirections(0,jvecind);
            jvec(1,0) = datavec[0].fDeformedDirections(1,jvecind);
            jvec(2,0) = datavec[0].fDeformedDirections(2,jvecind);
            
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<3; id++){
                for(int jd=0; jd<3; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }
            REAL prod = 0.;
            for(int id=0; id < 3;id++) prod += ivec(id,0)*jvecZ(id,0);
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
        ivec(0,0) = datavec[0].fDeformedDirections(0,ivecind);
        ivec(1,0) = datavec[0].fDeformedDirections(1,ivecind);
        ivec(2,0) = datavec[0].fDeformedDirections(2,ivecind);
        
        for (int jp=0; jp<phrp; jp++)
        {
            //dot product between  varphi and grad phi
            REAL prod = 0.;
            for(int iloc=0; iloc<3; iloc++)
            {
                prod += (ivec(iloc,0)*phiQ(ishapeind,0))*GradofPhiL2(iloc,jp);
            }
            
            REAL fact = weight*prod;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
            
        }
    }
    
    if(fReactionTerm)
    {
        REAL alpha = falpha;
        
        if(fReactionTermFunction)
        {
            TPZFNMatrix<3,STATE> resultMat;
            TPZManVector<STATE> res;
            fReactionTermFunction->Execute(datavec[1].x,res,resultMat);
            
            alpha = res[0];
        }
        
        for(int ip=0; ip<phrp; ip++)
        {
            for (int jp=0; jp<phrp; jp++)
            {
                REAL fact = (-1.)*weight*phip(ip,0)*phip(jp,0);
                // Matrix C
                ek(phrq+ip, phrq+jp) += alpha*fact;
            }
        }
    }

    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
}


///This method use piola contravariant mapping for nonlinear mappings
void TPZMatMixedPoisson3D::ContributeWithoutSecondIntegration(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    STATE force = fF;
    if(fForcingFunction) {
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(datavec[1].x,res);
        force = res[0];
    }
    
    TPZFNMatrix<9,REAL> PermTensor = fTensorK;
    TPZFNMatrix<9,REAL> InvPermTensor = fInvK;
    //int rtens = 2*fDim;
    if(fPermeabilityFunction)
    {

        TPZFNMatrix<3,STATE> resultMat;
        TPZManVector<STATE> res;
        fPermeabilityFunction->Execute(datavec[0].x,res,resultMat);
        
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
    TPZFNMatrix<40, REAL> divphi = datavec[0].divphi;
    
    
    int phrq, phrp;
    phrp = phip.Rows();
    phrq = datavec[0].fVecShapeIndex.NElements();
    
    //Calculate the matrix contribution for flux. Matrix A
    for(int iq=0; iq<phrq; iq++)
    {
        int ivecind = datavec[0].fVecShapeIndex[iq].first;
        int ishapeind = datavec[0].fVecShapeIndex[iq].second;
        TPZFNMatrix<3,REAL> ivec(3,1,0.);
        for(int id=0; id<fDim; id++){
            ivec(id,0) = datavec[0].fDeformedDirections(id,ivecind);
        }
        
        TPZFNMatrix<3,REAL> ivecZ(3,1,0.);
        TPZFNMatrix<3,REAL> jvecZ(3,1,0.);
        for (int jq=0; jq<phrq; jq++)
        {
            TPZFNMatrix<3,REAL> jvec(3,1,0.);
            int jvecind = datavec[0].fVecShapeIndex[jq].first;
            int jshapeind = datavec[0].fVecShapeIndex[jq].second;
            
            for(int id=0; id<fDim; id++){
                jvec(id,0) = datavec[0].fDeformedDirections(id,jvecind);
            }
            
            //dot product between Kinv[u]v
            jvecZ.Zero();
            for(int id=0; id<fDim; id++){
                for(int jd=0; jd<fDim; jd++){
                    jvecZ(id,0) += InvPermTensor(id,jd)*jvec(jd,0);
                }
            }

            REAL prod1 = ivec(0,0)*jvecZ(0,0) + ivec(1,0)*jvecZ(1,0) + ivec(2,0)*jvecZ(2,0);
            ek(iq,jq) += fvisc*weight*phiQ(ishapeind,0)*phiQ(jshapeind,0)*prod1;
        }
    }
    
    // Coupling terms between flux and pressure. Matrix B
    for(int iq=0; iq<phrq; iq++)
    {

        REAL divwq = divphi(iq,0)/datavec[0].detjac;
        
        for (int jp=0; jp<phrp; jp++) {
            
            REAL fact = (-1.)*weight*phip(jp,0)*divwq;
            // Matrix B
            ek(iq, phrq+jp) += fact;
            
            // Matrix B^T
            ek(phrq+jp,iq) += fact;
        }
    }
    
    //Calculate the matrix contribution for pressure. Matrix C: term: alpha*p
    if(fReactionTerm)
    {
        REAL alpha = falpha;
        if(fReactionTermFunction)
        {
            TPZFNMatrix<3,STATE> resultMat;
            TPZManVector<STATE> res;
            fReactionTermFunction->Execute(datavec[1].x,res,resultMat);
            
            alpha = res[0];
        }

        for(int ip=0; ip<phrp; ip++)
        {
            for (int jp=0; jp<phrp; jp++)
            {
                
                REAL fact = (-1.)*weight*phip(ip,0)*phip(jp,0);
                // Matrix C
                ek(phrq+ip, phrq+jp) += alpha*fact;
            }
        }
    }

    
    //termo fonte referente a equacao da pressao
    for(int ip=0; ip<phrp; ip++){
        ef(phrq+ip,0) += (-1.)*weight*force*phip(ip,0);
    }
    
    //
    //#ifdef PZ_LOG
    //    if(logdata.isDebugEnabled())
    //	{
    //        std::stringstream sout;
    //        sout<<"\n\n Matriz ek e vetor fk \n ";
    //        ek.Print("ekmph = ",sout,EMathematicaInput);
    //        ef.Print("efmph = ",sout,EMathematicaInput);
    //        LOGPZ_DEBUG(logdata,sout.str());
    //	}
    //#endif
    
}

void TPZMatMixedPoisson3D::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
#ifdef PZDEBUG
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
    if(!strcmp("Rhs",name.c_str()))            return 8;
    if(!strcmp("GradP",name.c_str()))          return 9;
    if(!strcmp("Divergence",name.c_str()))      return  10;
    if(!strcmp("ExactDiv",name.c_str()))        return  11;
    if(!strcmp("POrder",name.c_str()))        return  12;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZMatMixedPoisson3D::NSolutionVariables(int var){
    if(var == 1) return fDim;
    if(var == 2 || var == 12) return 1;
    if(var == 3) return 3;
    if(var == 4) return 3;
    if(var == 5) return 3;
    if(var == 6) return 1;
    if(var == 7) return fDim;
    if(var == 8 || var == 10 || var == 11) return 1;
    if(var == 9) return fDim;
    return TPZMaterial::NSolutionVariables(var);
}

// metodo para gerar vtk
void TPZMatMixedPoisson3D::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    
    TPZVec<STATE> SolP, SolQ;
    
   //  SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    
    if(var == 1){ //function (state variable Q)
        for (int ip = 0; ip<fDim; ip++)
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
    TPZFMatrix<STATE> flux(3,1);
    
    //Exact soluion
    if(var == 6){
        fExactSol->Execute(datavec[1].x, solExata,flux);
        Solout[0] = solExata[0];
        return;
    }//var6
    
    if(var == 7){
        fExactSol->Execute(datavec[0].x, solExata,flux);
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = flux(ip,0);
        }
        return;
    }//var7
    
    if(var == 8){
        REAL force = fF;
        if(fForcingFunction) {
            TPZManVector<STATE> res(1);
            fForcingFunction->Execute(datavec[1].x,res);
            force = res[0];
        }
        Solout[0] = force;
        
        return;
    }//var8
    
    if(var==9){
        
        TPZFNMatrix<660,STATE> GradofP;
        TPZAxesTools<STATE>::Axes2XYZ(datavec[1].dsol[0], GradofP, datavec[1].axes);
        //        int nc = GradofP.Cols();
        //        int nl = GradofP.Rows();
        
        for (int ip = 0; ip<fDim; ip++)
        {
            Solout[ip] = -1.0*GradofP(ip,0);
        }
        return;
    }
    
    if(var==10){
        REAL val = 0.;
        for(int i=0; i<fDim; i++){
            val += datavec[0].dsol[0](i,i);
        }
        Solout[0] = val;
        return;
    }
    
    if(var==11){
        fExactSol->Execute(datavec[0].x,solExata,flux);
        Solout[0]=flux(fDim,0);
        return;
    }
    
    if(var==12){
        Solout[0] = datavec[1].p;
        return;
    }
}

// metodo para computar erros
void TPZMatMixedPoisson3D::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( 3 /*this->NSolutionVariables(var)*/);
    // AQUI!!! //redefinicao feita  acima, antigamente mudava para 2, por exemplo, e nao ficava compativel com o resto que era 3
    
    if(var == 1){ //function (state variable Q)
        for (int ip = 0; ip<3; ip++)
        {
            Solout[ip] = data.sol[0][ip];
        }
        
        return;
    }
    
    if(var == 2){ //function (state variable p)
        
        TPZVec<STATE> SolP;
        SolP = data.sol[0];
        
        Solout[0] = SolP[0];
        return;
    }
    
    
}

#include "pzaxestools.h"
void TPZMatMixedPoisson3D::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes, int var,TPZVec<STATE> &Solout){
    
#ifndef STATE_COMPLEX
    Solout.Resize( this->NSolutionVariables( var ) );
    
    if(var == 1){
        int id;
        for(id=0 ; id<fDim; id++) {
            TPZFNMatrix<9,STATE> dsoldx;
            TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
            Solout[id] = dsoldx(id,0);//derivate
        }
        return;
    }
    if(var == 2) {
        Solout[0] = Sol[0];//function
        return;
    }//var == 2
    
#endif
    TPZMaterial::Solution(Sol, DSol, axes, var, Solout);
    
}//method

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



void TPZMatMixedPoisson3D::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                                  TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
                                  TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
    
    values.Resize(NEvalErrors());
    values.Fill(0.0);
    
    TPZManVector<STATE> sol(1),dsol(3,0.);
    Solution(u,dudx,axes,2,sol);
    Solution(u,dudx,axes,1,dsol);
    int id;
    //values[1] : eror em norma L2
    REAL  diff = fabs(sol[0] - u_exact[0]);
    values[1]  = diff*diff;
    //values[2] : erro em semi norma H1
    values[2] = 0.;
    for(id=0; id<fDim; id++) {
        DebugStop(); /// @omar:: need to be corrected
//        diff = fabs(dsol[id] - du_exact(id,0));
//        values[2]  += abs(fK)*diff*diff;
    }
    //values[0] : erro em norma H1 <=> norma Energia
    values[0]  = values[1]+values[2];
}

int TPZMatMixedPoisson3D::ClassId() const{
    return Hash("TPZMatMixedPoisson3D") ^ TPZMaterial::ClassId() << 1;
}


