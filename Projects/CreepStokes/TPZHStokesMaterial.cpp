
#include "TPZHStokesMaterial.h"

#include "pzaxestools.h"

#include <fstream>
using namespace std;



void TPZHStokesMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
    //2 = 1 Vel space + 1 Press space for datavecleft
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }

    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    
    
    data.fNeedsNormal = true;
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<3, STATE> normalM(fDimension,1,0.);
    for (int e=0; e<fDimension; e++) {
        normalM(e,0)=normal[e];
    }
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    TPZManVector<REAL, 3> tangent(fDimension,0.);
    TPZFNMatrix<3, STATE> tangentV(fDimension,1,0.);
    for(int i=0; i<fDimension; i++) tangent[i] = data.axes(0,i);
    for(int i=0; i<fDimension; i++) tangentV(i,0) = data.axes(0,i);
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].phi.Rows();
    nshapeP1 = datavecleft[pindex].phi.Rows();
    nshapeP2 = datavecright[pindex].phi.Rows();

//    if(nshapeP1 || nshapeP2) DebugStop();
    
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
      
        TPZFNMatrix<9, STATE> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.), phiVti(fDimension,1,0.);
        TPZFNMatrix<4, STATE> GradV1i(fDimension,fDimension,0.),
        GradV1it(fDimension,fDimension,0.),Du1i(fDimension,fDimension,0.),Du1ni(fDimension,1,0.),  Du1ti(fDimension,1,0.);
        REAL phiit = 0.;
        
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fDeformedDirections(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1i(e,f) = datavecleft[vindex].fDeformedDirections(e,ivec1)*dphiVx1(f,iphi1);
                //termo transposto:
                GradV1it(f,e) = datavecleft[vindex].fDeformedDirections(e,ivec1)*dphiVx1(f,iphi1);
            }
        }
        
        phiit = InnerVec(phiV1i,tangentV);
        for(int e = 0; e<fDimension; e++)
        {
            phiVti(e,0) = phiit*tangent[e];  //Aqui foi retirado um +=
        }
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1i(e,f)= (1./2.) * (GradV1i(e,f) + GradV1it(e,f));
            }
        }
        
        //Du1ni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du1ni(e,0) += Du1i(e,f)*normal[f] ;
                Du1ti(e,0) += Du1i(e,f)*tangent[f];
            }
        }
        
        
        TPZFNMatrix<9, STATE> GradV1nj(fDimension,1,0.);
        
        
        TPZFNMatrix<3, STATE> phiV1tti(fDimension,1,0.),normalM(fDimension,1,0.);
        for (int e=0; e<fDimension; e++) {
            normalM(e,0)=normal[e];
        }
        for (int e=0; e<fDimension; e++) {
            phiV1tti(e,0)=phiV1i(e,0)-InnerVec(phiV1i,normalM)*normal[e];
        }

        
        // K11 - (test V left) * (trial V left)
        for(int j1 = 0; j1 < nshapeV1; j1++)
        {
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            TPZFNMatrix<3, STATE> phiV1j(fDimension,1,0.),phiVtj(fDimension,1,0.), phiV1nj(1,1,0.);
            TPZFNMatrix<4, STATE> GradV1j(fDimension,fDimension,0.),GradV1jt(fDimension,fDimension,0.),
            Du1j(fDimension,fDimension,0.)
            ,Du1nj(fDimension,1,0.),Du1tj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++)
            {
                
                phiV1j(e,0)=datavecleft[vindex].fDeformedDirections(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                for (int f=0; f<fDimension; f++) {
                    GradV1j(e,f) = datavecleft[vindex].fDeformedDirections(e,jvec1)*dphiVx1(f,jphi1);
                    //termo transposto:
                    GradV1jt(f,e) = datavecleft[vindex].fDeformedDirections(e,jvec1)*dphiVx1(f,jphi1);
                    
                }
            }
            
            for(int e = 0; e<fDimension; e++)
            {
                phiVtj(e,0) = phiit*tangent[e];
            }
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1j(e,f) = (1./2.) * (GradV1j(e,f) + GradV1jt(e,f));
                }
            }
            
            //Du1nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du1nj(e,0) += Du1j(e,f)*normal[f];
                    Du1tj(e,0) += Du1j(e,f)*tangent[f];
                }
            }
            //Termo calculado pela tangente:
            
//            STATE fact = (-1.) * weight * 2.* fViscosity * (InnerVec(phiVti, Du1nj)+InnerVec(phiVtj,Du1ni));
//            ek(i1,j1) +=fact;

            //Termo calculado subtraindo a normal
            
            STATE fact = -2. * weight * InnerVec(phiV1tti,Du1nj);
            
            ek(i1,j1) +=fact;
            ek(j1,i1) +=fact;

            TPZFNMatrix<3, STATE> phiV1ttj(fDimension,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiV1ttj(e,0)=phiV1j(e,0)-InnerVec(phiV1j,normalM)*normal[e];
            }
            
            //Termo adionado, presente na formulacao
            STATE fact2 = fBeta * weight * InnerVec(phiV1tti,phiV1ttj);
            ek(i1,j1) +=fact2;
            
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * InnerVec(phiV1tti, phiV1ttj);
            ek(i1,j1) +=penalty;
            
            
        }
        
        
        // K13 and K31 - (trial V left) * (test V tangent right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            
            TPZFNMatrix<3, STATE> phiV2j(fDimension,1,0.),phiVtj(fDimension,1,0.), phiV2nj(1,1,0.);
            TPZFNMatrix<4, STATE> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),
            Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            TPZFNMatrix<3, STATE> phiV2ttj(fDimension,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiV2ttj(e,0)=phiV2j(e,0)-InnerVec(phiV2j,normalM)*normal[e];
            }
            
            TPZFNMatrix<3, STATE> phiVrj(fDimension,1,0.);
            for(int i=0; i<fDimension; i++)
            {
                phiVrj(i,0) = datavecright[vindex].phi(j2,0)*tangentV[i];
            }
            
            STATE fact1 = weight * 2. * fViscosity *InnerVec(Du1ni,phiVrj);
            STATE fact2 = -weight * fBeta * phiit * datavecright[vindex].phi(j2,0);
            
            ek(i1,j2+nshapeV1+nshapeP1) += (fact1+fact2);
            ek(j2+nshapeV1+nshapeP1,i1) += (fact1+fact2);
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * phiit * datavecright[vindex].phi(j2,0);
            ek(i1,j2+nshapeV1+nshapeP1) += -penalty;
            ek(j2+nshapeV1+nshapeP1,i1) += -penalty;
            
            
        }
    }
    
    //datavecright[vindex].phi.Print("phiRight = ",cout);
    
    for(int i2 = 0; i2 < nshapeV2; i2++ )
    {
        
        // K33 - (trial V right) * (test V right)
        for(int j2 = 0; j2 < nshapeV2; j2++)
        {
            
            STATE fact = weight * fBeta * datavecright[vindex].phi(i2,0) * datavecright[vindex].phi(j2,0);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
            
            //Penalidade:
            
            STATE penalty = fSigma * weight * fViscosity * datavecright[vindex].phi(i2,0) * datavecright[vindex].phi(j2,0);;
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) +=penalty;
            
        }
    }

    std::ofstream plotfileM("ekInterfaceH.txt");
    ek.Print("KintH = ",plotfileM,EMathematicaInput);

    
}



void TPZHStokesMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    //return;
    
    STATE rhsnorm = Norm(ef);
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    

    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }

    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    TPZFNMatrix<3, STATE> direction(fDimension,1,0.);
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    
    // Getting the linear combination or finite element approximations
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    
    int nshapeV, nshapeP;
    nshapeV = phiV.Rows();//datavec[vindex].fVecShapeIndex.NElements();
    nshapeP = datavec[pindex].phi.Rows();
    
    if((nshapeV != 0) && (bc.Type() < 10))
    {
        for(int i=0; i<fDimension; i++)
        {
            direction(0,0) = -datavec[vindex].axes(0,1);
            direction(1,0) = datavec[vindex].axes(0,0);
        }
    }
    else if(nshapeV != 0)
    {
        for(int i=0; i<fDimension; i++) direction(i,0) = datavec[vindex].axes(0,i);
    }
    //Adaptação para Hdiv
    int ekr= ek.Rows();
    
    //Vefifica se HDiv
    if(ekr!=nshapeP+nshapeV){
        DebugStop();
    }
    
    
    
    TPZFNMatrix<9, STATE> phiVi(fDimension,1,0.),phiVni(1,1,0.), phiVj(fDimension,1,0.),phiVnj(1,1,0.), phiPi(fDimension,1),phiPj(fDimension,1);
    
    TPZFNMatrix<3,STATE> v_2=bc.Val2();
    TPZFNMatrix<3,STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        case 10:
        {
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                TPZFMatrix<STATE> gradu;
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D = vbc[2];
            }
            
            if(fSpace==1){
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = direction[0] * v_2[0] + direction[1] * v_2[1];
                    
                    ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);
                    
                    for(int j = 0; j < nshapeV; j++){
                        
                        ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                        
                    }
                    
                }
                
                
                
            }else{
                DebugStop();
            }
            
            
            
        }
            break;
            
        case 1: //Neumann for continuous formulation
        case 11:
        {
            
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                bc.ForcingFunction()->Execute(datavec[pindex].x,vbc);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D  = vbc[2]*0.;
                
                
            }
            
            
            for(int i = 0; i < nshapeV; i++ )
            {
                
                //Adaptação para Hdiv
                
                STATE val = 0;
                for(int i=0; i<fDimension; i++) val += direction(i,0)*v_2(i,0);
                
                ef(i,0) += weight * val * datavec[vindex].phi(i,0);
                
            }
            
            
        }
            
            
            
            break;
            
        case 5: //Ponto pressao
        {
            p_D = bc.Val2()(0,0);
            
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                
                ef(i) += 1.0 * p_D * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    ek(i,j) += 1.0 * (phiP(i,0) * phiP(j,0));
                    
                }
                
            }
            
        }
            break;
            
            
//        case 5: // put a value on the diagonal of the pressure to mark a reference pressure
//            ek(0,0) += 1.;
//            break;
            
        default:
        {
            std::cout << "Boundary not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
    
    if(isnan(rhsnorm))
    {
        std::cout << "ef  has norm " << rhsnorm << std::endl;
    }
    
}


void TPZHStokesMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    DebugStop();
    
}

