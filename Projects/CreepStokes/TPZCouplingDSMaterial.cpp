/*
 *  TPZCouplingDSMaterial.cpp
 *  PZ
 *
 *  Created by Thiago Dias dos Santos on 12/01/2015.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZCouplingDSMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "pzmatwithmem.h"
#include "pzfmatrix.h"

//SetSpace
//#define IsHDivQ
//#define IsH1
//#define IsDGM


TPZCouplingDSMaterial::TPZCouplingDSMaterial() : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(){
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1;
    
}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::TPZCouplingDSMaterial(int matid, int dimension, REAL viscosity,REAL permeability, REAL theta) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(matid),fViscosity(viscosity),fTheta(theta),fDimension(dimension)
{
    // symmetric version
    //fTheta = -1;
    
    //fDim = 1;
    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=permeability;

}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::TPZCouplingDSMaterial(const TPZCouplingDSMaterial &mat) : TPZMatWithMem<TPZFMatrix<REAL>, TPZDiscontinuousGalerkin >(mat), fViscosity(mat.fViscosity), fTheta(mat.fTheta),fDimension(mat.fDimension)
{
       fk= mat.fk;
    
}

////////////////////////////////////////////////////////////////////

TPZCouplingDSMaterial::~TPZCouplingDSMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZCouplingDSMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("P", name.c_str()))  return 0;
    if (!strcmp("V", name.c_str()))  return 1;
    if (!strcmp("f", name.c_str()))         return 2;
    if (!strcmp("V_exact", name.c_str()))   return 3;
    if (!strcmp("P_exact", name.c_str()))   return 4;
//    if (!strcmp("V_exactBC", name.c_str()))   return 5;
   
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZCouplingDSMaterial::NSolutionVariables(int var) {
    
    switch(var) {
        
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return this->Dimension(); // Velocity, Vector
        case 2:
            return this->Dimension(); // f, Vector
        case 3:
            return this->Dimension(); // V_exact, Vector
        case 4:
            return this->Dimension(); // P_exact, Vector
//        case 5:
//            return this->Dimension(); // V_exactBC, Vector
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    
    //itapopo conferir esse metodo
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZManVector<REAL,3> v_h = datavec[vindex].sol[0];
    REAL p_h = datavec[pindex].sol[0][0];
    
   // TPZManVector<STATE> v_h = datavec[vindex].sol[0];
   // TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        
        case 0: //Pressure
        {
            Solout[0] = p_h;
        }
            break;
        
        case 1: //Velocity
        {
            Solout[0] = v_h[0]; // Vx
            Solout[1] = v_h[1]; // Vy
        }
            break;
        case 2: //f
        {
            TPZVec<double> f;
            if(this->HasForcingFunction()){
                this->ForcingFunction()->Execute(datavec[vindex].x, f);
            }
            Solout[0] = f[0]; // fx
            Solout[1] = f[1]; // fy
        }
            break;
        
        case 3: //v_exact
        {
            TPZVec<double> v;
            if(this->HasfForcingFunctionExact()){
//                this->ForcingFunctionExact()->Execute(datavec[vindex].x, v); // @omar::check it!
            }
            Solout[0] = v[0]; // vx
            Solout[1] = v[1]; // vy
        }
            break;
        
        case 4: //p_exact
        {
            TPZVec<double> p;
            if(this->HasfForcingFunctionExact()){
//                this->ForcingFunctionExactPressure()->Execute(datavec[pindex].x, p); // @omar::check it!
            }
            Solout[0] = p[0]; // px
            
        }
            break;
            
//        case 5: //v_exact
//        {
//            TPZVec<double> vbc;
//            if(this->HasffBCForcingFunction()){
//                this->ForcingFunctionBC()->Execute(datavec[vindex].x, vbc);
//            }
//            Solout[0] = vbc[0]; // vbcx
//            Solout[1] = vbc[1]; // vbcy
//        }
//            break;
            
            
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element
void TPZCouplingDSMaterial::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    
    //itapopo conferir esse método. Foi copiado do TPZDarcyFlow3D
    
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<STATE> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<STATE> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<STATE> Qaxes = datavec[ublock].axes;
    TPZFMatrix<STATE> QaxesT;
    TPZFMatrix<STATE> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<STATE> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<STATE> GradOfX;
    TPZFMatrix<STATE> GradOfXInverse;
    TPZFMatrix<STATE> VectorOnMaster;
    TPZFMatrix<STATE> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    if (HDivPiola == 1)
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fNormalVec(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fNormalVec(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fNormalVec(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            DivergenceofPhi(iq,0) =  (1.0/JacobianDet) * ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) +
                                                          dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) +
                                                          dphiuH1(2,ishapeindex)*VectorOnMaster(2,0) );
        }
    }
    else
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            /* Computing the divergence for constant jacobian elements */
            DivergenceofPhi(iq,0) =  datavec[ublock].fNormalVec(0,ivectorindex)*GradphiuH1(0,ishapeindex) +
            datavec[ublock].fNormalVec(1,ivectorindex)*GradphiuH1(1,ishapeindex) +
            datavec[ublock].fNormalVec(2,ivectorindex)*GradphiuH1(2,ishapeindex) ;
        }
    }
    
    return;
    
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Write(TPZStream &buf, int withclassid) {
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
 
    
}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::Read(TPZStream &buf, void *context) {
    
    TPZDiscontinuousGalerkin::Read(buf, context);

}

////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi){
    
   
    TPZFMatrix<STATE> &dphiV = dataV.dphix;
    
    const int dim = this->Dimension();
    
    GradPhi.clear();
    GradPhi.resize(dim);
    
    //for each shape
    for(int shape = 0; shape < dphiV.Rows(); shape++){
        
        TPZFMatrix<STATE> GPhi(dim,dim,0.);
        
        for(int i = 0; i < dim; i++){
            
            for(int j = 0; j < dim; j++){
                
                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
            
            }//j
        }//i
        
        GradPhi[shape] = GPhi;
        
    }//shape
    
}

// Contricucao dos elementos internos

void TPZCouplingDSMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
 
    return;
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavec[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavec[vindex]);
    }
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = datavec[vindex].normal;
    
    TPZManVector<REAL> n = datavec[0].normal;
    
       std::cout<<n<<std::endl;
    
        std::cout<<"____"<<std::endl;
    
//    TPZFNMatrix<4,REAL> nx(2,2),tx(2,2);
//    for (int i=0; i<2; i++) {
//        for (int j=0; j<2 ; j++) {
//            nx(i,j) = n[i]*n[j];
//            tx(i,j) = t[i]*t[j];
//        }
//    }
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    TPZVec<double> f;
    TPZFMatrix<STATE> phiVi(fDimension,1,0.0),phiVj(fDimension,1,0.0);
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<4> phiVti(1,1,0.);
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0) = phiV(iphi,0)*datavec[vindex].fNormalVec(e,ivec);
            
        }
        //phiVti(0,0) = normal[1];
        
        //phiVti(0,0) = phiVi(0,0)*(-normal[1])+phiVi(1,0)*(normal[0]);
        // matrix A - velocity * test-funtion velocity
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<4> phiVtj(1,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiVj(e,0) = phiV(jphi,0)*datavec[vindex].fNormalVec(e,jvec);
            }
            //phiVtj(0,0) = phiVj(0,0)*(-normal[1])+phiVj(1,0)*(normal[0]);
            
            //STATE val = phiVti(0,0)* phiVtj(0,0);
            
            STATE val = Inner(phiVi,phiVj);
            ek(i,j) += weight * fViscosity * pow(fk,-1./2.)* val ;
            
        }
        
        
    }

    
//    if(this->HasForcingFunction()){
//        this->ForcingFunction()->Execute(datavec[vindex].x, f);
//    }
//    
//    for (int i = 0; i < nshapeP; i++) {
//        
//        STATE factf= weight * phiP(i,0)*f[0];
//        ef(nshapeV+i,0) += factf;
//        
//    }
    

    
    //    std::cout<<ek<<std::endl;
    //
    //    std::cout<<ef<<std::endl;
    //    std::cout<<"____"<<std::endl;

    

}


void TPZCouplingDSMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    DebugStop();
    
}




////////////////////////////////////////////////////////////////////

void TPZCouplingDSMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){

    
    //this->ContributeInterface_t(data, datavecleft, datavecright, weight, ek, ef);
    
    //this->ContributeInterface_pn(data, datavecleft, datavecright, weight, ek, ef);
    
    //this->ContributeInterface_pn00(data, datavecleft, datavecright, weight, ek, ef);
    
    //return;

#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();

    
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &n = data.normal;
    TPZManVector<REAL,3> t(2,0.);
    t[0]=-n[1];
    t[1]=n[0];
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);

    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    TPZVec<double> f, v1, v2;
    TPZFMatrix<STATE> phiV1i(fDimension,1,0.0),phiV1j(fDimension,1,0.0);
    
    
    TPZFMatrix<STATE> phiV2i(fDimension,1,0.0),phiV2j(fDimension,1,0.0);
    
    for(int i2 = 0; i2 < nshapeV2; i2++ )
    {
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        TPZFNMatrix<4> phiV2ti(1,1,0.);
        for (int e=0; e<fDimension; e++) {
            phiV2i(e,0) = phiV2(iphi2,0)*datavecright[vindex].fNormalVec(e,ivec2);
            
        }

        
        phiV2ti(0,0) = phiV2i(0,0)*t[0]+phiV2i(1,0)*t[1];
        // matrix A - velocity * test-funtion velocity
        
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            
            TPZFNMatrix<4> phiV2tj(1,1,0.);
            for (int e=0; e<fDimension; e++) {
                phiV2j(e,0) = phiV2(jphi2,0)*datavecright[vindex].fNormalVec(e,jvec2);
            }
            phiV2tj(0,0) = phiV2j(0,0)*t[0]+phiV2j(1,0)*t[1];
            
            
            STATE val = phiV2ti(0,0) * phiV2tj(0,0) ;
            
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += weight * fViscosity * val;
            //ek(j2+nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += -weight * fViscosity * val ;

        }

        

        TPZFNMatrix<9> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.),GradV2ntj(1,1,0.);
        
        
//        // Tensor
//        for(int j2 = 0; j2 < nshapeV2; j2++){
//            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
//            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
//            TPZFNMatrix<9> GradV2nj(fDimension,1);
//            //TPZManVector<REAL,3> phiP1j(fDimension);
//            
//            GradV2nj.Zero();
//            
//            for (int e=0; e<fDimension; e++) {
//                for (int f=0; f<fDimension; f++) {
//                    GradV2nj(e,0) += datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*n[f];
//                }
//            }
//            
//            GradV2ntj(0,0)=GradV2nj(0,0)*t[0]+GradV2nj(1,0)*t[1];
//            
//            
//            STATE fact = (1.) * weight * fViscosity * phiV2ti(0,0) * GradV2ntj(0,0);
//            
//            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
//            ek(j2+nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
//            
//            
//            }

        
//        TPZVec<double> f;
//        
//        if(this->HasForcingFunction()){
//            this->ForcingFunction()->Execute(datavecright[vindex].x, f);
//        }
//        
//        TPZFNMatrix<9> GradV2_L(fDimension,1),GradV2ntj_L(1,1,0.);
//        //TPZManVector<REAL,3> phiP1j(fDimension);
//        
//        GradV2_L.Zero();
//        
//        for (int e=0; e<fDimension; e++) {
//            GradV2_L(e,0) = f[e];
//        }
//        
//        GradV2ntj_L(0,0)=GradV2_L(0,0)*t[0]+GradV2_L(1,0)*t[1];
//        
//        STATE fact = (-1.) * weight * fViscosity * phiV2ti(0,0) * GradV2ntj_L(0,0) ;
//        
//        ef(i2+nshapeV1+nshapeP1) += weight * fact;
//
//        
        
        }
    
    

// * pow(fk,-1./2.)
    
    
}




void TPZCouplingDSMaterial::ContributeInterface_pn00(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
    
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
#ifdef PZDEBUG
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
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    //Detjac
    REAL Detjac=fabs(data.detjac);
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    //TPZManVector<REAL,3> normalx(fDimension,phiP2.Cols());
    //TPZAxesTools<REAL>::Axes2XYZ(normal, normalx, data.axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        
        
        TPZFNMatrix<9> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.);
        
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
        }
        
        
        // K34 e K43- (trial V right) * (test P right)
        for(int j1 = 0; j1 < nshapeV1; j1++){
            
            TPZFNMatrix<9> phiV1j(fDimension,1),phiV1nj(1,1,0.);
            
            int jphi1 = datavecright[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecright[vindex].fVecShapeIndex[j1].first;
            
            for (int e=0; e<fDimension; e++) {
                
                phiV1j(e,0)=datavecright[vindex].fNormalVec(e,jvec1)*datavecright[vindex].phi(jphi1,0);
                phiV1nj(0,0)+=phiV1j(e,0)*normal[e];
                
            }
            
            STATE fact = (-1.) * weight * InnerVec(phiV1ni,phiV1nj);
            
            ek(i1,j1) += fact;
            //ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
        }
        
    }
    
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        TPZFNMatrix<9> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.);
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
        }
        
        // K34 e K43- (trial V right) * (test P right)
        for(int j2 = 0; j2 < nshapeV2; j2++){
            
            TPZFNMatrix<9> GradV2nj(fDimension,1),phiV2j(fDimension,1),phiV2nj(1,1,0.);
            
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            
            for (int e=0; e<fDimension; e++) {
                
                phiV2j(e,0)=datavecright[vindex].fNormalVec(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2nj(0,0)+=phiV2j(e,0)*normal[e];
                
            }
            
            STATE fact1 = (1.) * weight * InnerVec(phiV2ni,phiV2nj);
            
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact1;
            //ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
        }
        
    }
    
    
    
    
    
}



void TPZCouplingDSMaterial::ContributeInterface_pn(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    // Verificar que
    // os termos mistos devem estar sem viscosidade!
    
    
#ifdef PZDEBUG
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
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    // Setting forcing function
    /*STATE force = 0.;
     if(this->fForcingFunction) {
     TPZManVector<STATE> res(1);
     fForcingFunction->Execute(datavec[pindex].x,res);
     force = res[0];
     }*/
    
    //Gravity
    STATE rhoi = 900.; //itapopo
    STATE g = 9.81; //itapopo
    STATE force = rhoi*g;
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    //Detjac
    REAL Detjac=fabs(data.detjac);
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    //TPZManVector<REAL,3> normalx(fDimension,phiP2.Cols());
    //TPZAxesTools<REAL>::Axes2XYZ(normal, normalx, data.axes);
    
    TPZManVector<STATE> v2_h = datavecright[vindex].sol[0];
    TPZManVector<STATE> p2_h = datavecright[pindex].sol[0];
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        
        
        TPZFNMatrix<9> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.);
        
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1ni(e,0)+=datavecleft[vindex].fNormalVec(e,ivec1)*dphiVx1(f,iphi1)*normal[f];
            }
        }
        
        
        
        TPZFNMatrix<9> GradV1nj(fDimension,1,0.);
        
        
        // K12 e K21 - (trial V left) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            
            STATE fact = (-1.) * weight * Inner(phiV1ni,phiP1j);
            
            ek(i1,j1+nshapeV1) += fact;
            //ek(j1+nshapeV1,i1) += fact*fTheta;
            
        }
    
    }
    
    
//    for(int i2 = 0; i2 < nshapeV2; i2++ ){
//        
//
//        TPZFNMatrix<9> GradV2ni(fDimension,1),phiV2i(fDimension,1),phiV2ni(1,1,0.);
//        
//        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
//        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
//        
//        for (int e=0; e<fDimension; e++) {
//            
//            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
//            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
//            
//            for (int f=0; f<fDimension; f++) {
//                GradV2ni(e,0) += datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2)*normal[f];
//            }
//        }
//
//        
//        
//        // K33 - (trial V right) * (test V right)
//        for(int j2 = 0; j2 < nshapeV2; j2++){
//            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
//            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
//            TPZFNMatrix<9> GradV2nj(fDimension,1);
//            //TPZManVector<REAL,3> phiP1j(fDimension);
//            
//            GradV2nj.Zero();
//            
//            for (int e=0; e<fDimension; e++) {
//                for (int f=0; f<fDimension; f++) {
//                    GradV2nj(e,0) += datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*normal[f];
//                }
//            }
//            
//            STATE fact = (1.) * weight * fViscosity * InnerVec(phiV2i,GradV2nj);
//            
//            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
//            //ek(j2+nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
//            
//            
//        }
//        
//        // K34 e K43- (trial V right) * (test P right)
//        for(int j2 = 0; j2 < nshapeP2; j2++){
//            
//            TPZFNMatrix<9> phiP2j(1,1,0.);
//            phiP2j(0,0)=phiP2(j2,0);
//            
//            STATE fact = (-1.) * weight * InnerVec(phiV2ni,phiP2j);
//            
//            ek(i2+nshapeV1+nshapeP1,j2+2*nshapeV1+nshapeP1) += fact;
//            //ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
//        }
//        
//    }
//    

    
   
}


void TPZCouplingDSMaterial::ContributeInterface_pf(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left Darcy -> 1
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    
    // V - right Stokes -> 2
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left Darcy -> 1
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    
    // P - right Stokes -> 2
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    for(int i1 = 0; i1 < nshapeP1; i1++ )
    {
        
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            ek(i1+nshapeV1,j1+nshapeV1) += weight * (phiP1(i1,0) * phiP1(j1,0));
            
        }
        
    }
    
    for(int i2 = 0; i2 < nshapeP2; i2++ )
        {
    
        for(int j2 = 0; j2 < nshapeP2; j2++){
    
            ek(i2+nshapeV2+nshapeP1+nshapeV1,j2+nshapeV2+nshapeP1+nshapeV1) += weight * (phiP2(i2,0) * phiP2(j2,0));

            }
                
        }
 
    
}

void TPZCouplingDSMaterial::ContributeInterface_BJS(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    /* The first part is stisfied by construction of hdvi(Omega) */
    
    /* The second part  */
//    this->ContributeInterface_BJS_II(data, datavecleft, datavecright, weight, ek, ef);
    
    /* The third part  */
   // this->ContributeInterface_BJS_III(data, datavecleft, datavecright, weight, ek, ef);
    
    
    //
    //        TPZVec<double> v_2;
    //
    //        if(this->HasfForcingFunctionExact()){
    //            this->ForcingFunctionExact()->Execute(datavecright[vindex].x, v_2);
    //        }
    //
    //
    //            REAL vh2_n = v2_h[0];
    //            REAL v2_n = normal[0] * v_2[0] + normal[1] * v_2[1];
    //
    //            ef(i2+nshapeV1+nshapeP1,0) += weight * gBigNumber * (vh2_n - v2_n) * phiV2(i2,0);
    //
    //            for(int j2 = 0; j2 < nshapeV2; j2++){
    //
    //                ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += weight * gBigNumber * phiV2(j2,0) * phiV2(i2,0);
    //                
    //            }
    
    
    
    return;
    
}

void TPZCouplingDSMaterial::ContributeInterface_t(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left Darcy -> 1
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    
    // V - right Stokes -> 2
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left Darcy -> 1
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    
    // P - right Stokes -> 2
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    TPZManVector<STATE> v_h = datavecright[vindex].sol[0];
    TPZManVector<STATE> p_h = datavecright[pindex].sol[0];
    
    
    int sizek=ek.Rows();
    
    //Dirichlet
    

    
    for(int i2 = 0; i2 < nshapeV2; i2++ )
    {
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        TPZFNMatrix<9> GradV2ni(fDimension,1,0.),phiV2i(fDimension,1),phiV2ni(1,1,0.),phiV2ti(1,1,0.);
        GradV2ni.Zero();
        
        TPZFNMatrix<4> GradV2i(fDimension,fDimension,0.),GradV2it(fDimension,fDimension,0.),Du2i(fDimension,fDimension,0.),Du2ni(fDimension,1,0.);
        
        for (int e=0; e<fDimension; e++) {
            
            for (int f=0; f<fDimension; f++) {
                GradV2i(e,f) = datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2);
                //termo transposto:
                GradV2it(f,e) = datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2);
                
            }
        }
        
        //Du = 0.5(GradU+GradU^T)
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2i(e,f)= (1./2.) * (GradV2i(e,f) + GradV2it(e,f));
            }
        }
        
        //Duni
        for (int e=0; e<fDimension; e++) {
            for (int f=0; f<fDimension; f++) {
                Du2ni(e,0) += Du2i(e,f)*normal[f] ;
            }
        }
        
        
        for (int e=0; e<fDimension; e++) {
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
        }
        
        TPZManVector<REAL> n = data.normal;
        TPZManVector<REAL> t(2);
        t[0]=-n[1];
        t[1]=n[0];
        
        
        
        phiV2ti(0,0)= t[0] * phiV2i(0,0) + t[1] * phiV2i(1,0);
        TPZFNMatrix<9> phiV2tit(fDimension,1,0.);
        phiV2tit(0,0)=phiV2ti(0,0)*t[0];
        phiV2tit(1,0)=phiV2ti(0,0)*t[1];
        
//
//        REAL vh_t = v_h[1];
//        
//        REAL v_t = t[0] * v_2[0] + t[1] * v_2[1];
//        
//        TPZManVector<REAL> v_tt(2);
//        v_tt[0]=v_t*t[0];
//        v_tt[1]=v_t*t[1];
//        
//        TPZManVector<REAL> vh_tt(2);
//        vh_tt[0]=vh_t*t[0];
//        vh_tt[1]=vh_t*t[1];
//        
//        TPZFNMatrix<9> diffvt(fDimension,1,0.);
//        diffvt(0,0)=v_tt[0];
//        diffvt(1,0)=v_tt[1];
//        
//        STATE factf=(1.) * weight * fViscosity * InnerVec(diffvt,GradVni) ;
//        
        ef(i2,0) += 0.;
        
        
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            
            TPZFNMatrix<9> GradV2nj(fDimension,1),phiV2tj(1,1,0.),phiV2j(fDimension,1);
            
            for (int e=0; e<fDimension; e++) {
                phiV2j(e,0)=datavecright[vindex].fNormalVec(e,jvec2)*datavecright[vindex].phi(jphi2,0);
            }
            
            //std::cout<<phiVj<<std::endl;
            
            phiV2tj(0,0)= t[0] * phiV2j(0,0) + t[1] * phiV2j(1,0);
            
            
            
            TPZFNMatrix<4> GradV2j(fDimension,fDimension,0.),GradV2jt(fDimension,fDimension,0.),Du2j(fDimension,fDimension,0.),Du2nj(fDimension,1,0.);
            
            for (int e=0; e<fDimension; e++) {
                
                for (int f=0; f<fDimension; f++) {
                    GradV2j(e,f) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    //termo transposto:
                    GradV2jt(f,e) = datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2);
                    
                }
            }
            
            //Du = 0.5(GradU+GradU^T)
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2j(e,f)= (1./2.) * (GradV2j(e,f) + GradV2jt(e,f));
                }
            }
            
            //Du2nj
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    Du2nj(e,0) += Du2j(e,f)*normal[f] ;
                }
            }
            
            
            
            STATE fact =(-1.) * weight * 2. * fViscosity * InnerVec(phiV2tit, Du2nj) ;
            ek(i2+nshapeP1+nshapeV1,j2+nshapeP1+nshapeV1) += fact ;
            ek(j2+nshapeP1+nshapeV1,i2+nshapeP1+nshapeV1) += fTheta*fact;
            
        }
        
    }


}



void TPZCouplingDSMaterial::ContributeInterface_BJS_II(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left Darcy -> 1
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    
    // V - right Stokes -> 2
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left Darcy -> 1
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    
    // P - right Stokes -> 2
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    /* The second part  */
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        TPZFNMatrix<9> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.);
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fNormalVec(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];

        }
        
        
        /* Stokes Tensor part */
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            GradV2nj.Zero();
            
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradV2nj(e,0) += datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*normal[f];
                }
            }
            
            REAL fact = weight * fViscosity * InnerVec(phiV1i,GradV2nj);
            ek(i1,j2+nshapeV1+nshapeP1) += fact;
            
            
        }
        
        /* Darcy Pressure */
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            STATE fact = weight * Inner(phiV1ni,phiP1j);
            ek(i1,j1+nshapeV1) += fact;
            
        }
        
        
    }
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        /* Compute the phivec, phivec and gradient of velocity applied to normal vector */
        TPZFNMatrix<9> phiV2i(fDimension,1),phiV2ni(1,1,0.),GradV2ni(fDimension,1);
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV2ni(e,0) += datavecright[vindex].fNormalVec(e,ivec2)*dphiVx2(f,iphi2)*normal[f];
            }
        }
        
        
        /* Stokes Tensor part */
        for(int j2 = 0; j2 < nshapeV2; j2++){
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            TPZFNMatrix<9> GradV2nj(fDimension,1);
            //TPZManVector<REAL,3> phiP1j(fDimension);
            
            GradV2nj.Zero();
            
            for (int e=0; e<fDimension; e++) {
                for (int f=0; f<fDimension; f++) {
                    GradV2nj(e,0) += datavecright[vindex].fNormalVec(e,jvec2)*dphiVx2(f,jphi2)*normal[f];
                }
            }
            
            REAL fact = weight * fViscosity * InnerVec(phiV2i,GradV2nj);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
            
            
        }
        
        /* Stokes Pressure */
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            STATE fact = -1.0* weight * InnerVec(phiV2ni,phiP2j);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1+nshapeV2) += fact;
        }
        
        /* Darcy Pressure */
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            STATE fact = weight * Inner(phiV2ni,phiP1j);
            ek(i2+nshapeV1+nshapeP1,j1+nshapeV1) += fact;
            
        }
        
    }
    
}

void TPZCouplingDSMaterial::ContributeInterface_BJS_III(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
#ifdef PZDEBUG
    //2 = 1 Vel space + 1 Press space for datavecleft Darcy
    int nrefleft =  datavecleft.size();
    if (nrefleft != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
    
    //2 = 1 Vel space + 1 Press space for datavecright Stokes
    int nrefright =  datavecright.size();
    if (nrefright != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int vindex = this->VIndex();
    const int pindex = this->PIndex();
    
    if (datavecleft[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecleft[vindex]);
    }
    
    if (datavecright[vindex].fVecShapeIndex.size() == 0) {
        FillVecShapeIndex(datavecright[vindex]);
    }
    
    
    // Setting the phis
    // V - left Darcy -> 1
    TPZFMatrix<REAL> &phiV1 = datavecleft[vindex].phi;
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    
    // V - right Stokes -> 2
    TPZFMatrix<REAL> &phiV2 = datavecright[vindex].phi;
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left Darcy -> 1
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    
    // P - right Stokes -> 2
    TPZFMatrix<REAL> &phiP2 = datavecright[pindex].phi;
    TPZFMatrix<REAL> &dphiP2 = datavecright[pindex].dphix;
    
    //Normal
    TPZManVector<REAL,3> &normal = data.normal;
    
    //tangent
    TPZManVector<REAL,3> &tangent = normal;
    tangent[0] = normal[1];
    tangent[1] = -normal[0];
    
    
    TPZFNMatrix<220,REAL> dphiVx1(fDimension,dphiV1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV1, dphiVx1, datavecleft[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiVx2(fDimension,dphiV2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV2, dphiVx2, datavecright[vindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx1(fDimension,phiP1.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP1, dphiPx1, datavecleft[pindex].axes);
    
    
    TPZFNMatrix<220,REAL> dphiPx2(fDimension,phiP2.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP2, dphiPx2, datavecright[pindex].axes);
    
    
    int nshapeV1, nshapeV2, nshapeP1,nshapeP2;
    
    nshapeV1 = datavecleft[vindex].fVecShapeIndex.NElements();
    nshapeV2 = datavecright[vindex].fVecShapeIndex.NElements();
    
    nshapeP1 = phiP1.Rows();
    nshapeP2 = phiP2.Rows();
    
    REAL beta = 1.0;
    
    /* The third part  */
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        /* Compute the phivec, phivec and gradient of velocity applied to normal vector */
        TPZFNMatrix<9> phiV2i(fDimension,1),phiV2ti(1,1,0.);
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fNormalVec(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ti(0,0)+=phiV2i(e,0)*tangent[e];
            
        }
        
        for(int j1 = 0; j1 < nshapeV1; j1++){
            
            /* Compute the phivec, phivec and gradient of velocity applied to normal vector */
            TPZFNMatrix<9> phiV1j(fDimension,1),phiV1tj(1,1,0.);
            int jphi1 = datavecleft[vindex].fVecShapeIndex[j1].second;
            int jvec1 = datavecleft[vindex].fVecShapeIndex[j1].first;
            
            for (int e=0; e<fDimension; e++) {
                
                phiV1j(e,0)=datavecleft[vindex].fNormalVec(e,jvec1)*datavecleft[vindex].phi(jphi1,0);
                phiV1tj(0,0)+=phiV1j(e,0)*tangent[e];
                
            }
            
            REAL fact = - beta * weight * phiV1tj(0,0) * phiV2ti(0,0);
            ek(i2+nshapeV1+nshapeP1,j1) += fact;
            
        }
        
        for(int j2 = 0; j2 < nshapeV2; j2++){

            /* Compute the phivec, phivec and gradient of velocity applied to normal vector */
            TPZFNMatrix<9> phiV2j(fDimension,1),phiV2tj(1,1,0.);
            int jphi2 = datavecright[vindex].fVecShapeIndex[j2].second;
            int jvec2 = datavecright[vindex].fVecShapeIndex[j2].first;
            
            for (int e=0; e<fDimension; e++) {
                
                phiV2j(e,0)=datavecright[vindex].fNormalVec(e,jvec2)*datavecright[vindex].phi(jphi2,0);
                phiV2tj(0,0)+=phiV2j(e,0)*tangent[e];
                
            }
            
            REAL fact =  beta * weight * phiV2tj(0,0) * phiV2ti(0,0);
            ek(i2+nshapeV1+nshapeP1,j2+nshapeV1+nshapeP1) += fact;
        }
        
    }
    
}

void TPZCouplingDSMaterial::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
   
    DebugStop();
    
}



////////////////////////////////////////////////////////////////////

STATE TPZCouplingDSMaterial::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two tensors

    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }

    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}


////////////////////////////////////////////////////////////////////

STATE TPZCouplingDSMaterial::InnerVec(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two vectors
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0;
    
    for(int j = 0; j < S.Cols(); j++){
        for(int i = 0; i < S.Rows(); i++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



////////////////////////////////////////////////////////////////////

STATE TPZCouplingDSMaterial::Tr( TPZFMatrix<REAL> &GradU ){
 
#ifdef DEBUG
    if( GradU.Rows() != GradU.Cols() ) {
        DebugStop();
    }
#endif
    
    STATE Val = 0.;
    
    for(int i = 0; i < GradU.Rows(); i++){
        Val += GradU(i,i);
    }
    
    return Val;
}


/// transform a H1 data structure to a vector data structure
void TPZCouplingDSMaterial::FillVecShapeIndex(TPZMaterialData &data)
{
    data.fNormalVec.Resize(fDimension,fDimension);
    data.fNormalVec.Identity();
    data.fVecShapeIndex.Resize(fDimension*data.phi.Rows());
    for (int d=0; d<fDimension; d++) {
        for (int i=0; i<data.phi.Rows(); i++) {
            data.fVecShapeIndex[i*fDimension+d].first = d;
            data.fVecShapeIndex[i*fDimension+d].second = i;
        }
    }
}



void TPZCouplingDSMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors)
{
 
    DebugStop();
}