/*
 *  TPZDarcyPMaterial.cpp
 *  PZ
 *
 *  Created by Pablo G. S. Carvalho on 08/09/2016.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZDarcyPMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"
#include "TPZMatWithMem.h"
#include "pzfmatrix.h"


TPZDarcyPMaterial::TPZDarcyPMaterial() : TPZMatWithMem<TPZFMatrix<STATE>, TPZMaterial >(){

    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    fk=1.;
    
}

////////////////////////////////////////////////////////////////////

TPZDarcyPMaterial::TPZDarcyPMaterial(int matid, int dimension, int space, STATE viscosity, STATE permeability, STATE theta) : TPZMatWithMem<TPZFMatrix<STATE>, TPZMaterial >(matid),fDimension(dimension),fSpace(space),fViscosity(viscosity),fk(permeability),fTheta(theta)
{

    TPZFNMatrix<3,STATE> Vl(1,1,0.);
    this->SetDefaultMem(Vl);
    
}

////////////////////////////////////////////////////////////////////

TPZDarcyPMaterial::TPZDarcyPMaterial(const TPZDarcyPMaterial &mat) : TPZMatWithMem<TPZFMatrix<STATE>, TPZMaterial >(mat),fDimension(mat.fDimension),fSpace(mat.fSpace),fViscosity(mat.fViscosity),fk(mat.fk), fTheta(mat.fTheta)
{
    
    
}

////////////////////////////////////////////////////////////////////

TPZDarcyPMaterial::~TPZDarcyPMaterial(){
    
    
}


////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::FillDataRequirements(TPZMaterialData &data)
{
    TPZMaterial::FillDataRequirements(data);
    data.fNeedsSol = true;
}

////////////////////////////////////////////////////////////////////


void TPZDarcyPMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsNormal = true;
        datavec[idata].fNeedsSol = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.fNeedsNormal = true;
}


////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZDarcyPMaterial::VariableIndex(const std::string &name) {
    
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

int TPZDarcyPMaterial::NSolutionVariables(int var) {
    
    switch(var) {
            
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return this->Dimension(); // Velocity, Vector
        case 2:
            return 1; // f, Scalar
        case 3:
            return this->Dimension(); // V_exact, Vector
        case 4:
            return 1; // P_exact, Scalar

        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) {
    
    
    //itapopo conferir esse metodo
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZManVector<STATE,3> v_h = datavec[vindex].sol[0];
    STATE p_h = datavec[pindex].sol[0][0];
    
    TPZFNMatrix<9,STATE> gradu(2,1);

    
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
            TPZVec<STATE> f;
            if(this->HasForcingFunction()){
                this->ForcingFunction()->Execute(datavec[vindex].x, f, gradu);
            }
            Solout[0] = f[0]; // fx
        }
            break;
            
        case 3: //v_exact
        {
            TPZVec<STATE> sol;
            if(this->HasExactSol()){
                this->fExactSol->Execute(datavec[vindex].x, sol, gradu); // @omar::check it!
            }
            Solout[0] = sol[0]; // vx
            Solout[1] = sol[1]; // vy
        }
            break;
            
        case 4: //p_exact
        {
            TPZVec<STATE> sol;
            if(this->HasExactSol()){
                this->fExactSol->Execute(datavec[pindex].x, sol, gradu); // @omar::check it!
            }
            Solout[0] = sol[2]; // px
            
        }
            break;

           
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
}

////////////////////////////////////////////////////////////////////

// Divergence on deformed element
void TPZDarcyPMaterial::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
{
    
    //itapopo conferir esse método. Foi copiado do TPZDarcyFlow3D
    
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFMatrix<REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFMatrix<REAL> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFMatrix<REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFMatrix<REAL> Qaxes = datavec[ublock].axes;
    TPZFMatrix<REAL> QaxesT;
    TPZFMatrix<REAL> Jacobian = datavec[ublock].jacobian;
    TPZFMatrix<REAL> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFMatrix<REAL> GradOfX;
    TPZFMatrix<REAL> GradOfXInverse;
    TPZFMatrix<REAL> VectorOnMaster;
    TPZFMatrix<REAL> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fDeformedDirections(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fDeformedDirections(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fDeformedDirections(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            DivergenceofPhi(iq,0) =  (1.0/JacobianDet) * ( dphiuH1(0,ishapeindex)*VectorOnMaster(0,0) +
                                                          dphiuH1(1,ishapeindex)*VectorOnMaster(1,0) +
                                                          dphiuH1(2,ishapeindex)*VectorOnMaster(2,0) );
        }
    }
    
    return;
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::Write(TPZStream &buf, int withclassid) const{
    
    TPZMaterial::Write(buf, withclassid);
    
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::Read(TPZStream &buf, void *context) {
    
    TPZMaterial::Read(buf, context);
    
}

// Divergence on master element
void TPZDarcyPMaterial::ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU)
{
    int ublock = 0;
    
    // Getting test and basis functions
    TPZFNMatrix<100,REAL> phiuH1         = datavec[ublock].phi;   // For H1  test functions Q
    TPZFNMatrix<300,REAL> dphiuH1       = datavec[ublock].dphi; // Derivative For H1  test functions
    TPZFNMatrix<300,REAL> dphiuH1axes   = datavec[ublock].dphix; // Derivative For H1  test functions
    TPZFNMatrix<9,STATE> gradu = datavec[ublock].dsol[0];
    TPZFNMatrix<9,STATE> graduMaster;
    gradu.Transpose();
    
    TPZFNMatrix<660> GradphiuH1;
    TPZAxesTools<REAL>::Axes2XYZ(dphiuH1axes, GradphiuH1, datavec[ublock].axes);
    
    int nphiuHdiv = datavec[ublock].fVecShapeIndex.NElements();
    
    DivergenceofPhi.Resize(nphiuHdiv,1);
    
    REAL JacobianDet = datavec[ublock].detjac;
    
    TPZFNMatrix<9,REAL> Qaxes = datavec[ublock].axes;
    TPZFNMatrix<9,REAL> QaxesT;
    TPZFNMatrix<9,REAL> Jacobian = datavec[ublock].jacobian;
    TPZFNMatrix<9,REAL> JacobianInverse = datavec[ublock].jacinv;
    
    TPZFNMatrix<9,REAL> GradOfX;
    TPZFNMatrix<9,REAL> GradOfXInverse;
    TPZFNMatrix<9,REAL> VectorOnMaster;
    TPZFNMatrix<9,REAL> VectorOnXYZ(3,1,0.0);
    Qaxes.Transpose(&QaxesT);
    QaxesT.Multiply(Jacobian, GradOfX);
    JacobianInverse.Multiply(Qaxes, GradOfXInverse);
    
    int ivectorindex = 0;
    int ishapeindex = 0;
    
    {
        for (int iq = 0; iq < nphiuHdiv; iq++)
        {
            ivectorindex = datavec[ublock].fVecShapeIndex[iq].first;
            ishapeindex = datavec[ublock].fVecShapeIndex[iq].second;
            
            VectorOnXYZ(0,0) = datavec[ublock].fDeformedDirections(0,ivectorindex);
            VectorOnXYZ(1,0) = datavec[ublock].fDeformedDirections(1,ivectorindex);
            VectorOnXYZ(2,0) = datavec[ublock].fDeformedDirections(2,ivectorindex);
            
            GradOfXInverse.Multiply(VectorOnXYZ, VectorOnMaster);
            VectorOnMaster *= JacobianDet;
            
            /* Contravariant Piola mapping preserves the divergence */
            REAL dot = 0.0;
            for (int i = 0;  i < fDimension; i++) {
                dot += dphiuH1(i,ishapeindex)*VectorOnMaster(i,0);
            }
            DivergenceofPhi(iq,0) = (1.0/JacobianDet) * dot;
        }
        
    }
    
    return;
    
}

////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi){
    
    
    TPZFMatrix<REAL> &dphiV = dataV.dphix;
    
    const int dim = this->Dimension();
    
    GradPhi.clear();
    GradPhi.resize(dim);
    
    //for each shape
    for(int shape = 0; shape < dphiV.Rows(); shape++){
        
        TPZFMatrix<REAL> GPhi(dim,dim,0.);
        
        for(int i = 0; i < dim; i++){
            
            for(int j = 0; j < dim; j++){
                
                GPhi(i,j) = dphiV(j,shape);// itapopo H1 ??
                
            }//j
        }//i
        
        GradPhi[shape] = GPhi;
        
    }//shape
    
}

// Contricucao dos elementos internos

void TPZDarcyPMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){
    
    
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
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    TPZVec<STATE> f(fDimension);
    for (int e=0; e<fDimension; e++) {
        f[e] = 0.;
    }
    
    TPZFMatrix<STATE> phiVi(fDimension,1,0.0),phiVj(fDimension,1,0.0);
    
    TPZFNMatrix<100,STATE> divphi;
    STATE divu;
    this->ComputeDivergenceOnMaster(datavec, divphi, divu);
    
    
    for(int i = 0; i < nshapeV; i++ )
    {
        int iphi = datavec[vindex].fVecShapeIndex[i].second;
        int ivec = datavec[vindex].fVecShapeIndex[i].first;
        TPZFNMatrix<4> GradVi(fDimension,fDimension);
        for (int e=0; e<fDimension; e++) {
            phiVi(e,0) = phiV(iphi,0)*datavec[vindex].fDeformedDirections(e,ivec);
            for (int f=0; f<fDimension; f++) {
                GradVi(e,f) = datavec[vindex].fDeformedDirections(e,ivec)*dphiVx(f,iphi);
                
            }
            
        }
        
        // matrix A - velocity * test-funtion velocity
        for(int j = 0; j < nshapeV; j++){
            int jphi = datavec[vindex].fVecShapeIndex[j].second;
            int jvec = datavec[vindex].fVecShapeIndex[j].first;
            
            TPZFNMatrix<4> GradVj(fDimension,fDimension);
            for (int e=0; e<fDimension; e++) {
                phiVj(e,0) = phiV(jphi,0)*datavec[vindex].fDeformedDirections(e,jvec);
                
            }
            
            STATE val = InnerVec(phiVi, phiVj);
            ek(i,j) += weight * (1./fk) * val ;
            
        }
        
        
        
        // matrix B and Bt - pressure * test-funtion velocity and symmetry
        for (int j = 0; j < nshapeP; j++) {
            
            TPZManVector<REAL,3> GradPj(fDimension);
            for (int e=0; e<fDimension; e++) {
                GradPj[e] = dphiPx(e,j);
            }
            
            STATE fact;
            {
                fact  = weight * phiP(j,0) * divphi(i,0); ///p*div(U)
                
            }
            
            
            // colocar vectoriais vezes pressao
            // Matrix B
            ek(i, nshapeV+j) += -fact;
            
            // colocar pressao vezes vectoriais
            // Matrix B^T
            ek(nshapeV+j,i) += fact;
            
            
        }
        
        
    }
    

    if(this->HasForcingFunction()){
        
        TPZFMatrix<STATE> gradu;
        
        this->ForcingFunction()->Execute(datavec[vindex].x, f, gradu);
    }
    else{
        f[0] = 0.0;
    }
    
    for (int i = 0; i < nshapeP; i++) {
        
        STATE factf= weight * phiP(i,0)*f[0];
        ef(nshapeV+i,0) += factf;
        
    }
 
    
}


void TPZDarcyPMaterial::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    
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
    
    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &phiV = datavec[vindex].phi;
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    
    // Getting the linear combination or finite element approximations
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    //Adaptação para Hdiv
//    int ekr= ek.Rows();
  
    //Vefifica se HDiv
//    if(ekr!=nshapeP+nshapeV){
//        nshapeV=nshapeV/2;
//    }
    
    
    int gy=v_h.size();
    
    
    TPZFNMatrix<9> phiVi(fDimension,1), phiVj(fDimension,1), phiPi(fDimension,1),phiPj(fDimension,1);
    
    TPZFMatrix<STATE> v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    switch (bc.Type()) {
        case 0: //Dirichlet for continuous formulation
        {
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                TPZFNMatrix<9,STATE> gradu(2,1);
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D = vbc[2];
            }
            
            if(fSpace==1){
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    
                    //Adaptação para Hdiv
                    
                    TPZManVector<REAL> n = datavec[0].normal;
                    
                    REAL vh_n = v_h[0];
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                    
                    ef(i,0) += -weight * gBigNumber * (vh_n - v_n) * phiV(i,0);
                    
                    for(int j = 0; j < nshapeV; j++){
                        
                        ek(i,j) += weight * gBigNumber * phiV(j,0) * phiV(i,0);
                        
                    }
                    
                }
                
            }else{
                
                for(int i = 0; i < nshapeV; i++ )
                {
                    int iphi = datavec[vindex].fVecShapeIndex[i].second;
                    int ivec = datavec[vindex].fVecShapeIndex[i].first;
                    
                    for (int e=0; e<fDimension; e++) {
                        phiVi(e,0)=datavec[vindex].fDeformedDirections(e,ivec)*phiV(iphi,0);
                    }
                    
                    
                    //Adaptação para Hdiv
                    
                    STATE factef=0.0;
                    for(int is=0; is<gy ; is++){
                        factef += -1.0*(v_h[is] - v_2(is,0)) * phiVi(is,0);
                    }
                    
                    ef(i,0) += weight * gBigNumber * factef;
                    
                    
                    for(int j = 0; j < nshapeV; j++){
                        int jphi = datavec[vindex].fVecShapeIndex[j].second;
                        int jvec = datavec[vindex].fVecShapeIndex[j].first;
                        
                        
                        
                        for (int e=0; e<fDimension; e++) {
                            phiVj(e,0)=datavec[vindex].fDeformedDirections(e,jvec)*phiV(jphi,0);
                        }
                        
                        //Adaptação para Hdiv
                        
                        STATE factek=0.0;
                        for(int is=0; is<gy ; is++){
                            factek += phiVj(is,0) * phiVi(is,0);
                        }
                        
                        ek(i,j) += weight * gBigNumber * factek;
                        
                    }
                    
                }
            }
            
        }
            break;
            
        case 1: //Neumann for continuous formulation
        {
            
            if(bc.HasForcingFunction())
            {
                TPZManVector<STATE> vbc(3);
                TPZFNMatrix<9,STATE> gradu(2,1);
                bc.ForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
                v_2(0,0) = vbc[0];
                v_2(1,0) = vbc[1];
                p_D = vbc[2];
            }
            
            bool strong_Neumann = false; // @omar:: why enforced Neumann data results in a wrong solution?
            if(strong_Neumann){
                
                TPZManVector<REAL,3> n = datavec[0].normal;
                
                REAL Vn = n[0]*v_2(0,0) + n[1]*v_2(1,0);
                int nshapeVn = phiV.Rows();
                for(int i = 0; i < nshapeVn; i++ )
                {
                    
                    ef(i) += gBigNumber * Vn * phiV(i,0);
                    
                    for(int j = 0; j < nshapeVn; j++){
                        
                        ek(i,j) += gBigNumber * (phiV(i,0) * phiV(j,0));
                        
                    }
                    
                }
                
            }
            else{
                
                for(int i = 0; i < nshapeP; i++ )
                {
                    
                    TPZManVector<REAL> n = datavec[0].normal;
                    
                    REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
                    
                    STATE factf=(-1.) * weight * v_n * phiP(i,0) ;
                    
                    ef(i+nshapeV,0) += fTheta*factf ;
                    
                    
                }
            }
            
        }
            
            break;
            
        case 5: //Ponto pressao
        {
            p_D = bc.Val2()(0,0);
            
            
            for(int i = 0; i < nshapeP; i++ )
            {
                
                ef(i+phiV.Rows()) += 1.0 * p_D * phiP(i,0);
                
                for(int j = 0; j < nshapeP; j++){
                    
                    ek(i+phiV.Rows(),j+phiV.Rows()) += 1.0 * (phiP(i,0) * phiP(j,0));
                    
                }
                
            }
            
        }
            break;
            
            
            
        default:
        {
            std::cout << "Boundary not implemented " << std::endl;
            DebugStop();
        }
            break;
    }
    
    
    
}




////////////////////////////////////////////////////////////////////

void TPZDarcyPMaterial::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavecleft, std::map<int, TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    
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
    
    
    // Setting the phis
    // V - left
    TPZFMatrix<REAL> &dphiV1 = datavecleft[vindex].dphix;
    // V - right
    TPZFMatrix<REAL> &dphiV2 = datavecright[vindex].dphix;
    
    // P - left
    TPZFMatrix<REAL> &phiP1 = datavecleft[pindex].phi;
    TPZFMatrix<REAL> &dphiP1 = datavecleft[pindex].dphix;
    // P - right
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
    
    
    
    for(int i1 = 0; i1 < nshapeV1; i1++ )
    {
        int iphi1 = datavecleft[vindex].fVecShapeIndex[i1].second;
        int ivec1 = datavecleft[vindex].fVecShapeIndex[i1].first;
        
        
        
        TPZFNMatrix<9> GradV1ni(fDimension,1,0.),phiV1i(fDimension,1),phiV1ni(1,1,0.);
        
        
        for (int e=0; e<fDimension; e++) {
            
            phiV1i(e,0)=datavecleft[vindex].fDeformedDirections(e,ivec1)*datavecleft[vindex].phi(iphi1,0);
            phiV1ni(0,0)+=phiV1i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV1ni(e,0)+=datavecleft[vindex].fDeformedDirections(e,ivec1)*dphiVx1(f,iphi1)*normal[f];
            }
        }
        
        TPZFNMatrix<9> GradV1nj(fDimension,1,0.);
        
        
        // K12 e K21 - (trial V left) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            
            REAL fact = (1./2.) * weight * Inner(phiV1ni,phiP1j);
            
            ek(i1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i1) += fact*fTheta;
            
        }
        
        
        
        // K14 e K41 - (trial V left) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            REAL fact = (1./2.) * weight * InnerVec(phiV1ni,phiP2j);
            
            ek(i1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i1) += fact*fTheta;
            
        }
        
    }
    
    
    for(int i2 = 0; i2 < nshapeV2; i2++ ){
        
        TPZFNMatrix<9> GradV2ni(fDimension,1,0.),phiV2i(fDimension,1),phiV2ni(1,1,0.);
        
        int iphi2 = datavecright[vindex].fVecShapeIndex[i2].second;
        int ivec2 = datavecright[vindex].fVecShapeIndex[i2].first;
        
        for (int e=0; e<fDimension; e++) {
            
            phiV2i(e,0)=datavecright[vindex].fDeformedDirections(e,ivec2)*datavecright[vindex].phi(iphi2,0);
            phiV2ni(0,0)+=phiV2i(e,0)*normal[e];
            
            for (int f=0; f<fDimension; f++) {
                GradV2ni(e,0) += datavecright[vindex].fDeformedDirections(e,ivec2)*dphiVx2(f,iphi2)*normal[f];
            }
        }
        
        
        // K32 - (trial V right) * (test P left)
        for(int j1 = 0; j1 < nshapeP1; j1++){
            
            TPZFNMatrix<9> phiP1j(1,1,0.);
            phiP1j(0,0)=phiP1(j1,0);
            
            REAL fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP1j);
            
            ek(i2+nshapeV1+nshapeP1,j1+nshapeV1) += fact;
            ek(j1+nshapeV1,i2+nshapeV1+nshapeP1) += fact*fTheta;
            
        }
        
        
        // K34 - (trial V right) * (test P right)
        for(int j2 = 0; j2 < nshapeP2; j2++){
            
            TPZFNMatrix<9> phiP2j(1,1,0.);
            phiP2j(0,0)=phiP2(j2,0);
            
            REAL fact = (-1./2.) * weight * InnerVec(phiV2ni,phiP2j);
            
            ek(i2+nshapeV1+nshapeP1,j2+2*nshapeV1+nshapeP1) += fact;
            ek(j2+2*nshapeV1+nshapeP1,i2+nshapeV1+nshapeP1) += fact*fTheta;
        }
        
    }
    
    
}


void TPZDarcyPMaterial::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
    
    
    
#ifdef IsH1
    //Caso H1 -> return
    return;
#endif
    
    
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

    
    // Setting the phis
    // V
    TPZFMatrix<REAL> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<REAL> &phiP = datavec[pindex].phi;
    TPZFMatrix<REAL> &dphiP = datavec[pindex].dphix;
    
    TPZFNMatrix<220,REAL> dphiVx(fDimension,dphiV.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiV, dphiVx, datavec[vindex].axes);
    
    TPZFNMatrix<220,REAL> dphiPx(fDimension,phiP.Cols());
    TPZAxesTools<REAL>::Axes2XYZ(dphiP, dphiPx, datavec[pindex].axes);
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = datavec[vindex].fVecShapeIndex.NElements();
    
    TPZManVector<STATE> v_h = datavec[vindex].sol[0];
    TPZManVector<STATE> p_h = datavec[pindex].sol[0];
    
    //Dirichlet
    
    TPZFMatrix<STATE> v_2=bc.Val2();
    TPZFMatrix<STATE> v_1=bc.Val1();
    STATE p_D = bc.Val1()(0,0);
    
    if(bc.HasForcingFunction())
    {
        TPZManVector<STATE> vbc(3);
        TPZFNMatrix<9,STATE> gradu(2,1);
        bc.ForcingFunction()->Execute(datavec[vindex].x,vbc,gradu);
        v_2(0,0) = vbc[0];
        v_2(1,0) = vbc[1];
        p_D=vbc[2];
        
    }
    
    for(int i = 0; i < nshapeP; i++ )
    {
        
        TPZManVector<REAL> n = data.normal;
        
        REAL v_n = n[0] * v_2[0] + n[1] * v_2[1];
        
        REAL v_t = n[1] * v_2[0] + n[0] * v_2[1];
        
        if(fSpace==1){
            
            STATE factf=(-1.) * weight * v_t * phiP(i,0) ;
            
            ef(i+nshapeV,0) += fTheta*factf*0. ;
            
        }
        
        if(fSpace==3){
            
            STATE factf2=(-1.) * weight * v_n  * phiP(i,0) ;
            
            ef(i+nshapeV,0) += fTheta*factf2*0.;
            
        }
        
    }
    
    
}



////////////////////////////////////////////////////////////////////
template <typename TVar>
TVar TPZDarcyPMaterial::Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    
    //inner product of two tensors
    
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    TVar Val = 0;
    
    for(int i = 0; i < S.Cols(); i++){
        for(int j = 0; j < S.Cols(); j++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}


////////////////////////////////////////////////////////////////////
template <typename TVar>
TVar TPZDarcyPMaterial::InnerVec(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T){
    
    //inner product of two vectors
    
    
#ifdef DEBUG
    if( S.Rows() != S.Cols() || T.Cols() != T.Rows() || S.Rows() != T.Rows() ) {
        DebugStop();
    }
#endif
    
    TVar Val = 0;
    
    for(int j = 0; j < S.Cols(); j++){
        for(int i = 0; i < S.Rows(); i++){
            Val += S(i,j)*T(i,j);
        }
    }
    
    return Val;
    
}



////////////////////////////////////////////////////////////////////

STATE TPZDarcyPMaterial::Tr( TPZFMatrix<REAL> &GradU ){
    
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
void TPZDarcyPMaterial::FillVecShapeIndex(TPZMaterialData &data)
{
    data.fDeformedDirections.Resize(fDimension,fDimension);
    data.fDeformedDirections.Identity();
    data.fVecShapeIndex.Resize(fDimension*data.phi.Rows());
    for (int d=0; d<fDimension; d++) {
        for (int i=0; i<data.phi.Rows(); i++) {
            data.fVecShapeIndex[i*fDimension+d].first = d;
            data.fVecShapeIndex[i*fDimension+d].second = i;
        }
    }
}



void TPZDarcyPMaterial::Errors(TPZVec<TPZMaterialData> &data, TPZVec<REAL> &errors)
{

    TPZVec<STATE> u_exact(1);
    TPZFMatrix<STATE> du_exact(3,1);
    fExactSol->Execute(data[0].x, u_exact, du_exact);
    errors.Resize(NEvalErrors());
    errors.Fill(0.0);
    TPZManVector<STATE> Velocity, Pressure;
    Velocity.Fill(0.0);
    Pressure.Fill(0.0);
    
    this->Solution(data,VariableIndex("V"), Velocity);
    this->Solution(data,VariableIndex("P"), Pressure);
    
    int vindex = this->VIndex();
    int pindex = this->PIndex();
    
    TPZFMatrix<REAL> dudx(Dimension(),Dimension());
    TPZFMatrix<STATE> &dsol = data[vindex].dsol[0];
    TPZFMatrix<STATE> &dsolp = data[pindex].dsol[0];
    //std::cout<<dsol<<std::endl;
    
    //Adaptação feita para Hdiv
    dsol.Resize(Dimension(),Dimension());
    
    TPZFNMatrix<2,STATE> dsolxy(2,2), dsolxyp(2,1);
    TPZAxesTools<STATE>::Axes2XYZ(dsol, dsolxy, data[vindex].axes);
    TPZAxesTools<STATE>::Axes2XYZ(dsolp, dsolxyp, data[pindex].axes);
    
    
    
    int shift = 3;
    // velocity
    
    //values[2] : erro norma L2
    REAL diff, diffp;
    errors[1] = 0.;
    for(int i=0; i<Dimension(); i++) {
        diff = Velocity[i] - u_exact[i];
        errors[1]  += diff*diff;
    }
    
    ////////////////////////////////////////////////// H1 / GD
    
    if(fSpace==2){
        
        //values[2] : erro em semi norma H1
        errors[2] = 0.;
        TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
        for(int i=0; i<Dimension(); i++) {
            for(int j=0; j<Dimension(); j++) {
                S(i,j) = dsolxy(i,j) - du_exact(i,j);
            }
        }
        
        diff = Inner(S, S);
        errors[2]  += diff;
        
        //values[0] : erro em norma H1 <=> norma Energia
        errors[0]  = errors[1]+errors[2];
        
    }
    
    if(fSpace==3){
        
        //values[2] : erro em semi norma H1
        errors[2] = 0.;
        TPZFMatrix<STATE> S(Dimension(),Dimension(),0.0);
        for(int i=0; i<Dimension(); i++) {
            for(int j=0; j<Dimension(); j++) {
                S(i,j) = dsolxy(i,j) - du_exact(i,j);
            }
        }
        
        diff = Inner(S, S);
        errors[2]  += diff;
        
        //values[0] : erro em norma H1 <=> norma Energia
        errors[0]  = errors[1]+errors[2];
        
    }
    
    
    ////////////////////////////////////////////////// H1 / GD
    
    // pressure
    
    /// values[1] : eror em norma L2
    diffp = Pressure[0] - u_exact[2];
    errors[shift+1]  = diffp*diffp;
    
    // pressure gradient error ....
    
    errors[shift+2] = 0.;
    TPZFMatrix<STATE> Sp(Dimension(),1,0.0);
    for(int i=0; i<Dimension(); i++) {
        Sp(i,0) = dsolxyp(i,0) - du_exact(2,i);
    }
    
    diffp = InnerVec(Sp, Sp);
    errors[shift+2]  += diffp;
    
    //values[0] : erro em norma H1 <=> norma Energia
    errors[shift]  = errors[1+shift]+errors[2+shift];
    
    ////////////////////////////////////////////////// HDIV
    
    if(fSpace==1){
        /// erro norma HDiv
        
        STATE Div_exact=0., Div=0.;
        for(int i=0; i<Dimension(); i++) {
            Div_exact+=du_exact(i,i);
            Div+=dsolxy(i,i);
        }
        
        diff = Div-Div_exact;
        
        errors[2]  = diff*diff;
        
        errors[0]  = errors[1]+errors[2];
        
    }
    
    ////////////////////////////////////////////////// HDIV
    
}
