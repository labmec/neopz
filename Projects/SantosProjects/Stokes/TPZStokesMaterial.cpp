/*
 *  TPZStokesMaterial.cpp
 *  PZ
 *
 *  Created by Thiago Dias dos Santos on 12/01/2015.
 *  Copyright 2015 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZStokesMaterial.h"
#include "pzbndcond.h"
#include "pzaxestools.h"


TPZStokesMaterial::TPZStokesMaterial() : TPZDiscontinuousGalerkin(){
    
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(int matid) : TPZDiscontinuousGalerkin(matid){
    
    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::TPZStokesMaterial(const TPZStokesMaterial &mat) : TPZDiscontinuousGalerkin(mat){

    
}

////////////////////////////////////////////////////////////////////

TPZStokesMaterial::~TPZStokesMaterial(){
    
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec)
{
    int ndata = datavec.size();
    for (int idata=0; idata < ndata ; idata++) {
        datavec[idata].SetAllRequirements(false);
        datavec[idata].fNeedsSol = true;
        datavec[idata].fNeedsNormal = true;
    }
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Print(std::ostream &out) {
    out << "\t Base class print:\n";
    out << " name of material : " << this->Name() << "\n";
    TPZMaterial::Print(out);
}

////////////////////////////////////////////////////////////////////

int TPZStokesMaterial::VariableIndex(const std::string &name) {
    
    if (!strcmp("Pressure", name.c_str())) return this->PIndex();
    if (!strcmp("Velocity", name.c_str())) return this->VIndex();
   
    std::cout  << " Var index not implemented " << std::endl;
    DebugStop();
    return 0;
}

////////////////////////////////////////////////////////////////////

int TPZStokesMaterial::NSolutionVariables(int var) {
    
    switch(var) {
        
        case 0:
            return 1; // Pressure, Scalar
        case 1:
            return this->Dimension(); // Velocity, Vector
    
        default:
        {
            std::cout  << " Var index not implemented " << std::endl;
            DebugStop();
        }
    }
    return 0;
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout) {
    
    
    //itapopo conferir esse metodo
    
    int Vblock = this->VIndex();
    int Pblock = this->PIndex();
    
    TPZManVector<REAL,3> V = datavec[Vblock].sol[0];
    
    REAL P = datavec[Pblock].sol[0][0];
    
    Solout.Resize(this->NSolutionVariables(var));
    
    switch(var) {
        
        case 0: //Pressure
        {
            Solout[0] = P;
        }
            break;
        
        case 1: //Velocity
        {
            Solout[0] = V[0]; // Vx
            Solout[1] = V[1]; // Vy
            Solout[2] = V[2]; // Vy
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
void TPZStokesMaterial::ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi)
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

void TPZStokesMaterial::Write(TPZStream &buf, int withclassid) const{
    
    TPZDiscontinuousGalerkin::Write(buf, withclassid);
 
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::Read(TPZStream &buf, void *context) {
    
    TPZDiscontinuousGalerkin::Read(buf, context);

}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<STATE> > &GradPhi){
    
   
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

////////////////////////////////////////////////////////////////////
void TPZStokesMaterial::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){

#ifdef DEBUG
    //2 = 1 Vel space + 1 Press space
    int nref =  datavec.size();
    if (nref != 2 ) {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    const int pindex = this->PIndex();
    const int vindex = this->VIndex();
    
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
    TPZFMatrix<STATE> &phiV = datavec[vindex].phi;
    TPZFMatrix<STATE> &dphiV = datavec[vindex].dphix;
    // P
    TPZFMatrix<STATE> &phiP = datavec[pindex].phi;
    TPZFMatrix<STATE> &dphiP = datavec[pindex].dphix;
    
    int nshapeV, nshapeP;
    nshapeP = phiP.Rows();
    nshapeV = phiV.Rows(); //datavec[0].fVecShapeIndex.NElements();
    
    TPZVec<TPZFMatrix<STATE> > GradU;
    this->FillGradPhi(datavec[vindex], GradU);
    
    const STATE Visc = 1.; //itapopo
    
    // Integral value - Matrix A and B
    for(int i = 0; i < nshapeV; i++){
        
        // matrix A - gradV
        for(int j = 0; j < nshapeV; j++){
          
            //itapopo verificar termo simétrico
            ek(i,j) += 2. * weight * Visc * Inner( GradU[j], GradU[i] ) ; ///Visc*(GradU+GradU^T):GradPhi
            
        }//j
        
        // matrix B - pressure and velocity
        for (int j = 0; j < nshapeP; j++) {
            
            STATE fact = (-1.) * weight * phiP(j,0) * Tr( GradU[i] ); ///p*div(U)
            
            // Matrix B
            ek(i, nshapeV+j) += fact;
            
            // Matrix B^T
            ek(nshapeV+j,i) += fact;
        }//j
        
        // force vector
        ef(i,0) += (-1.)*weight*force*phiV(i,0);//itapopo conferir termo e sinal

    }//i
    
    
}

////////////////////////////////////////////////////////////////////

void TPZStokesMaterial::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    const int pindex = this->PIndex();
    const int vindex = this->VIndex();
    
    
    TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
    TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    TPZManVector<REAL,3> &normal = data.normal;
    
    TPZFNMatrix<660> dphiL, dphiR;
    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
    
    int &LeftPOrder=dataleft.p;
    int &RightPOrder=dataright.p;
    
    REAL &faceSize=data.HSize;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    int il,jl,ir,jr,id;
    
    //Convection term
    REAL ConvNormal = 0.;
    for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id] * normal[id];
    if(ConvNormal > 0.) {
        for(il=0; il<nrowl; il++) {
            for(jl=0; jl<nrowl; jl++) {
                ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
            }
        }
        for(ir=0; ir<nrowr; ir++) {
            for(jl=0; jl<nrowl; jl++) {
                ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl);
            }
        }
    } else {
        for(ir=0; ir<nrowr; ir++) {
            for(jr=0; jr<nrowr; jr++) {
                ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr);
            }
        }
        for(il=0; il<nrowl; il++) {
            for(jr=0; jr<nrowr; jr++) {
                ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr);
            }
        }
    }
    
    if(IsZero(fK)) return;
    //diffusion term
    STATE leftK, rightK;
    leftK  = this->fK;
    rightK = this->fK;
    
    // 1) phi_I_left, phi_J_left
    for(il=0; il<nrowl; il++) {
        REAL dphiLinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiLinormal += dphiL(id,il)*normal[id];
        }
        for(jl=0; jl<nrowl; jl++) {
            REAL dphiLjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiLjnormal += dphiL(id,jl)*normal[id];
            }
            ek(il,jl) += (STATE)(weight * ( this->fSymmetry * (0.5)*dphiLinormal*phiL(jl,0)-(0.5)*dphiLjnormal*phiL(il,0))) * leftK;
        }
    }
    
    // 2) phi_I_right, phi_J_right
    for(ir=0; ir<nrowr; ir++) {
        REAL dphiRinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiRinormal += dphiR(id,ir)*normal[id];
        }
        for(jr=0; jr<nrowr; jr++) {
            REAL dphiRjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiRjnormal += dphiR(id,jr)*normal[id];
            }
            ek(ir+nrowl,jr+nrowl) += (STATE)(weight * (this->fSymmetry * ((-0.5) * dphiRinormal * phiR(jr) ) + (0.5) * dphiRjnormal * phiR(ir))) * rightK;
        }
    }
    
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        REAL dphiLinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiLinormal += dphiL(id,il)*normal[id];
        }
        for(jr=0; jr<nrowr; jr++) {
            REAL dphiRjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiRjnormal += dphiR(id,jr)*normal[id];
            }
            ek(il,jr+nrowl) += (STATE)weight * ((STATE)fSymmetry * ((STATE)((-0.5) * dphiLinormal * phiR(jr)) * leftK ) - (STATE)((0.5) * dphiRjnormal * phiL(il))* rightK );
        }
    }
    
    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        REAL dphiRinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiRinormal += dphiR(id,ir)*normal[id];
        }
        for(jl=0; jl<nrowl; jl++) {
            REAL dphiLjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiLjnormal += dphiL(id,jl)*normal[id];
            }
            ek(ir+nrowl,jl) += (STATE)weight * (
                                                (STATE)(fSymmetry * (0.5) * dphiRinormal * phiL(jl)) * rightK + (STATE)((0.5) * dphiLjnormal * phiR(ir)) * leftK
                                                );
        }
    }
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
    
    if (this->fPenaltyConstant == 0.) return;
    
    leftK  = this->fK;
    rightK = this->fK;
    
    
    
    //penalty = <A p^2>/h
    REAL penalty = fPenaltyConstant * (0.5 * (abs(leftK)*LeftPOrder*LeftPOrder + abs(rightK)*RightPOrder*RightPOrder)) / faceSize;
    
    if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
        
        // 1) left i / left j
        for(il=0; il<nrowl; il++) {
            for(jl=0; jl<nrowl; jl++) {
                ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
            }
        }
        
        // 2) right i / right j
        for(ir=0; ir<nrowr; ir++) {
            for(jr=0; jr<nrowr; jr++) {
                ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
            }
        }
        
        // 3) left i / right j
        for(il=0; il<nrowl; il++) {
            for(jr=0; jr<nrowr; jr++) {
                ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
            }
        }
        
        // 4) right i / left j
        for(ir=0; ir<nrowr; ir++) {
            for(jl=0; jl<nrowl; jl++) {
                ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
            }
        }
        
    }
    
    if (this->fPenaltyType == EFluxPenalty || this->fPenaltyType == EBoth){
        
        REAL NormalFlux_i = 0.;
        REAL NormalFlux_j = 0.;
        
        // 1) left i / left j
        for(il=0; il<nrowl; il++) {
            NormalFlux_i = 0.;
            for(id=0; id<fDim; id++) {
                NormalFlux_i += dphiL(id,il)*normal[id];
            }
            for(jl=0; jl<nrowl; jl++) {
                NormalFlux_j = 0.;
                for(id=0; id<fDim; id++) {
                    NormalFlux_j += dphiL(id,jl)*normal[id];
                }
                ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
            }
        }
        
        // 2) right i / right j
        for(ir=0; ir<nrowr; ir++) {
            NormalFlux_i = 0.;
            for(id=0; id<fDim; id++) {
                NormalFlux_i += dphiR(id,ir)*normal[id];
            }
            for(jr=0; jr<nrowr; jr++) {
                NormalFlux_j = 0.;
                for(id=0; id<fDim; id++) {
                    NormalFlux_j += dphiR(id,jr)*normal[id];
                }      
                ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
            }
        }
        
        // 3) left i / right j
        for(il=0; il<nrowl; il++) {
            NormalFlux_i = 0.;
            for(id=0; id<fDim; id++) {
                NormalFlux_i += dphiL(id,il)*normal[id];
            }
            for(jr=0; jr<nrowr; jr++) {
                NormalFlux_j = 0.;
                for(id=0; id<fDim; id++) {
                    NormalFlux_j += dphiR(id,jr)*normal[id];
                }      
                ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
            }
        }
        
        // 4) right i / left j
        for(ir=0; ir<nrowr; ir++) {
            NormalFlux_i = 0.;
            for(id=0; id<fDim; id++) {
                NormalFlux_i += dphiR(id,ir)*normal[id];
            }
            for(jl=0; jl<nrowl; jl++) {
                NormalFlux_j = 0.;
                for(id=0; id<fDim; id++) {
                    NormalFlux_j += dphiL(id,jl)*normal[id];
                }
                ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
            }
        }
        
    }
    
    
}

////////////////////////////////////////////////////////////////////

STATE TPZStokesMaterial::Inner(TPZFMatrix<STATE> &S, TPZFMatrix<STATE> &T){
    
    //inner product of two tensors
    
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

STATE TPZStokesMaterial::Tr( TPZFMatrix<REAL> &GradU ){
 
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
