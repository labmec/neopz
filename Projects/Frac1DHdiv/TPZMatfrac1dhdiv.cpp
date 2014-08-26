//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzlog.h"
#include "TPZMatfrac1dhdiv.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(): TPZDiscontinuousGalerkin()
{
    fDim = -1;
    fTheta = 0.0;
    fDeltaT = 0.0;
    fmu = 0.;
}

TPZMatfrac1dhdiv::TPZMatfrac1dhdiv(int matid, int dim, REAL viscosity): TPZDiscontinuousGalerkin(matid)
{
    // Two-dimensional analysis
    
    if(dim != 1){
        DebugStop();
    }
    
    fDim = dim;
    fTheta = 0.0;
    fDeltaT = 0.0;
    fmu = viscosity;
}


TPZMatfrac1dhdiv::~TPZMatfrac1dhdiv()
{
    
}

int TPZMatfrac1dhdiv::Dimension() const {return fDim;};

int TPZMatfrac1dhdiv::NStateVariables() {return 2;}

void TPZMatfrac1dhdiv::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
    out << "Base Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}


void TPZMatfrac1dhdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    
#ifdef DEBUG
    int nref =  datavec.size();
    if (nref != 2 )
    {
        std::cout << " Erro. The size of the datavec is different from 2 \n";
        DebugStop();
    }
#endif
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiP =  datavec[1].phi;
    
    TPZFMatrix<REAL> &dphiP =   datavec[1].dphix;
    
    // number of test functions for each state variable
    int phrQ, phrP;
    phrQ = phiQ.Rows();
    phrP = phiP.Rows();
    
    // blocks
    int FirstQ  = 0;
    int FirstP  = phrQ+FirstQ;
    
    //  Getting and computing another required data
    
    // Derivadas na ordem = {Q,P,DPDX}
    MFad sol_q(datavec[0].sol[0][0],0);
    MFad sol_p(datavec[1].sol[0][0],1);
    MFad dsol_p(datavec[1].dsol[0](0,0),2);
    
    //  Contribution of domain integrals for Jacobian matrix and Residual vector
    //  n + 1 time step
    if(gState == ECurrentState)
    {
        REAL cte = 1./this->fmu;
        for (int iq = 0; iq < phrQ; iq++) {
            MFad res00 = cte * sol_q * phiQ(iq,0);
            MFad res01 = dsol_p * phiQ(iq,0);
            ef(FirstQ+iq,0) += weight * ( res00.val() + res01.val() );
            for (int jq = 0; jq < phrQ; jq++) {
                ek(FirstQ+iq,FirstQ+jq) += weight * (res00.dx(0) * phiQ(jq));
            }
            for (int jp = 0; jp < phrP; jp++) {
                ek(FirstQ+iq,FirstP+jp) += weight * (res01.dx(2) * dphiP(0,jp));
            }
        }
        for (int ip = 0; ip < phrP ; ip++) {
            MFad res10 = sol_q * dphiP(0,ip);
            ef(FirstP+ip,0) += weight * (res10.val());
            for (int jq = 0; jq < phrQ; jq++) {
                ek(FirstP+ip,FirstQ+jq) += weight * (res10.dx(0) * phiQ(jq,0));
            }
        }
    }
    //  End of contribution of domain integrals for Jacobian matrix and Residual vector
    
    
    //  n time step
    //  This values are constant in Newton iteration
    if(gState == ELastState)
    {
        // nothing yet. Will change when put dwdt
    }
    
    
}


void TPZMatfrac1dhdiv::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    //  n time step
    //  This values are constant in Newton iteration
    if(gState == ELastState)
    {
        return;
    }
    
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiQR = dataright[0].phi;
    TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
    TPZFMatrix<REAL> &phiPR = dataright[1].phi;
    
    // number of test functions for each state variable
    int phrQL, phrQR, phrPL, phrPR;
    phrQL = phiQL.Rows();
    phrPL = phiPL.Rows();
    phrQR = phiQR.Rows();
    phrPR = phiPR.Rows();
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    // Derivadas na ordem = {QLeft,PLeft,DPDXLeft,QRight,PRight,DPDXRight}
    MFad sol_qL(dataleft[0].sol[0][0],0);
    MFad sol_pL(dataleft[1].sol[0][0],1);
    MFad dsol_pL(dataleft[1].dsol[0](0,0),2);
    MFad sol_qR(dataright[0].sol[0][0],3);
    MFad sol_pR(dataright[1].sol[0][0],4);
    MFad dsol_pR(dataright[1].dsol[0](0,0),5);
    
    // blocks
    int FirstQL  = 0;
    int FirstPL  = phrQL+FirstQL;

    int FirstQR  = 0;
    int FirstPR  = phrQR+FirstQR;
    
    int ileftblock  = FirstPL + phrPL;
    int jrightblock = FirstPR + phrPR;
    
    //  Contribution of contour integrals for Jacobian matrix and Residual vector
    
    for (int iql = 0; iql < phrQL; iql++) {
        MFad res01 = (sol_pR-sol_pL) * phiQL(iql,0);
        ef(FirstQL+iql,0) += weight * (res01.val());
        for (int jpl = 0; jpl < phrPL; jpl++)
        {
            ek(FirstQL+iql,FirstPL+jpl) += weight * (res01.dx(1) * phiPL(jpl,0));
        }
        for (int jpr = 0; jpr < phrPR; jpr++)
        {
            ek(FirstQL+iql,jrightblock + FirstPR+jpr) += weight * (res01.dx(3) * phiPR(jpr,0));
        }
    }
    
    
    for (int jpl = 0; jpl < phrPL; jpl++) {
        MFad res10 = (phiPL(jpl,0) * sol_qL);
        for (int iql = 0; iql < phrQL; iql++) {
            ek(jpl,iql) = 0.0; //Complete the transpose term;
        }
    }
    
    
    
    // End of contribution of contour integrals for Jacobian matrix and Residual vector
    
    
    
    
}

void TPZMatfrac1dhdiv::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    // n time step
    if(gState == ELastState)
    {
        return;
    }
    
    int nref =  dataleft.size();
    if (nref != 2) {
        std::cout << " Error:: datavec size must to be equal to 2 \n" << std::endl;
        DebugStop();
    }
    if (bc.Val2().Rows() != 2)
    {
        std::cout << " Error:: This material need boundary conditions for qn, p (fluid pressure).\n" << std::endl;
        std::cout << " give me one matrix with this form Val2(2,1).\n" << std::endl;
        DebugStop();
    }
    
    if (bc.Val1().Rows() != 2)
    {
        std::cout << " Error:: This material need boundary conditions for qn, p (fluid pressure) .\n" << std::endl;
        DebugStop();
    }
    
    
    // n+1 time step
    if(gState == ECurrentState)
    {
        
        
        
    }
    
    
    switch (bc.Type()) {
        case 0 :    // Flux bc: qn = (+ or -)g (Dirichlet)
        {
            ApplyQnD(data, dataleft, weight, ek, ef, bc);
        }
        case 1 :    // Pressure bc: PN = (+)k (Neumann)
        {
            ApplyPN(data, dataleft, weight, ek, ef, bc);
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
    
    
    
}

void TPZMatfrac1dhdiv::ApplyQnD(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    
    TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
    TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
    int PRowsleft = phiPL.Rows();
    
    int iRightInterfaceBlock = QRowsleft + PRowsleft;
    int jRightInterfaceBlock = QRowsleft + PRowsleft;
    
    int FirstQL = 0;
    int FirstPL = QRowsleft + FirstQL;
    
    TPZManVector<REAL,3> &normal = data.normal;
    REAL n1 = normal[0];
    REAL n2 = normal[1];

}


void TPZMatfrac1dhdiv::ApplyPN        (TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<> &ek,TPZFMatrix<> &ef,TPZBndCond &bc)
{
    
    TPZFMatrix<REAL> &phiQL     = dataleft[0].phi;
    TPZFMatrix<REAL> &phiPL     = dataleft[1].phi;
    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
    int PRowsleft = phiPL.Rows();
    
    int iRightInterfaceBlock = QRowsleft + PRowsleft;
    int jRightInterfaceBlock = QRowsleft + PRowsleft;
    
    int FirstQL = 0;
    int FirstPL = QRowsleft + FirstQL;
    
    TPZManVector<REAL,3> &normal = data.normal;
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    
}

/** Returns the variable index associated with the name */
int TPZMatfrac1dhdiv::VariableIndex(const std::string &name){
    if(!strcmp("BulkVelocity",name.c_str()))        return  1;
    if(!strcmp("Pressure",name.c_str()))    return  2;
    
    return TPZMaterial::VariableIndex(name);
}

int TPZMatfrac1dhdiv::NSolutionVariables(int var){
    if(var == 1) return 2;
    if(var == 2) return 1;
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMatfrac1dhdiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    TPZVec<STATE> SolQ, SolP;
    SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 1){ //function (state variable Q)
        Solout[0] = SolQ[0];
        Solout[1] = SolQ[1];
        return;
    }
    
    if(var == 2){
        Solout[0] = SolP[0];//function (state variable p)
        return;
    }
}


void TPZMatfrac1dhdiv::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
    int nref = datavec.size();
    for(int i = 0; i<nref; i++ )
    {
        datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNeighborSol = true;
        datavec[i].fNeedsNeighborCenter = false;
        datavec[i].fNeedsNormal = true;
    }
    
}

void TPZMatfrac1dhdiv::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;        
    }
}
