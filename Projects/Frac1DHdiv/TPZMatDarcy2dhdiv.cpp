//
//  TPZMatDarcy2dhdiv.cpp
//  PZ
//
//  Created by omar and nathan on 9/3/14.
//
//

#include "TPZMatDarcy2dhdiv.h"

#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"

#include <iostream>

#ifdef PZ_LOG
static PZLogger logger("pz.multiphase");
#endif

#ifdef PZ_LOG
static PZLogger logdata("pz.material.multiphase.data");
#endif

TPZMatDarcy2dhdiv::TPZMatDarcy2dhdiv(): TPZMaterial()
{
    fDim = 2;
    fNotContribute = false;
}

TPZMatDarcy2dhdiv::TPZMatDarcy2dhdiv(int matid): TPZMaterial(matid)
{
    fDim = 2;
    fNotContribute = false;
}


TPZMatDarcy2dhdiv::~TPZMatDarcy2dhdiv()
{
}

int TPZMatDarcy2dhdiv::Dimension() const {return fDim;};


void TPZMatDarcy2dhdiv::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
    out << "Base Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}


// Contribute methods

void TPZMatDarcy2dhdiv::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
  DebugStop(); // Should be datavec
}


void TPZMatDarcy2dhdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
  
  if(fNotContribute) return;
    
#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 2 )
    {
        std::cout << " Error. datavec size is different from 2 \n";
        DebugStop();
    }
#endif
    
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiQ     =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiP     =  datavec[1].phi;
    
    TPZFMatrix<REAL> &dphiP     =  datavec[1].dphix;
    
    // number of test functions for each state variable
    int phrQ, phrP;
    phrQ = datavec[0].fVecShapeIndex.NElements();
    phrP = phiP.Rows();
    
    // blocks
    int FirstQ  = 0;
    int FirstP  = phrQ + FirstQ;
    
    
    
    //  Getting and computing another required data
    REAL TimeStep = fData->TimeStep();
    REAL Theta = fData->Theta();
    TPZFMatrix<STATE> Kabsolute = fData->K();
    TPZFMatrix<STATE> Kinverse = fData->Kinv();

    TPZManVector<STATE,3> sol_q =    datavec[0].sol[0];
    TPZManVector<STATE,3> sol_p =    datavec[1].sol[0];
    
    TPZFMatrix<STATE> dsol_q =datavec[0].dsol[0];
    TPZFMatrix<STATE> dsol_p =datavec[1].dsol[0];
    REAL Pressure = sol_p[0];
    
    
    REAL rockporosity, oildensity, oilviscosity;
    REAL drockporositydp, doildensitydp, doilviscositydp;
    fData->Porosity(Pressure, rockporosity, drockporositydp);
    fData->Density(Pressure, oildensity, doildensitydp);
    fData->Viscosity(Pressure, oilviscosity, doilviscositydp);
    
    REAL bulklambda, dbulklambdadp;
    bulklambda = (oildensity/oilviscosity);
    dbulklambdadp = /*-1.0*((oildensity*doilviscositydp)/(oilviscosity*oilviscosity)) + */ (doildensitydp/oilviscosity);
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of domain integrals for Residual Vector
    //  n time step
    if(fData->IsLastState())
    {
        
        //  Second Block (Equation Two) Bulk flux  equation
        for(int ip=0; ip < phrP; ip++)
        {
            
            // Integrate[W*(d(\phi*(rho)/dt)), Omega_{e}] (Equation Two)
            REAL Integrating = phiP(ip,0) * rockporosity * (oildensity);
            ef(ip + FirstP) += (1.0) * weight * Integrating;
            
            
            //  Second Block (Equation Two) Bulk flux  equation
            // Integrate[dot(grad(W),q), Omega_{e}] (Equation Two)
            //  Compute grad(W)
            TPZManVector<STATE> dphip(2,0);
            dphip[0] = dphiP(0,ip)*datavec[1].axes(0,0)+dphiP(1,ip)*datavec[1].axes(1,0);
            dphip[1] = dphiP(0,ip)*datavec[1].axes(0,1)+dphiP(1,ip)*datavec[1].axes(1,1);
            
            REAL dotprod =
            (dphip[0]) * (sol_q[0]) +
            (dphip[1]) * (sol_q[1]);
            
            ef(ip + FirstP) += (1.0-Theta) * (TimeStep) * weight * dotprod;

        }
        
        return;
        
    }
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  End of contribution of domain integrals for Residual Vector
    
    
    
        //  ////////////////////////// Jacobian Matrix and Residual Vector ///////////////////////////////////
        //  Contribution of domain integrals for Jacobian matrix and Residual Vector
        // n+1 time step
    
        //  First Block (Equation One) constitutive law
        // Integrate[(Kinv/bulklambda)*dot(v,v), Omega_{e} ]  (Equation One)
        REAL OneOverLambda = 1.0/bulklambda;
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[0].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[0].fVecShapeIndex[iq].second;
                
            REAL Kqdotv =
            (Kinverse(0,0)*sol_q[0]+Kinverse(0,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(0,ivectorindex)) +
            (Kinverse(1,0)*sol_q[0]+Kinverse(1,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(1,ivectorindex)) ;

            ef(iq + FirstQ) +=  OneOverLambda * weight * Kqdotv;    //  dot(K q,v)

            
            for (int jq=0; jq<phrQ; jq++)
            {
                int jvectorindex    = datavec[0].fVecShapeIndex[jq].first;
                int jshapeindex     = datavec[0].fVecShapeIndex[jq].second;

                REAL vec1 = (Kinverse(0,0)*datavec[0].fDeformedDirections(0,ivectorindex)+Kinverse(0,1)*datavec[0].fDeformedDirections(1,ivectorindex));
                REAL vec2 = (Kinverse(1,0)*datavec[0].fDeformedDirections(0,ivectorindex)+Kinverse(1,1)*datavec[0].fDeformedDirections(1,ivectorindex));

                REAL Kvdotv =
                (phiQ(ishapeindex,0) * vec1) * (phiQ(jshapeindex,0)*datavec[0].fDeformedDirections(0,jvectorindex)) +
                (phiQ(ishapeindex,0) * vec2) * (phiQ(jshapeindex,0)*datavec[0].fDeformedDirections(1,jvectorindex)) ;    //  dot(K vj,vi)

                ek(iq + FirstQ,jq + FirstQ) += weight * OneOverLambda * Kvdotv;
            }
            
            //  First Block (Equation One) constitutive law
            // Integrate[(d(1/bulklambdal)/dP)*dot(q,v), Omega_{e} ]    (Equation One)
            for (int jp=0; jp<phrP; jp++)
            {
                ek(iq + FirstQ,jp + FirstP) += (-1.0) * weight * dbulklambdadp  * OneOverLambda * OneOverLambda * phiP(jp,0) * Kqdotv;
            }
            
            
            //  First Block (Equation One) constitutive law
            //  Integrate [ K dot(v,grad(P)) , Omega_{e}]   (Equation One)
            //  Compute grad(P)
            TPZManVector<STATE> dsolp(2,0);
            dsolp[0] = dsol_p(0,0)*datavec[1].axes(0,0)+dsol_p(1,0)*datavec[1].axes(1,0);
            dsolp[1] = dsol_p(0,0)*datavec[1].axes(0,1)+dsol_p(1,0)*datavec[1].axes(1,1);

            REAL gradPdotv   =   (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(0,ivectorindex))*(dsolp[0]) +
                                 (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(1,ivectorindex))*(dsolp[1]);
            
            ef(iq + FirstQ) += weight * (gradPdotv);
            
            //  First Block (Equation One) constitutive law
            //  Integrate [ K dot(v,grad(P)) , Omega_{e}]   (Equation One)
            for (int jp=0; jp<phrP; jp++)
            {
                //  Compute grad(W)
                TPZManVector<STATE> dphip(2,0);
                dphip[0] = dphiP(0,jp)*datavec[1].axes(0,0)+dphiP(1,jp)*datavec[1].axes(1,0);
                dphip[1] = dphiP(0,jp)*datavec[1].axes(0,1)+dphiP(1,jp)*datavec[1].axes(1,1);
                
                REAL gradPhiPdotv   =   (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(0,ivectorindex))*(dphip[0]) +
                                        (phiQ(ishapeindex,0)*datavec[0].fDeformedDirections(1,ivectorindex))*(dphip[1]);
                
                ek(iq + FirstQ,jp + FirstP) += weight * gradPhiPdotv;
                
            }
            
        }
        

        //  Second Block (Equation Two) Bulk flux  equation
        for(int ip=0; ip < phrP; ip++)
        {
            
            // Integrate[W*(d(\phi*(rho)/dt)), Omega_{e}] (Equation Two)
            REAL Integrating = phiP(ip,0) * rockporosity * (oildensity);
            ef(ip + FirstP) += (-1.0) * weight * Integrating;
            
            // d(porosity)/dPalpha and d(oildensity)/dPalpha
            for (int jp=0; jp < phrP; jp++)
            {
                REAL Integrating = phiP(ip,0) * (drockporositydp * oildensity + rockporosity * doildensitydp) * phiP(jp,0) ;
                ek(ip + FirstP,jp + FirstP) +=  (-1.0) * weight * Integrating;
            }
            
            
            //  Second Block (Equation Two) Bulk flux  equation
            // Integrate[dot(grad(W),q), Omega_{e}] (Equation Two)
            //  Compute grad(W)
            TPZManVector<STATE> dphip(2,0);
            dphip[0] = dphiP(0,ip)*datavec[1].axes(0,0)+dphiP(1,ip)*datavec[1].axes(1,0);
            dphip[1] = dphiP(0,ip)*datavec[1].axes(0,1)+dphiP(1,ip)*datavec[1].axes(1,1);
            
            REAL gradwdotq =
            (dphip[0]) * (sol_q[0]) +
            (dphip[1]) * (sol_q[1]);
            
            ef(ip + FirstP) += (Theta) * (TimeStep) * weight * gradwdotq;
            
            // d(q)/dQalpha
            for (int jq=0; jq<phrQ; jq++)
            {
                
                int jvectorindex    = datavec[0].fVecShapeIndex[jq].first;
                int jshapeindex     = datavec[0].fVecShapeIndex[jq].second;
                
                REAL gradwdotv =
                (dphip[0]) * (phiQ(jshapeindex,0)*datavec[0].fDeformedDirections(0,jvectorindex)) +
                (dphip[1]) * (phiQ(jshapeindex,0)*datavec[0].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstP,jq + FirstQ) += (Theta) * (TimeStep) * weight * gradwdotv;
            }
            
        }
        
        //  ////////////////////////// Jacobian Matrix and Residual Vector ///////////////////////////////////
        //  End of contribution of domain integrals for Jacobian matrix and Residual Vector
    
    

    
}

//  Residual vector contribution
void TPZMatDarcy2dhdiv::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
  TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.);
  this->Contribute(datavec, weight, ekfake, ef);
}


void TPZMatDarcy2dhdiv::ContributeInterface(TPZVec<TPZMaterialData> &datavec,TPZVec<TPZMaterialData> &dataleftvec,TPZVec<TPZMaterialData> &datarightvec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(fNotContribute) return;
  
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiQR = dataright[0].phi;
    
    TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
    TPZFMatrix<REAL> &phiPR = dataright[1].phi;
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    //  REAL n3 = normal[2];
    
    TPZManVector<STATE,3> sol_qL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_qR =dataright[0].sol[0];
    
    TPZManVector<STATE,3> sol_pL =dataleft[1].sol[0];
    TPZManVector<STATE,3> sol_pR =dataright[1].sol[0];
    
    
    //  Getting Q solution for left and right side
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL qxR = sol_qR[0];
    REAL qyR = sol_qR[1];
    REAL qnL = (qxL*n1) + (qyL*n2);
    REAL qnR = (qxR*n1) + (qyR*n2);

    //  Getting P solution for left and right side
    REAL PressureL = sol_pL[0];
    REAL PressureR = sol_pR[0];

    //  Getting another required data
    REAL TimeStep = fData->TimeStep();
    REAL Theta = fData->Theta();
    
    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
    int QRowsRight = dataright[0].fVecShapeIndex.NElements();
    
    int PRowsleft = phiPL.Rows();
    int PRowsRight = phiPR.Rows();
    
    int iRightInterfaceBlock = QRowsleft + PRowsleft;
    int jRightInterfaceBlock = QRowsleft + PRowsleft;
    
    int FirstQL = 0;
    int FirstPL = QRowsleft + FirstQL;
    
    int FirstQR = 0;
    int FirstPR = QRowsRight + FirstQR;

    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of contour integrals for Residual Vector
    //  Time step n
    
    if(fData->IsLastState())
    {
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            
            ef(ip + FirstPL) += (-1.0) * (1.0-Theta) * (TimeStep) * weight * qnL * phiPL(ip,0);
            
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Right-Right Part
        for (int ip=0; ip < PRowsRight; ip++)
        {
            
            ef(ip + FirstPR + iRightInterfaceBlock) += (1.0) * (1.0-Theta) * (TimeStep) * weight * qnR * phiPR(ip,0);

        }
        
        return;
        
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  End of contribution of contour integrals for Residual Vector
    
    
        //  ////////////////////////// Jacobian Matrix and Residual Vector ///////////////////////////////////
        //  Contribution of contour integrals for Jacobian matrix and Residual Vector
        //  Time step n+1
        
        
        //  First Block (Equation One) constitutive law
        //  Integrate[L dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
        
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;
            
            REAL vnL   = (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex))*(n1) +
                         (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex))*(n2);
            
            ef(iq + FirstQL) += (-1.0) * weight * PressureL * vnL;
            
            for (int jp=0; jp < PRowsleft; jp++)
            {
                ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * phiPL(jp,0) * vnL;
            }

        }
        
        for (int iq=0; iq < QRowsRight; iq++)
        {
            int iRvectorindex       = dataright[0].fVecShapeIndex[iq].first;
            int iRshapeindex        = dataright[0].fVecShapeIndex[iq].second;
            
            REAL vnR   = (phiQR(iRshapeindex,0)*dataright[0].fDeformedDirections(0,iRvectorindex))*(n1) +
                         (phiQR(iRshapeindex,0)*dataright[0].fDeformedDirections(1,iRvectorindex))*(n2);
            
            ef(iq + iRightInterfaceBlock + FirstQR) += (1.0) * weight * PressureR * vnR;
            
            for (int jp=0; jp < PRowsRight; jp++)
            {
                ek(iq + FirstQR + iRightInterfaceBlock,jp + FirstPR + jRightInterfaceBlock) +=  (1.0) * weight * phiPR(jp,0) * vnR;
                
            }
            
        }
        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
         
            ef(ip + FirstPL) += (-1.0) * (Theta) * (TimeStep) * weight * qnL * phiPL(ip,0);
            
            for (int jq=0; jq<QRowsleft; jq++)
            {
                
                int jvectorindex    = dataleft[0].fVecShapeIndex[jq].first;
                int jshapeindex     = dataleft[0].fVecShapeIndex[jq].second;
                
                REAL vnL =
                (n1) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(0,jvectorindex)) +
                (n2) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstPL,jq + FirstQL) += (-1.0) * (Theta) * (TimeStep) * weight * phiPL(ip,0) * vnL;
                
            }
            
        }
        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Right-Right Part
        for (int ip=0; ip < PRowsRight; ip++)
        {
            
            ef(ip + FirstPR + iRightInterfaceBlock) += (1.0) * (Theta) * (TimeStep) * weight * qnR * phiPR(ip,0);
            
            for (int jq=0; jq<QRowsRight; jq++)
            {
                
                int jvectorindex    = dataright[0].fVecShapeIndex[jq].first;
                int jshapeindex     = dataright[0].fVecShapeIndex[jq].second;
                
                REAL vnR =
                (n1) * (phiQR(jshapeindex,0)*dataright[0].fDeformedDirections(0,jvectorindex)) +
                (n2) * (phiQR(jshapeindex,0)*dataright[0].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstPR + iRightInterfaceBlock,jq + FirstQR + jRightInterfaceBlock) += (1.0) * (Theta) * (TimeStep) * weight * phiPR(ip,0) * vnR;
                
            }
            
        }

    
    //  ////////////////////////// Jacobian Matrix and Residual Vector ///////////////////////////////////
    //  End of contribution of countour integrals for Jacobian matrix and Residual Vector

    
    
}



void TPZMatDarcy2dhdiv::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
  
  TPZFMatrix<STATE> ekfake(ef.Rows(),ef.Rows(),0.);
  this->ContributeInterface(data, dataleft, dataright, weight, ekfake, ef);
  
}

void TPZMatDarcy2dhdiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    

    return; // This method is called but shouldn t do anything
    
}

void TPZMatDarcy2dhdiv::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    DebugStop();
}

void TPZMatDarcy2dhdiv::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    if(fNotContribute) return;
  
    int nref =  dataleft.size();
    if (nref != 2) {
        std::cout << " Error:: datavec size must to be equal to 4 \n" << std::endl;
        DebugStop();
    }
    
    if (bc.Val2().Rows() != 3){
        std::cout << " Error:: This material need boundary conditions for qx, qy, p (pore pressure) and s (Saturation) .\n" << std::endl;
        std::cout << " give me one matrix with this form Val2(3,1).\n" << std::endl;
        DebugStop();
    }
    
    if (bc.Val1().Rows() != 3){
        std::cout << " Error:: This material need boundary conditions for ux, uy, qx, qy, p (pore pressure) and s (Saturation) .\n" << std::endl;
        DebugStop();
    }
    
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Regular Controur integrals
    
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    //  REAL n3 = normal[2];
    
    TPZManVector<STATE,3> sol_qL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[1].sol[0];
    
    
    //  Getting Q solution for left and right side
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL qnL = (qxL*n1) + (qyL*n2);
    
    //  Getting another required data
    REAL TimeStep = fData->TimeStep();
    REAL Theta = fData->Theta();
    
    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
    int PRowsleft = phiPL.Rows();
    
    
    int FirstQL = 0;
    int FirstPL = QRowsleft + FirstQL;
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of contour integrals for Residual Vector
    //  Time step n
    
    if(fData->IsLastState())
    {
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            ef(ip + FirstPL) += (-1.0) * (1.0-Theta) * (TimeStep) * weight * phiPL(ip,0) * qnL;
        }
        
        return;
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  End of contribution of contour integrals for Residual Vector
    
    
    
    
//    if(!fData->IsLastState())
//    {
    
        //  ////////////////////////// Jacobian Matrix and Residual Vector ///////////////////////////////////
        //  Contribution of contour integrals for Jacobian matrix and Residual Vector
        //  Time step n+1
        
        
//    //  First Block (Equation One) constitutive law
//    //  Integrate[L dot(v, n), Gamme_{e}]  (Equation One) Left-Left part
//        
//    for (int iq=0; iq < QRowsleft; iq++)
//    {
//        
//        int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
//        int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;
//        
//        REAL vnL   =   (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex))*(n1) +
//                       (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex))*(n2);
//        
//        ef(iq + FirstQL) += (-1.0) * weight * vnL * PressureL;
//        
//        for (int jp=0; jp < PRowsleft; jp++)
//        {
//            ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * vnL * phiPL(jp,0);
//        }
//        
//        
//    }
    
    
    //  Second Block (Equation Two) Bulk flux  equation
    // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
    for (int ip=0; ip < PRowsleft; ip++)
    {
        
        ef(ip + FirstPL) += (-1.0) * (Theta) * (TimeStep) * weight * qnL * phiPL(ip,0);
        
        for (int jq=0; jq<QRowsleft; jq++)
        {
            int jvectorindex    = dataleft[0].fVecShapeIndex[jq].first;
            int jshapeindex     = dataleft[0].fVecShapeIndex[jq].second;
            
            REAL vnL =
            (n1) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(0,jvectorindex)) +
            (n2) * (phiQL(jshapeindex,0)*dataleft[0].fDeformedDirections(1,jvectorindex)) ;
            
            ek(ip + FirstPL,jq + FirstQL) += (-1.0) * (Theta) * (TimeStep) * weight * phiPL(ip,0) * vnL;
            
        }
    }
    

    
//    }
    //  Regular Controur integrals
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    switch (bc.Type()) {
        case 0 :    // Dirichlet BC  bc: Qn
        {
            ApplyQnD(data,dataleft,weight,ek,ef,bc);
        }
            break;
            
        case 1 :    // Neumann BC  bc: P
        {
            ApplyPN(data,dataleft,weight,ek,ef,bc);
        }
            break;
            
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
    
}

void TPZMatDarcy2dhdiv::ApplyQnD       (TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
//    TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    
    TPZManVector<STATE,3> sol_qL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[1].sol[0];
    
    
    //  Getting Q solution for left and right side
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL qnL = (qxL*n1) + (qyL*n2);
//    REAL PressureL = sol_pL[0];
    
    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
//    int PRowsleft = phiPL.Rows();
    
    int FirstQL = 0;
//    int FirstPL = QRowsleft + FirstQL;
    
    
    STATE v2[3];
    v2[0] = bc.Val2()(0,0);    //  qx
    v2[1] = bc.Val2()(1,0);    //  qy
    v2[2] = bc.Val2()(2,0); //  Pressure
    REAL qN = (v2[0]*n1 + v2[1]*n2);    // Normal Flux
    
    
//    //  First Block (Equation One) constitutive law
//    //  Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
//    for (int iq=0; iq < QRowsleft; iq++)
//    {
//        int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
//        int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;
//        
//        REAL vnL   =   (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex))*(n1) +
//        (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex))*(n2);
//        
//        ef(iq + FirstQL) += (-1.0) * weight * vnL * PressureL;
//        
//        for (int jp=0; jp < PRowsleft; jp++)
//        {
//            ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * vnL * phiPL(jp,0);
//        }
//        
//    }
    
    for(int iq=0; iq < QRowsleft; iq++)
    {
        int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
        int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;
        
        REAL vni    =   (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex)*n1)+(phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex)*n2);
        ef(iq + FirstQL) += weight * ( (gBigNumber * ( qnL - qN ) * vni ) );
        
        for (int jq=0; jq < QRowsleft; jq++)
        {
            int jLvectorindex       = dataleft[0].fVecShapeIndex[jq].first;
            int jLshapeindex        = dataleft[0].fVecShapeIndex[jq].second;
            
            REAL vnj    =   (phiQL(jLshapeindex,0)*dataleft[0].fDeformedDirections(0,jLvectorindex)*n1)+(phiQL(jLshapeindex,0)*dataleft[0].fDeformedDirections(1,jLvectorindex)*n2);
            ek(iq + FirstQL,jq + FirstQL) += weight * ( (gBigNumber * ( vnj ) * vni ) );
        }
    }
}


void TPZMatDarcy2dhdiv::ApplyPN        (TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    
    TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    
    TPZManVector<STATE,3> sol_qL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[1].sol[0];

    
    //  Getting P solution for left and right side
    REAL PressureL = sol_pL[0];
    

    
    int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
    
    int PRowsleft = phiPL.Rows();
    
    
    int FirstQL = 0;
    int FirstPL = QRowsleft + FirstQL;
    
    
    STATE v2[3];
    v2[0] = bc.Val2()(0,0);    //  qx
    v2[1] = bc.Val2()(1,0);    //  qy
    v2[2] = bc.Val2()(2,0); //  Pressure
    
    //  First Block (Equation One) constitutive law
    //  Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
    for (int iq=0; iq < QRowsleft; iq++)
    {
        int iLvectorindex       = dataleft[0].fVecShapeIndex[iq].first;
        int iLshapeindex        = dataleft[0].fVecShapeIndex[iq].second;

        REAL vnL   =   (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(0,iLvectorindex))*(n1) +
                       (phiQL(iLshapeindex,0)*dataleft[0].fDeformedDirections(1,iLvectorindex))*(n2);
        
        ef(iq + FirstQL) += (1.0) * weight * vnL * (v2[2]-PressureL);
        
        for (int jp=0; jp < PRowsleft; jp++)
        {
            ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * vnL * phiPL(jp,0);
        }
        
    }
    

    
}

void TPZMatDarcy2dhdiv::SetNotContribute(bool setNotCont){
  fNotContribute = setNotCont;
}

/** Returns the variable index associated with the name */
int TPZMatDarcy2dhdiv::VariableIndex(const std::string &name){
    if(!strcmp("MassVelocity",name.c_str()))        return  1;
    if(!strcmp("Pressure",name.c_str()))    return  2;
    if(!strcmp("MassVelocityExact",name.c_str()))        return  3;
    if(!strcmp("PressureExact",name.c_str()))    return  4;

    
    return TPZMaterial::VariableIndex(name);
}

int TPZMatDarcy2dhdiv::NSolutionVariables(int var){
    if(var == 1) return 2;
    if(var == 2) return 1;
    if(var == 3) return 2;
    if(var == 4) return 1;
    
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMatDarcy2dhdiv::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    TPZVec<STATE> SolQ, SolP;
    SolQ = datavec[0].sol[0];
    SolP = datavec[1].sol[0];
    
    if(var == 1){ //function (state variable Q)
        Solout[0] = SolQ[0];
        Solout[1] = SolQ[1];
        return;
    }
    
    if(var == 2){//function (state variable p)
        Solout[0] = SolP[0];
        return;
    }

    if(var == 3){ //function (state variable Q anal)
        TPZVec<STATE> solExact(1);
        TPZFMatrix<STATE> flux(1,1);
        fTimedependentFunctionExact->Execute(datavec[1].x, fData->Time(), solExact,flux);
        Solout[0] = flux(0,0);
        Solout[1] = 0.0;
        return;
    }
    
    if(var == 4){//function (state variable p anal)
        TPZVec<STATE> solExact(1);
        TPZFMatrix<STATE> flux(1,1);
        fTimedependentFunctionExact->Execute(datavec[1].x, fData->Time(), solExact,flux);
        Solout[0] = solExact[0];
        return;
    }
    
}


void TPZMatDarcy2dhdiv::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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

void TPZMatDarcy2dhdiv::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;
    }
}
