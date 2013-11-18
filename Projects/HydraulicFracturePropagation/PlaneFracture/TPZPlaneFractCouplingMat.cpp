//
//  pznlfluidstructure2d.cpp
//  PZ
//
//  Created by Cesar Lucci on 05 nov 2013.
//
//

#include "TPZPlaneFractCouplingMat.h"
#include "TPZPlaneFractureData.h"
#include "pzbndcond.h"

TPZPlaneFractCouplingMat::EState TPZPlaneFractCouplingMat::gState = ECurrentState;


TPZPlaneFractCouplingMat::TPZPlaneFractCouplingMat() : TPZElast3Dnlinear()
{
    
}

TPZPlaneFractCouplingMat::TPZPlaneFractCouplingMat(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                                                   STATE preStressXX, STATE preStressYY, STATE preStressZZ, STATE visc) :
                                                   TPZElast3Dnlinear(nummat, E, poisson, force, preStressXX, preStressYY, preStressZZ)
{
    fVisc = visc;
}

TPZPlaneFractCouplingMat::~TPZPlaneFractCouplingMat()
{

}

void TPZPlaneFractCouplingMat::Contribute(TPZVec<TPZMaterialData> &datavec,
                                          STATE weight,
                                          TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef)
{
    if(gState == ELastState)
    {
        return;
    }
    else
    {
        TPZElast3Dnlinear::Contribute(datavec[0], weight, ek, ef);
        ContributePressure(datavec, weight, ek, ef);
    }
}

void TPZPlaneFractCouplingMat::ContributePressure(TPZVec<TPZMaterialData> &datavec,
                                                  REAL weight,
                                                  TPZFMatrix<REAL> &ek,
                                                  TPZFMatrix<REAL> &ef)
{
    if(!datavec[1].phi) return;
    
    TPZFMatrix<REAL>  & phi_p = datavec[1].phi;
    TPZFMatrix<REAL>  & dphi_p = datavec[1].dphix;
    TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
    TPZFMatrix<REAL> & dsol_p = datavec[1].dsol[0];
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
    
    int nsolu = datavec[0].sol.NElements();
    REAL w = 0.;
    
    for(int s = 0; s < nsolu; s++)
    {
        TPZManVector<REAL,3> sol_u = datavec[0].sol[s];
        REAL uy = sol_u[1];
        w += 2.*uy;
    }
    
    int phipCols = phi_p.Rows();
    int phiuCols = phi_u.Cols();
    
    REAL visc = this->fVisc;
    REAL deltaT = globFractInput3DData.actDeltaT();
    
	if(gState == ECurrentState) //current state (n+1): Matrix stiffnes
    {
        REAL actQl = globFractInput3DData.QlFVl(datavec[1].gelElId, sol_p[0]);
        REAL actdQldp = globFractInput3DData.dQlFVl(datavec[1].gelElId, sol_p[0]);
        
        for(int in = 0; in < phipCols; in++)
        {
            //----Residuo----
            //termo (wˆ3/(12*mi))*gradP * gradVp
            ef(phiuCols+in,0) += (-1.) * weight * (w*w*w/(12.*visc)) * dsol_p(0,0) * dphi_p(0,in);
            
            //termo w/deltaT * Vp
            ef(phiuCols+in,0) += (-1.) * weight * w/deltaT * phi_p(in,0);
            
            //termo 2Ql * Vp
            ef(phiuCols+in,0) += (-1.) * weight * (2.*actQl) * phi_p(in,0);
            
            
            //------Matriz tangente-----
            for(int jn = 0; jn < phiuCols; jn++)
            {
                //termo D[ (wˆ3/(12*mi))*gradP * gradVp , w ]
                ek(phiuCols+in, jn) += (+1.) * weight * ( 3.*w*w/(12.*visc) * (2.*phi_u(1,jn)) ) * dsol_p(0,0) * dphi_p(0,in);
                
                //termo D[ w/deltaT * Vp , w ]
                ek(phiuCols+in, jn) += (+1.) * weight * ( 2./deltaT * phi_u(1,jn) ) * phi_p(in,0);
            }
            for(int jn = 0; jn < phipCols; jn++)
            {
                //termo D[ (wˆ3/(12*mi))*gradP * gradVp , p ]
                ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (w*w*w/(12.*visc)) * dphi_p(0,in) * dphi_p(0,jn);
                
                //termo D[ 2Ql * Vp , p]
                ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (2.*actdQldp) * phi_p(in,0) * phi_p(jn,0);
            }
        }
    }
    
    //Last state (n): Matrix mass
	if(gState == ELastState)
    {
        for(int phip = 0; phip < phipCols; phip++)
        {
            //termo w/deltaT * Vp
            ef(phiuCols+phip,0) += (+1.) * weight * w/deltaT * phi_p(phip,0);
        }
    }
}

void TPZPlaneFractCouplingMat::ContributeBC(TPZVec<TPZMaterialData> &datavec,
                                            STATE weight,
                                            TPZFMatrix<STATE> &ek,
                                            TPZFMatrix<STATE> &ef,
                                            TPZBndCond &bc)
{
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype!=datavec[0].EVecShape && datavec[0].phi.Cols()!= 0)
    {
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
	
	switch (bc.Type())
    {
        case 0:// Dirichlet condition only on the elasticity equation
        {
            // Calculate the matrix contribution for pressure
            ApplyDirichlet_U(datavec, weight, ek,ef,bc);
            break;
        }
        case 1:// Neumann condition in  both elasticity and pressure equations
        {
            ApplyNeumann_U(datavec, weight,ek,ef,bc);
            ContributePressure(datavec, weight, ek, ef);
            break;
		}
        case 2:// Mixed condition only on the elasticity equation
        {
            ApplyMixed_U(datavec, weight, ek,ef,bc);
            break;
        }
        case 3:// Neumann condition only on the pressure equation pressure
        {
            ApplyNeumann_P(datavec, weight, ek,ef,bc);
            break;
        }
        default:
        {
            DebugStop();
        }
    }
}

void TPZPlaneFractCouplingMat::Solution(TPZVec<TPZMaterialData> &datavec,
                                        int var,
                                        TPZVec<STATE> &Solout)
{
    if(var < 19)
    {
        TPZElast3Dnlinear::Solution(datavec[0], var, Solout);
    }
    else
    {
        DebugStop();
    }
}

void TPZPlaneFractCouplingMat::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec,
                                                STATE weight,
                                                TPZFMatrix<> &ek,
                                                TPZFMatrix<> &ef,
                                                TPZBndCond &bc)
{
    if(gState == ELastState)
    {
        return;
    }
    else
    {
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
    }
}

void TPZPlaneFractCouplingMat::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec,
                                              STATE weight,
                                              TPZFMatrix<> &ek,
                                              TPZFMatrix<> &ef,
                                              TPZBndCond &bc)
{
    if(gState == ELastState)
    {
        return;
    }
    else
    {
        TPZFMatrix<REAL> &phi_u = datavec[0].phi;
        TPZFMatrix<REAL> &phi_p = datavec[1].phi;
        TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
        
        bc.Val2().Zero();
        bc.Val2()(1,0) = sol_p[0];
        
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
        
        for (int in = 0; in < phi_u.Cols(); in++)
        {
            for(int jp = 0; jp < phi_p.Rows(); jp++)
            {
                ek(in,phi_u.Cols()+jp) += weight * phi_u(0,in)*phi_p(jp,0)*datavec[0].normal[0] +
                                          weight * phi_u(1,in)*phi_p(jp,0)*datavec[0].normal[1] +
                                          weight * phi_u(2,in)*phi_p(jp,0)*datavec[0].normal[2];
            }
        }
    }
}

void TPZPlaneFractCouplingMat::ApplyMixed_U(TPZVec<TPZMaterialData> &datavec,
                                            STATE weight,
                                            TPZFMatrix<> &ek,
                                            TPZFMatrix<> &ef,
                                            TPZBndCond &bc)
{
    if(gState == ELastState)
    {
        return;
    }
    else
    {
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
    }
}

void TPZPlaneFractCouplingMat::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec,
                                              STATE weight,
                                              TPZFMatrix<> &ek,
                                              TPZFMatrix<> &ef,
                                              TPZBndCond &bc)
{
    if(gState == ELastState)
    {
        return;
    }
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> & phi_p = datavec[1].phi;
    int  phrp = phi_p.Rows();
    
    REAL Qinj = bc.Val2()(2,0);
    for(int in = 0; in < phrp; in++)
    {
        ef(in+c_u,0) += (-1.) * weight * Qinj * phi_p(in,0);
    }
}

void TPZPlaneFractCouplingMat::FillDataRequirements(TPZVec<TPZMaterialData> &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i < nref; i++)
	{
		datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = false;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
}

void TPZPlaneFractCouplingMat::FillBoundaryConditionDataRequirement(int type, TPZVec<TPZMaterialData> &datavec)
{
    int nref = datavec.size();
	for(int i = 0; i < nref; i++)
	{
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
	}
}
