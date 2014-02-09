//
//  TPZPlaneFractCouplingMat.cpp
//  PZ
//
//  Created by Cesar Lucci on 11/01/13.
//
//

#include "TPZPlaneFractCouplingMat.h"
#include "TPZPlaneFractureData.h"
#include "pzmat2dlin.h"
#include "pzbndcond.h"

TPZPlaneFractCouplingMat::EState TPZPlaneFractCouplingMat::gState = EActualState;


TPZPlaneFractCouplingMat::TPZPlaneFractCouplingMat() : TPZElast3Dnlinear()
{
    fVisc = 0.;
    fCl = 0.;
    fPe = 0.;
    fgradPref = 0.;
    fvsp = 0.;

}

TPZPlaneFractCouplingMat::TPZPlaneFractCouplingMat(int nummat, STATE E, STATE poisson, TPZVec<STATE> &force,
                                                   STATE preStressXX, STATE preStressYY, STATE preStressZZ,
                                                   STATE visc,
                                                   STATE Cl,
                                                   STATE Pe,
                                                   STATE gradPref,
                                                   STATE vsp) :
                                                   TPZElast3Dnlinear(nummat, E, poisson, force, preStressXX, preStressYY, preStressZZ)
{
    fVisc = visc;
    fCl = Cl;
    fPe = Pe;
    fgradPref = gradPref;
    fvsp = vsp;
}

TPZPlaneFractCouplingMat::~TPZPlaneFractCouplingMat()
{
    fVisc = 0.;
    fCl = 0.;
    fPe = 0.;
    fgradPref = 0.;
    fvsp = 0.;
}

int TPZPlaneFractCouplingMat::NSolutionVariables(int var)
{
    if(var == 19)//Pressure
    {
        return 1;
    }
    else//Elasticity
    {
        return TPZElast3Dnlinear::NSolutionVariables(var);
    }
}

int TPZPlaneFractCouplingMat::VariableIndex(const std::string &name)
{
    if(!strcmp("Pressure",name.c_str())) return 19;
    else//Elasticity
    {
        return TPZElast3Dnlinear::VariableIndex(name);
    }
}

void TPZPlaneFractCouplingMat::Contribute(TPZVec<TPZMaterialData> &datavec,
                                          STATE weight,
                                          TPZFMatrix<STATE> &ek,
                                          TPZFMatrix<STATE> &ef)
{
    if(gState == EPastState)
    {
        return;
    }
    else
    {
        TPZElast3Dnlinear::Contribute(datavec[0], weight, ek, ef);
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
    
#ifdef DEBUG
    if(dsol_p.Rows() != 2 || dphi_p.Rows() != 2)
    {
        DebugStop();
    }
#endif
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
    
    REAL w = 0.;
    {
        int nsolu = datavec[0].sol.NElements();
        for(int s = 0; s < nsolu; s++)
        {
            TPZManVector<REAL,3> sol_u = datavec[0].sol[s];
            REAL uy = sol_u[1];
            w += 2.*uy;
        }
    }
    
    int phipRows = phi_p.Rows();
    int phiuCols = phi_u.Cols();
    
    REAL visc = this->fVisc;
    REAL deltaT = globTimeControl.actDeltaT();
    
	if(gState == EActualState) //current state (n+1): Matrix stiffnes
    {
        REAL actQl = globLeakoffStorage.QlFVl(datavec[1].gelElId, sol_p[0], deltaT, fCl, fPe, fgradPref, fvsp);
        REAL actdQldp = globLeakoffStorage.dQlFVl(datavec[1].gelElId, sol_p[0], deltaT, fCl, fPe, fgradPref, fvsp);
        
        for(int in = 0; in < phipRows; in++)
        {
            //----Residuo----
            //termo (wˆ3/(12*mi)) * gradP * gradVp
            ef(phiuCols+in,0) += (-1.) * weight * (w*w*w/(12.*visc)) * dsol_p(0,0) * dphi_p(0,in);
            ef(phiuCols+in,0) += (-1.) * weight * (w*w*w/(12.*visc)) * dsol_p(1,0) * dphi_p(1,in);
            
            //termo w/deltaT * Vp
            ef(phiuCols+in,0) += (-1.) * weight * w/deltaT * phi_p(in,0);
            
            //termo 2Ql * Vp
            ef(phiuCols+in,0) += (-1.) * weight * (2.*actQl) * phi_p(in,0);
            
            
            //------Matriz tangente-----
            for(int jn = 0; jn < phiuCols; jn++)
            {
                //termo D[ (wˆ3/(12*mi)) * gradP * gradVp , w ]
                ek(phiuCols+in, jn) += (+1.) * weight * ( 3.*w*w/(12.*visc) * (2.*phi_u(1,jn)) ) * dsol_p(0,0) * dphi_p(0,in);
                ek(phiuCols+in, jn) += (+1.) * weight * ( 3.*w*w/(12.*visc) * (2.*phi_u(1,jn)) ) * dsol_p(1,0) * dphi_p(1,in);
                
                //termo D[ w/deltaT * Vp , w ]
                ek(phiuCols+in, jn) += (+1.) * weight * ( 2./deltaT * phi_u(1,jn) ) * phi_p(in,0);
            }
            for(int jn = 0; jn < phipRows; jn++)
            {
                //termo D[ (wˆ3/(12*mi))*gradP * gradVp , p ]
                ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (w*w*w/(12.*visc)) * dphi_p(0,in) * dphi_p(0,jn);
                ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (w*w*w/(12.*visc)) * dphi_p(1,in) * dphi_p(1,jn);
                
                //termo D[ 2Ql * Vp , p]
                ek(phiuCols+in, phiuCols+jn) += (+1.) * weight * (2.*actdQldp) * phi_p(in,0) * phi_p(jn,0);
            }
        }
    }
	else if(gState == EPastState)//Past state (n): Matrix mass
    {
        for(int in = 0; in < phipRows; in++)
        {
            //termo w/deltaT * Vp
            ef(phiuCols+in,0) += (+1.) * weight * w/deltaT * phi_p(in,0);
        }
    }
#ifdef DEBUG
    else
    {
        DebugStop();//mas hein???
    }
#endif
}

void TPZPlaneFractCouplingMat::ContributeBC(TPZVec<TPZMaterialData> &datavec,
                                            STATE weight,
                                            TPZFMatrix<STATE> &ek,
                                            TPZFMatrix<STATE> &ef,
                                            TPZBndCond &bc)
{
    TPZMaterialData::MShapeFunctionType shapetype = datavec[0].fShapeType;
    if(shapetype != datavec[0].EVecShape && datavec[0].phi.Cols() != 0)
    {
        std::cout << " The space to elasticity problem must be reduced space.\n";
		DebugStop();
    }
	
	switch (bc.Type())
    {
        case 0:// Dirichlet condition only on the elasticity equation
        {
            // Calculate the matrix contribution for pressure
            ApplyDirichlet_U(datavec, weight, ek, ef, bc);
            break;
        }
        case 1:// Neumann condition in  both elasticity and pressure equations
        {
            ApplyNeumann_U(datavec, weight, ek, ef, bc);
            ContributePressure(datavec, weight, ek, ef);
            break;
		}
        case 2:
        {   //Estao querendo utilizar condicao mista!!!
            DebugStop();
            break;
        }
        case 3:// Block directional condition only on the elasticity equation
        {
            ApplyBlockedDir_U(datavec, weight, ek, ef, bc);
            break;
        }
        case 4:// Neumann condition only on the pressure equation pressure
        {
            ApplyNeumann_P(datavec, weight, ek, ef, bc);
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
#ifdef DEBUG
    else
    {
        DebugStop();
    }
#endif
}

void TPZPlaneFractCouplingMat::ApplyDirichlet_U(TPZVec<TPZMaterialData> &datavec,
                                                STATE weight,
                                                TPZFMatrix<> &ek,
                                                TPZFMatrix<> &ef,
                                                TPZBndCond &bc)
{
    if(gState == EPastState)
    {
        return;
    }
    else
    {
        DebugStop();//Validar!!!
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
    }
}

void TPZPlaneFractCouplingMat::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec,
                                              STATE weight,
                                              TPZFMatrix<> &ek,
                                              TPZFMatrix<> &ef,
                                              TPZBndCond &bc)
{
    if(gState == EPastState)
    {
        return;
    }
    else
    {
        TPZManVector<REAL,3> sol_p = datavec[1].sol[0];
        
        bc.Val2().Zero();
        bc.Val2()(1,0) = sol_p[0];
        
        //Elastica
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
        
        //Pressao
        TPZFMatrix<REAL> &phi_u = datavec[0].phi;
        TPZFMatrix<REAL> &phi_p = datavec[1].phi;
        for (int in = 0; in < phi_u.Cols(); in++)
        {
            for(int jp = 0; jp < phi_p.Rows(); jp++)
            {
                ek(in,phi_u.Cols()+jp) -= weight * ( phi_u(1,in)*phi_p(jp,0) );
            }
        }
    }
}

void TPZPlaneFractCouplingMat::ApplyBlockedDir_U(TPZVec<TPZMaterialData> &datavec,
                                                 STATE weight,
                                                 TPZFMatrix<> &ek,
                                                 TPZFMatrix<> &ef,
                                                 TPZBndCond &bc)
{
    if(gState == EPastState)
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
    if(gState == EPastState)
    {
        return;
    }
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> & phi_p = datavec[1].phi;
    int  phrp = phi_p.Rows();
    
    REAL Qinj = bc.Val2()(0,0);
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
