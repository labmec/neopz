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
    fLastElastFunction = new TPZLastElastFunction();
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
    this->fLastElastFunction = new TPZLastElastFunction();
    this->fVisc = visc;
    this->fCl = Cl;
    this->fPe = Pe;
    this->fgradPref = gradPref;
    this->fvsp = vsp;
}

TPZPlaneFractCouplingMat::~TPZPlaneFractCouplingMat()
{
    this->fLastElastFunction = NULL;
    this->fVisc = 0.;
    this->fCl = 0.;
    this->fPe = 0.;
    this->fgradPref = 0.;
    this->fvsp = 0.;
}

void TPZPlaneFractCouplingMat::SetLastElastCMesh(TPZCompMesh * LastElastCMesh)
{
    fLastElastFunction->SetLastElastCMesh(LastElastCMesh);
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
                                                  TPZFMatrix<STATE> &ek,
                                                  TPZFMatrix<STATE> &ef)
{
    if(!datavec[1].phi) return;
    
    TPZFMatrix<REAL>  & phi_p = datavec[1].phi;
    TPZFMatrix<REAL>  & dphi_p = datavec[1].dphix;
    TPZManVector<STATE,3> sol_p = datavec[1].sol[0];
    TPZFMatrix<STATE> & dsol_p = datavec[1].dsol[0];
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
    
    int phipRows = phi_p.Rows();
    int phiuCols = phi_u.Cols();
    
    REAL deltaT = globTimeControl.actDeltaT();
    
	if(gState == EActualState) //current state (n+1): Matrix stiffnes
    {
        TPZManVector<STATE,3> sol_u = datavec[0].sol[0];
        REAL uy = sol_u[1];
        REAL w = 2. * uy;

        REAL actQl = 0.;
        REAL actdQldp = 0.;
        if(w > 0.)
        {
            REAL pfrac = sol_p[0];
            actQl = globLeakoffStorage.QlFVl(datavec[1].gelElId, pfrac, deltaT, fCl, fPe, fgradPref, fvsp);
            actdQldp = globLeakoffStorage.dQlFVl(datavec[1].gelElId, pfrac, deltaT, fCl, fPe, fgradPref, fvsp);
        }
        
        REAL visc = this->fVisc;
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
            REAL uy;
            this->fLastElastFunction->Execute(datavec[1].x,uy);
            REAL w = 2.*uy;
            ef(phiuCols+in,0) += (+1.) * weight * w/deltaT * phi_p(in,0);
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
        case 5:// Dirichlet condition only on the pressure equation pressure
        {
            ApplyDirichlet_P(datavec, weight, ek, ef, bc);
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
    std::cout << "\n\n\nEntrei no TPZPlaneFractCouplingMat::Solution\n\n\n";
    DebugStop();
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
                                                TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef,
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

void TPZPlaneFractCouplingMat::ApplyNeumann_U(TPZVec<TPZMaterialData> &datavec,
                                              STATE weight,
                                              TPZFMatrix<STATE> &ek,
                                              TPZFMatrix<STATE> &ef,
                                              TPZBndCond &bc)
{
    if(gState == EPastState)
    {
        return;
    }
    else
    {
        TPZManVector<STATE,3> sol_p = datavec[1].sol[0];
        
        bc.Val2().Zero();
        bc.Val2()(1,0) = sol_p[0];
        
        /**
         * D(u,alphaElast) calculo da matriz tangente.
         */
        //Elastica
        TPZElast3Dnlinear::ContributeBC(datavec[0], weight, ek, ef, bc);
        
        /**
         * D(u,alphaPressure) calculo da matriz tangente.
         */
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
                                                 TPZFMatrix<STATE> &ek,
                                                 TPZFMatrix<STATE> &ef,
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

void TPZPlaneFractCouplingMat::ApplyDirichlet_P(TPZVec<TPZMaterialData> &datavec,
                                                STATE weight,
                                                TPZFMatrix<STATE> &ek,
                                                TPZFMatrix<STATE> &ef,
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
    
    TPZManVector<STATE> sol_pMat = datavec[1].sol[0];
    REAL solpress = sol_pMat[0];
    REAL imposedPress = bc.Val2()(0,0);
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    for(int in = 0; in < phrp; in++)
    {
        ef(in+c_u,0) += (+1.) * weight * BIGNUMBER * imposedPress * phi_p(in,0);
        ef(in+c_u,0) += (-1.) * weight * BIGNUMBER * solpress * phi_p(in,0);
        
        for(int jn = 0; jn < phrp; jn++)
        {
            ek(in+c_u,jn+c_u) += weight * BIGNUMBER * (phi_p(in,0) * phi_p(jn,0));
        }
    }
}

void TPZPlaneFractCouplingMat::ApplyNeumann_P(TPZVec<TPZMaterialData> &datavec,
                                              STATE weight,
                                              TPZFMatrix<STATE> &ek,
                                              TPZFMatrix<STATE> &ef,
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
    
//    {//Newmann original (funcao constante)
//        REAL Qinj1wing_Hbullet = bc.Val2()(0,0);
//        for(int in = 0; in < phrp; in++)
//        {
//            ef(in+c_u,0) += (-1.) * weight * Qinj1wing_Hbullet * phi_p(in,0);
//        }
//    }

    {//Newmann suave por funcao parabola
        REAL Qinj1wing_hBullet = bc.Val2()(0,0);
        REAL z = datavec[1].x[2];
        REAL zf = globLayerStruct.DownBulletDepth();
        REAL hBullet = globLayerStruct.HBullet();
        REAL val = (6.*Qinj1wing_hBullet*(z - zf)*(hBullet - z + zf))/(hBullet*hBullet);

        for(int in = 0; in < phrp; in++)
        {
            ef(in+c_u,0) += (-1.) * weight * val * phi_p(in,0);
        }
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


//=====================================================================

/*
TPZPlaneFractBulletMat::TPZPlaneFractBulletMat()
{
    this->fDiameter = 0.;
    this->fvisc = 0.;
    this->fQinj_hbullet = 0.;
}

TPZPlaneFractBulletMat::TPZPlaneFractBulletMat(int nummat,
                                               REAL Diameter,
                                               REAL visc,
                                               REAL Qinj_hbullet) : TPZMaterial(nummat)
{
    this->fDiameter = Diameter;
    this->fvisc = visc;
    this->fQinj_hbullet = Qinj_hbullet;
}

TPZPlaneFractBulletMat::~TPZPlaneFractBulletMat()
{
    
}

void TPZPlaneFractBulletMat::Contribute(TPZVec<TPZMaterialData> &datavec,
                                        STATE weight,
                                        TPZFMatrix<STATE> &ek,
                                        TPZFMatrix<STATE> &ef)
{
    if(TPZPlaneFractCouplingMat::IsPastState())
    {
        return;
    }
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	int c_u = phi_u.Cols();
    
    TPZFMatrix<REAL> & phi_p = datavec[1].phi;
    int  phrp = phi_p.Rows();
    
    TPZFMatrix<REAL> & dphi_p = datavec[1].dphix;
    
    TPZFMatrix<REAL> & dsol_p = datavec[1].dsol[0];
    REAL dpdx = dsol_p(0,0);
    
    const REAL D = this->fDiameter;
    const REAL piD4_128mi = M_PI * D*D*D*D / (128. * this->fvisc);
    
    for(int in = 0; in < phrp; in++)
    {
        ef(in+c_u,0) += (-1.) * weight * ( dphi_p(in,0) * piD4_128mi * dpdx  +  this->fQinj_hbullet * phi_p(in,0) );
        
        for(int jn = 0; jn < phrp; jn++)
        {
            ek(in+c_u,jn+c_u) += (+1.) * weight * ( piD4_128mi + dphi_p(in,0) * dphi_p(jn,0) );
        }
    }
}

void TPZPlaneFractBulletMat::ContributeInterface(TPZVec<TPZMaterialData> &datavec,
                                                 TPZVec<TPZMaterialData> &dataleftvec,
                                                 TPZVec<TPZMaterialData> &datarightvec,
                                                 REAL weight,
                                                 TPZFMatrix<STATE> &ek,
                                                 TPZFMatrix<STATE> &ef)
{
    if(TPZPlaneFractCouplingMat::IsPastState())
    {
        return;
    }
    
    TPZFMatrix<REAL> & phi_u = datavec[0].phi;
	int nShapeU = phi_u.Cols();
    
    TPZFMatrix<REAL> & phi_Left = dataleftvec[1].phi;
    TPZManVector<REAL> p_LeftVec = dataleftvec[1].sol[0];
    REAL p_Left = p_LeftVec[0];
    
    TPZFMatrix<REAL> & phi_Right = datarightvec[0].phi;
    TPZManVector<REAL> p_RightVec = datarightvec[0].sol[0];
    REAL p_Right = p_RightVec[0];
    
    const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
    
    int nShapeLeft = phi_Left.Rows();
    int nShapeRight = phi_Right.Rows();
    
    for(int iL = 0; iL < nShapeLeft; iL++)
    {
        ef(iL,0) += (-1.) * weight * BIGNUMBER * ( (p_Left - p_Right) * phi_Left(iL,0) );
        
        for(int jL = 0; jL < nShapeLeft; jL++)
        {
            ek(iL,jL) += (+1.) * weight * BIGNUMBER * ( phi_Left(jL,0) * phi_Left(iL,0) );
        }
        for(int jR = 0; jR < nShapeRight; jR++)
        {
            ek(iL,nShapeLeft+jR) += (-1.) * weight * BIGNUMBER * ( phi_Left(iL,0) * phi_Right(jR,0) );
        }
    }
    for(int iR = 0; iR < nShapeRight; iR++)
    {
        ef(nShapeLeft + iR,0) += (-1.) * weight * BIGNUMBER * ( (p_Right - p_Left) * phi_Right(iR) );
        
        for(int jL = 0; jL < nShapeLeft; jL++)
        {
            ek(nShapeLeft+iR,jL) += (-1.) * weight * BIGNUMBER * ( phi_Left(jL,0) * phi_Right(iR,0) );
        }
        for(int jR = 0; jR < nShapeRight; jR++)
        {
            ek(nShapeLeft+iR,nShapeLeft+jR) += (+1.) * weight * BIGNUMBER * ( phi_Right(iR,0) * phi_Right(jR,0) );
        }
    }
}
    

void TPZPlaneFractBulletMat::FillDataRequirementsInterface(TPZVec<TPZMaterialData> &datavec)
{
    int nref = datavec.size();
	for(int i = 0; i < nref; i++)
	{
		datavec[i].SetAllRequirements(false);
        datavec[i].fNeedsSol = true;
	}
}
*/
