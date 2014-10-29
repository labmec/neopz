//
//  pzelasticSest2D.cpp
//  PZ
//
//  Created by Diogo Cecilio on 9/23/14.
//
//

#include "pzelasticSest2D.h"

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D():TPZMatElasticity2D()
{
}
/**
 * @brief Creates an elastic material with:
 * @param id material id
 * @param E elasticity modulus
 * @param nu poisson coefficient
 * @param fx forcing function \f$ -x = fx \f$
 * @param fy forcing function \f$ -y = fy \f$
 * @param plainstress \f$ plainstress = 1 \f$ indicates use of plainstress
 */
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id, REAL E, REAL nu, REAL fx, REAL fy, int plainstress):TPZMatElasticity2D(id,E,nu,fx,fy,plainstress)
{
 
    
}

TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(int id):TPZMatElasticity2D(id)
{
}

/** @brief Copies the data of one TPZElasticityMaterial object to another */
TPZElasticityMaterialSest2D::TPZElasticityMaterialSest2D(const TPZMatElasticity2D &copy):TPZMatElasticity2D(copy)
{
}

/// Compute post processed values
void TPZElasticityMaterialSest2D::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout)
{
    enum EVar {sigx=5,sigy=6,sigz=12};
    int numbersol = datavec[0].dsol.size();
    int ipos = 0;
    if (fPostProcIndex < numbersol) {
        ipos = fPostProcIndex;
    }
    
    TPZVec<STATE> &Sol = datavec[0].sol[ipos];
    TPZFMatrix<STATE> &DSol = datavec[0].dsol[ipos];
    TPZFMatrix<REAL> &axes = datavec[0].axes;
    TPZFNMatrix<4,STATE> DSolxy(2,2);
    
    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL epsz = 0.;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL TauXY;
    
    // dudx - dudy
    DSolxy(0,0) = DSol(0,0)*axes(0,0)+DSol(1,0)*axes(1,0);
    DSolxy(1,0) = DSol(0,0)*axes(0,1)+DSol(1,0)*axes(1,1);
    // dvdx - dvdy
    DSolxy(0,1) = DSol(0,1)*axes(0,0)+DSol(1,1)*axes(1,0);
    DSolxy(1,1) = DSol(0,1)*axes(0,1)+DSol(1,1)*axes(1,1);
    
    epsx = DSolxy(0,0);// du/dx
    epsy = DSolxy(1,1);// dv/dy
    epsxy = 0.5*(DSolxy(1,0)+DSolxy(0,1));
    epsz = datavec[1].sol[ipos][0];
    
    REAL lambda = GetLambda();
    REAL mu = GetMU();
    if (this->fPlaneStress) {
        DebugStop();
    }
    TauXY = 2*mu*epsxy+fPreStressXY;
#ifdef DEBUG
    REAL TauXY2 = fE*epsxy/(1.+fnu)+fPreStressXY;
#ifdef REALfloat
    if (fabs(TauXY-TauXY2) > 1.e-10) {
        DebugStop();
    }
#else
    if (fabs(TauXY-TauXY2) > 1.e-6) {
        DebugStop();
    }
#endif
#endif
    if (this->fPlaneStress){
        DebugStop();
    }
    else
    {
        SigX = fE/((1.-2.*fnu)*(1.+fnu))*((1.-fnu)*epsx+fnu*(epsy+epsz))+fPreStressXX;
        SigY = fE/((1.-2.*fnu)*(1.+fnu))*(fnu*epsx+(1.-fnu)*(epsy+epsz))+fPreStressYY;
        SigZ = fPreStressZZ+lambda*(epsx+epsy+epsz)+2.*mu*epsz;
    }

    switch (var) {
        case sigx:
            Solout[0] = SigX;
            break;
        case sigy:
            Solout[0] = SigY;
            break;
        case sigz:
            Solout[0] = SigZ;
            break;
        default:
            DebugStop();
            break;
    }
    
}




