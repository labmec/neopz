//
//  TPZReynoldsFlow.cpp
//  PZ
//
//  Created by Cesar Lucci on 06/02/13.
//
//

#include "TPZReynoldsFlow.h"
#include "pzbndcond.h"

TPZReynoldsFlow::TPZReynoldsFlow() :
TPZRegisterClassId(&TPZReynoldsFlow::ClassId),
TPZMaterial()
{
    f_visc = 0.;
    f_deltaT = 0.;
    f_staticPotential = 0.;
    f_nplus1Computation = false;
}

TPZReynoldsFlow::TPZReynoldsFlow(int matId, REAL visc, REAL deltaT, REAL staticPotential) :
TPZRegisterClassId(&TPZReynoldsFlow::ClassId),
TPZMaterial(matId)
{
    f_visc = visc;
    f_deltaT = deltaT;
    f_staticPotential = staticPotential;
    f_nplus1Computation = false;
}

TPZReynoldsFlow::TPZReynoldsFlow(const TPZReynoldsFlow &cp) :
TPZRegisterClassId(&TPZReynoldsFlow::ClassId),
TPZMaterial(cp)
{
    f_visc = cp.f_visc;
    f_deltaT = cp.f_deltaT;
    f_staticPotential = cp.f_staticPotential;
    f_nplus1Computation = cp.f_nplus1Computation;
}


TPZReynoldsFlow::~TPZReynoldsFlow()
{

}


int TPZReynoldsFlow::Dimension() const
{
    return 2;
}


int TPZReynoldsFlow::NStateVariables() const
{
    return 1;
}

void TPZReynoldsFlow::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    if(f_nplus1Computation == false)//estamos no passo n
    {
        REAL simmetryy = 2.;
        STATE wn = simmetryy;
        wn *= data.sol[1][1];//data.sol = {{p},{ux,uy,uz}} e estamos interessados em uy
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            ef(i,0) += (STATE)(weight * (1./f_deltaT) * data.phi(i,0)) * wn;
        }
    }
    else//estamos no passo n+1
    {
        REAL simmetryy = 2.;
        STATE wnplus1 = (STATE)simmetryy * data.sol[1][1];//data.sol = {{p},{ux,uy,uz}} e estamos interessados em uy
        STATE w3 = wnplus1*wnplus1*wnplus1;
        REAL carterGAMMA = 1.;//TODO : Outra cmesh ou estrutura de dados? Lembre-se que valria com o tempo e no espaco!!!
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            ef(i,0) += (STATE)(weight*data.phi(i,0)) * ((STATE)(carterGAMMA * f_staticPotential) - wnplus1*(STATE)(1./f_deltaT));
        }
        for(int i = 0; i < data.phi.Rows(); i++)
        {
            for(int j = 0; j < data.phi.Rows(); j++)
            {
                ek(i,j) += w3 * (STATE)(weight * (1./(12.*f_visc)) * (data.dphix(0,j)*data.dphix(0,i) + data.dphix(1,j)*data.dphix(1,i)));
                ek(i,j) += (STATE)(weight*carterGAMMA*data.phi(j,0)*data.phi(i,0));
            }
        }
    }
}


void TPZReynoldsFlow::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    if(f_nplus1Computation == false)
    {
        return;
    }
    
    if(bc.Type() != 1)
    {
        DebugStop();
    }
    
    for(int i = 0; i < data.phi.Rows(); i++)
    {
        ef(i,0) += (STATE)(weight* data.phi(i,0)) * bc.Val2()(0,0);
    }
}


TPZMaterial * TPZReynoldsFlow::NewMaterial()
{
    return new TPZReynoldsFlow(*this);
}

void TPZReynoldsFlow::Write(TPZStream &buf, int withclassid) const
{
    DebugStop();
}

void TPZReynoldsFlow::Read(TPZStream &buf, void *context)
{
    DebugStop();
}

int TPZReynoldsFlow::ClassId() const{
    return Hash("TPZReynoldsFlow") ^ TPZMaterial::ClassId() << 1;
}
