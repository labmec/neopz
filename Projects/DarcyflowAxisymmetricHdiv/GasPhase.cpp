//
//  GasPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "GasPhase.h"


GasPhase::GasPhase() : ReducedPVT()
{

    /** @brief Specific gravity of Gas fraction $ */
    fGas_gamma = 0.7;
    
    /** @brief Carbon dioxide content yCO2 molar fraction $ */
    fyCO2 = 0.0;
    
    /** @brief Acid sulfhidric content yH2S molar fraction $ */
    fyH2S = 0.0;
    
    /** @brief Nitrogen content yN2 molar fraction $ */
    fyN2 = 0.0;
    
    /** @brief Pressure at standard contidions [Pa] - 14.696 [Psi] $ */
    fPstd = 101325.353;
    
    /** @brief Temperature at standard contidions [K] - 60 [F] $ */
    fTstd = 288.706;
    
    /** @brief air mass density at standard contidions [kg/m3] - 0.0715112288 [lb/feet3]  $ */
    fRho_air_std = 1.1455;

}

GasPhase::~GasPhase()
{
    
}


/** @brief Density - kg/m3  $\rho$ */
void GasPhase::Density(TPZManVector<REAL> &rho, TPZManVector<REAL> state_vars)
{
    rho[0] = 0.0;
    rho[2] = 0.0;
}

/** @brief viscosity - Pa s  $\mu$ */
void GasPhase::Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars)
{
    mu[0] = 0.0;
    mu[2] = 0.0;
}

/** @brief Compressibility - 1/pa $c$ */
void GasPhase::Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars)
{
    c[0] = 0.0;
    c[2] = 0.0;
}

/** @brief Computes the pseudo critic pressure of Gas $ */
REAL GasPhase::Ppc(){
    REAL PpcHC, Ppct;
    PpcHC = 706. - 51.7 * fGas_gamma -11.1*(fGas_gamma*fGas_gamma);
    Ppct= (PpcHC - 80.0*fyCO2 + 130.0*fyH2S - 250.0*fyN2) * (6894.76);
    return Ppct; /* Pascal */
}

/** @brief Computes the pseudo critic temperature of Gas $ */
REAL GasPhase::Tpc(){
    REAL TpcHC, Tpct;
    TpcHC = 187.0 + 330.0*fGas_gamma -71.5*(fGas_gamma*fGas_gamma);
    Tpct= (TpcHC+440.0*fyCO2+600.0*fyH2S-170.0*fyN2) * (5.0/9.0);
    return Tpct; /* K degrees */
}

/** @brief Computes the compressibility factor using Beggs and Brill correlation (1974)  $ */
void GasPhase::Z(TPZManVector<REAL> &z, TPZManVector<REAL> state_vars){

    /* Getting the pressure gas phase */
    REAL Pg = state_vars[1];
    
    /* Computing Tr value */
    REAL Tpc= 0.0;
    REAL Tr =0.0;
    REAL T = this->GetTRes()/this->GetTRef();
    Tpc = this->Tpc()/GetTRef();
    Tr = T/Tpc;
    
    /* Computing Pr value */
    REAL Ppc= 0.0;
    REAL Pr =0.0;

    Ppc = this->Ppc()/(GetPRef());
    Pr = Pg/Ppc;
    
    if ((1.2 && Tr && 2.4) ||  (0.0 < Pr < 13.0) ) {
        std::cout<< " Tr or Pr put pf this conditions = 1.2 < Tr < 2.4;  0 < Pr < 13"<< std::endl;
        std::cout<< " Pr =" << Pr << std::endl;
        std::cout<< " Tr =" << Tr << std::endl;
        DebugStop();
    }
    
    /* Computing z factor and derivatives */
    
    REAL a,a1=1.39,a2=0.92,a3=0.36,a4=0.101;
    a = a1 * (sqrt(Tr - a2)) - a3 * Tr - a4;
    
    REAL b,b1=0.62,b2=0.23,b3=0.066,b4=0.86,b5=0.037,b6=0.32,b7=9.0,dbdPr;
    b = (b1 - b2*Tr)*Pr + (Pr*Pr)*((b3/(Tr - b4)) - b5) + (b6*pow(Pr,6.0))*(pow(10.0,b7*(1.0 - Tr)));
    
    dbdPr = b1 + 3.0 * pow(2.0,1.0 + b7*(1.0 - Tr)) * pow(5.0,b7*(1.0 - Tr)) * b6 * pow(5.0,Pr) - b2*Tr +
    2.0 * Pr *(-b5 + (b3/(-b4 + Tr)));
    
    REAL c,c1=0.132,c2=0.32;
    c = c1 - c2 * log10(Tr);
    
    REAL d,d1=0.3106,d2=0.49,d3=0.1824;
    d = pow(10.0,(d1 - d2 * Tr  + d3 * Tr*Tr));
    
    REAL zv = 0.0;
    REAL dzdPr = 0.0;
    
    zv = a + (1.0 - a) / (exp(b)) + c * pow(Pr,d);
    dzdPr = c*d*pow(Pr,d-1.0) - (1.0 - a)*exp(-b)*dbdPr;
    
    z[0] = zv;
    z[2] = dzdPr;
    
}
