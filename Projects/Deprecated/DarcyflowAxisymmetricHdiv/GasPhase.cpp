//
//  GasPhase.cpp
//  PZ
//
//  Created by Omar on 9/30/15.
//
//

#include "GasPhase.h"


GasPhase::GasPhase() : Phase()
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
    TPZManVector<REAL> bg(5,0.0);
    REAL rho_g, drho_gdPg;
    this->Bg(bg, state_vars);
    REAL rho_gas_std;
    rho_gas_std = fGas_gamma * fRho_air_std;

    rho_g = (rho_gas_std / bg[0]);
    drho_gdPg = - (rho_gas_std * bg[2] / (bg[0]*bg[0]));
    
    rho[0] = rho_g / GetRho();
    rho[2] = drho_gdPg / GetRho();
}

/** @brief viscosity - Pa s  $\mu$ */
void GasPhase::Viscosity(TPZManVector<REAL> &mu, TPZManVector<REAL> state_vars)
{
    /* Getting the pressure gas phase */
    REAL Pg = state_vars[1] + fPstd / GetPRef();
    
    REAL a0,a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13,a14,a15;
    
    REAL mu_std;
    mu_std = this->Viscosity_std();
    
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
    
    a0 = -2.4621182;
    a1 = 2.97054714;
    a2 = -0.286264054;
    a3 = 0.00805420522;
    a4 = 2.80860949;
    a5 = -3.49803305;
    a6 = 0.36037302;
    a7 = -0.0104432413;
    a8 = -0.793385684;
    a9 = 1.39643306;
    a10 = -0.149144925;
    a11 = 0.00441015512;
    a12 = 0.0839387178;
    a13 = -0.186408848;
    a14 = 0.0203367881;
    a15 = -0.000609579263;

    REAL b,b0,b1,b2,b3;
    REAL db,db0,db1,db2,db3;
    REAL mu_gas, dmudP;
    
    b0 = a0 + a1 * Pr + a2 * std::pow(Pr,2.0) + a3 * std::pow(Pr,3.0);
    db0 = a1 + 2.0 * a2 * Pr + 3.0 * a3 * std::pow(Pr,2.0);
    
    b1 = a4 + a5 * Pr + a6 * std::pow(Pr,2.0) + a7 * std::pow(Pr,3.0);
    db1 = a5 + 2.0 * a6 * Pr + 3.0 * a7 * std::pow(Pr,2.0);
    
    b2 = a8 + a9 * Pr + a10 * std::pow(Pr,2.0) + a11 * std::pow(Pr,3.0);
    db2 = a9 + 2.0 * a10 * Pr + 3.0 * a11 * std::pow(Pr,2.0);
    
    b3 = a12 + a13 * Pr + a14 * std::pow(Pr,2.0) + a15 * std::pow(Pr,3.0);
    db3 = a13 + 2.0 * a14 * Pr + 3.0 * a15 * Pr;
    
    b = b0 + b1 * Tr + b2 * std::pow(Tr,3.0) + b3 * std::pow(Tr,3.0);
    db = db0 + db1 * Tr + db2 * std::pow(Tr,3.0) + db3 * std::pow(Tr,3.0);
    
    mu_gas = (mu_std / Tr) * exp(b) * (0.001);
    dmudP = mu_gas * db * (0.001);
    
    mu[0] = mu_gas / GetMu();
    mu[2] = dmudP / GetMu();
}

/** @brief Compressibility - 1/pa $c$ */
void GasPhase::Compressibility(TPZManVector<REAL> &c, TPZManVector<REAL> state_vars)
{
    c[0] = 0.0;
    c[2] = 0.0;
}

/** @brief Kr - $k_{r}$ */
void GasPhase::Kr(TPZManVector<REAL> &kr, TPZManVector<REAL> state_vars){

    REAL So = state_vars[2];
    REAL Sw = state_vars[3];
    
    kr[0] = 1-So-Sw;
    kr[1] = 0.0;
    kr[2] = 0.0;
    kr[3] = -1.0;
    kr[4] = -1.0;
    
//    if (fIsNonlinearKrQ) {
//        kr[0] = Se*Se;
//        kr[1] = 0.0;
//        kr[2] = 0.0;
//        kr[3] = 2.0*Se*(1.0)/(1.0-Swr-Sor);
//        kr[4] = 0.0;
//    }
//    else{
//        
//        kr[0] = Se;
//        kr[1] = 0.0;
//        kr[2] = 0.0;
//        kr[3] = 1.0*(1.0)/(1.0-Swr-Sor);
//        kr[4] = 0.0;
//        
//    }
    
}

/** @brief Pc - $P_{c}$ */
void GasPhase::Pc(TPZManVector<REAL> &pc, TPZManVector<REAL> state_vars){
    
//    REAL So = state_vars[2];
//    REAL Sw = state_vars[3];
    
    pc[0] = 0.0;
    pc[1] = 0.0;
    pc[2] = 0.0;
    pc[3] = 0.0;
    pc[4] = 0.0;
    
}

/** @brief Computes the pseudo critic pressure of Gas $ */
REAL GasPhase::Ppc(){
    REAL PpcHC, Ppct;
    PpcHC = 706.0 - 51.7 * fGas_gamma - 11.1 * (fGas_gamma*fGas_gamma);
    Ppct = (PpcHC - 80.0*fyCO2 + 130.0*fyH2S - 250.0*fyN2) * (6894.76);
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
    REAL Pg = fabs(state_vars[1] + fPstd / GetPRef()); // Negative values generate complex numbers
    
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
    
    if (!(1.2 < Tr && Tr < 2.4) ||  !(0.0 < Pr < 13.0) ) {
        std::cout<< " Tr or Pr put pf this conditions = 1.2 < Tr < 2.4;  0 < Pr < 13"<< std::endl;
        std::cout<< " Pr =" << Pr << std::endl;
        std::cout<< " Tr =" << Tr << std::endl;
        DebugStop();
    }
    
    /* Computing z factor and derivatives */
    
    REAL a,a1=1.39,a2=0.92,a3=0.36,a4=0.101;
    a = a1 * (sqrt(Tr - a2)) - a3 * Tr - a4;
    
    REAL b,b1=0.62,b2=0.23,b3=0.066,b4=0.86,b5=0.037,b6=0.32,b7=9.0,dbdPr;
    b = (b1 - b2*Tr)*Pr + (Pr*Pr)*((b3/(Tr - b4)) - b5) + (b6*std::pow(Pr,6.0))*(std::pow(10.0,b7*(1.0 - Tr)));
    
    dbdPr = b1 + 3.0 * std::pow(2.0,1.0 + b7*(1.0 - Tr)) * std::pow(5.0,b7*(1.0 - Tr)) * b6 * std::pow(5.0,Pr) - b2*Tr +
    2.0 * Pr *(-b5 + (b3/(-b4 + Tr)));
    
    REAL c,c1=0.132,c2=0.32;
    c = c1 - c2 * log10(Tr);
    
    REAL d,d1=0.3106,d2=0.49,d3=0.1824;
    d = std::pow(10.0,(d1 - d2 * Tr  + d3 * Tr*Tr));
    
    REAL zv = 0.0;
    REAL dzdPr = 0.0;
    
    zv = a + (1.0 - a) / (exp(b)) + c * std::pow(Pr,d);
    dzdPr = c*d*std::pow(Pr,d-1.0) - (1.0 - a)*exp(-b)*dbdPr;
    
//    if (sign)
//    {
//        z[0] = zv;
//        z[2] = dzdPr*(1.0/Ppc);
//        return;
//    }
    
    z[0] = zv;
    z[2] = dzdPr*(1.0/Ppc);
    return;
    
}

/** @brief Computes the gas formation volume factor  $ */
void GasPhase::Bg(TPZManVector<REAL> &bg, TPZManVector<REAL> state_vars){
    
    /* Getting the pressure gas phase */
    REAL Pg = (state_vars[1] + fPstd / GetPRef());
    
    TPZManVector<REAL> z(5,0.0);
    this->Z(z, state_vars);
    
    REAL bgv, dbgdP;
    
    bgv = z[0] * ((fPstd/GetPRef())*GetTRes())/(fTstd*Pg) ;
    dbgdP = z[2] * ((fPstd/GetPRef())*GetTRes())/(fTstd*Pg) - z[0] * ((fPstd/GetPRef())*GetTRes())/(fTstd*Pg*Pg) ;
    
    bg[0] = bgv;
    bg[2] = dbgdP;
    
}

/** @brief viscosity - Pa s  $\mu$ */
REAL GasPhase::Viscosity_std(){
    
    REAL logdg, DeltaUg1N2, DeltaUg1CO2, DeltaUg1H2S, MUg1HC, MUg1;
    REAL Tres_F = (GetTRes() - 273.15) * (9.0/5.0) + 32.0;
    
    logdg = log(fGas_gamma) / log(10.0);
    DeltaUg1N2 = fyN2 * (0.00848 *logdg + 0.00959) ;
    DeltaUg1CO2 = fyCO2 * (0.00908 * logdg + 0.00624);
    DeltaUg1H2S = fyH2S * (0.00848 * logdg + 0.00373) ;
    MUg1HC = (0.00001709 - 0.000002062 * fGas_gamma) * Tres_F + 0.008188 - 0.00615 * logdg ;
    MUg1 = MUg1HC + DeltaUg1N2 + DeltaUg1CO2 + DeltaUg1H2S;
    
    return MUg1;
    
}
