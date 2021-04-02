//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzmultiphase.h"
#include "pzlog.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"
#include <math.h> 
#include <iostream>

#ifdef PZ_LOG
static TPZLogger logger("pz.multiphase");
#endif

#ifdef PZ_LOG
static TPZLogger logdata("pz.material.multiphase.data");
#endif

TPZMultiphase::TPZMultiphase(): 
TPZRegisterClassId(&TPZMultiphase::ClassId),
TPZMaterial()
{
    fDim = 2;
    fTheta = 1.0;
    fGamma = 1.0;
    fDeltaT = 1.0;
    fmatId = 0;
    fLref = 1.0;
    fPref= 1.0;
    fKref = 1.0;
    fRhoref = 1.0;
    fEtaref = 1.0;
    fxi = 1.0;    
}

TPZMultiphase::TPZMultiphase(int matid, int dim):
TPZRegisterClassId(&TPZMultiphase::ClassId),
TPZMaterial(matid)
{
    // Two-dimensional analysis
    
    if(dim<0 || dim >2){
        DebugStop();
    }
    
    fDim = dim;
    fTheta = 1.0;
    fGamma = 1.0;
    fDeltaT = 1.0;
    fmatId = matid;
    fLref = 1.0;
    fPref= 1.0;
    fKref = 1.0;
    fRhoref = 1.0;
    fEtaref = 1.0;
    fxi = 1.0;
}


TPZMultiphase::~TPZMultiphase()
{
}

int TPZMultiphase::Dimension() const {return fDim;};

int TPZMultiphase::MatId() {return fmatId;}

int TPZMultiphase::NStateVariables() const {return 8;}

void TPZMultiphase::Print(std::ostream &out) {
    out << "name of material : " << Name() << "\n";
    out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
    out << "Base Class properties :";
    TPZMaterial::Print(out);
    out << "\n";
}



// Data set

/**
 * @brief Lame First Parameter.
 * \f$ lamelamda \f$
 */
STATE TPZMultiphase::LameLambda()
{
    REAL lamelambda = 1.0;
    return lamelambda;
}

/**
 * @brief Undrained Lame First Parameter.
 * \f$ lamelamdaU \f$
 */
STATE TPZMultiphase::LameLambdaU()
{
    REAL lamelambdaU = 1.0;
    return lamelambdaU;
}

/**
 * @brief Lame Second Parameter.
 * \f$ lamemu \f$
 */
STATE TPZMultiphase::LameMu()
{
    REAL lameMu = 0.5;
    return lameMu;
}

/**
 * @brief Biot parameter Parameter.
 * \f$ lamemu \f$
 */
STATE TPZMultiphase::BiotAlpha()
{
    REAL alpha = 0.0;
    return alpha;
}

/**
 * @brief //Se o 1/M coeficiente poroelastico de armazenamento a volume constante.
 * \f$ Se \f$
 */
STATE TPZMultiphase::Se()
{
    REAL Se = 0.0;
    return Se;
}


/** Capilar pressure \f$ pc = pc( Sw ) \f$ */
void TPZMultiphase::CapillaryPressure(REAL Sw, REAL &pc, REAL &DpcDSo){
//   Bentsen and Anli Capillary Pressure model
    REAL PsiToPa = 6894.76;
    REAL Pct = 1.0, Pcs = 0.5, Swi = 0.0;
    
//     pc = 0.0;
//     DpcDSo = 0.0;    
    
     pc = (Pct - Pcs * log ((0.000000001 + Sw - Swi)/(1.0 - Swi)))*PsiToPa/this->fPref;
     DpcDSo = -(Pcs/(0.000000001 + Sw - Swi))*PsiToPa/this->fPref;
}

/** Oil relative permeability \f$ Kro = 1 - Sw \f$ */
void TPZMultiphase::Kro(REAL Sw, REAL &Kro, REAL &dKroDSw){
    
    
//  // Zero model
//  Kro = 0.0;
//  dKroDSw = 0.0;
    
    // linear model
    Kro = 1.0-Sw;
    dKroDSw = -1.0;
    
    //  // Non-linear model
    //  Kro = (1-Sw)*(1-Sw);
    //  dKroDSw = 2.0*(1-Sw)*(-1.0);
    
    //  // Non-linear model
    //  Kro = (1-Sw)*(1-Sw)*(1-Sw);
    //  dKroDSw = 3.0*(1-Sw)*(1-Sw)*(-1.0);
}


/** Water relative permeability \f$ Krw = Sw \f$ */
void TPZMultiphase::Krw(REAL Sw, REAL &Krw, REAL &dKrwDSw){
    
//  // Zero model
//  Krw = 1;
//  dKrwDSw = 0.0;
    
    // linear model
    Krw = Sw;
    dKrwDSw = 1.0;
    
    //  // Non-linear model
    //  Krw = (Sw)*(Sw);
    //  dKrwDSw = 2.0*(Sw)*(1.0);
    
    //  // Non-linear model
    //  Krw = (Sw)*(Sw)*(Sw);
    //  dKrwDSw = 3.0*(Sw)*(Sw)*(1.0);
    
}

/** Water saturation maximum value of the fractional flow product function \f$ S* = Sw* \f$ */
void TPZMultiphase::SWaterstar(REAL &Swstar, REAL &Po, REAL &Sw)
{
    
    REAL waterdensity,oildensity;
    REAL dwaterdensitydp,doildensitydp;
    REAL waterviscosity, oilviscosity;
    REAL dwaterviscositydp, doilviscositydp;
//    REAL SUser=0.98; 
    
    RhoWater(Po,waterdensity,dwaterdensitydp);
    RhoOil(Po,oildensity,doildensitydp);
    WaterViscosity(Po,waterviscosity,dwaterviscositydp);
    OilViscosity(Po,oilviscosity,doilviscositydp);
    
    // Zero model
    Swstar = Sw;
    
    // linear model
    Swstar = (waterviscosity*oildensity)/(oilviscosity*waterdensity+waterviscosity*oildensity);
    
//     if(Swstar <= SUser)
//     {
//         Swstar = SUser;
//     }
    
//      // Non-linear model
//     Swstar = (sqrt(oilviscosity*waterviscosity*waterdensity*oildensity)-(waterviscosity*oildensity))/
//               (oilviscosity*waterdensity-waterviscosity*oildensity);

    
    //  // Non-linear model cubic it is a large simbolic expression, so we got a problem here!!!!!
    //  Krw = (Sw)*(Sw)*(Sw);
    //  dKrwDSw = 3.0*(Sw)*(Sw)*(1.0);
    
}

/** Rock porosity \f$ Phi = Phi( P ) \f$ */
void TPZMultiphase::Porosity(REAL po, REAL &poros, REAL &dPorosDpo){
    const REAL comp = (0.0e-10)*(fPref);
    const REAL pref = (1.0e6)/(fPref);
    const REAL porosRef = 0.1;
    poros = porosRef*exp(comp*(po-pref));
    dPorosDpo = comp*porosRef*exp(comp*(po-pref));
}


/** Oil density  \f$ RhoOil = RhoOil( P ) \f$ */
void TPZMultiphase::RhoOil(REAL po, REAL &RhoOil, REAL &dRhoOilDpo){
    const REAL Oilcomp = (0.0e-8)*(fPref);
    const REAL pref = (1.0e6)/(fPref);
    RhoOil = RhoOilSC()*exp(Oilcomp*(po-pref))/(fRhoref);
    dRhoOilDpo = Oilcomp*RhoOilSC()*exp(Oilcomp*(po-pref))/(fRhoref);
}

/** Water density  \f$ RhoWater = RhoWater( P ) \f$ */
void TPZMultiphase::RhoWater(REAL po, REAL &RhoWater, REAL &dRhoWaterDpo){
    const REAL Watercomp = (0.0e-9)*(fPref);
    const REAL pref = (1.0e6)/(fPref);
    RhoWater = RhoWaterSC()*exp(Watercomp*(po-pref))/(fRhoref);
    dRhoWaterDpo = Watercomp*RhoWaterSC()*exp(Watercomp*(po-pref))/(fRhoref);
}


/** Oil viscosity  \f$ OilViscosity = OilViscosity( P ) \f$ */
void TPZMultiphase::OilViscosity(REAL po, REAL &OilViscosity, REAL &dOilViscosityDpo){
    const REAL OilViscRef = (1.0e-3)/(fEtaref);
    OilViscosity = OilViscRef;
    dOilViscosityDpo = 0;
}

/** Water viscosity  \f$ WaterViscosity = WaterViscosity( P ) \f$ */
void TPZMultiphase::WaterViscosity(REAL po, REAL &WaterViscosity, REAL &dWaterViscosityDpo){
    const REAL WaterViscRef = (1.0e-3)/(fEtaref);
    WaterViscosity = WaterViscRef;
    dWaterViscosityDpo = 0;
}

/** Oil mobility. \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$  */
void TPZMultiphase::OilLabmda(REAL &OilLabmda, REAL Po, REAL Sw, REAL &dOilLabmdaDPo, REAL &dOilLabmdaDSw){
    
    REAL krOil,Oilviscosity,OilDensity;
    REAL dKroDSw,DOilviscosityDp,DOilDensityDp;
    
    Kro(Sw, krOil, dKroDSw);
    OilViscosity(Po, Oilviscosity,DOilviscosityDp);
    RhoOil(Po, OilDensity,DOilDensityDp);
    
    OilLabmda = ((krOil)/(Oilviscosity))*(OilDensity);
    dOilLabmdaDPo = (DOilDensityDp/Oilviscosity)*(krOil);
    dOilLabmdaDSw = (OilDensity/Oilviscosity)*(dKroDSw);
}


/** Water mobility. \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::WaterLabmda(REAL &WaterLabmda, REAL Pw, REAL Sw, REAL &dWaterLabmdaDPw, REAL &dWaterLabmdaDSw){
    
    REAL krWater,Waterviscosity,WaterDensity;
    REAL dKrwDSw,DWaterviscosityDp,DWaterDensityDp;
    
    Krw(Sw, krWater, dKrwDSw);
    WaterViscosity(Pw, Waterviscosity,DWaterviscosityDp);
    RhoWater(Pw, WaterDensity,DWaterDensityDp);
    
    WaterLabmda = ((krWater)/(Waterviscosity))*(WaterDensity);
    dWaterLabmdaDPw = (DWaterDensityDp/Waterviscosity)*(krWater);
    dWaterLabmdaDSw = (WaterDensity/Waterviscosity)*(dKrwDSw);
}



/** Bulk mobility. \f$ \lambda = \lambda( pw , Sw ) \f$  */
void TPZMultiphase::Labmda(REAL &Labmda, REAL Pw, REAL Sw, REAL &dLabmdaDPw, REAL &dLabmdaDSw){
    
    REAL OilLabmdav;
    REAL dOilLabmdaDpo,dOilLabmdaDSw;
    
    REAL WaterLabmdav;
    REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
    
    OilLabmda(OilLabmdav, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
    WaterLabmda(WaterLabmdav, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);
    
    Labmda = OilLabmdav + WaterLabmdav;
    dLabmdaDPw = dOilLabmdaDpo+dWaterLabmdaDpo;
    dLabmdaDSw = dOilLabmdaDSw+dWaterLabmdaDSw;
}

/** Oil fractional flux. \f$ f_{Oil} = f_{Oil}( pw , Sw ) \f$  */
void TPZMultiphase::fOil(REAL &fOil, REAL Pw, REAL Sw, REAL &dfOilDPw, REAL &dfOilDSw){
    
    REAL OilLabmdab;
    REAL dOilLabmdaDpo,dOilLabmdaDSw;
    
    //  REAL WaterLabmdav;
    //  REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
    
    REAL Lambdab;
    REAL dLabmdaDPw, dLabmdaDSw;
    
    OilLabmda(OilLabmdab, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
    //  WaterLabmda(WaterLabmdav, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);
    Labmda(Lambdab,Pw,Sw,dLabmdaDPw,dLabmdaDSw);
    
    fOil = ((OilLabmdab)/(Lambdab));
    dfOilDPw = ((dOilLabmdaDpo)/(Lambdab))-((OilLabmdab)/(Lambdab*Lambdab))*dLabmdaDPw;
    dfOilDSw = ((dOilLabmdaDSw)/(Lambdab))-((OilLabmdab)/(Lambdab*Lambdab))*dLabmdaDSw;
    
}


/** Water fractional flux. \f$ f_{Water} = f_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::fWater(REAL &fWater, REAL Pw, REAL Sw, REAL &dfWaterDPw, REAL &dfWaterDSw){
    
    //  REAL OilLabmdab;
    //  REAL dOilLabmdaDpo,dOilLabmdaDSw;
    
    REAL WaterLabmdab;
    REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
    
    REAL Lambdab;
    REAL dLabmdaDPw, dLabmdaDSw;
    
    //  OilLabmda(OilLabmdab, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
    WaterLabmda(WaterLabmdab, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);
    Labmda(Lambdab,Pw,Sw,dLabmdaDPw,dLabmdaDSw);
    
    fWater = ((WaterLabmdab)/(Lambdab));
    dfWaterDPw = ((dWaterLabmdaDpo)/(Lambdab))-((WaterLabmdab)/(Lambdab*Lambdab))*dLabmdaDPw;
    dfWaterDSw = ((dWaterLabmdaDSw)/(Lambdab))-((WaterLabmdab)/(Lambdab*Lambdab))*dLabmdaDSw;
    
}

/**
* Fractional product function, water decrease direction (dSw/dt < 0).
* \f$ fw*  \f$
*/
void TPZMultiphase::fWaterstar(REAL &fWstar, REAL Pw, REAL Sw, REAL &dfWstarDPw, REAL &dfWstarDSw)
{
    REAL fwater, foil;
    REAL dfwaterdP, dfoildP;
    REAL dfwaterdSw, dfoildSw; 
    
    REAL Swcutoff;
    SWaterstar(Swcutoff,Pw,Sw);      
    
    if(Sw <= Swcutoff)
    {
        fOil(foil,Pw,Sw,dfoildP,dfoildSw);
        fWater(fwater,Pw,Sw,dfwaterdP,dfwaterdSw);
        fWstar = foil*fwater;
        dfWstarDPw = (dfoildP)*fwater + (dfwaterdP)*foil;
        dfWstarDSw = (dfoildSw)*fwater + (dfwaterdSw)*foil;        
    }
    else
    {
        fOil(foil,Pw,Swcutoff,dfoildP,dfoildSw);
        fWater(fwater,Pw,Swcutoff,dfwaterdP,dfwaterdSw);
        fWstar = foil*fwater;
        dfWstarDPw = 0.0*((dfoildP)*fwater + (dfwaterdP)*foil);
        dfWstarDSw = 0.0*((dfoildSw)*fwater + (dfwaterdSw)*foil);        
    }
    
    
}

/**
* Fractional product function, oil decrease direction (dSw/dt > 0).
* \f$ fo*  \f$
*/
void TPZMultiphase::fOilstar(REAL &fOstar, REAL Pw, REAL Sw, REAL &dfOstarDPw, REAL &dfOstarDSw)
{
    REAL fwater, foil;
    REAL dfwaterdP, dfoildP;
    REAL dfwaterdSw, dfoildSw;
    REAL Swcutoff;   
    SWaterstar(Swcutoff,Pw,Sw);
    
    
    if(Sw >= Swcutoff)
    {
        fOil(foil,Pw,Sw,dfoildP,dfoildSw);
        fWater(fwater,Pw,Sw,dfwaterdP,dfwaterdSw);
        fOstar = foil*fwater;
        dfOstarDPw = (dfoildP)*fwater + (dfwaterdP)*foil;
        dfOstarDSw = (dfoildSw)*fwater + (dfwaterdSw)*foil;        
    }
    else
    {
        fOil(foil,Pw,Swcutoff,dfoildP,dfoildSw);
        fWater(fwater,Pw,Swcutoff,dfwaterdP,dfwaterdSw);
        fOstar = foil*fwater;
        dfOstarDPw = 0.0*((dfoildP)*fwater + (dfwaterdP)*foil);
        dfOstarDSw = 0.0*((dfoildSw)*fwater + (dfwaterdSw)*foil);        
    }    
}

/**
*  Fractional product upwind function, Gdotn > 0 means water decrease (dSw/dt < 0), Gdotn < 0 means water increase (dSw/dt > 0).
* \f$ f*  \f$
*/
void TPZMultiphase::fstar(REAL &fStar, REAL Pw, REAL Sw, REAL Gdotn, REAL &dfstarDPw, REAL &dfstarDSw)
{
    REAL fwaterStar, foilStar;
    REAL dfwaterStardP, dfoilStardP;
    REAL dfwaterStardSw, dfoilStardSw;

    if(Gdotn > 0.0)
    {
            fWaterstar(fwaterStar,Pw,Sw,dfwaterStardP,dfwaterStardSw);
            fStar = fwaterStar;
            dfstarDPw = dfwaterStardP;
            dfstarDSw = dfwaterStardSw;
    }
    else
    {
            fOilstar(foilStar,Pw,Sw,dfoilStardP,dfoilStardSw);
            fStar = foilStar;
            dfstarDPw = dfoilStardP;
            dfstarDSw = dfoilStardSw;        
//             if(fabs(Gdotn) <= 1.0e-14)
//             {
//                 fStar = 0.0;
//                 dfstarDPw = 0.0;
//                 dfstarDSw = 0.0;
//             }

    }
    
}


// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////


/** Capilar pressure \f$ pc = pc( Sw ) \f$ */
void TPZMultiphase::CapillaryPressure(BFadREAL So, BFadREAL &pc){
    pc = 0.0;
}

/** Oil relative permeability \f$ Kro = 1 - Sw \f$ */
void TPZMultiphase::Kro(BFadREAL Sw, BFadREAL &Kro){
    Kro = 1.0-Sw.val();
}


/** Water relative permeability \f$ Krw = Sw \f$ */
void TPZMultiphase::Krw(BFadREAL Sw, BFadREAL &Krw){
    Krw = Sw;
}


/** Rock porosity \f$ Phi = Phi( P ) \f$ */
void TPZMultiphase::Porosity(BFadREAL po, BFadREAL &poros){
    const REAL comp = 1.0e-10;
    const REAL pref = 1.0e6;
    const REAL porosRef = 0.3;
    poros = porosRef*exp(comp*((po.val())-pref));
}


/** Oil density  \f$ RhoOil = RhoOil( P ) \f$ */
void TPZMultiphase::RhoOil(BFadREAL po, BFadREAL &RhoOil){
    const REAL Oilcomp = 0.0e-8;
    const REAL pref = 1.0e6;
    RhoOil = RhoOilSC()*exp(Oilcomp*((po.val())-pref));
}

/** Water density  \f$ RhoWater = RhoWater( P )  \f$ */
void TPZMultiphase::RhoWater(BFadREAL po, BFadREAL &RhoWater){
    const REAL Watercomp = 0.0e-9;
    const REAL pref = 1.0e6;
    RhoWater = RhoWaterSC()*exp(Watercomp*((po.val())-pref));
}


/** Oil viscosity  \f$ OilViscosity = OilViscosity( P ) \f$ */
void TPZMultiphase::OilViscosity(BFadREAL po, BFadREAL &OilViscosity){
    const REAL OilViscRef = 1.0e-3;
    OilViscosity = OilViscRef;
}

/** Water viscosity  \f$ WaterViscosity = WaterViscosity( P ) \f$ */
void TPZMultiphase::WaterViscosity(BFadREAL po, BFadREAL &WaterViscosity){
    const REAL WaterViscRef = 1.0e-3;
    WaterViscosity = WaterViscRef;
}

/** Oil mobility. \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$  */
void TPZMultiphase::OilLabmda(BFadREAL OilLabmda, BFadREAL Po, BFadREAL &Sw){
    
    BFadREAL krOil,Oilviscosity,OilDensity;
    
    Kro(Sw, krOil);
    OilViscosity(Po,Oilviscosity);
    RhoOil(Po, OilDensity);
    
    OilLabmda = ((krOil)/(Oilviscosity))*(OilDensity);
    
}


/** Water mobility. \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::WaterLabmda(BFadREAL WaterLabmda, BFadREAL Pw, BFadREAL &Sw){
    
    BFadREAL krWater,Waterviscosity,WaterDensity;
    
    Krw(Sw, krWater);
    WaterViscosity(Pw,Waterviscosity);
    RhoWater(Pw,WaterDensity);
    
    WaterLabmda = ((krWater)/(Waterviscosity))*(WaterDensity);
    
}


/** Bulk mobility. \f$ \lambda = \lambda( pw , Sw ) \f$  */
void TPZMultiphase::Labmda(BFadREAL Labmda, BFadREAL Pw, BFadREAL &Sw){
    
    BFadREAL OilLabmdaf, WaterLabmdaf;
    
    OilLabmda(OilLabmdaf, Pw, Sw);
    WaterLabmda(WaterLabmdaf, Pw, Sw);
    
    Labmda = OilLabmdaf + WaterLabmdaf;
    
}


/** Oil fractional flux. \f$ f_{Oil} = f_{Oil}( pw , Sw ) \f$  */
void TPZMultiphase::fOil(BFadREAL fOil, BFadREAL Pw, BFadREAL &Sw)
{
}


/** Water fractional flux. \f$ f_{Water} = f_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::fWater(BFadREAL fWater, BFadREAL Pw, BFadREAL &Sw)
{
    
}


// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////


// Constant data

/** Oil density at standar conditions - kg/m3 */
REAL TPZMultiphase::RhoOilSC(){
    return 1000.0;
}

/** Water density at standar conditions - kg/m3 */
REAL TPZMultiphase::RhoWaterSC(){
    return 1000.0;
}

/** Gravity */
TPZFMatrix<REAL> TPZMultiphase::Gravity(){
    TPZFMatrix<REAL> gravity(3,1,0.0);
    gravity(0,0) =  0.0;
    gravity(1,0) =  0.0;
    gravity(2,0) =  0.0;
    gravity *= (1.0/(fPref/(fLref*fRhoref)));
    return gravity;
}

/** Permeabilidade absoluta */
void TPZMultiphase::K(TPZFMatrix<REAL> &Kab){
    Kab.Resize(3,3);
    Kab.Zero();
    Kab(0,0) = 1.0e-13;
    Kab(1,1) = 1.0e-13;
    Kab(2,2) = 0.0e-13;
    Kab *= 1.0/(fKref);
    
}

/** Inverse Permeabilidade absoluta */
TPZFMatrix<REAL>  TPZMultiphase::Kinv(TPZFMatrix<REAL> &Kab){
    TPZFMatrix<REAL> Kinverse(3,3,0.0);
    
    REAL Constant = (-1.0*Kab(0,1)*Kab(1,0)+Kab(0,0)*Kab(1,1));
    
    Kinverse(0,0) =     1.0*Kab(1,1)/(Constant);
    Kinverse(0,1) = -   1.0*Kab(0,1)/(Constant);
    Kinverse(1,0) = -   1.0*Kab(1,0)/(Constant);
    Kinverse(1,1) =     1.0*Kab(0,0)/(Constant);
    return Kinverse;
}


/** @brief Absolute permeability. */
void TPZMultiphase::LoadKMap(std::string MaptoRead)
{
    
    std::string FileName;
    std::string stringTemp;
    FileName = MaptoRead;
    
    // Definitions
    int64_t numelements=0;
    int64_t numKData=0;
    
    //  Scanning for total Number geometric elements
    int64_t NumEntitiestoRead;
    TPZStack <std::string> SentinelString;
    {
        
        // reading a general mesh information by filter
        std::ifstream read (FileName.c_str());
        std::string FlagString;
        //      int flag = -1;
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            if(str == "KABSOLUTE" || str == "KABSOLUTE\r")
            {
                SentinelString.push_back(str);
            }
            
        }
        
        FlagString = "EndReading";
        SentinelString.push_back(FlagString);
    }
    
    NumEntitiestoRead = SentinelString.size();
    
    TPZStack <int> GeneralData(NumEntitiestoRead,0);
    TPZStack <int> DataToProcess(NumEntitiestoRead,-1);
    
    {
        // reading a general map information by filter
        std::ifstream read (FileName.c_str());
        std::string FlagString;
        int64_t cont = 0;
        
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            
            // Reading General Data
            if(str == SentinelString[cont])
            {
                FlagString = str;
            }
            if(SentinelString[cont] == "" || SentinelString[cont] == "\r")
            {
                cont++;
            }
            
            if(SentinelString[cont] == "EndReading")
            {
                break;
            }
            
            if( (str != "" || str != "\r") && FlagString == SentinelString[cont])
            {
                // Data scaning
                while (read) {
                    char buftemp[1024];
                    read.getline(buftemp, 1024);
                    std::string strtemp(buftemp);
                    GeneralData[cont]++;
                    if(strtemp == "" || strtemp == "\r")
                    {
                        FlagString = "";
                        GeneralData[cont]--;
                        cont++;
                        break;
                    }
                }
                
            }
            
            
        }
    }
    
    
    for (int64_t i = 0 ; i < NumEntitiestoRead; i++ )
    {
        
        if(SentinelString[i] == "KABSOLUTE" || SentinelString[i] == "KABSOLUTE\r")
        {
            numKData=GeneralData[i];
            DataToProcess[i]=0;
        }
    }
    
    numelements=numKData;
    TPZFMatrix<REAL> Kabsolute(3,3);
    Kabsolute.Resize(3,3);
    Kabsolute.Zero();
    //  TPZStack<TPZFMatrix<REAL> > KabsoluteMap(numelements,0);
    TPZStack<TPZFMatrix<REAL> > KabsoluteMap(100,0);
    
    int64_t elementId = 0;
    int64_t ContOfKs = 0;
    
    REAL kxx , kxy, kxz;
    REAL kyx , kyy, kyz;
    REAL kzx , kzy, kzz;
    
    
    
    {
        
        // reading a general mesh information by filter
        std::ifstream read (FileName.c_str());
        std::string FlagString;
        int64_t cont = 0;
        //      int dim = 0;
        //      int flag = 0;
        while(read)
        {
            char buf[1024];
            read.getline(buf, 1024);
            std::string str(buf);
            std::string strtemp="InitialState";
            
            // Reading General Data
            if(str == SentinelString[cont])
            {
                FlagString = str;
            }
            
            if(SentinelString[cont] == "" || SentinelString[cont] == "\r")
            {
                cont++;
            }
            
            if(SentinelString[cont] == "EndReading")
            {
                break;
            }
            
            if( (str != "" || str != "\r" )&& FlagString == SentinelString[cont])
            {
                // Data scaning
                while (read)
                {
                    
                    switch (DataToProcess[cont])
                    {
                        case 0:
                        {
                            //"KABSOLUTE"
                            if (GeneralData[cont] != 0)
                            {
                                read >> elementId;
                                read >> kxx;
                                read >> kxy;
                                read >> kxz;
                                read >> kyx;
                                read >> kyy;
                                read >> kyz;
                                read >> kzx;
                                read >> kzy;
                                read >> kzz;
                                
                                Kabsolute(0,0)=(1.0/fKref)*kxx;
                                Kabsolute(0,1)=(1.0/fKref)*kxy;
                                Kabsolute(0,2)=(1.0/fKref)*kxz;
                                Kabsolute(1,0)=(1.0/fKref)*kyx;
                                Kabsolute(1,1)=(1.0/fKref)*kyy;
                                Kabsolute(1,2)=(1.0/fKref)*kyz;
                                Kabsolute(2,0)=(1.0/fKref)*kzx;
                                Kabsolute(2,1)=(1.0/fKref)*kzy;
                                Kabsolute(2,2)=(1.0/fKref)*kzz;
                                
                                KabsoluteMap[elementId]=Kabsolute;
                                
                                ContOfKs++;
                            }
                            if(ContOfKs == numKData)
                            {
                                strtemp = "";
                            }
                        }
                            break;
                        default:
                        {
                            strtemp = "";
                        }
                            break;
                    }
                    
                    if(strtemp == "" || strtemp == "\r")
                    {
                        FlagString = "";
                        cont++;
                        break;
                    }
                }
                
            }
            
            
        }
    }
    
    
    SetKMap(KabsoluteMap);
    
}



// Contribute methods

void TPZMultiphase::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
    DebugStop();
}


void TPZMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{

#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 5 )
    {
        std::cout << " Error. datavec size is different from 5 \n";
        DebugStop();
    }
#endif
    
    
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU     =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiQ     =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiP     =  datavec[2].phi;
    TPZFMatrix<REAL>  &phiS     =  datavec[3].phi;
//    TPZFMatrix<REAL>  &phiQG    =  datavec[4].phi;
    
    TPZFMatrix<REAL> &dphiU     =  datavec[0].dphix;
    //    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP     =  datavec[2].dphix;
    TPZFMatrix<REAL> &dphiS     =  datavec[3].dphix;
    
    // number of test functions for each state variable
    int phrU, phrQ, phrP, phrS, phrQG;
    phrU = phiU.Rows();
    phrQ = datavec[1].fVecShapeIndex.NElements();
    phrP = phiP.Rows();
    phrS = phiS.Rows();
    phrQG = datavec[4].fVecShapeIndex.NElements();
    
    // blocks
    int UStateVar = 2;
    int FirstU  = 0;
    int FirstQ  = phrU * UStateVar + FirstU;
    int FirstP  = phrQ + FirstQ;
    int FirstS  = phrP + FirstP;
//    int FirstQG = FirstS;
    
    //  Getting and computing another required data
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    REAL Gamma = this->fGamma;
    int GeoID = datavec[1].gelElId;
    TPZFMatrix<REAL> Kabsolute;
    TPZFMatrix<REAL> Kinverse;
    TPZFMatrix<REAL> Gfield;
    if (fYorN)
    {
        Kabsolute=this->fKabsoluteMap[GeoID];
    }
    else
    {
        this->K(Kabsolute);
    }
    
    Kinverse=this->Kinv(Kabsolute);
    Gfield = this->Gravity();
    
    TPZManVector<STATE,3> sol_u =    datavec[0].sol[0];
    TPZManVector<STATE,3> sol_q =    datavec[1].sol[0];
    TPZManVector<STATE,3> sol_p =    datavec[2].sol[0];
    TPZManVector<STATE,3> sol_s =    datavec[3].sol[0];
    TPZManVector<STATE,3> sol_qg =   datavec[4].sol[0];

    TPZFMatrix<STATE> dsol_u = datavec[0].dsol[0];
    TPZFMatrix<STATE> dsol_q =datavec[1].dsol[0];
    TPZFMatrix<STATE> dsol_p =datavec[2].dsol[0];
    TPZFMatrix<STATE> dsol_s =datavec[3].dsol[0];
    
    REAL LambdaL, LambdaLU, MuL;
    REAL Balpha, Sestr;
    
    REAL rockporosity, oildensity, waterdensity;
    REAL drockporositydp, doildensitydp, dwaterdensitydp;
    
    REAL oilviscosity, waterviscosity;
    REAL doilviscositydp, dwaterviscositydp;
    
    REAL bulklambda, oillambda, waterlambda;
    REAL dbulklambdadp, doillambdadp, dwaterlambdadp;
    REAL dbulklambdads, doillambdads, dwaterlambdads;
    
    REAL bulkfoil, bulkfwater;
    REAL dbulkfoildp, dbulkfwaterdp;
    REAL dbulkfoilds, dbulkfwaterds;
    
    // Functions computed at point x_{k} for each integration point
    LambdaL     = this->LameLambda();
    LambdaLU    = this->LameLambdaU();
    MuL         = this->LameMu();
    Balpha      = this->BiotAlpha();
    Sestr       = this->Se();

    REAL Pressure = sol_p[0];
    
    int VecPos= 0;
    this->Porosity(sol_p[VecPos], rockporosity, drockporositydp);
    this->RhoOil(sol_p[VecPos], oildensity, doildensitydp);
    this->RhoWater(sol_p[VecPos], waterdensity, dwaterdensitydp);
    this->OilViscosity(sol_p[VecPos], oilviscosity, doilviscositydp);
    this->WaterViscosity(sol_p[VecPos], waterviscosity, dwaterviscositydp);
    this->OilLabmda(oillambda, sol_p[VecPos], sol_s[VecPos], doillambdadp, doillambdads);
    this->WaterLabmda(waterlambda, sol_p[VecPos], sol_s[VecPos], dwaterlambdadp, dwaterlambdads);
    this->Labmda(bulklambda, sol_p[VecPos], sol_s[VecPos], dbulklambdadp, dbulklambdads);
    this->fOil(bulkfoil, sol_p[VecPos], sol_s[VecPos], dbulkfoildp, dbulkfoilds);
    this->fWater(bulkfwater, sol_p[VecPos], sol_s[VecPos], dbulkfwaterdp, dbulkfwaterds);
    
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    //  Contribution of domain integrals for Jacobian matrix
    // n+1 time step
    if(gState == ECurrentState)
    {

        //  Elasticity Block (Equation for elasticity )
        //	Elastic equation
		//	Linear strain operator
		//	Ke Matrix
        TPZFMatrix<REAL>	du(2,2);
		for(int iu = 0; iu < phrU; iu++ )
		{
			//	Derivative for Vx
			du(0,0) = dphiU(0,iu)*datavec[0].axes(0,0)+dphiU(1,iu)*datavec[0].axes(1,0);
			//	Derivative for Vy
			du(1,0) = dphiU(0,iu)*datavec[0].axes(0,1)+dphiU(1,iu)*datavec[0].axes(1,1);

			for(int ju = 0; ju < phrU; ju++)
			{
				//	Derivative for Ux
				du(0,1) = dphiU(0,ju)*datavec[0].axes(0,0)+dphiU(1,ju)*datavec[0].axes(1,0);
				//	Derivative for Uy
				du(1,1) = dphiU(0,ju)*datavec[0].axes(0,1)+dphiU(1,ju)*datavec[0].axes(1,1);
				
				if (this->fPlaneStress == 1)
				{
					/* Plain stress state */
					ek(2*iu + FirstU, 2*ju + FirstU)	     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*du(0,1)		+ (2*MuL)*du(1,0)*du(1,1));
					
					ek(2*iu + FirstU, 2*ju+1 + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*du(1,1)			+ (2*MuL)*du(1,0)*du(0,1));
					
					ek(2*iu+1 + FirstU, 2*ju + FirstU)       += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*du(0,1)			+ (2*MuL)*du(0,0)*du(1,1));
					
					ek(2*iu+1 + FirstU, 2*ju+1 + FirstU)     += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*du(1,1)		+ (2*MuL)*du(0,0)*du(0,1));
				}
				else
				{                   
					/* Plain Strain State */
					ek(2*iu + FirstU,2*ju + FirstU)         += weight*	((LambdaL + 2*MuL)*du(0,0)*du(0,1)	+ (MuL)*du(1,0)*du(1,1));
					
					ek(2*iu + FirstU,2*ju+1 + FirstU)       += weight*	(LambdaL*du(0,0)*du(1,1)			+ (MuL)*du(1,0)*du(0,1));
					
					ek(2*iu+1 + FirstU,2*ju + FirstU)       += weight*	(LambdaL*du(1,0)*du(0,1)			+ (MuL)*du(0,0)*du(1,1));
					
					ek(2*iu+1 + FirstU,2*ju+1 + FirstU)     += weight*	((LambdaL + 2*MuL)*du(1,0)*du(1,1)	+ (MuL)*du(0,0)*du(0,1));
                    
				}
			}
		}
        
		//	Matrix Qc
		//	Coupling matrix for Elastic equation
		for(int in = 0; in < phrU; in++ )
		{
			du(0,0) = dphiU(0,in)*datavec[0].axes(0,0)+dphiU(1,in)*datavec[0].axes(1,0);
			du(1,0) = dphiU(0,in)*datavec[0].axes(0,1)+dphiU(1,in)*datavec[0].axes(1,1);
			
			for(int jp = 0; jp < phrP; jp++)
			{
				ek(2*in + FirstU,jp + FirstP)    += (-1.0)*Balpha*weight*(phiP(jp,0)*du(0,0));
				ek(2*in+1 + FirstU,jp + FirstP)  += (-1.0)*Balpha*weight*(phiP(jp,0)*du(1,0));
			}
		}
        

        REAL SaturationAtnplusOne = sol_s[0];
        //  First Block (Equation One) constitutive law
        // Integrate[(Kinv/bulklambda)*dot(v,v), Omega_{e} ]  (Equation One)
        REAL OneOverLambda = 1.0/bulklambda;
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            for (int jq=0; jq<phrQ; jq++)
            {
                
                int jvectorindex    = datavec[1].fVecShapeIndex[jq].first;
                int jshapeindex     = datavec[1].fVecShapeIndex[jq].second;
                
                if (fnewWS)
                {    
                    REAL vec1 = (Kinverse(0,0)*datavec[1].fDeformedDirections(0,ivectorindex)+Kinverse(0,1)*datavec[1].fDeformedDirections(1,ivectorindex));
                    REAL vec2 = (Kinverse(1,0)*datavec[1].fDeformedDirections(0,ivectorindex)+Kinverse(1,1)*datavec[1].fDeformedDirections(1,ivectorindex));
                    
                    REAL dotprod =
                    (phiQ(ishapeindex,0)*vec1) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(0,jvectorindex)) +
                    (phiQ(ishapeindex,0)*vec2) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(1,jvectorindex)) ;  //  dot(K q,v)
                    
                    ek(iq + FirstQ,jq + FirstQ) += weight * OneOverLambda * dotprod;
                }
                else 
                {
                    REAL dotprod =
                    (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(0,jvectorindex)) +
                    (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(1,jvectorindex)) ; //  dot(q,v)
                    
                    ek(iq + FirstQ,jq + FirstQ) += weight * OneOverLambda * dotprod;
                }

            }
            
        }
        
        
        //  First Block (Equation One) constitutive law
        // Integrate[(d(1/bulklambdal)/dS)*dot(q,v), Omega_{e} ]    (Equation One)
        /*REAL OneOverLambda = 1/bulklambda;*/
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            for (int jsat=0; jsat<phrS; jsat++)
            {
                if (fnewWS)
                {
                    REAL dotprod =
                    (Kinverse(0,0)*sol_q[0]+Kinverse(0,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                    (Kinverse(1,0)*sol_q[0]+Kinverse(1,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                    
                    ek(iq + FirstQ,jsat + FirstS) -= weight * dbulklambdads  * OneOverLambda * OneOverLambda * phiS(jsat,0) * dotprod;
                }
                else
                {
                    REAL dotprod =
                    (sol_q[0]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                    (sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) ;
                    
                    ek(iq + FirstQ,jsat + FirstS) -= weight * dbulklambdads  * OneOverLambda * OneOverLambda * phiS(jsat,0) * dotprod;
                }

            }

        }
        
        
        //  First Block (Equation One) constitutive law
        // Integrate[(d(1/bulklambdal)/dP)*dot(q,v), Omega_{e} ]    (Equation One)
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            for (int jp=0; jp<phrP; jp++)
            {
                
                if (fnewWS)
                {
                    REAL dotprod =
                    (Kinverse(0,0)*sol_q[0]+Kinverse(0,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                    (Kinverse(1,0)*sol_q[0]+Kinverse(1,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) ;
                    
                    ek(iq + FirstQ,jp + FirstP) -= weight * dbulklambdadp  * OneOverLambda * OneOverLambda * phiP(jp,0) * dotprod;
                }
                else
                {
                    REAL dotprod =
                    (sol_q[0]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                    (sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) ;
                    
                    ek(iq + FirstQ,jp + FirstP) -= weight * dbulklambdadp  * OneOverLambda * OneOverLambda * phiP(jp,0) * dotprod;
                }
            }
        }        
        
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K dot(v,grad(P)) , Omega_{e}]   (Equation One)
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            for (int jp=0; jp<phrP; jp++)
            {
                
                
                //  Compute grad(W)
                TPZManVector<STATE> dsolp(2,0);
                dsolp[0] = dphiP(0,jp)*datavec[1].axes(0,0)+dphiP(1,jp)*datavec[1].axes(1,0);
                dsolp[1] = dphiP(0,jp)*datavec[1].axes(0,1)+dphiP(1,jp)*datavec[1].axes(1,1);
                
                if (fnewWS)
                {
                    
                    REAL e1e1   =   (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex))*(dsolp[0]);
                    REAL e2e2   =   (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex))*(dsolp[1]);
                    
                    ek(iq + FirstQ,jp + FirstP) += weight * ( e1e1 + e2e2 );
                    
                }
                else 
                {
                    REAL e1e1   =   (Kabsolute(0,0)*(dsolp[0])+
                                     Kabsolute(0,1)*(dsolp[1]))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                    
                    REAL e2e2   =   (Kabsolute(1,0)*(dsolp[0])+
                                     Kabsolute(1,1)*(dsolp[1]))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                    
                    ek(iq + FirstQ,jp + FirstP) += weight * ( e1e1 + e2e2 );
                }
            }
            
        }
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K (dot(f1*rho1*grad(g*z),v)+dot(f1*rho1*grad(g*z),v)) , Omega_{e}]  (Equation One)
        //  dS/dPalpha;
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (Gfield(0,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                REAL e2e2   =   (Gfield(1,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jp=0; jp<phrP; jp++)
                {
                    ek(iq + FirstQ,jp + FirstP) -= weight * ( ( bulkfwater * dwaterdensitydp + bulkfoil * doildensitydp ) * phiP(jp,0) * (e1e1 + e2e2));
                }
                
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(Gfield(0,0))+
                                 Kabsolute(0,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                
                REAL e2e2   =   (Kabsolute(1,0)*(Gfield(0,0))+
                                 Kabsolute(1,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jp=0; jp<phrP; jp++)
                {
                    ek(iq + FirstQ,jp + FirstP) -= weight * ( ( bulkfwater * dwaterdensitydp + bulkfoil * doildensitydp ) * phiP(jp,0) * (e1e1 + e2e2));
                }
            }
            
        }
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K (dot(f1*rho1*grad(g*z),v)+dot(f1*rho1*grad(g*z),v)) , Omega_{e}]  (Equation One)
        //  dP/dPalpha;
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (Gfield(0,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                REAL e2e2   =   (Gfield(1,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jp=0; jp<phrP; jp++)
                {
                    ek(iq + FirstQ,jp + FirstP) -= weight * ( ( dbulkfwaterdp * waterdensity + dbulkfoildp * oildensity ) * phiP(jp,0) * (e1e1 + e2e2));
                }
                
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(Gfield(0,0))+
                                 Kabsolute(0,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                
                REAL e2e2   =   (Kabsolute(1,0)*(Gfield(0,0))+
                                 Kabsolute(1,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jp=0; jp<phrP; jp++)
                {
                    ek(iq + FirstQ,jp + FirstP) -= weight * ( ( dbulkfwaterdp * waterdensity + dbulkfoildp * oildensity ) * phiP(jp,0) * (e1e1 + e2e2));
                }
            }
            
        }        
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K (dot(f1*rho1*grad(g*z),v)+dot(f1*rho1*grad(g*z),v)) , Omega_{e}]  (Equation One)
        //  dS/dSalpha;
        for(int iq=0; iq < phrQ; iq++)  
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (Gfield(0,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                REAL e2e2   =   (Gfield(1,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jsat=0; jsat < phrS; jsat++)
                {
                    ek(iq+ FirstQ,jsat + FirstS) -= weight * ( ( dbulkfwaterds * waterdensity + dbulkfoilds * oildensity ) * phiS(jsat,0) * (e1e1 + e2e2));
                }
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(Gfield(0,0))+
                                 Kabsolute(0,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                
                REAL e2e2   =   (Kabsolute(1,0)*(Gfield(0,0))+
                                 Kabsolute(1,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                for (int jsat=0; jsat<phrS; jsat++)
                {
                    ek(iq + FirstQ,jsat + FirstS) -= weight * ( ( dbulkfwaterds * waterdensity + dbulkfoilds * oildensity ) * phiS(jsat,0) * (e1e1 + e2e2));
                }
            }
            
        }
        
        //  Poroelastic Contribution
        REAL divphiU, divU, dsolU[2][2];
        
        dsolU[0][0] = dsol_u(0,0)*datavec[0].axes(0,0)+dsol_u(1,0)*datavec[0].axes(1,0); // dUx/dx
        dsolU[1][0] = dsol_u(0,0)*datavec[0].axes(0,1)+dsol_u(1,0)*datavec[0].axes(1,1); // dUx/dy

        dsolU[0][1] = dsol_u(0,1)*datavec[0].axes(0,0)+dsol_u(1,1)*datavec[0].axes(1,0); // dUy/dx
        dsolU[1][1] = dsol_u(0,1)*datavec[0].axes(0,1)+dsol_u(1,1)*datavec[0].axes(1,1); // dUy/dy
        divU = dsolU[0][0]+dsolU[1][1]+0.0;
        
        REAL PhiStar = Balpha * divU + Sestr * Pressure;       
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(porosity)/dPalpha
        for(int ip=0; ip < phrP; ip++)
        {
            for (int ju = 0; ju < phrU; ju++)
            {
                TPZFMatrix<REAL>    du(2,2);               

                du(0,0) = dphiU(0,ju)*datavec[0].axes(0,0); // dPhiU/dx dalphax
                du(1,0) = dphiU(0,ju)*datavec[0].axes(0,1); // dPhiU/dy dalphax
                
                du(0,1) = dphiU(1,ju)*datavec[0].axes(1,0); // dPhiU/dx dalphay
                du(1,1) = dphiU(1,ju)*datavec[0].axes(1,1); // dPhiU/dy dalphay        

                divphiU = du(0,0)+du(1,0)+0.0;

                REAL dPhiStardUx = Balpha * (du(0,0)+du(1,0));
                REAL dPhiStardUy = Balpha * (du(0,1)+du(1,1));

                REAL Integratingx = phiP(ip,0) * dPhiStardUx * (waterdensity * (SaturationAtnplusOne) + oildensity * (1 - SaturationAtnplusOne));
                REAL Integratingy = phiP(ip,0) * dPhiStardUy * (waterdensity * (SaturationAtnplusOne) + oildensity * (1 - SaturationAtnplusOne));
                
                ek(ip + FirstP,2*ju + FirstU)       += (-1.0) * weight * Integratingx;
                ek(ip + FirstP,2*ju + 1 + FirstU)   += (-1.0) * weight * Integratingy;
            }
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(porosity)/dPalpha
        for(int ip=0; ip < phrP; ip++)
        {
            for (int jp=0; jp < phrP; jp++)
            {
                REAL dPhiStardP = Sestr * phiP(jp,0);                  
                REAL Integrating = phiP(ip,0) * dPhiStardP * (waterdensity * (SaturationAtnplusOne) + oildensity * (1 - SaturationAtnplusOne));
                ek(ip + FirstP,jp + FirstP) += (-1.0) * weight * Integrating;
            }
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(porosity)/dPalpha
        for(int ip=0; ip < phrP; ip++)
        {
            for (int jp=0; jp < phrP; jp++)
            {
                REAL Integrating = phiP(ip,0) * PhiStar * (dwaterdensitydp * (SaturationAtnplusOne) + doildensitydp * (1 - SaturationAtnplusOne)) * phiP(jp,0);
                ek(ip + FirstP,jp + FirstP) -=  weight * Integrating;
            }
        }        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(porosity)/dPalpha
        for(int ip=0; ip < phrP; ip++)
        {
            for (int jp=0; jp < phrP; jp++)
            {
                REAL Integrating = phiP(ip,0) * drockporositydp * phiP(jp,0) * (waterdensity * (SaturationAtnplusOne) + oildensity * (1 - SaturationAtnplusOne));
                ek(ip + FirstP,jp + FirstP) +=  (-1.0) * weight * Integrating;
            }
        }
        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(porosity)/dPalpha
        for(int ip=0; ip < phrP; ip++)
        {
            for (int jp=0; jp < phrP; jp++)
            {
                REAL Integrating = phiP(ip,0) * rockporosity * (dwaterdensitydp * (SaturationAtnplusOne) + doildensitydp * (1 - SaturationAtnplusOne)) * phiP(jp,0);
                ek(ip + FirstP,jp + FirstP) -=  weight * Integrating;
            }
        }
        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(S1)/dSalpha
        for(int ip=0; ip<phrP; ip++)
        {
            for (int jsat=0; jsat<phrS; jsat++)
            {
                REAL Integrating = phiP(ip,0) * PhiStar * (waterdensity - oildensity) * phiS(jsat,0);
                ek(ip + FirstP,jsat + FirstS) -=  weight * Integrating;
            }
            
        }        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        // d(S1)/dSalpha
        for(int ip=0; ip<phrP; ip++)
        {
            for (int jsat=0; jsat<phrS; jsat++)
            {
                REAL Integrating = phiP(ip,0) * rockporosity * (waterdensity - oildensity) * phiS(jsat,0);
                ek(ip + FirstP,jsat + FirstS) -=  weight * Integrating;
            }
            
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[dot(grad(w),v), Omega_{e}] (Equation Two)
        for(int ip=0; ip<phrP; ip++)
        {
            //  Compute grad(W)
            TPZManVector<STATE> dsolp(2,0);
            dsolp[0] = dphiP(0,ip)*datavec[2].axes(0,0)+dphiP(1,ip)*datavec[2].axes(1,0);
            dsolp[1] = dphiP(0,ip)*datavec[2].axes(0,1)+dphiP(1,ip)*datavec[2].axes(1,1);
            
            for (int jq=0; jq<phrQ; jq++)
            {
                
                int jvectorindex    = datavec[1].fVecShapeIndex[jq].first;
                int jshapeindex     = datavec[1].fVecShapeIndex[jq].second;
                
                REAL dotprod =
                (dsolp[0]) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(0,jvectorindex)) +
                (dsolp[1]) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstP,jq + FirstQ) += (Gamma) * (TimeStep) * weight * dotprod;
            }
            
        }
        
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            
            for (int jsat=0; jsat<phrS; jsat++)
            {
                ek(isat + FirstS,jsat+ FirstS) +=  (PhiStar * waterdensity) * weight * phiS(isat,0) * phiS(jsat,0);
            }
            
        }        
        
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            
            for (int jsat=0; jsat<phrS; jsat++)
            {
                ek(isat + FirstS,jsat+ FirstS) +=  (rockporosity * waterdensity) * weight * phiS(isat,0) * phiS(jsat,0);
            }
            
        }
        
        //  Poroelastic contribution        
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            for (int ju=0; ju<phrU; ju++)
            {
                TPZFMatrix<REAL>    du(2,2);          

                du(0,0) = dphiU(0,ju)*datavec[0].axes(0,0); // dPhiU/dx dalphax
                du(1,0) = dphiU(0,ju)*datavec[0].axes(0,1); // dPhiU/dy dalphax

                du(0,1) = dphiU(1,ju)*datavec[0].axes(1,0); // dPhiU/dx dalphay
                du(1,1) = dphiU(1,ju)*datavec[0].axes(1,1); // dPhiU/dy dalphay        

                divphiU = du(0,0)+du(1,0)+0.0;

                REAL dPhiStardUx = Balpha * (du(0,0)+du(1,0));
                REAL dPhiStardUy = Balpha * (du(0,1)+du(1,1));

                ek(isat + FirstS,2*ju + FirstU)     += weight * (dPhiStardUx * waterdensity) * SaturationAtnplusOne * phiS(isat,0);
                ek(isat + FirstS,2*ju + 1 + FirstU) += weight * (dPhiStardUy * waterdensity) * SaturationAtnplusOne * phiS(isat,0);                
            }
        }           
        
        //  Poroelastic contribution
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            for (int jp=0; jp<phrP; jp++)
            {
                REAL dPhiStardP =  Sestr * phiP(jp,0);                   
                ek(isat + FirstS,jp + FirstP) += weight * (dPhiStardP * waterdensity + PhiStar * dwaterdensitydp * phiP(jp,0)) * SaturationAtnplusOne * phiS(isat,0);
            }
        }

        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            
            for (int jp=0; jp<phrP; jp++)
            {
                ek(isat + FirstS,jp + FirstP) += weight * (drockporositydp * waterdensity + rockporosity * dwaterdensitydp) * SaturationAtnplusOne * phiS(isat,0) * phiP(jp,0);
            }
            
        }
        
        //  Third Vector Block (Equation three)
        // Integrate[dot(d(f1(S1)/ds)q,grad(L)), Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            
            //  Compute grad(L)
            TPZManVector<STATE> Gradphis(2,0);
            Gradphis[0] = dphiS(0,isat)*datavec[3].axes(0,0)+dphiS(1,isat)*datavec[3].axes(1,0);
            Gradphis[1] = dphiS(0,isat)*datavec[3].axes(0,1)+dphiS(1,isat)*datavec[3].axes(1,1);
            
            
            REAL dotprod =
            (Gradphis[0]) * (sol_q[0]) +
            (Gradphis[1]) * (sol_q[1]);
            
            for (int jsat=0; jsat<phrS; jsat++)
            {
                ek(isat + FirstS,jsat + FirstS) -= (Theta) * (TimeStep) * weight * dbulkfwaterds * phiS(jsat,0) * dotprod;   
            }
        }
        
        //  Third Vector Block (Equation three)
        // Integrate[dot(f1(S1)v,grad(L)), Omega_{e}]   (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            //  Compute grad(L)
            TPZManVector<STATE> Gradphis(2,0);
            Gradphis[0] = dphiS(0,isat)*datavec[3].axes(0,0)+dphiS(1,isat)*datavec[3].axes(1,0);
            Gradphis[1] = dphiS(0,isat)*datavec[3].axes(0,1)+dphiS(1,isat)*datavec[3].axes(1,1);
            
            for (int jq=0; jq<phrQ; jq++)
            {
                
                int jvectorindex    = datavec[1].fVecShapeIndex[jq].first;
                int jshapeindex     = datavec[1].fVecShapeIndex[jq].second;
                REAL dotprod =
                (Gradphis[0]) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(0,jvectorindex)) +
                (Gradphis[1]) * (phiQ(jshapeindex,0)*datavec[1].fDeformedDirections(1,jvectorindex));
                
                ek(isat + FirstS,jq + FirstQ) -= (Theta) * (TimeStep) * weight * bulkfwater * dotprod;
                
            }
            
        }
        
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
        //  End of contribution of domain integrals for Jacobian matrix
        
    }
    
    
    //  n time step
    //  This values are constant in Newton iteration
    if(gState == ELastState)
    {
        
    }
    
    
    
    // Compute Residual contribution ef for domain integrals
    this->Contribute(datavec, weight, ef);
    
    
}

//  Residual vector contribution
void TPZMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    // Getting weight functions
    TPZFMatrix<REAL>  &phiU =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiQ =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiP =  datavec[2].phi;
    TPZFMatrix<REAL>  &phiS =  datavec[3].phi;
//    TPZFMatrix<REAL>  &phiQG    =  datavec[4].phi;
    
    TPZFMatrix<REAL> &dphiU = datavec[0].dphix;
    //    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[2].dphix;
    TPZFMatrix<REAL> &dphiS = datavec[3].dphix;// This a null value since S is constant by element.
    
    
    // number of test functions for each state variable
    int phrU, phrQ, phrP, phrS, phrQG;
    phrU = phiU.Rows();
    phrQ = datavec[1].fVecShapeIndex.NElements(); //    phiQ.Rows();
    phrP = phiP.Rows();
    phrS = phiS.Rows();
    phrQG   = datavec[4].fVecShapeIndex.NElements();
    
    
    // blocks
    int UStateVar = 2;
    int FirstU  = 0;
    int FirstQ  = phrU * UStateVar + FirstU;
    int FirstP  = phrQ + FirstQ;
    int FirstS  = phrP + FirstP;
    
    //  Getting and computing another required data
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    REAL Gamma = this->fGamma;
    int GeoID = datavec[1].gelElId;
    TPZFMatrix<REAL> Kabsolute;
    TPZFMatrix<REAL> Kinverse;
    TPZFMatrix<REAL> Gfield;
    if (fYorN)
    {
        Kabsolute=this->fKabsoluteMap[GeoID];
    }
    else
    {
        this->K(Kabsolute);
    }
    
    Kinverse=this->Kinv(Kabsolute);
    Gfield = this->Gravity();
    TPZManVector<STATE,3> sol_u =datavec[0].sol[0];
    TPZManVector<STATE,3> sol_q =datavec[1].sol[0];
    TPZManVector<STATE,3> sol_p =datavec[2].sol[0];
    TPZManVector<STATE,3> sol_s =datavec[3].sol[0];
    TPZManVector<STATE,3> sol_qg =datavec[4].sol[0];
    
    TPZFMatrix<STATE> dsol_u = datavec[0].dsol[0];
    TPZFMatrix<STATE> dsol_q =datavec[1].dsol[0];
    TPZFMatrix<STATE> dsol_p =datavec[2].dsol[0];
    TPZFMatrix<STATE> dsol_s =datavec[3].dsol[0];
    
    TPZFMatrix<> axesQ, axesP;
    axesQ=datavec[1].axes;
    
    REAL Pressure = sol_p[0];
    
    REAL LambdaL, LambdaLU, MuL;
    REAL Balpha, Sestr;    
    
    REAL rockporosity, oildensity, waterdensity;
    REAL drockporositydp, doildensitydp, dwaterdensitydp;
    
    REAL oilviscosity, waterviscosity;
    REAL doilviscositydp, dwaterviscositydp;
    
    REAL bulklambda, oillambda, waterlambda;
    REAL dbulklambdadp, doillambdadp, dwaterlambdadp;
    REAL dbulklambdads, doillambdads, dwaterlambdads;
    
    REAL bulkfoil, bulkfwater;
    REAL dbulkfoildp, dbulkfwaterdp;
    REAL dbulkfoilds, dbulkfwaterds;
    
    // Functions computed at point x_{k} for each integration point
    LambdaL     = this->LameLambda();
    LambdaLU    = this->LameLambdaU();
    MuL         = this->LameMu();
    Balpha      = this->BiotAlpha();
    Sestr       = this->Se();
    
    int VecPos= 0;
    //    REAL PressureRef = 1.0e6;
    this->Porosity(sol_p[VecPos], rockporosity, drockporositydp);
    this->RhoOil(sol_p[VecPos], oildensity, doildensitydp);
    this->RhoWater(sol_p[VecPos], waterdensity, dwaterdensitydp);
    this->OilViscosity(sol_p[VecPos], oilviscosity, doilviscositydp);
    this->WaterViscosity(sol_p[VecPos], waterviscosity, dwaterviscositydp);
    this->OilLabmda(oillambda, sol_p[VecPos], sol_s[VecPos], doillambdadp, doillambdads);
    this->WaterLabmda(waterlambda, sol_p[VecPos], sol_s[VecPos], dwaterlambdadp, dwaterlambdads);
    this->Labmda(bulklambda, sol_p[VecPos], sol_s[VecPos], dbulklambdadp, dbulklambdads);
    this->fOil(bulkfoil, sol_p[VecPos], sol_s[VecPos], dbulkfoildp, dbulkfoilds);
    this->fWater(bulkfwater, sol_p[VecPos], sol_s[VecPos], dbulkfwaterdp, dbulkfwaterds);
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of domain integrals for Residual Vector
    
    // n+1 time step
    if(gState == ECurrentState)
    {
        //  Elastic equation
        //  Linear strain operator
        //  Ke Matrix
        TPZFMatrix<REAL>    du(2,2);        
        TPZFMatrix<REAL>    dux(2,2);
        TPZFMatrix<REAL>    duy(2,2);        
        // Required check out of this implementation
        //  Derivative for Ux
        dux(0,1) = dsol_u(0,0)*datavec[0].axes(0,0)+dsol_u(1,0)*datavec[0].axes(1,0); // dUx/dx
        dux(1,1) = dsol_u(0,0)*datavec[0].axes(0,1)+dsol_u(1,0)*datavec[0].axes(1,1); // dUx/dy
        
        //  Derivative for Uy 
        duy(0,1) = dsol_u(0,1)*datavec[0].axes(0,0)+dsol_u(1,1)*datavec[0].axes(1,0); // dUy/dx
        duy(1,1) = dsol_u(0,1)*datavec[0].axes(0,1)+dsol_u(1,1)*datavec[0].axes(1,1); // dUy/dy
        
        for(int iu = 0; iu < phrU; iu++ )
        {
            //  Derivative for Vx
            du(0,0) = dphiU(0,iu)*datavec[0].axes(0,0)+dphiU(1,iu)*datavec[0].axes(1,0);
            //  Derivative for Vy
            du(1,0) = dphiU(0,iu)*datavec[0].axes(0,1)+dphiU(1,iu)*datavec[0].axes(1,1);
            
            //  Fu Vector Force right hand term  check the gravity term
//             ef(2*iu + FirstU)     +=    weight*Gfield(0,0)*phiU(iu, 0);
//             ef(2*iu+1 + FirstU)   +=    weight*Gfield(1,0)*phiU(iu, 0);
            
                if (fPlaneStress == 1)
                {
                    /* Plain stress state */
                    ef(2*iu + FirstU)           += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(0,0)*dux(0,1)      + (2*MuL)*du(1,0)*dux(1,1));
                    
                    ef(2*iu + FirstU)           += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(0,0)*duy(1,1)         + (2*MuL)*du(1,0)*duy(0,1));
                    
                    ef(2*iu+1 + FirstU)         += weight*((2*(MuL)*(LambdaL)/(LambdaL+2*MuL))*du(1,0)*dux(0,1)         + (2*MuL)*du(0,0)*dux(1,1));
                    
                    ef(2*iu+1 + FirstU)         += weight*((4*(MuL)*(LambdaL+MuL)/(LambdaL+2*MuL))*du(1,0)*duy(1,1)     + (2*MuL)*du(0,0)*duy(0,1));
                }
                else
                {
                    
//                     /* Plain Strain State */
//                     ek(2*iu + FirstU,2*ju + FirstU)         += weight*  ((LambdaL + 2*MuL)*du(0,0)*du(0,1)  + (MuL)*du(1,0)*du(1,1));
//                     
//                     ek(2*iu + FirstU,2*ju+1 + FirstU)       += weight*  (LambdaL*du(0,0)*du(1,1)            + (MuL)*du(1,0)*du(0,1));
//                     
//                     ek(2*iu+1 + FirstU,2*ju + FirstU)       += weight*  (LambdaL*du(1,0)*du(0,1)            + (MuL)*du(0,0)*du(1,1));
//                     
//                     ek(2*iu+1 + FirstU,2*ju+1 + FirstU)     += weight*  ((LambdaL + 2*MuL)*du(1,0)*du(1,1)  + (MuL)*du(0,0)*du(0,1));
                    
                    /* Plain Strain State */
                    ef(2*iu + FirstU)           += weight*  ((LambdaL + 2*MuL)*du(0,0)*dux(0,1)  + (MuL)*du(1,0)*(dux(1,1)));
                    
                    ef(2*iu + FirstU)           += weight*  (LambdaL*du(0,0)*duy(1,1)            + (MuL)*du(1,0)*(duy(0,1)));
                     
                    ef(2*iu+1 + FirstU)         += weight*  (LambdaL*du(1,0)*dux(0,1)            + (MuL)*du(0,0)*(dux(1,1)));
                    
                    ef(2*iu+1 + FirstU)         += weight*  ((LambdaL + 2*MuL)*du(1,0)*duy(1,1)  + (MuL)*du(0,0)*(duy(0,1)));
                }
        }
        
        //  Matrix Qc
        //  Coupling matrix for Elastic equation
        for(int in = 0; in < phrU; in++ )
        {
            du(0,0) = dphiU(0,in)*datavec[0].axes(0,0)+dphiU(1,in)*datavec[0].axes(1,0);
            du(1,0) = dphiU(0,in)*datavec[0].axes(0,1)+dphiU(1,in)*datavec[0].axes(1,1);
            
            ef(2*in + FirstU)    += (-1.0)*Balpha*weight*(Pressure*du(0,0));
            ef(2*in+1 + FirstU)  += (-1.0)*Balpha*weight*(Pressure*du(1,0));
        }        
        
        //  This block was verified
        REAL SaturationAtnplusOne = sol_s[0];
        //  First Block (Equation One) constitutive law
        //  Integrate[(viscosity/density)*dot(q,v), Omega_{e}]   (Equation One)
        //  REAL ViscOverdensity = waterviscosity/waterdensity;
        REAL OneOverLambda = 1.0/bulklambda;
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                REAL dotprod =
                (Kinverse(0,0)*sol_q[0]+Kinverse(0,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                (Kinverse(1,0)*sol_q[0]+Kinverse(1,1)*sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) ;
                
                ef(iq + FirstQ) +=  OneOverLambda * weight * dotprod;
            }
            else
            {
                REAL dotprod =
                (sol_q[0]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex)) +
                (sol_q[1]) * (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex)) ;
                
                ef(iq + FirstQ) +=  OneOverLambda * weight * dotprod;
            }

        }
        
        
        //  This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K dot(v,grad(P)) , Omega_{e}]   (Equation One)
        //  Compute grad(P)
        TPZManVector<STATE> dsolp(2,0);
        dsolp[0] = dsol_p(0,0)*datavec[2].axes(0,0)+dsol_p(1,0)*datavec[2].axes(1,0);
        dsolp[1] = dsol_p(0,0)*datavec[2].axes(0,1)+dsol_p(1,0)*datavec[2].axes(1,1);
        
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                REAL e1e1   =   (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex))*(dsolp[0]);
                REAL e2e2   =   (phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex))*(dsolp[1]);
                
                ef(iq + FirstQ) += weight * (e1e1 + e2e2);
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(dsolp[0])+
                                 Kabsolute(0,1)*(dsolp[1]))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                
                REAL e2e2   =   (Kabsolute(1,0)*(dsolp[0])+
                                 Kabsolute(1,1)*(dsolp[1]))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));

                ef(iq + FirstQ) += weight * (e1e1 + e2e2);
            }
        }
        
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate [ K (dot(f1*rho1*grad(g*z),v)+dot(f1*rho1*grad(g*z),v)) , Omega_{e}]  (Equation One)
        
        for(int iq=0; iq<phrQ; iq++)
        {
            
            int ivectorindex        = datavec[1].fVecShapeIndex[iq].first;
            int ishapeindex         = datavec[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (Gfield(0,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                REAL e2e2   =   (Gfield(1,0))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                ef(iq + FirstQ) -= weight * ( ( bulkfwater * waterdensity + bulkfoil * oildensity )* (e1e1 + e2e2));
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(Gfield(0,0))+
                                 Kabsolute(0,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(0,ivectorindex));
                
                REAL e2e2   =   (Kabsolute(1,0)*(Gfield(0,0))+
                                 Kabsolute(1,1)*(Gfield(1,0)))*(phiQ(ishapeindex,0)*datavec[1].fDeformedDirections(1,ivectorindex));
                
                ef(iq + FirstQ) -= weight * ((bulkfwater * waterdensity + bulkfoil * oildensity)* (e1e1 + e2e2));
            }
            
        }
        
        //  Poroelastic Contribution
        REAL divU, dsolU[2][2];
        dsolU[0][0] = dsol_u(0,0)*datavec[0].axes(0,0)+dsol_u(1,0)*datavec[0].axes(1,0); // dUx/dx
        dsolU[1][0] = dsol_u(0,0)*datavec[0].axes(0,1)+dsol_u(1,0)*datavec[0].axes(1,1); // dUx/dy

        dsolU[0][1] = dsol_u(0,1)*datavec[0].axes(0,0)+dsol_u(1,1)*datavec[0].axes(1,0); // dUy/dx
        dsolU[1][1] = dsol_u(0,1)*datavec[0].axes(0,1)+dsol_u(1,1)*datavec[0].axes(1,1); // dUy/dy

        divU = dsolU[0][0]+dsolU[1][1]+0.0;
        REAL PhiStar = Balpha * divU + Sestr * Pressure;
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        for(int ip=0; ip<phrP; ip++)
        {
            // Here rockporosity is the poroelastic contribution
            REAL Integrating = phiP(ip,0) * PhiStar * (waterdensity * SaturationAtnplusOne + oildensity * (1 - SaturationAtnplusOne));
            ef(ip + FirstP) += (-1.0) * weight * Integrating;
        }
       
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        for(int ip=0; ip<phrP; ip++)
        {
            // Here rockporosity is the poroelastic contribution
            REAL Integrating = phiP(ip,0) * rockporosity * (waterdensity * SaturationAtnplusOne + oildensity * (1 - SaturationAtnplusOne));
            ef(ip + FirstP) += (-1.0) * weight * Integrating;
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[dot(grad(W),q), Omega_{e}] (Equation Two)
        //      This block was verified
        for(int ip=0; ip<phrP; ip++)
        {
            //  Compute grad(W)
            TPZManVector<STATE> dsolp(2,0);
            dsolp[0] = dphiP(0,ip)*datavec[2].axes(0,0)+dphiP(1,ip)*datavec[2].axes(1,0);
            dsolp[1] = dphiP(0,ip)*datavec[2].axes(0,1)+dphiP(1,ip)*datavec[2].axes(1,1);
            
            REAL dotprod =
            (dsolp[0]) * (sol_q[0]) +
            (dsolp[1]) * (sol_q[1]);
            
            ef(ip + FirstP) += (Gamma) * (TimeStep) * weight * dotprod;
        }
        
        //  Poroelastic Contribution
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            // Here rockporosity is the poroelastic contribution            
            ef(isat + FirstS) += weight * (PhiStar * waterdensity) * phiS(isat,0) * SaturationAtnplusOne;
        }        
        
        
        //  This block was verified
        //  Third Vector Block (Equation three)
        //  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            // Here rockporosity is the poroelastic contribution            
            ef(isat + FirstS) += weight * (rockporosity * waterdensity) * phiS(isat,0) * SaturationAtnplusOne;
        }
        
        
        //  Third Vector Block (Equation three)
        //  Integrate[dot(f1(S1)q,grad(L)), Omega_{e}]   (Equation three)
        //  std::cout << "phrS:   " << phrS << "   \n" << std::endl;
        for(int isat=0; isat<phrS; isat++)
        {
            //  Compute grad(L)
            TPZManVector<STATE> Gradphis(2,0);
            Gradphis[0] = dphiS(0,isat)*datavec[3].axes(0,0)+dphiS(1,isat)*datavec[3].axes(1,0);
            Gradphis[1] = dphiS(0,isat)*datavec[3].axes(0,1)+dphiS(1,isat)*datavec[3].axes(1,1);

            REAL dotprod =
            (Gradphis[0]) * (sol_q[0]) +
            (Gradphis[1]) * (sol_q[1]);

            ef(isat + FirstS) -= (Theta) * (TimeStep) * weight * bulkfwater * dotprod;
        }
        
    }
    
    
    // n time step
    if(gState == ELastState)
    {
        REAL SaturationAtnTimeStep = sol_s[0]; //   Gettin Saturation at n time step
        
        //  Poroelastic contribution
        REAL divU, dsolU[2][2];
        dsolU[0][0] = dsol_u(0,0)*datavec[0].axes(0,0)+dsol_u(1,0)*datavec[0].axes(1,0); // dUx/dx
        dsolU[1][0] = dsol_u(0,0)*datavec[0].axes(0,1)+dsol_u(1,0)*datavec[0].axes(1,1); // dUx/dy

        dsolU[0][1] = dsol_u(0,1)*datavec[0].axes(0,0)+dsol_u(1,1)*datavec[0].axes(1,0); // dUy/dx
        dsolU[1][1] = dsol_u(0,1)*datavec[0].axes(0,1)+dsol_u(1,1)*datavec[0].axes(1,1); // dUy/dy

        divU = dsolU[0][0]+dsolU[1][1]+0.0;
        REAL PhiStar = Balpha * divU + Sestr * Pressure;
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        for(int ip = 0; ip < phrP; ip++)
        {
            // Here rockporosity is the poroelastic contribution            
            REAL Integrating = phiP(ip,0) * PhiStar * (waterdensity * SaturationAtnTimeStep + oildensity * (1 - SaturationAtnTimeStep));
            ef(ip + FirstP) += weight * Integrating;
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[W*(d(\phi*(rho1 * S1 + rho2 * S2)/dt)), Omega_{e}] (Equation Two)
        for(int ip = 0; ip < phrP; ip++)
        {
            // Here rockporosity is the poroelastic contribution            
            REAL Integrating = phiP(ip,0) * rockporosity * (waterdensity * SaturationAtnTimeStep + oildensity * (1 - SaturationAtnTimeStep));
            ef(ip + FirstP) += weight * Integrating;
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[dot(grad(w),q), Omega_{e}] (Equation Two)
        //      This block was verified
        for(int ip=0; ip<phrP; ip++)
        {
            //  Compute grad(W)
            TPZManVector<STATE> dsolp(2,0);
            dsolp[0] = dphiP(0,ip)*datavec[2].axes(0,0)+dphiP(1,ip)*datavec[2].axes(1,0);
            dsolp[1] = dphiP(0,ip)*datavec[2].axes(0,1)+dphiP(1,ip)*datavec[2].axes(1,1);
            
            REAL dotprod =
            (dsolp[0]) * (sol_q[0]) +
            (dsolp[1]) * (sol_q[1]);
            
            ef(ip + FirstP) += (1-Gamma) * (TimeStep) * weight * dotprod;
        }
        
        
        //  Poroelastic Contribution
        //  Third Vector Block (Equation three)
        // Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            // Here rockporosity is the poroelastic contribution            
            ef(isat + FirstS) += (-1.0) * weight * (PhiStar * waterdensity) * phiS(isat,0) * SaturationAtnTimeStep;
        }
        
        //  Third Vector Block (Equation three)
        //  (-1.0) * Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            // Here rockporosity is the poroelastic contribution            
            ef(isat + FirstS) += (-1.0) * weight * (rockporosity * waterdensity) * phiS(isat,0) * SaturationAtnTimeStep;
        }
        
        //  Third Vector Block (Equation three)
        //  Integrate[dot(f1(S1)q,grad(L)), Omega_{e}]   (Equation three)
        for(int isat=0; isat<phrS; isat++)
        {
            //  Compute grad(L)
            TPZManVector<STATE> Gradphis(2,0);
            Gradphis[0] = dphiS(0,isat)*datavec[3].axes(0,0)+dphiS(1,isat)*datavec[3].axes(1,0);
            Gradphis[1] = dphiS(0,isat)*datavec[3].axes(0,1)+dphiS(1,isat)*datavec[3].axes(1,1);
            
            REAL dotprod =
            (Gradphis[0]) * (sol_q[0]) +
            (Gradphis[1]) * (sol_q[1]) ;
            
            ef(isat + FirstS) -= (1-Theta) * (TimeStep) * weight * bulkfwater * dotprod;
        }
    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  End of contribution of domain integrals for Residual Vector
    
    
}




void TPZMultiphase::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
    DebugStop();
}

void TPZMultiphase::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
#ifdef PZDEBUG
    for(int i=0; i<5; i++)
    {
        if(dataleft.find(i) == dataleft.end()) DebugStop();
        if(dataright.find(i) == dataright.end()) DebugStop();
    }
#endif
    TPZFMatrix<REAL> &phiUL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiUR = dataright[0].phi;
    
    TPZFMatrix<REAL> &phiQL = dataleft[1].phi;
    TPZFMatrix<REAL> &phiQR = dataright[1].phi;
    
    TPZFMatrix<REAL> &phiPL = dataleft[2].phi;
    TPZFMatrix<REAL> &phiPR = dataright[2].phi;
    
    TPZFMatrix<REAL> &phiSL = dataleft[3].phi;
    TPZFMatrix<REAL> &phiSR = dataright[3].phi;
    
    TPZFMatrix<REAL> &phiQGL = dataleft[4].phi;
//    TPZFMatrix<REAL> &phiQGR = dataright[4].phi;
    
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    //  REAL n3 = normal[2];

    TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
    TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
    TPZManVector<STATE,3> sol_qgL =dataleft[4].sol[0];

    TPZManVector<STATE,3> sol_uR =dataright[0].sol[0];
    TPZManVector<STATE,3> sol_qR =dataright[1].sol[0];
    TPZManVector<STATE,3> sol_pR =dataright[2].sol[0];
    TPZManVector<STATE,3> sol_sR =dataright[3].sol[0];
    TPZManVector<STATE,3> sol_qgR =dataright[4].sol[0];
    
    //  Getting Q solution for left and right side
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL qxR = sol_qR[0];
    REAL qyR = sol_qR[1];
    REAL dotqnL = (qxL*n1) + (qyL*n2);
    REAL dotqnR = (qxR*n1) + (qyR*n2);
    
    //  Getting QG solution for left and right side
//    REAL qgxL = sol_qgL[0];
//    REAL qgyL = sol_qgL[1];
//    REAL qgxR = sol_qgR[0];
//    REAL qgyR = sol_qgR[1];
//    REAL dotqgnL = (qgxL*n1) + (qgyL*n2);
//    REAL dotqgnR = (qgxR*n1) + (qgyR*n2);    
    
    //  Getting S solution for left and right side
    REAL SaturationL    =   sol_sL[0];
    REAL SaturationR    =   sol_sR[0];
    
    //  Getting another required data
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    REAL Gamma = this->fGamma;
    this->fnewWS=true;
    
    // Getting Harmonic mean of permeabilities
    TPZFMatrix<REAL> KabsoluteLeft;
    TPZFMatrix<REAL> KabsoluteRight;
    TPZFMatrix<REAL> Gfield;
    TPZFMatrix<REAL> Kmean(3,3,0);
    
    int GeoIDLeft = dataleft[1].gelElId;
    int GeoIDRight = dataright[1].gelElId;
    
    REAL rockporosityl, oildensityl, waterdensityl;
    REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;
    REAL rockporosityr, oildensityr, waterdensityr;
    REAL drockporositydpr, doildensitydpr, dwaterdensitydpr;
    
    REAL oilviscosityl, waterviscosityl;
    REAL doilviscositydpl, dwaterviscositydpl;
    REAL oilviscosityr, waterviscosityr;
    REAL doilviscositydpr, dwaterviscositydpr;
    
    REAL bulklambdal, oillambdal, waterlambdal;
    REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
    REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;
    REAL bulklambdar, oillambdar, waterlambdar;
    REAL dbulklambdadpr, doillambdadpr, dwaterlambdadpr;
    REAL dbulklambdadsr, doillambdadsr, dwaterlambdadsr;
    
    REAL bulkfoill, bulkfwaterl;
    REAL dbulkfoildpl, dbulkfwaterdpl;
    REAL dbulkfoildsl, dbulkfwaterdsl;
    
    REAL bulkfoilr, bulkfwaterr;
    REAL dbulkfoildpr, dbulkfwaterdpr;
    REAL dbulkfoildsr, dbulkfwaterdsr;
    
    REAL bulkfStarl;
    REAL dbulkfStardpl;
    REAL dbulkfStardsl;
    
    REAL bulkfStarr;
    REAL dbulkfStardpr;
    REAL dbulkfStardsr;     
    
    REAL Pcl;
    REAL dPcdSl;
    
    REAL Pcr;
    REAL dPcdSr;    
    
    // Functions computed at point x_{k} for each integration point
    int VecPos= 0;
    //  REAL PressureRef = 1.0e6;
    this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
    this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
    this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
    this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
    this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
    this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
    this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
    this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
    this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
    this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);
    this->CapillaryPressure(sol_sL[VecPos],Pcl,dPcdSl);
    
    this->Porosity(sol_pR[VecPos], rockporosityr, drockporositydpr);
    this->RhoOil(sol_pR[VecPos], oildensityr, doildensitydpr);
    this->RhoWater(sol_pR[VecPos], waterdensityr, dwaterdensitydpr);
    this->OilViscosity(sol_pR[VecPos], oilviscosityr, doilviscositydpr);
    this->WaterViscosity(sol_pR[VecPos], waterviscosityr, dwaterviscositydpr);
    this->OilLabmda(oillambdar, sol_pR[VecPos], sol_sR[VecPos], doillambdadpr, doillambdadsr);
    this->WaterLabmda(waterlambdar, sol_pR[VecPos], sol_sR[VecPos], dwaterlambdadpr, dwaterlambdadsr);
    this->Labmda(bulklambdar, sol_pR[VecPos], sol_sR[VecPos], dbulklambdadpr, dbulklambdadsr);
    this->fOil(bulkfoilr, sol_pR[VecPos], sol_sR[VecPos], dbulkfoildpr, dbulkfoildsr);
    this->fWater(bulkfwaterr, sol_pR[VecPos], sol_sR[VecPos], dbulkfwaterdpr, dbulkfwaterdsr);
    this->CapillaryPressure(sol_sR[VecPos],Pcr,dPcdSr);    
    
    if (fYorN)
    {
        
        KabsoluteLeft=fKabsoluteMap[GeoIDLeft];
        KabsoluteRight=fKabsoluteMap[GeoIDRight];
        
    }
    else
    {
        this->K(KabsoluteLeft);
        this->K(KabsoluteRight);
    }
    
    Gfield = this->Gravity();
    REAL Gravitydotnl =  Gfield(0,0)*(n1) + Gfield(1,0)*(n2);
    REAL Gravitydotnr =  Gfield(0,0)*(n1) + Gfield(1,0)*(n2);    
    this->fstar(bulkfStarl,sol_pL[VecPos], sol_sL[VecPos],1.0*Gravitydotnl,dbulkfStardpl,dbulkfStardsl);
    this->fstar(bulkfStarr,sol_pR[VecPos], sol_sR[VecPos],-1.0*Gravitydotnr,dbulkfStardpr,dbulkfStardsr);
    
    // Min of fstar at Gamma
    REAL bulkfstar = 0.0;
    REAL dbulkfstardp = 0.0;
    REAL dbulkfstards = 0.0;    
    if(fabs(bulkfStarl) >= fabs(bulkfStarr))
    {
        bulkfstar = bulkfStarr;
        dbulkfstardp = dbulkfStardpr;
        dbulkfstards = dbulkfStardsr;
    }
    else
    {
        bulkfstar = bulkfStarl;
        dbulkfstardp = dbulkfStardpl;
        dbulkfstards = dbulkfStardsl;        
    }      
    
    
    for (int in=0; in < 3; in++) {
        for (int jn=0; jn < 3; jn++) {
            if (KabsoluteLeft(in,jn)+KabsoluteRight(in,jn)==0) {
                Kmean(in,jn)= 0.0;
            }
            else
            {
                Kmean(in,jn)= (2.0*KabsoluteLeft(in,jn)*KabsoluteRight(in,jn))/(KabsoluteLeft(in,jn)+KabsoluteRight(in,jn));
            }
            
        }
    }
    
    TPZFMatrix<STATE> KGL(3,3,0.0),KGR(3,3,0.0), KG(3,3,0.0);
    KGL(0,0) = KabsoluteLeft(0,0)*Gfield(0,0)+KabsoluteLeft(0,1)*Gfield(1,0);
    KGL(1,0) = KabsoluteLeft(1,0)*Gfield(0,0)+KabsoluteLeft(1,1)*Gfield(1,0);        
    KGR(0,0) = KabsoluteRight(0,0)*Gfield(0,0)+KabsoluteRight(0,1)*Gfield(1,0);
    KGR(1,0) = KabsoluteRight(1,0)*Gfield(0,0)+KabsoluteRight(1,1)*Gfield(1,0);
    
    KG(0,0) = Kmean(0,0)*Gfield(0,0)+Kmean(0,1)*Gfield(1,0);
    KG(1,0) = Kmean(1,0)*Gfield(0,0)+Kmean(1,1)*Gfield(1,0);         
    
//    REAL GravityFluxL   =   (KabsoluteLeft(0,0)*Gfield(0,0) + KabsoluteLeft(0,1)*Gfield(1,0))*(n1)+
//    (KabsoluteLeft(1,0)*Gfield(0,0) + KabsoluteLeft(1,1)*Gfield(1,0))*(n2);
    
//    REAL GravityFluxR   =   (KabsoluteRight(0,0)*Gfield(0,0) + KabsoluteRight(0,1)*Gfield(1,0))*(n1)+
//    (KabsoluteRight(1,0)*Gfield(0,0) + KabsoluteRight(1,1)*Gfield(1,0))*(n2);
    
    
    int URowsleft = phiUL.Rows();
    int URowsRight = phiUR.Rows();
    
    int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
    int QRowsRight = dataright[1].fVecShapeIndex.NElements();
    
    int PRowsleft = phiPL.Rows();
    int PRowsRight = phiPR.Rows();
    
    int SRowsleft = phiSL.Rows();
    int SRowsRight = phiSR.Rows();
    
    int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
//    int QGRowsRight = dataright[4].fVecShapeIndex.NElements();
    
    int UStateVar = 2;    
    
    int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
    int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
    
    int FirstUL = 0;
    int FirstQL = URowsleft * UStateVar + FirstUL;
    int FirstPL = QRowsleft + FirstQL;
    int FirstSL = PRowsleft + FirstPL;
    int FirstQGL = SRowsleft + FirstSL;
    
    int FirstUR = 0;
    int FirstQR = URowsRight * UStateVar + FirstUR;
    int FirstPR = QRowsRight + FirstQR;
    int FirstSR = PRowsRight + FirstPR;
//    int FirstQGR = SRowsRight + FirstSR;
    
    
    if(gState == ECurrentState)
    {
        
        //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
        //  Contribution of contour integrals for Jacobian matrix

        //  First Block (Equation One) constitutive law
        //  Integrate[L dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
        
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            
            for (int jp=0; jp < PRowsleft; jp++)
            {

                if (fnewWS)
                {
                    REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
                    REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
                    ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * (e1e1 + e2e2 ) * phiPL(jp,0);
                }
                else
                {
                    REAL e1e1   =   (Kmean(0,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                     Kmean(0,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n1);
                    
                    REAL e2e2   =   (Kmean(1,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                     Kmean(1,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n2);
                    ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * (e1e1 + e2e2 ) * phiPL(jp,0);
                }
                
            }
            
        }
        
        //  First Block (Equation One) constitutive law
        // Integrate[L dot(K v, n), Gamme_{e}]  (Equation One) Right-Right Part
        for (int iq=0; iq < QRowsRight; iq++)
        {
            int iRvectorindex       = dataright[1].fVecShapeIndex[iq].first;
            int iRshapeindex        = dataright[1].fVecShapeIndex[iq].second;
            
            for (int jp=0; jp < PRowsRight; jp++)
            {

                if (fnewWS)
                {
                    
                    REAL e1e1   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))*(n1);
                    REAL e2e2   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex))*(n2);
                    
                    ek(iq + FirstQR + iRightInterfaceBlock,jp + FirstPR + jRightInterfaceBlock) +=  (1.0) * weight * (e1e1 + e2e2 ) * phiPR(jp,0) ;
                }
                else
                {
                    REAL e1e1   =   (Kmean(0,0)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))+
                                     Kmean(0,1)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex)))*(n1);
                    
                    REAL e2e2   =   (Kmean(1,0)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))+
                                     Kmean(1,1)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex)))*(n2);
                    
                    ek(iq + FirstQR + iRightInterfaceBlock,jp + FirstPR + jRightInterfaceBlock) +=  (1.0) * weight * (e1e1 + e2e2 ) * phiPR(jp,0) ;
                }
                
                
                
            }
            
        }
        
        REAL dSwPcdSL = SaturationL * dPcdSl + Pcl;
        REAL dSwPcdSR = SaturationR * dPcdSr + Pcr;
        
        //This block was verified
        //  First Block (Equation One) constitutive law Capillary Pressure
        // Integrate[Sw Pc dot( v, n), Gamme_{e}]  (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            for (int jsat=0; jsat < SRowsleft; jsat++)
            {
		REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
		REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
		ek(iq + FirstQL, jsat + FirstSL) += (-1.0) * (-1.0) * weight * dSwPcdSL * (e1e1 + e2e2 ) * phiSL(jsat,0);
            }            
        }
        
        //This block was verified
        //  First Block (Equation One) constitutive law Capillary Pressure
        // Integrate[Sw Pc dot(v, n), Gamme_{e}]  (Equation One) Right-Right Part
        for (int iq=0; iq < QRowsRight; iq++)
        {
            int iRvectorindex       = dataright[1].fVecShapeIndex[iq].first;
            int iRshapeindex        = dataright[1].fVecShapeIndex[iq].second;
	    
            for (int jsat=0; jsat < SRowsleft; jsat++)
            {
		REAL e1e1   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))*(n1);
		REAL e2e2   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex))*(n2);
		ek(iq + FirstQR + iRightInterfaceBlock,jsat + FirstSR + jRightInterfaceBlock) += (-1.0) * (1.0) * weight * dSwPcdSR * (e1e1 + e2e2 ) * phiSR(jsat,0) ;
            }
            
        }        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            
            for (int jq=0; jq<QRowsleft; jq++)
            {
                
                int jvectorindex    = dataleft[1].fVecShapeIndex[jq].first;
                int jshapeindex     = dataleft[1].fVecShapeIndex[jq].second;
                
                REAL dotprod =
                (n1) * (phiQL(jshapeindex,0)*dataleft[1].fDeformedDirections(0,jvectorindex)) +
                (n2) * (phiQL(jshapeindex,0)*dataleft[1].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstPL,jq + FirstQL) += (-1.0) * (Gamma) * (TimeStep) * weight * dotprod * phiPL(ip,0);
                
            }
            
        }
        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamma_{e}]    (Equation Two) Right-Right Part
        for (int ip=0; ip < PRowsRight; ip++)
        {
            
            for (int jq=0; jq<QRowsRight; jq++)
            {
                
                int jvectorindex    = dataright[1].fVecShapeIndex[jq].first;
                int jshapeindex     = dataright[1].fVecShapeIndex[jq].second;
                
                REAL dotprod =
                (n1) * (phiQR(jshapeindex,0)*dataright[1].fDeformedDirections(0,jvectorindex)) +
                (n2) * (phiQR(jshapeindex,0)*dataright[1].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstPR + iRightInterfaceBlock,jq + FirstQR + jRightInterfaceBlock) += (1.0) * (Gamma) * (TimeStep) * weight * dotprod * phiPR(ip,0);
                
            }
            
        }
        
        
        
        //  Upwind scheme
        //  Third Vector Block (Equation three) Saturation  equation
        
        REAL UpwindSaturation = 0.0;
        
        if (dotqnL > 0.0)
        {
            UpwindSaturation = bulkfwaterl;
            
            //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Left-Left Part
            for(int isat=0; isat<SRowsleft; isat++)
            {
                for(int jsat=0; jsat<SRowsleft; jsat++)
                {
                    ek(isat + FirstSL ,jsat + FirstSL)
                    += weight * (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
                }
            }
            
            //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Right-Left Part
            for(int isat=0; isat<SRowsRight; isat++)
            {
                
                for(int jsat=0; jsat<SRowsleft; jsat++)
                {
                    ek(isat + FirstSR + iRightInterfaceBlock,jsat + FirstSL)
                    -= weight * (Theta) * (TimeStep) * phiSR(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
                }
            }
            
        }
        
        else
            
        {
            UpwindSaturation = bulkfwaterr;
            
            //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Left-Right Part
            for(int isat=0; isat<SRowsleft; isat++)
            {
                for(int jsat=0; jsat<SRowsRight; jsat++)
                {
                    ek(isat + FirstSL,jsat +  FirstSR + jRightInterfaceBlock)
                    += weight * (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsr * phiSR(jsat,0) * dotqnR;
                }
            }
            //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Right-Right Part
            for(int isat=0; isat<SRowsRight; isat++)
            {
                
                for(int jsat=0; jsat<SRowsRight; jsat++)
                {
                    ek(isat + FirstSR + iRightInterfaceBlock,jsat + FirstSR + jRightInterfaceBlock)
                    -= weight * (Theta) * (TimeStep) * phiSR(isat,0) * dbulkfwaterdsr * phiSR(jsat,0) * dotqnR;
                }
            }
            
        }
        
        
        //  Third Vector Block (Equation three) Saturation  equation
        //  Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]  (Equation three) Left-Left Part
        for (int isat=0; isat < SRowsleft; isat++) {
            
            for (int jq=0; jq < QRowsleft; jq++)
            {
                int jLvectorindex       = dataleft[1].fVecShapeIndex[jq].first;
                int jLshapeindex        = dataleft[1].fVecShapeIndex[jq].second;
                
                REAL dotprodL =
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(0,jLvectorindex)) * (n1) +
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(1,jLvectorindex)) * (n2) ;//+
                //              (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(2,jLvectorindex)) * (n3) ;    //  dot(q,n)    left
                
                ek(isat + FirstSL ,jq + FirstQL) += weight * (Theta) * (TimeStep) * phiSL(isat,0) * UpwindSaturation * dotprodL;
                
            }
        }
        
        //  Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]  (Equation three) Right-Left Part
        for (int isat=0; isat < SRowsRight; isat++)
        {
            for (int jq=0; jq < QRowsleft; jq++)
            {
                int jLvectorindex       = dataleft[1].fVecShapeIndex[jq].first;
                int jLshapeindex        = dataleft[1].fVecShapeIndex[jq].second;
                
                REAL dotprodL =
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(0,jLvectorindex)) * (n1) +
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(1,jLvectorindex)) * (n2) ;//+
                //              (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(2,jLvectorindex)) * (n3) ;    //  dot(q,n)    left
                
                ek(isat+ FirstSR + iRightInterfaceBlock,jq + FirstQL)  -= weight * (Theta) * (TimeStep) * phiSR(isat,0) * UpwindSaturation * dotprodL;
                
            }
        }

        REAL QGgstarL = bulklambdal * (waterdensityl - oildensityl) * (KGL(0,0)*n1 + KGL(1,0)*n2);
        REAL QGgstarR = bulklambdar * (waterdensityr - oildensityr) * (KGR(0,0)*n1 + KGR(1,0)*n2);        
        
        // It is needed to implement the complete derivative of P
        REAL dQGgstarLdP = ((bulklambdal * (dwaterdensitydpl - doildensitydpl)) + 
                           (dbulklambdadpl * (waterdensityl - oildensityl))) * (KGL(0,0)*n1 + KGL(1,0)*n2);
                           
        REAL dQGgstarRdP = ((bulklambdar * (dwaterdensitydpr - doildensitydpr)) + 
                           (dbulklambdadpr * (waterdensityr - oildensityr))) * (KGR(0,0)*n1 + KGR(1,0)*n2);
        
        REAL dQGgstarLdS = (dbulklambdadsl) * (waterdensityl - oildensityl) * (KGL(0,0)*n1 + KGL(1,0)*n2);
        REAL dQGgstarRdS = (dbulklambdadsr) * (waterdensityr - oildensityr) * (KGR(0,0)*n1 + KGR(1,0)*n2);      
        
        REAL QGstar = 0.0;
        REAL dQGstardP = 0.0;
        REAL dQGstardS = 0.0;        
        
//         if(Gravitydotnl > 0.0 )
//         {
        if(fabs(1.0*bulkfStarl*QGgstarL) < fabs(1.0*bulkfStarr*QGgstarR))
        {           
            QGstar = 1.0 * bulkfStarl * QGgstarL;
            dQGstardS = 1.0*(dbulkfStardsl * QGgstarL + bulkfStarl * dQGgstarLdS);            
            dQGstardP = 1.0*(dbulkfStardpl * QGgstarL + bulkfStarl * dQGgstarLdP);
            
            // Gravitational segregation scheme
            //  Four Block (Equation Four) gravitational flux constitutive law
            // Integrate[L dot(K v, n), Gamma_{e}]  (Equation One) Left-Left part
            for (int iqg=0; iqg < QGRowsleft; iqg++)
            {
                int iLvectorindex       = dataleft[4].fVecShapeIndex[iqg].first;
                int iLshapeindex        = dataleft[4].fVecShapeIndex[iqg].second;
                REAL e1e1i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(0,iLvectorindex))*(n1);
                REAL e2e2i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(1,iLvectorindex))*(n2);
                
                for(int jsat=0; jsat<SRowsleft; jsat++)
                {            
                    // This degree of freedom is equal for Left and Right part
                    ek(iqg + FirstQGL, jsat + FirstSL) 
                    -= weight * (dQGstardS * phiSL(jsat,0)) * (e1e1i + e2e2i);
                }
                
                for(int jp=0; jp<PRowsleft; jp++)
                {
                    // This degree of freedom is equal for Left and Right part
                    ek(iqg + FirstQGL, jp + FirstPL) 
                    -= weight * (dQGstardP * phiPL(jp,0)) * (e1e1i + e2e2i);
                }                
                
            }                          
            
        }
        else
        {
            QGstar = 1.0 *bulkfStarr*QGgstarR;
            dQGstardS = 1.0*(dbulkfStardsr * QGgstarR + bulkfStarr * dQGgstarRdS);            
            dQGstardP = 1.0*(dbulkfStardpr * QGgstarR + bulkfStarr * dQGgstarRdP);
            
            // Gravitational segregation scheme
            //  Four Block (Equation Four) gravitational flux constitutive law
            // Integrate[L dot(K v, n), Gamma_{e}]  (Equation One) Left-Left part
            for (int iqg=0; iqg < QGRowsleft; iqg++)
            {
                int iLvectorindex       = dataleft[4].fVecShapeIndex[iqg].first;
                int iLshapeindex        = dataleft[4].fVecShapeIndex[iqg].second;
                REAL e1e1i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(0,iLvectorindex))*(n1);
                REAL e2e2i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(1,iLvectorindex))*(n2);
                
                for(int jsat=0; jsat<SRowsRight; jsat++)
                {            
                    // This degree of freedom is equal for Left and Right part
                    ek(iqg + FirstQGL, jRightInterfaceBlock + jsat + FirstSR) 
                   -= weight * (dQGstardS * phiSR(jsat,0)) * (e1e1i + e2e2i);
                }
                
                for(int jp=0; jp<PRowsRight; jp++)
                {            
                    // This degree of freedom is equal for Left and Right part
                    ek(iqg + FirstQGL, jRightInterfaceBlock + jp + FirstPR ) 
                    -= weight * (dQGstardP * phiPR(jp,0)) * (e1e1i + e2e2i);
                }                
                
            }            

        }
        
        
        // Gravitational segregation scheme
        //  Four Block (Equation Four) gravitational flux constitutive law
        // Integrate[L dot(K v, n), Gamma_{e}]  (Equation One) Left-Left part
        for (int iqg=0; iqg < QGRowsleft; iqg++)
        {
            int iLvectorindex       = dataleft[4].fVecShapeIndex[iqg].first;
            int iLshapeindex        = dataleft[4].fVecShapeIndex[iqg].second;
            REAL e1e1i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(0,iLvectorindex))*(n1);
            REAL e2e2i   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(1,iLvectorindex))*(n2);
            
            for (int jqg=0; jqg < QGRowsleft; jqg++)
            {
                int jLvectorindex       = dataleft[4].fVecShapeIndex[jqg].first;
                int jLshapeindex        = dataleft[4].fVecShapeIndex[jqg].second;
                REAL e1e1j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(0,jLvectorindex))*(n1);
                REAL e2e2j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(1,jLvectorindex))*(n2);
            
                // This degree of freedom is equal for Left and Right part
                ek(iqg + FirstQGL, jqg + FirstQGL) 
                += weight * (e1e1j + e2e2j) * (e1e1i + e2e2i);
            }
            
        }
        
     
        // Gravitational segregation scheme
        // (Theta) * deltat * Integrate[L*dot(qg,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat < SRowsleft; isat++)
        {
            for (int jqg=0; jqg < QGRowsleft; jqg++)
            {
                int jLvectorindex       = dataleft[4].fVecShapeIndex[jqg].first;
                int jLshapeindex        = dataleft[4].fVecShapeIndex[jqg].second;
                REAL e1e1j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(0,jLvectorindex))*(n1);
                REAL e2e2j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(1,jLvectorindex))*(n2);
                
                REAL ResidualPart   =   (Theta) * (TimeStep) * ( e1e1j + e2e2j );
                ek(isat + FirstSL,jqg + FirstQGL) 
                += weight * phiSL(isat,0) * ResidualPart;
            }
        }
        
        // (Theta) * deltat * Integrate[L* dot(qg,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat < SRowsRight; isat++)
        {
            for (int jqg=0; jqg < QGRowsleft; jqg++)
            {
                int jLvectorindex       = dataleft[4].fVecShapeIndex[jqg].first;
                int jLshapeindex        = dataleft[4].fVecShapeIndex[jqg].second;
                REAL e1e1j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(0,jLvectorindex))*(n1);
                REAL e2e2j   =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(1,jLvectorindex))*(n2);
                
                REAL ResidualPart   =   (Theta) * (TimeStep) * ( e1e1j + e2e2j );
                ek(isat + FirstSR + iRightInterfaceBlock,jqg + FirstQGL) 
                -= weight * phiSR(isat,0) * ResidualPart;
            }
        }   
 
        
    }
    
    
    //  ////////////////////////// Jacobian Matrix ///////////////////////////////////
    //  End of contribution of countour integrals for Jacobian matrix
    
    
    
    //  n time step
    //  This values are constant in Newton iteration
    if(gState == ELastState)
    {
        
    }
    
    //  if (fnewWS) {
    // Compute Residual contribution ef for contour integrals
    this->ContributeInterface(data, dataleft, dataright, weight, ef);
    //  }
    
    
}



void TPZMultiphase::ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
#ifdef PZDEBUG
    for(int i=0; i<5; i++)
    {
        if(dataleft.find(i) == dataleft.end()) DebugStop();
        if(dataright.find(i) == dataright.end()) DebugStop();
    }
#endif

    TPZFMatrix<REAL> &phiUL = dataleft[0].phi;
    TPZFMatrix<REAL> &phiUR = dataright[0].phi;    

    TPZFMatrix<REAL> &phiQL = dataleft[1].phi;
    TPZFMatrix<REAL> &phiQR = dataright[1].phi;

    TPZFMatrix<REAL> &phiPL = dataleft[2].phi;
    TPZFMatrix<REAL> &phiPR = dataright[2].phi;

    TPZFMatrix<REAL> &phiSL = dataleft[3].phi;
    TPZFMatrix<REAL> &phiSR = dataright[3].phi;

    TPZFMatrix<REAL> &phiQGL = dataleft[4].phi;
//    TPZFMatrix<REAL> &phiQGR = dataright[4].phi;      
    
    TPZManVector<REAL,3> &normal = data.normal;
    
    
    REAL n1 = normal[0];
    REAL n2 = normal[1];
    
    TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
    TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
    TPZManVector<STATE,3> sol_qgL =dataleft[4].sol[0];
    
    TPZManVector<STATE,3> sol_uR =dataright[0].sol[0];
    TPZManVector<STATE,3> sol_qR =dataright[1].sol[0];
    TPZManVector<STATE,3> sol_pR =dataright[2].sol[0];
    TPZManVector<STATE,3> sol_sR =dataright[3].sol[0];
    TPZManVector<STATE,3> sol_qgR    =dataright[4].sol[0];
    
    
    //  Getting Q solution for left and right side
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL qxR = sol_qR[0];
    REAL qyR = sol_qR[1];
    //  REAL qz = sol_qR[2];
    
    //  Getting Qg solution for left and right side
    REAL qgxL = sol_qgL[0];
    REAL qgyL = sol_qgL[1];
    REAL qgxR = sol_qgR[0];
    REAL qgyR = sol_qgR[1];    
    
    //  Getting P solution for left and right side
    REAL PseudoPressureL    =   sol_pL[0];
    REAL PseudoPressureR    =   sol_pR[0];
    
    //  Getting S solution for left and right side
    REAL SaturationL    =   sol_sL[0];
    REAL SaturationR    =   sol_sR[0];
    
    REAL dotqnL = (qxL*n1) + (qyL*n2);
    REAL dotqnR = (qxR*n1) + (qyR*n2);
    
    REAL dotqgnL = (qgxL*n1) + (qgyL*n2);
    REAL dotqgnR = (qgxR*n1) + (qgyR*n2);     
    
    //  Getting another required data
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    REAL Gamma = this->fGamma;
    this->fnewWS=true;
    
    // Getting Harmonic mean of permeabilities
    TPZFMatrix<REAL> KabsoluteLeft;
    TPZFMatrix<REAL> KabsoluteRight;
    TPZFMatrix<REAL> Gfield;
    TPZFMatrix<REAL> Kmean(3,3,0);
    
    int GeoIDLeft = dataleft[1].gelElId;
    int GeoIDRight = dataright[1].gelElId;
    
    REAL rockporosityl, oildensityl, waterdensityl;
    REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;
    REAL rockporosityr, oildensityr, waterdensityr;
    REAL drockporositydpr, doildensitydpr, dwaterdensitydpr;
    
    REAL oilviscosityl, waterviscosityl;
    REAL doilviscositydpl, dwaterviscositydpl;
    REAL oilviscosityr, waterviscosityr;
    REAL doilviscositydpr, dwaterviscositydpr;
    
    REAL bulklambdal, oillambdal, waterlambdal;
    REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
    REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;
    REAL bulklambdar, oillambdar, waterlambdar;
    REAL dbulklambdadpr, doillambdadpr, dwaterlambdadpr;
    REAL dbulklambdadsr, doillambdadsr, dwaterlambdadsr;
    
    
    REAL bulkfoill, bulkfwaterl;
    REAL dbulkfoildpl, dbulkfwaterdpl;
    REAL dbulkfoildsl, dbulkfwaterdsl;
    
    REAL bulkfoilr, bulkfwaterr;
    REAL dbulkfoildpr, dbulkfwaterdpr;
    REAL dbulkfoildsr, dbulkfwaterdsr;
    
    REAL bulkfStarl;
    REAL dbulkfStardpl;
    REAL dbulkfStardsl;
    
    REAL bulkfStarr;
    REAL dbulkfStardpr;
    REAL dbulkfStardsr;
    
    REAL Pcl;
    REAL dPcdSl;
    
    REAL Pcr;
    REAL dPcdSr;     
    
    // Functions computed at point x_{k} for each integration point
    int VecPos= 0;
    //  REAL PressureRef = 1.0e6;
    this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
    this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
    this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
    this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
    this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
    this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
    this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
    this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
    this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
    this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);
    this->CapillaryPressure(sol_sL[VecPos],Pcl,dPcdSl);
    
    this->Porosity(sol_pR[VecPos], rockporosityr, drockporositydpr);
    this->RhoOil(sol_pR[VecPos], oildensityr, doildensitydpr);
    this->RhoWater(sol_pR[VecPos], waterdensityr, dwaterdensitydpr);
    this->OilViscosity(sol_pR[VecPos], oilviscosityr, doilviscositydpr);
    this->WaterViscosity(sol_pR[VecPos], waterviscosityr, dwaterviscositydpr);
    this->OilLabmda(oillambdar, sol_pR[VecPos], sol_sR[VecPos], doillambdadpr, doillambdadsr);
    this->WaterLabmda(waterlambdar, sol_pR[VecPos], sol_sR[VecPos], dwaterlambdadpr, dwaterlambdadsr);
    this->Labmda(bulklambdar, sol_pR[VecPos], sol_sR[VecPos], dbulklambdadpr, dbulklambdadsr);
    this->fOil(bulkfoilr, sol_pR[VecPos], sol_sR[VecPos], dbulkfoildpr, dbulkfoildsr);
    this->fWater(bulkfwaterr, sol_pR[VecPos], sol_sR[VecPos], dbulkfwaterdpr, dbulkfwaterdsr);
    this->CapillaryPressure(sol_sR[VecPos],Pcr,dPcdSr);    
    
    if (fYorN)
    {
        
        KabsoluteLeft=fKabsoluteMap[GeoIDLeft];
        
        KabsoluteRight=fKabsoluteMap[GeoIDRight];
    }
    else
    {
        this->K(KabsoluteLeft);
        this->K(KabsoluteRight);
    }
    
    Gfield = this->Gravity();
    REAL Gravitydotnl =  Gfield(0,0)*(n1)  + Gfield(1,0)*(n2);
    REAL Gravitydotnr =  Gfield(0,0)*(n1)  + Gfield(1,0)*(n2);    
    this->fstar(bulkfStarl,sol_pL[VecPos], sol_sL[VecPos],1.0*Gravitydotnl,dbulkfStardpl,dbulkfStardsl);
    this->fstar(bulkfStarr,sol_pR[VecPos], sol_sR[VecPos],-1.0*Gravitydotnr,dbulkfStardpr,dbulkfStardsr);
    
    // Min of fstar at Gamma
    REAL bulkfstar = 0.0;
    REAL dbulkfstardp = 0.0;
    REAL dbulkfstards = 0.0;    
    if(fabs(bulkfStarl) >= fabs(bulkfStarr))
    {
        bulkfstar = bulkfStarr;
        dbulkfstardp = dbulkfStardpr;
        dbulkfstards = dbulkfStardsr;
    }
    else
    {
        bulkfstar = bulkfStarl;
        dbulkfstardp = dbulkfStardpl;
        dbulkfstards = dbulkfStardsl;        
    }      
    
    
    for (int in=0; in < 3; in++) {
        for (int jn=0; jn < 3; jn++) {
            if (KabsoluteLeft(in,jn)+KabsoluteRight(in,jn)<=0) {
                Kmean(in,jn)= 0.0;
            }
            else
            {
                Kmean(in,jn)= (2.0*KabsoluteLeft(in,jn)*KabsoluteRight(in,jn))/(KabsoluteLeft(in,jn)+KabsoluteRight(in,jn));
            }
            
        }
    }

    TPZFMatrix<STATE> KGL(3,3,0.0),KGR(3,3,0.0), KG(3,3,0.0);
    KGL(0,0) = KabsoluteLeft(0,0)*Gfield(0,0)+KabsoluteLeft(0,1)*Gfield(1,0);
    KGL(1,0) = KabsoluteLeft(1,0)*Gfield(0,0)+KabsoluteLeft(1,1)*Gfield(1,0);        
    KGR(0,0) = KabsoluteRight(0,0)*Gfield(0,0)+KabsoluteRight(0,1)*Gfield(1,0);
    KGR(1,0) = KabsoluteRight(1,0)*Gfield(0,0)+KabsoluteRight(1,1)*Gfield(1,0);
    
    KG(0,0) = Kmean(0,0)*Gfield(0,0)+Kmean(0,1)*Gfield(1,0);
    KG(1,0) = Kmean(1,0)*Gfield(0,0)+Kmean(1,1)*Gfield(1,0);    
    
//    REAL GravityFluxL   =   (KabsoluteLeft(0,0)*Gfield(0,0) + KabsoluteLeft(0,1)*Gfield(1,0))*(n1)+
//    (KabsoluteLeft(1,0)*Gfield(0,0) + KabsoluteLeft(1,1)*Gfield(1,0))*(n2);
//    
//    REAL GravityFluxR   =   (KabsoluteRight(0,0)*Gfield(0,0) + KabsoluteRight(0,1)*Gfield(1,0))*(n1)+
//    (KabsoluteRight(1,0)*Gfield(0,0) + KabsoluteRight(1,1)*Gfield(1,0))*(n2);

    int URowsleft = phiUL.Rows();
    int URowsRight = phiUR.Rows();    
    int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
    int QRowsRight = dataright[1].fVecShapeIndex.NElements();
    int PRowsleft = phiPL.Rows();
    int PRowsRight = phiPR.Rows();
    int SRowsleft = phiSL.Rows();
    int SRowsRight = phiSR.Rows();
    
    int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
//    int QGRowsRight = dataright[4].fVecShapeIndex.NElements();    
    
    int UStateVar = 2;    
    
    int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//    int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;

    int FirstUL = 0;
    int FirstQL = URowsleft * UStateVar + FirstUL;
    int FirstPL = QRowsleft + FirstQL;
    int FirstSL = PRowsleft + FirstPL;
    int FirstQGL = SRowsleft + FirstSL;
    
    int FirstUR = 0;
    int FirstQR = URowsRight * UStateVar + FirstUR;
    int FirstPR = QRowsRight + FirstQR;
    int FirstSR = PRowsRight + FirstPR;
//    int FirstQGR = SRowsRight + FirstSR;
    
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  Contribution of contour integrals for Residual Vector
    
    // n+1 time step
    if(gState == ECurrentState)
    {
        
        
        //This block was verified
        //  First Block (Equation One) constitutive law
        // Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
                REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
                
                ef(iq + FirstQL) += (-1.0) * weight * (e1e1 + e2e2 ) * PseudoPressureL;
            }
            else
            {
                REAL e1e1   =   (Kmean(0,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                 Kmean(0,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n1);
                
                REAL e2e2   =   (Kmean(1,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                 Kmean(1,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n2);
                
                ef(iq + FirstQL) += (-1.0) * weight * (e1e1 + e2e2 ) * PseudoPressureL;
            }
            
            
        }
        
        
        
        //This block was verified
        //  First Block (Equation One) constitutive law
        // Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Right-Right Part
        for (int iq=0; iq < QRowsRight; iq++)
        {
            
            int iRvectorindex       = dataright[1].fVecShapeIndex[iq].first;
            int iRshapeindex        = dataright[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))*(n1);
                REAL e2e2   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex))*(n2);
                
                
                ef(iq + iRightInterfaceBlock + FirstQR) += (1.0) * weight * (e1e1 + e2e2 ) * PseudoPressureR;
            }
            else
            {
                REAL e1e1   =   (Kmean(0,0)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))+
                                 Kmean(0,1)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex)))*(n1);
                
                REAL e2e2   =   (Kmean(1,0)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))+
                                 Kmean(1,1)*(phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex)))*(n2);
                
                ef(iq + iRightInterfaceBlock + FirstQR) += (1.0) * weight * (e1e1 + e2e2 ) * PseudoPressureR;
            }
            
        }
        
        
        REAL SwPcL = SaturationL * Pcl;
        REAL SwPcR = SaturationR * Pcr;
	
        //This block was verified
        //  First Block (Equation One) constitutive law Capillary Pressure
        // Integrate[Sw Pc dot( v, n), Gamme_{e}]  (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
	    REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
	    REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
	    ef(iq + FirstQL) += (-1.0) * (-1.0) * weight * (e1e1 + e2e2 ) * SwPcL;

        }
        
        //This block was verified
        //  First Block (Equation One) constitutive law Capillary Pressure
        // Integrate[Sw Pc dot(v, n), Gamme_{e}]  (Equation One) Right-Right Part
        for (int iq=0; iq < QRowsRight; iq++)
        {
            int iRvectorindex       = dataright[1].fVecShapeIndex[iq].first;
            int iRshapeindex        = dataright[1].fVecShapeIndex[iq].second;    
	    REAL e1e1   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(0,iRvectorindex))*(n1);
	    REAL e2e2   =   (phiQR(iRshapeindex,0)*dataright[1].fDeformedDirections(1,iRvectorindex))*(n2);

	    ef(iq + iRightInterfaceBlock + FirstQR) += (-1.0) * (1.0) * weight * (e1e1 + e2e2 ) * SwPcR;
        }        
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(q, n), Gamme_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            ef(ip + FirstPL) += (-1.0) * (Gamma) * (TimeStep) * weight * dotqnL * phiPL(ip,0);
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(q, n), Gamme_{e}]    (Equation Two) Right-Right Part
        for (int ip=0; ip < PRowsRight; ip++)
        {
            ef(ip + FirstPR + iRightInterfaceBlock) += (1.0) * (Gamma) * (TimeStep) * weight * dotqnR * phiPR(ip,0);
        }
        
        
        REAL UpwindSaturation = 0.0;
        
        if (dotqnL > 0.0)
        {
            UpwindSaturation = bulkfwaterl;
            
        }
        else
        {
            UpwindSaturation = bulkfwaterr;
            
        }
        
        //  Third Vector Block (Equation three)
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat < SRowsleft; isat++)
        {
            REAL ResidualPart   =   (Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
            ef(isat + FirstSL) += weight * phiSL(isat,0) * ResidualPart;
        }
        
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat < SRowsRight; isat++)
        {
            REAL ResidualPart   =   (Theta) * (TimeStep) * ( UpwindSaturation * dotqnR );
            ef(isat + iRightInterfaceBlock + FirstSR) -= weight * phiSR(isat,0) * ResidualPart;
        }

        REAL QGgstarL = bulklambdal * (waterdensityl - oildensityl) * (KGL(0,0)*n1 + KGL(1,0)*n2);
        REAL QGgstarR = bulklambdar * (waterdensityr - oildensityr) * (KGR(0,0)*n1 + KGR(1,0)*n2);        
        REAL QGstar = 0.0;
        
//         if(Gravitydotnl > 0.0 )
//         {
        if(fabs(1.0*bulkfStarl*QGgstarL) < fabs(1.0*bulkfStarr*QGgstarR))
        {    
            QGstar = 1.0*bulkfStarl*QGgstarL;
        }
        else
        {
            QGstar = 1.0*bulkfStarr*QGgstarR;            
        }
        
// #ifdef PZ_LOG
//             if(logdata.isDebugEnabled())
//             {
//                 std::stringstream sout;
//                 sout <<  "SLeft =" << sol_sL[VecPos] << " bulkfStarl =" << bulkfStarl << " Gravitydotnl =" << Gravitydotnl << std::endl;
//                 sout <<  "SRight =" << sol_sR[VecPos] << " bulkfStarr =" << bulkfStarr << " Gravitydotnr =" << -1.0*Gravitydotnr << std::endl;
//                 sout <<  "QGstarL =" << 1.0*bulkfStarl*QGgstarL << " Gravitydotnl =" << Gravitydotnl << std::endl;
//                 sout <<  "QGstarR =" << 1.0*bulkfStarr*QGgstarR << " Gravitydotnl =" << -1.0*Gravitydotnr << std::endl;             
//                 sout <<  "QGstar =" << QGstar << std::endl;                 
//                 LOGPZ_DEBUG(logdata,sout.str());
//             }
// #endif

        // Gravitational segregation scheme
        //  Four Block (Equation Four) gravitational flux constitutive law
        // Integrate[L dot(K v, n), Gamma{e}]  (Equation One) Left-Left part
        for (int iqg=0; iqg < QGRowsleft; iqg++)
        {
            int iLvectorindex       = dataleft[4].fVecShapeIndex[iqg].first;
            int iLshapeindex        = dataleft[4].fVecShapeIndex[iqg].second;
            REAL e1e1   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(0,iLvectorindex))*(n1);
            REAL e2e2   =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(1,iLvectorindex))*(n2);
            // This degree of freedom is the same for Left and Right part
            ef(iqg + FirstQGL) += weight * (dotqgnL-QGstar) * (e1e1 + e2e2);
            
        }
        
        // Gravitational segregation scheme
        // (Theta) * deltat * Integrate[L*dot(qg,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat < SRowsleft; isat++)
        {
            REAL ResidualPart   =   (Theta) * (TimeStep) * ( dotqgnL );
            ef(isat + FirstSL) += weight * phiSL(isat,0) * ResidualPart;
        }
        
        // (Theta) * deltat * Integrate[L* dot(qg,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat < SRowsRight; isat++)
        {
            REAL ResidualPart   =   (Theta) * (TimeStep) * ( dotqgnL );
            ef(isat + iRightInterfaceBlock + FirstSR) -= weight * phiSR(isat,0) * ResidualPart;
        }        
        
//         //  Four Block (Equation Four) gravitational flux constitutive law
//         // Integrate[L dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
//         for (int iqg=0; iqg < QGRowsleft; iqg++)
//         {
//             
//             int iLvectorindex       = dataleft[3].fVecShapeIndex[iqg].first;
//             int iLshapeindex        = dataleft[3].fVecShapeIndex[iqg].second;
//                 
//                 if (fnewWS)
//                 {
//                     REAL e1e1   =   (phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(0,iLvectorindex))*(n1);
//                     REAL e2e2   =   (phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(1,iLvectorindex))*(n2);
//                     ef(iqg+QRowsleft+PRowsleft+SRowsleft) 
//                     += (-1.0) * weight * (e1e1 + e2e2 ) * (waterdensityl - oildensityl);
//                 }
//                 else
//                 {
//                     REAL e1e1   =   (Kmean(0,0)*(phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(0,iLvectorindex))+
//                                      Kmean(0,1)*(phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(1,iLvectorindex)))*(n1);
//                     
//                     REAL e2e2   =   (Kmean(1,0)*(phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(0,iLvectorindex))+
//                                      Kmean(1,1)*(phiQGL(iLshapeindex,0)*dataleft[3].fDeformedDirections(1,iLvectorindex)))*(n2);
//                     ef(iqg+QRowsleft+PRowsleft+SRowsleft) 
//                     += (-1.0) * weight * (e1e1 + e2e2 ) * (waterdensityl - oildensityl);
//                 }
//             
//         }
//         
//         //  Four Block (Equation Four) gravitational flux constitutive law
//         // Integrate[L dot(K v, n), Gamme_{e}]  (Equation One) Right-Right Part
//         for (int iqg=0; iqg < QGRowsRight; iqg++)
//         {
//             int iRvectorindex       = dataright[3].fVecShapeIndex[iqg].first;
//             int iRshapeindex        = dataright[3].fVecShapeIndex[iqg].second;
// 
//                 if (fnewWS)
//                 {
//                     
//                     REAL e1e1   =   (phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(0,iRvectorindex))*(n1);
//                     REAL e2e2   =   (phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(1,iRvectorindex))*(n2);
//                     
//                     ef(iqg + iRightInterfaceBlock + QGRowsRight + PRowsRight + SRowsRight) 
//                     +=  (1.0) * weight * (e1e1 + e2e2 ) * (waterdensityr - oildensityr);
//                 }
//                 else
//                 {
//                     REAL e1e1   =   (Kmean(0,0)*(phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(0,iRvectorindex))+
//                                      Kmean(0,1)*(phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(1,iRvectorindex)))*(n1);
//                     
//                     REAL e2e2   =   (Kmean(1,0)*(phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(0,iRvectorindex))+
//                                      Kmean(1,1)*(phiQGR(iRshapeindex,0)*dataright[3].fDeformedDirections(1,iRvectorindex)))*(n2);
//                     
//                     ef(iqg + iRightInterfaceBlock + QGRowsRight + PRowsRight + SRowsRight) 
//                     +=  (1.0) * weight * (e1e1 + e2e2 ) * (waterdensityr - oildensityr);
//                 }
//             
//         }        
//         
        
        
    }
    
    
    
    
    // n time step
    if(gState == ELastState)
    {
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(q, n), Gamme_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            ef(ip + FirstPL) -= (1-Gamma) * (TimeStep) * weight * dotqnL * phiPL(ip,0);
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(q, n), Gamme_{e}]    (Equation Two) Right-Right Part
        for (int ip=0; ip < PRowsRight; ip++)
        {
            ef(ip + iRightInterfaceBlock + FirstPR) += (1-Gamma) * (TimeStep) * weight * dotqnR * phiPR(ip,0);
        }
        
        REAL dotqnL = (qxL*n1) + (qyL*n2);
        //      REAL dotqnR = (qxR*n1) + (qyR*n2);
        REAL UpwindSaturation = 0.0;
        
        if (dotqnL > 0.0)
        {
            UpwindSaturation = bulkfwaterl;
            
        }
        else
        {
            UpwindSaturation = bulkfwaterr;
            
        }
        
        
        //  Third Vector Block (Equation three)
        // (1-Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
            ef(isat + FirstSL) += weight * phiSL(isat,0) * ResidualPart;
        }
        
        // (1-Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat<SRowsRight; isat++)
        {
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
            ef(isat + iRightInterfaceBlock + FirstSR) -= weight * phiSR(isat,0) * ResidualPart;
        }
        
        // Gravitational segregation scheme
        // (Theta) * deltat * Integrate[L*dot(qg,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( dotqgnL );
            ef(isat + FirstSL) += weight * phiSL(isat,0) * ResidualPart;
        }
        
        // (Theta) * deltat * Integrate[L* dot(qg,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat<SRowsRight; isat++)
        {
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( dotqgnR);
            ef(isat + iRightInterfaceBlock + FirstSR) -= weight * phiSR(isat,0) * ResidualPart;
        }         

    }
    
    //  ////////////////////////// Residual Vector ///////////////////////////////////
    //  End of contribution of contour integrals for Residual Vector
    
}



void TPZMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}

void TPZMultiphase::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
    DebugStop();
}


void TPZMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
    
#ifdef PZDEBUG
    int nref =  datavec.size();
    if (nref != 5 ) {
        std::cout << " Error.!! datavec size not equal to 4 \n";
        DebugStop();
    }
#endif
    
    
    
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    DebugStop();
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
    DebugStop();
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
#ifdef PZDEBUG
    for(int i=0; i<5; i++)
    {
        if(dataleft.find(i) == dataleft.end()) DebugStop();
    }
#endif

    int nref =  dataleft.size();
    if (nref != 5) {
        std::cout << " Error:: datavec size must to be equal to 5 \n" << std::endl;
        DebugStop();
    }
    
    if (bc.Val2().Rows() != 6){
        std::cout << " Error:: This material need boundary conditions for qx, qy, p (pore pressure) and s (Saturation) .\n" << std::endl;
        std::cout << " give me one matrix with this form Val2(3,1).\n" << std::endl;
        DebugStop();
    }
    
    if (bc.Val1().Rows() != 6){
        std::cout << " Error:: This material need boundary conditions for ux, uy, qx, qy, p (pore pressure) and s (Saturation) .\n" << std::endl;
        DebugStop();
    }

    TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
    TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
    TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
    TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
    TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

    int URowsleft = phiUL.Rows();
    int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
    //      int QRowsleft1 = phiQL.Rows();
    int PRowsleft = phiPL.Rows();
    int SRowsleft = phiSL.Rows();
    int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
    
    int UStateVar = 2;     
    
//    int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//    int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
    
    int FirstUL = 0;
    int FirstQL = URowsleft * UStateVar + FirstUL;
    int FirstPL = QRowsleft + FirstQL;
    int FirstSL = PRowsleft + FirstPL;
    int FirstQGL = SRowsleft + FirstSL;

    TPZManVector<REAL,3> &normal = data.normal;
    REAL n1 = normal[0];
    REAL n2 = normal[1];

    TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
    TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
    TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
    TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
    
    //  Getting Q solution for left and right side
//    REAL uxL = sol_uL[0];
//    REAL uyL = sol_uL[1];
    REAL qxL = sol_qL[0];
    REAL qyL = sol_qL[1];
    REAL dotqnL = (qxL*n1) + (qyL*n2);
    
    //  Getting P solution for left side
    REAL PseudoPressureL    =   sol_pL[0];
    
    //  Getting S solution for left side
    //      REAL SaturationL    =   sol_sL[0];
    
    //  Getting another required data
    
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    REAL Gamma = this->fGamma;
    this->fnewWS=true;
    
    // Getting Harmonic mean of permeabilities
    int GeoIDLeft = dataleft[1].gelElId;
    TPZFMatrix<REAL> Kabsolute;
    TPZFMatrix<REAL> Kinverse;
    TPZFMatrix<REAL> Gfield;
    
    if (fYorN)
    {
        Kabsolute=this->fKabsoluteMap[GeoIDLeft];
    }
    else
    {
        this->K(Kabsolute);
    }
    
    Kinverse = this->Kinv(Kabsolute);
    Gfield = this->Gravity();
    
    REAL GravityFluxL   =   (Kabsolute(0,0)*Gfield(0,0) + Kabsolute(0,1)*Gfield(1,0))*(n1)+
                            (Kabsolute(1,0)*Gfield(0,0) + Kabsolute(1,1)*Gfield(1,0))*(n2);
    
    REAL rockporosityl, oildensityl, waterdensityl;
    REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;
    
    REAL oilviscosityl, waterviscosityl;
    REAL doilviscositydpl, dwaterviscositydpl;
    
    REAL bulklambdal, oillambdal, waterlambdal;
    REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
    REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;
    
    REAL bulkfoill, bulkfwaterl;
    REAL dbulkfoildpl, dbulkfwaterdpl;
    REAL dbulkfoildsl, dbulkfwaterdsl;
    
    
    // Functions computed at point x_{k} for each integration point
    int VecPos= 0;
    
    this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
    this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
    this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
    this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
    this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
    this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
    this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
    this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
    this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
    this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);
    
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    //  Regular Controur integrals
    
    // n+1 time step
    if(gState == ECurrentState)
    {

        //First Block (Equation One) constitutive law
        //Integrate[L dot(K v, n), Gamma_{e}]   (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            for (int jp=0; jp < PRowsleft; jp++)
            {
                int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
                int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
                
                if (fnewWS)
                {
                    REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
                    REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
                    
                    ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * (e1e1 + e2e2 ) * phiPL(jp,0);
                    
                }
                else
                {
                    REAL e1e1   =   (Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                     Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n1);
                    
                    REAL e2e2   =   (Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                     Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n2);
                    
                    ek(iq + FirstQL, jp + FirstPL) += (-1.0) * weight * (e1e1 + e2e2 ) * phiPL(jp,0);
                    
                }
            }
        }
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(v, n), Gamme_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            
            for (int jq=0; jq< QRowsleft; jq++)
            {
                
                int jvectorindex    = dataleft[1].fVecShapeIndex[jq].first;
                int jshapeindex     = dataleft[1].fVecShapeIndex[jq].second;
                
                REAL dotprod =
                (n1) * (phiQL(jshapeindex,0)*dataleft[1].fDeformedDirections(0,jvectorindex)) +
                (n2) * (phiQL(jshapeindex,0)*dataleft[1].fDeformedDirections(1,jvectorindex)) ;
                
                ek(ip + FirstPL,jq + FirstQL) += (-1.0) * (Gamma) * (TimeStep) * weight * dotprod * phiPL(ip,0);
                
            }
            
        }
        
        
        //      This block was verified
        //  First Block (Equation One) constitutive law
        // Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
                REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
                
                ef(iq + FirstQL) += (-1.0) * weight * (e1e1 + e2e2) * PseudoPressureL;
                
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n1);
                
                REAL e2e2   =   (Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n2);
                
                
                ef(iq + FirstQL) += (-1.0) * weight * (e1e1 + e2e2) * PseudoPressureL;
                
            }
            
        }
        


        for (int ip=0; ip < PRowsleft; ip++)
        {
            ef(ip+FirstPL) +=  (-1.0) * (Gamma) * (TimeStep) * weight * dotqnL * phiPL(ip,0);
        }
        
        // Null gravitational flux on boundaries        
        for(int iqg=0; iqg < QGRowsleft; iqg++)
        {
            int iLvectorindex       = dataleft[4].fVecShapeIndex[iqg].first;
            int iLshapeindex        = dataleft[4].fVecShapeIndex[iqg].second;
            
            REAL vni    =   (phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(0,iLvectorindex)*n1)+(phiQGL(iLshapeindex,0)*dataleft[4].fDeformedDirections(1,iLvectorindex)*n2);
            
            for (int jqg=0; jqg < QGRowsleft; jqg++)
            {
                int jLvectorindex       = dataleft[4].fVecShapeIndex[jqg].first;
                int jLshapeindex        = dataleft[4].fVecShapeIndex[jqg].second;
                
                REAL vnj    =   (phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(0,jLvectorindex)*n1)+(phiQGL(jLshapeindex,0)*dataleft[4].fDeformedDirections(1,jLvectorindex)*n2);

                ek(iqg + FirstQGL,jqg + FirstQGL) += weight * (gBigNumber * ( vnj ) * vni );
            }
        }        
        
         
     }
    
    
    // n time step
    if(gState == ELastState)
    {
        
        //  Second Block (Equation Two) Bulk flux  equation
        // Integrate[L dot(q, n), Gamme_{e}]    (Equation Two) Left-Left Part
        for (int ip=0; ip < PRowsleft; ip++)
        {
            ef(ip + FirstPL) -=  (1-Gamma) * (TimeStep) * weight * dotqnL * phiPL(ip,0);
        }
        
        REAL UpwindFSaturation = 0.0;
        REAL signofG = 1.0;
        
        if (dotqnL > 0.0)
        {
            UpwindFSaturation = bulkfwaterl;
            
        }
        else
        {
            UpwindFSaturation = 0.0;
            
        }
        
        //  Third Vector Block (Equation three)
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Left-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( UpwindFSaturation * bulkfoill * bulklambdal * (waterdensityl - oildensityl) * GravityFluxL );
            ef(isat + FirstSL) += signofG * (-1.0) * weight * phiSL(isat,0) * ResidualPart;
        }
        
    }
    
    
    //  Regular Controur integrals
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    STATE v2[6];
    v2[0] = bc.Val2()(0,0);         //  ux
    v2[1] = bc.Val2()(1,0);         //  uy
    v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
    v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
    v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
    v2[5] = bc.Val2()(5,0);         //  Saturation
//    REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux
    
    
    switch (bc.Type()) {
        case 1 :    // inflow BC  bc: Ux, Uy, Qn, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;
            
        case 2 :    // inflow BC  bc: Tx, Ty, Qn, Sin
        {
            if(gState == ECurrentState)
            {
                ApplySigmaN(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            

        }
        break;
            
        case 3 :    // inflow BC  bc: Ux, Qn, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;            

        case 4 :    // inflow BC  bc: Uy, Qn, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;
        
        case 5 :    // outflow BC  bc: Ux, Uy, Qn, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break;
            
        case 6 :    // outflow BC  bc: Tx, Ty, Qn, Sout
        {
            if(gState == ECurrentState)
            {
                ApplySigmaN(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            

        }
        break;
            
        case 7 :    // outflow BC  bc: Ux, Qn, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break;            

        case 8 :    // outflow BC  bc: Uy, Qn, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyQnD(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break;          
        
        case 9 :    // inflow BC  bc: Ux, Uy, PN, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;
            
        case 10 :    // inflow BC  bc: Tx, Ty, PN, Sin
        {
            if(gState == ECurrentState)
            {
                ApplySigmaN(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            

        }
        break;
            
        case 11 :    // inflow BC  bc: Ux, PN, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;            

        case 12 :    // inflow BC  bc: Uy, PN, Sin
        {
            if(gState == ECurrentState)
            {
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySin(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySin(data,dataleft,weight,ef,bc);
            }            
        }
        break;
        
        case 13 :    // outflow BC  bc: Ux, Uy, PN, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break;
            
        case 14 :    // outflow BC  bc: Tx, Ty, PN, Sout
        {
            if(gState == ECurrentState)
            {
                ApplySigmaN(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            

        }
        break;
            
        case 15 :    // outflow BC  bc: Ux, PN, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUxD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break;            

        case 16 :    // outflow BC  bc: Uy, PN, Sout
        {
            if(gState == ECurrentState)
            {
                ApplyUyD(data,dataleft,weight,ek,ef,bc);
                ApplyPN(data,dataleft,weight,ek,ef,bc);
                ApplySout(data,dataleft,weight,ek,ef,bc);
            }
            
            if(gState == ELastState)
            {
                ApplySout(data,dataleft,weight,ef,bc);
            }            
        }
        break; 
        
        default: std::cout << "This BC doesn't exist." << std::endl;
        {
            
            DebugStop();
        }
            break;
    }
      
}

    void TPZMultiphase::ApplyUxD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {

        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
//        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
//        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;

        int URowsleft = phiUL.Rows();
//        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        //      int QRowsleft1 = phiQL.Rows();
//        int PRowsleft = phiPL.Rows();
//        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();

//        int UStateVar = 2;     

//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   

//        int FirstUL = 0;
//        int FirstQL = URowsleft * UStateVar + FirstUL;
//        int FirstPL = QRowsleft + FirstQL;
//        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;


        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];

        //  Getting Q solution for left and right side
        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation    
        
                //  Dirichlet condition for each state variable
                //  Elasticity Equation
                for(int iu = 0 ; iu < URowsleft; iu++)
                {
                    //  Contribution for load Vector
                    ef(2*iu)      += this->fxi * gBigNumber * weight * (uxL - v2[0]) * phiUL(iu,0);   // X displacement Value      
//                     ef(2*iu+1)    += gBigNumber * weight * (uyL - v2[1]) * phiUL(iu,0);   // y displacement Value 
                    
                    for(int ju = 0 ; ju < URowsleft; ju++) 
                    {
                        //  Contribution for Stiffness Matrix
                        ek(2*iu,2*ju)       += this->fxi * gBigNumber * weight * phiUL(iu,0) * phiUL(ju,0);  // X displacement
//                         ek(2*iu+1,2*ju+1)   += gBigNumber * weight * phiUL(iu,0) * phiUL(ju,0);  // Y displacement
                    }
                }
        
        
        
    }
    
    void TPZMultiphase::ApplyUyD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {

        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
//        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
//        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;

        int URowsleft = phiUL.Rows();
//        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        //      int QRowsleft1 = phiQL.Rows();
//        int PRowsleft = phiPL.Rows();
//        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();

//        int UStateVar = 2;     

//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   

//        int FirstUL = 0;
//        int FirstQL = URowsleft * UStateVar + FirstUL;
//        int FirstPL = QRowsleft + FirstQL;
//        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;


        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];

        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
        REAL uyL = sol_uL[1];
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation    
        
                //  Dirichlet condition for each state variable
                //  Elasticity Equation
                for(int iu = 0 ; iu < URowsleft; iu++)
                {
                    //  Contribution for load Vector
//                    ef(2*iu)      += gBigNumber * weight * (uxL - v2[0]) * phiUL(iu,0);   // X displacement Value      
                     ef(2*iu+1)    += this->fxi * gBigNumber * weight * (uyL - v2[1]) * phiUL(iu,0);   // y displacement Value 
                    
                    for(int ju = 0 ; ju < URowsleft; ju++) 
                    {
                        //  Contribution for Stiffness Matrix
//                        ek(2*iu,2*ju)       += gBigNumber * weight * phiUL(iu,0) * phiUL(ju,0);  // X displacement
                         ek(2*iu+1,2*ju+1)   += this->fxi * gBigNumber * weight * phiUL(iu,0) * phiUL(ju,0);  // Y displacement
                    }
                }
        
        
    }
    
    void TPZMultiphase::ApplySigmaN    (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
//        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
//        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;

        int URowsleft = phiUL.Rows();
//        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
//        //      int QRowsleft1 = phiQL.Rows();
//        int PRowsleft = phiPL.Rows();
//        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
//
//        int UStateVar = 2;     

//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   

//        int FirstUL = 0;
//        int FirstQL = URowsleft * UStateVar + FirstUL;
//        int FirstPL = QRowsleft + FirstQL;
//        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;


        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];

        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation    
        
                //  Dirichlet condition for each state variable
                //  Elasticity Equation
                for(int iu = 0 ; iu < URowsleft; iu++)
                {
                    //  Contribution for load Vector
                    ef(2*iu)      +=  weight * ( v2[0]) * phiUL(iu,0);   // X Traction Value      
                    ef(2*iu+1)    +=  weight * ( v2[1]) * phiUL(iu,0);   // y Traction Value 
                    
                }        
        
    }
    
    void TPZMultiphase::ApplyQnD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
//        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
//        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
//        int PRowsleft = phiPL.Rows();
//        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
//        int FirstPL = QRowsleft + FirstQL;
//        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        
        //  Getting Q solution for left and right side
        REAL qxL = sol_qL[0];
        REAL qyL = sol_qL[1];
        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux        
        
        
        //  Phil's Hint
        for(int iq=0; iq < QRowsleft; iq++)
        {
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            
            REAL vni    =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex)*n1)+(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)*n2);
            // The coefficient 0.0001 is required for balance the residual contribution
            ef(iq + FirstQL) += weight * ( this->fxi * (gBigNumber * ( dotqnL - qN ) * vni ) );
            //                  ef(iq) += weight * ( (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN  ) * vni ) );
            //                  ef(iq) += weight * ( (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * vni ) );
            
            for (int jq=0; jq < QRowsleft; jq++)
            {
                int jLvectorindex       = dataleft[1].fVecShapeIndex[jq].first;
                int jLshapeindex        = dataleft[1].fVecShapeIndex[jq].second;
                
                REAL vnj    =   (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(0,jLvectorindex)*n1)+(phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(1,jLvectorindex)*n2);
                ek(iq + FirstQL,jq + FirstQL) += weight * ( this->fxi * (gBigNumber * ( vnj ) * vni ) );
                //                      ek(iq,jq) += weight * ( 2.0 * (gBigNumber * ( dotqnL - qN ) * ( vnj ) * vni ) );
                //                      ek(iq,jq) += weight * ( 4.0 * (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * ( vnj ) * vni ) );
            }
        }        
    }

    
    void TPZMultiphase::ApplyPN        (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
//        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
//        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
//        int PRowsleft = phiPL.Rows();
//        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
//        int FirstPL = QRowsleft + FirstQL;
//        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
        TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
        
        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
//        REAL qxL = sol_qL[0];
//        REAL qyL = sol_qL[1];
//        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        //  Getting P solution for left side
//        REAL PseudoPressureL    =   sol_pL[0];
        
        //  Getting another required data
        
//        REAL TimeStep = this->fDeltaT;
//        REAL Theta = this->fTheta;
//        REAL Gamma = this->fGamma;
        this->fnewWS=true;
        
        // Getting Harmonic mean of permeabilities
        int GeoIDLeft = dataleft[1].gelElId;
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;        
        
        if (fYorN)
        {
            Kabsolute=this->fKabsoluteMap[GeoIDLeft];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse = this->Kinv(Kabsolute);
        Gfield = this->Gravity();    
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
//        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux         
        
        
        //  This block was verified
        //  First Block (Equation One) constitutive law
        //  Integrate[P dot(K v, n), Gamme_{e}]  (Equation One) Left-Left part
        for (int iq=0; iq < QRowsleft; iq++)
        {
            
            int iLvectorindex       = dataleft[1].fVecShapeIndex[iq].first;
            int iLshapeindex        = dataleft[1].fVecShapeIndex[iq].second;
            
            if (fnewWS)
            {
                
                REAL e1e1   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))*(n1);
                REAL e2e2   =   (phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex))*(n2);
                
                ef(iq + FirstQL) += (1.0) * weight * (e1e1 + e2e2 ) * (v2[4]);
                
            }
            else
            {
                REAL e1e1   =   (Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n1);
                
                REAL e2e2   =   (Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(0,iLvectorindex))+
                                Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[1].fDeformedDirections(1,iLvectorindex)))*(n2);
                
                ef(iq + FirstQL) += (1.0) * weight * (e1e1 + e2e2 ) * (v2[4]);
                
            }
        }        
        
    }
    
    void TPZMultiphase::ApplySin       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        int PRowsleft = phiPL.Rows();
        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
        int FirstPL = QRowsleft + FirstQL;
        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
        TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
        
        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
//        REAL qxL = sol_qL[0];
//        REAL qyL = sol_qL[1];
//        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        //  Getting P solution for left side
//        REAL PseudoPressureL    =   sol_pL[0];
        
        //  Getting another required data
        
        REAL TimeStep = this->fDeltaT;
        REAL Theta = this->fTheta;
//        REAL Gamma = this->fGamma;
        this->fnewWS=true;
        
        // Getting Harmonic mean of permeabilities
        int GeoIDLeft = dataleft[1].gelElId;
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;        
        
        if (fYorN)
        {
            Kabsolute=this->fKabsoluteMap[GeoIDLeft];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse = this->Kinv(Kabsolute);
        Gfield = this->Gravity();

//        REAL GravityFluxL   =   (Kabsolute(0,0)*Gfield(0,0) + Kabsolute(0,1)*Gfield(1,0))*(n1)+
//                                (Kabsolute(1,0)*Gfield(0,0) + Kabsolute(1,1)*Gfield(1,0))*(n2);

        REAL rockporosityl, oildensityl, waterdensityl;
        REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;

        REAL oilviscosityl, waterviscosityl;
        REAL doilviscositydpl, dwaterviscositydpl;

        REAL bulklambdal, oillambdal, waterlambdal;
        REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
        REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;

        REAL bulkfoill, bulkfwaterl;
        REAL dbulkfoildpl, dbulkfwaterdpl;
        REAL dbulkfoildsl, dbulkfwaterdsl;


        // Functions computed at point x_{k} for each integration point
        int VecPos= 0;

        this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
        this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
        this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
        this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
        this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
        this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
        this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
        this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
        this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
        this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux             
        
                //  Upwind scheme
                //  Third Vector Block (Equation three) Saturation  equation
                REAL UpwindSaturation = 0.0;
                
//              if (dotqnL > 0.0)
//              {
//                    this->fWater(bulkfwaterl, sol_pL[VecPos], v2[5], dbulkfwaterdpl, dbulkfwaterdsl);
//                    UpwindSaturation = bulkfwaterl;
//                  
//                  //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Bc-Left Part
//                  for(int isat=0; isat<SRowsleft; isat++)
//                  {
//                      for(int jsat=0; jsat<SRowsleft; jsat++)
//                      {
//                          ek(isat+QRowsleft+PRowsleft,jsat+QRowsleft+PRowsleft) -= (-1.0) * weight *
//                          (Theta) * (TimeStep) * phiSL(isat,0)* dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
//                      }
//                  }
//                  
//                  
//              }
//              
//              else
//              {
                    this->fWater(bulkfwaterl, sol_pL[VecPos], v2[5], dbulkfwaterdpl, dbulkfwaterdsl);
                    UpwindSaturation = bulkfwaterl;		    
                    
//              }
             
             
//              //  Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]  (Equation three) Right-Left Part
//              for (int isat=0; isat < SRowsleft; isat++) {
//                  
//                  for (int jq=0; jq<QRowsleft; jq++)
//                  {
//                      int jLvectorindex       = dataleft[1].fVecShapeIndex[jq].first;
//                      int jLshapeindex        = dataleft[1].fVecShapeIndex[jq].second;
//                      
//                      REAL dotprodL =
//                      (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(0,jLvectorindex)) * (n1) +
//                      (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(1,jLvectorindex)) * (n2) ;//+
//                      
//                      ek(isat+QRowsleft+PRowsleft,jq) -= (-1.0) * weight * (Theta) * (TimeStep) * phiSL(isat,0) * (UpwindSaturation) * dotprodL;
//                      
//                  }
//                  
//              }
                
                // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
                for(int isat=0; isat < SRowsleft; isat++)
                {
                    REAL ResidualPart   =   (Theta) * (TimeStep) * ( (UpwindSaturation) * qN );
                    ef(isat + FirstSL) -= (-1.0) * weight * phiSL(isat,0) * ResidualPart;
                }        
        
    }
    
    void TPZMultiphase::ApplySout      (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
    
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        int PRowsleft = phiPL.Rows();
        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
        int FirstPL = QRowsleft + FirstQL;
        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
        TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
        
        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
        REAL qxL = sol_qL[0];
        REAL qyL = sol_qL[1];
        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        //  Getting P solution for left side
//        REAL PseudoPressureL    =   sol_pL[0];
        
        //  Getting another required data
        
        REAL TimeStep = this->fDeltaT;
        REAL Theta = this->fTheta;
//        REAL Gamma = this->fGamma;
        this->fnewWS=true;
        
        // Getting Harmonic mean of permeabilities
        int GeoIDLeft = dataleft[1].gelElId;
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;        
        
        if (fYorN)
        {
            Kabsolute=this->fKabsoluteMap[GeoIDLeft];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse = this->Kinv(Kabsolute);
        Gfield = this->Gravity();

//        REAL GravityFluxL   =   (Kabsolute(0,0)*Gfield(0,0) + Kabsolute(0,1)*Gfield(1,0))*(n1)+
//                                (Kabsolute(1,0)*Gfield(0,0) + Kabsolute(1,1)*Gfield(1,0))*(n2);

        REAL rockporosityl, oildensityl, waterdensityl;
        REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;

        REAL oilviscosityl, waterviscosityl;
        REAL doilviscositydpl, dwaterviscositydpl;

        REAL bulklambdal, oillambdal, waterlambdal;
        REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
        REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;

        REAL bulkfoill, bulkfwaterl;
        REAL dbulkfoildpl, dbulkfwaterdpl;
        REAL dbulkfoildsl, dbulkfwaterdsl;


        // Functions computed at point x_{k} for each integration point
        int VecPos= 0;

        this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
        this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
        this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
        this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
        this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
        this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
        this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
        this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
        this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
        this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
//        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux             
        
        
        //  Upwind scheme
        //  Third Vector Block (Equation three) Saturation  equation
        REAL UpwindSaturation = 0.0;
        
        if (dotqnL > 0.0)
        {
            
            UpwindSaturation = bulkfwaterl;
            //  Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}] (Equation three) Bc-Left Part
            for(int isat=0; isat<SRowsleft; isat++)
            {
                for(int jsat=0; jsat<SRowsleft; jsat++)
                {
                    ek(isat + FirstSL,jsat + FirstSL) -= (-1.0) * weight
                    * (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
                }
            }
        }
        
        else
            
        {
            UpwindSaturation = bulkfwaterl;
            if (dotqnL < 0.0 && fabs(dotqnL) > 1.0e-10) { std::cout << "Boundary condition error: inflow detected in outflow boundary condition: dotqnL = " << dotqnL << "\n";}
        }
        
        //  Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]  (Equation three) Right-Left Part
        for (int isat=0; isat < SRowsleft; isat++) {
            
            for (int jq=0; jq<QRowsleft; jq++)
            {
                int jLvectorindex       = dataleft[1].fVecShapeIndex[jq].first;
                int jLshapeindex        = dataleft[1].fVecShapeIndex[jq].second;
                
                
                REAL dotprodL =
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(0,jLvectorindex)) * (n1) +
                (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(1,jLvectorindex)) * (n2) ;//+
                //              (phiQL(jLshapeindex,0)*dataleft[1].fDeformedDirections(2,jLvectorindex)) * (n3) ;    //  dot(q,n)    left
                
                ek(isat + FirstSL,jq + FirstQL) -= (-1.0) * weight * (Theta) * (TimeStep) * phiSL(isat,0) * (UpwindSaturation) * dotprodL;
                
            }
        }
        
        
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {
            REAL ResidualPart   =   (Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
            ef(isat + FirstSL) -= (-1.0) * weight * phiSL(isat,0) * ResidualPart;
        }        
        
    }
    
    void TPZMultiphase::ApplySin       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        int PRowsleft = phiPL.Rows();
        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
        int FirstPL = QRowsleft + FirstQL;
        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
        TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
        
        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
//        REAL qxL = sol_qL[0];
//        REAL qyL = sol_qL[1];
//        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        //  Getting P solution for left side
//        REAL PseudoPressureL    =   sol_pL[0];
        
        //  Getting another required data
        
        REAL TimeStep = this->fDeltaT;
        REAL Theta = this->fTheta;
//        REAL Gamma = this->fGamma;
        this->fnewWS=true;
        
        // Getting Harmonic mean of permeabilities
        int GeoIDLeft = dataleft[1].gelElId;
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;        
        
        if (fYorN)
        {
            Kabsolute=this->fKabsoluteMap[GeoIDLeft];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse = this->Kinv(Kabsolute);
        Gfield = this->Gravity();

//        REAL GravityFluxL   =   (Kabsolute(0,0)*Gfield(0,0) + Kabsolute(0,1)*Gfield(1,0))*(n1)+
//                                (Kabsolute(1,0)*Gfield(0,0) + Kabsolute(1,1)*Gfield(1,0))*(n2);

        REAL rockporosityl, oildensityl, waterdensityl;
        REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;

        REAL oilviscosityl, waterviscosityl;
        REAL doilviscositydpl, dwaterviscositydpl;

        REAL bulklambdal, oillambdal, waterlambdal;
        REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
        REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;

        REAL bulkfoill, bulkfwaterl;
        REAL dbulkfoildpl, dbulkfwaterdpl;
        REAL dbulkfoildsl, dbulkfwaterdsl;


        // Functions computed at point x_{k} for each integration point
        int VecPos= 0;

        this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
        this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
        this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
        this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
        this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
        this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
        this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
        this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
        this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
        this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux         
        
        REAL UpwindSaturation = 0.0;
        
//              
//              if (dotqnL > 0.0)
//              {
//                  this->fWater(bulkfwaterl, sol_pL[VecPos], v2[3], dbulkfwaterdpl, dbulkfwaterdsl);
//                  UpwindSaturation = bulkfwaterl;
//              }
//              
//              else
//              {
            this->fWater(bulkfwaterl, sol_pL[VecPos], v2[5], dbulkfwaterdpl, dbulkfwaterdsl);
            UpwindSaturation = bulkfwaterl;
            
//                  
//              }
//              
        
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {
            
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( (UpwindSaturation) * qN );
            ef(isat + FirstSL) -= (-1.0) * weight * phiSL(isat,0) * ResidualPart;
        }        
        
    }
    
    void TPZMultiphase::ApplySout      (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
    {
        
        TPZFMatrix<REAL> &phiUL     = dataleft[0].phi;    
//        TPZFMatrix<REAL> &phiQL     = dataleft[1].phi;
        TPZFMatrix<REAL> &phiPL     = dataleft[2].phi;
        TPZFMatrix<REAL> &phiSL     = dataleft[3].phi;
//        TPZFMatrix<REAL> &phiQGL    = dataleft[4].phi;    

        int URowsleft = phiUL.Rows();
        int QRowsleft = dataleft[1].fVecShapeIndex.NElements();
        int PRowsleft = phiPL.Rows();
        int SRowsleft = phiSL.Rows();
//        int QGRowsleft = dataleft[4].fVecShapeIndex.NElements();
        
        int UStateVar = 2;     
        
//        int iRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;
//        int jRightInterfaceBlock = URowsleft * UStateVar + QRowsleft + PRowsleft + SRowsleft + QGRowsleft;   
        
        int FirstUL = 0;
        int FirstQL = URowsleft * UStateVar + FirstUL;
        int FirstPL = QRowsleft + FirstQL;
        int FirstSL = PRowsleft + FirstPL;
//        int FirstQGL = SRowsleft + FirstSL;

        TPZManVector<REAL,3> &normal = data.normal;
        REAL n1 = normal[0];
        REAL n2 = normal[1];

        TPZManVector<STATE,3> sol_uL =dataleft[0].sol[0];
        TPZManVector<STATE,3> sol_qL =dataleft[1].sol[0];
        TPZManVector<STATE,3> sol_pL =dataleft[2].sol[0];
        TPZManVector<STATE,3> sol_sL =dataleft[3].sol[0];
        
        //  Getting Q solution for left and right side
//        REAL uxL = sol_uL[0];
//        REAL uyL = sol_uL[1];
        REAL qxL = sol_qL[0];
        REAL qyL = sol_qL[1];
        REAL dotqnL = (qxL*n1) + (qyL*n2);
        
        //  Getting P solution for left side
//        REAL PseudoPressureL    =   sol_pL[0];
        
        //  Getting another required data
        
        REAL TimeStep = this->fDeltaT;
        REAL Theta = this->fTheta;
//        REAL Gamma = this->fGamma;
        this->fnewWS=true;
        
        // Getting Harmonic mean of permeabilities
        int GeoIDLeft = dataleft[1].gelElId;
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;        
        
        if (fYorN)
        {
            Kabsolute=this->fKabsoluteMap[GeoIDLeft];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse = this->Kinv(Kabsolute);
        Gfield = this->Gravity();

//        REAL GravityFluxL   =   (Kabsolute(0,0)*Gfield(0,0) + Kabsolute(0,1)*Gfield(1,0))*(n1)+
//                                (Kabsolute(1,0)*Gfield(0,0) + Kabsolute(1,1)*Gfield(1,0))*(n2);

        REAL rockporosityl, oildensityl, waterdensityl;
        REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;

        REAL oilviscosityl, waterviscosityl;
        REAL doilviscositydpl, dwaterviscositydpl;

        REAL bulklambdal, oillambdal, waterlambdal;
        REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
        REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;

        REAL bulkfoill, bulkfwaterl;
        REAL dbulkfoildpl, dbulkfwaterdpl;
        REAL dbulkfoildsl, dbulkfwaterdsl;


        // Functions computed at point x_{k} for each integration point
        int VecPos= 0;

        this->Porosity(sol_pL[VecPos], rockporosityl, drockporositydpl);
        this->RhoOil(sol_pL[VecPos], oildensityl, doildensitydpl);
        this->RhoWater(sol_pL[VecPos], waterdensityl, dwaterdensitydpl);
        this->OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
        this->WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);
        this->OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
        this->WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
        this->Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
        this->fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
        this->fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        
        
        STATE v2[6];
        v2[0] = bc.Val2()(0,0);         //  ux
        v2[1] = bc.Val2()(1,0);         //  uy
        v2[2] = bc.Val2()(2,0)/this->fQref;    //  qx
        v2[3] = bc.Val2()(3,0)/this->fQref;    //  qy
        v2[4] = bc.Val2()(4,0)/this->fPref; //  Pressure
        v2[5] = bc.Val2()(5,0);         //  Saturation
//        REAL qN = (v2[2]*n1 + v2[3]*n2);    // Normal Flux       
        
        REAL UpwindSaturation = 0.0;
        
        if (dotqnL > 0.0)
        {
            UpwindSaturation = bulkfwaterl;
        }
        else
        {
            UpwindSaturation = bulkfwaterl;
            if (dotqnL < 0.0 && fabs(dotqnL) > 1.0e-12) {std::cout << "Boundary condition error: inflow detected in outflow boundary condition: dotqnL = " << dotqnL << "\n";}
        }
        
        // (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
        for(int isat=0; isat<SRowsleft; isat++)
        {   
            REAL ResidualPart   =   (1-Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
            ef(isat + FirstSL) -= (-1.0) * weight * phiSL(isat,0) * ResidualPart;
        }         
        
    }    

/** Returns the variable index associated with the name */
int TPZMultiphase::VariableIndex(const std::string &name){
    if(!strcmp("BulkVelocity",name.c_str()))        return  1;
    if(!strcmp("WaterVelocity",name.c_str()))        return  2;
    if(!strcmp("OilVelocity",name.c_str()))        return  3;
    if(!strcmp("WeightedPressure",name.c_str()))    return  4;
    if(!strcmp("WaterSaturation",name.c_str()))    return  5;
    if(!strcmp("OilSaturation",name.c_str()))    return  6;
    if(!strcmp("WaterDensity",name.c_str()))    return  7;
    if(!strcmp("OilDensity",name.c_str()))    return  8;
    if(!strcmp("RockPorosity",name.c_str()))    return  9;
    if(!strcmp("GravityVelocity",name.c_str()))        return  10;
    if(!strcmp("Kabsolute",name.c_str()))    return  11;
    if(!strcmp("SwExact",name.c_str()))    return  12;
    if(!strcmp("Displacement",name.c_str()))    return  13;
    if(!strcmp("SigmaX",name.c_str()))                      return  14;
    if(!strcmp("SigmaY",name.c_str()))                      return  15;    
    
    return TPZMaterial::VariableIndex(name);
}

int TPZMultiphase::NSolutionVariables(int var){
    if(var == 1) return 2;
    if(var == 2) return 2;
    if(var == 3) return 2;
    if(var == 4) return 1;
    if(var == 5) return 1;
    if(var == 6) return 1;
    if(var == 7) return 1;
    if(var == 8) return 1;
    if(var == 9) return 1;
    if(var == 10) return 2;
    if(var == 11) return 3;
    if(var == 12) return 1;
    if(var == 13) return 3;  
    if(var == 14) return 1;
    if(var == 15) return 1;     
    
    return TPZMaterial::NSolutionVariables(var);
}

void TPZMultiphase::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
    
    Solout.Resize( this->NSolutionVariables(var));
    TPZVec<STATE> SolU, SolQ, SolP, SolS, SolQG, SolSExact(1);
    TPZFMatrix<STATE> SolSExactD;
    SolU = datavec[0].sol[0];    
    SolQ = datavec[1].sol[0];
    SolP = datavec[2].sol[0];
    SolS = datavec[3].sol[0];
    SolQG   = datavec[4].sol[0];
    
    TPZFNMatrix <6,STATE> dSolU;    
    dSolU = datavec[0].dsol[0];
    TPZFNMatrix <9> axesU;
    axesU = datavec[0].axes;
    
    REAL epsx;
    REAL epsy;
    REAL epsxy;
    REAL SigX;
    REAL SigY;
    REAL SigZ;
    REAL Tau, DSolxy[2][2];
    REAL divu;  
    
    if(var == 1){ //function (state variable Q)
        Solout[0] = SolQ[0];
        Solout[1] = SolQ[1];
        return;
    }
    
    if(var == 4){
        Solout[0] = SolP[0];//function (state variable p)
        return;
    }
    
    if(var == 5){
        Solout[0] = SolS[0];//function (state variable S)
        return;
    }
    
    if(var == 6){
        Solout[0] = 1.0-SolS[0];//function (state variable S)
        return;
    }   
    
    if(var == 10){ //function (state variable QG)
        Solout[0] = SolQG[0];
        Solout[1] = SolQG[1];
        return;
    }
    
    if(var == 11){ //function (state variable QG)
        
        TPZFMatrix<REAL> Kabsolute;
        TPZFMatrix<REAL> Kinverse;
        TPZFMatrix<REAL> Gfield;
        if (fYorN)
        {
            int iel  = datavec[0].gelElId;
            Kabsolute=this->fKabsoluteMap[iel];
        }
        else
        {
            this->K(Kabsolute);
        }
        
        Kinverse=this->Kinv(Kabsolute);
        Gfield = this->Gravity();

        Solout[0] = Kabsolute(0,0);
        Solout[1] = Kabsolute(1,1);
        Solout[2] = Kabsolute(2,2);
        
        return;
    }
    
    
    if(var == 12){ // ExactSolution
		fTimedependentFunctionExact->Execute(datavec[2].x, fTime, SolSExact,SolSExactD);
		Solout[0] = SolSExact[0];
        return;
    }
    
    if(var == 13){ //function (state variable U)
        Solout[0] = SolU[0];
        Solout[1] = SolU[1];
        Solout[2] = 0.0;        
        return;
    }    
    

                        
    DSolxy[0][0] = dSolU(0,0)*axesU(0,0)+dSolU(1,0)*axesU(1,0); // dUx/dx
    DSolxy[1][0] = dSolU(0,0)*axesU(0,1)+dSolU(1,0)*axesU(1,1); // dUx/dy
    
    DSolxy[0][1] = dSolU(0,1)*axesU(0,0)+dSolU(1,1)*axesU(1,0); // dUy/dx
    DSolxy[1][1] = dSolU(0,1)*axesU(0,1)+dSolU(1,1)*axesU(1,1); // dUy/dy
    
    divu = DSolxy[0][0]+DSolxy[1][1]+0.0;   
    
    REAL lambda,lambdau, mu, Balpha, Sestr;
        // Functions computed at point x_{k} for each integration point
    lambda     = this->LameLambda();
    lambdau    = this->LameLambdaU();
    mu         = this->LameMu();
    Balpha      = this->BiotAlpha();
    Sestr       = this->Se();
    
    epsx = DSolxy[0][0];// du/dx
    epsy = DSolxy[1][1];// dv/dy
    epsxy = 0.5*(DSolxy[1][0]+DSolxy[0][1]);
    REAL C11 = 4*(mu)*(lambda+mu)/(lambda+2*mu);
    REAL C22 = 2*(mu)*(lambda)/(lambda+2*mu);
    
    if (this->fPlaneStress==1)
    {
        SigX = C11*epsx+C22*epsy;
        SigY = C11*epsy+C22*epsx;
        SigZ = 0.0;
        Tau = 2.0*mu*epsxy;
    }
    else
    {
        SigX = ((lambda + 2*mu)*(epsx) + (lambda)*epsy - Balpha*SolP[0]);
        SigY = ((lambda + 2*mu)*(epsy) + (lambda)*epsx - Balpha*SolP[0]);    
        SigZ = lambda*divu - Balpha*SolP[0];
        Tau = 2.0*mu*epsxy;        
    }
    
    if(var == 14){ // (SigmaX)
        Solout[0] = SigX;        
        return;
    }       

    if(var == 15){ // (SigmaY)
        Solout[0] = SigY;        
        return;
    }           
    
}


void TPZMultiphase::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
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

void TPZMultiphase::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
    for(int i = 0; i<nref; i++)
    {
        datavec[i].fNeedsSol = true;
        datavec[i].fNeedsNormal = true;
        datavec[i].fNeedsNeighborSol = true;        
    }
}

int TPZMultiphase::ClassId() const{
    return Hash("TPZMultiphase") ^ TPZMaterial::ClassId() << 1;
}
