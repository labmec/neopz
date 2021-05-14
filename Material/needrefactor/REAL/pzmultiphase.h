//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_pzmultiphase_h
#define PZ_pzmultiphase_h

#include "TPZMaterial.h"
#include "fad.h"
#include "pzstack.h"
#include <iostream>
#include <fstream>
#include <string>
/**
 * @ingroup material
 * @author Omar Duran
 * @since 19/08/2013
 * @brief Material to solve a 2d multiphase transport problem by multiphysics simulation
 * @brief Here is used L, Hdiv ... spaces and first order upwind scheme
 */


class TPZMultiphase : public TPZMaterial {
    
protected:
    
    /** @brief Problem dimension */
    int fDim;
    
    /** @brief Material id */
    int fmatId; 
    
    /** @brief Definition of constants */
    REAL ff;
    
    /** @brief Big number balance */
    REAL fxi;    
    
    /** @brief State: Stiffness or Mass Matrix Calculations */
    enum EState { ELastState = 0, ECurrentState = 1 };
    EState gState;
    

    typedef TFad<2, REAL> BFadREAL; 
    

public:
    
    bool fnewWS;
    
    TPZMultiphase();
    
    TPZMultiphase(int matid, int dim);
    
    virtual ~TPZMultiphase();
    
    /** @brief copy constructor */
    TPZMultiphase(const TPZMultiphase &copy);
    
    TPZMultiphase &operator=(const TPZMultiphase &copy);
    
    virtual void Print(std::ostream & out) override;
    
    virtual std::string Name()  override { return "TPZMultiphase"; }
    
    virtual int Dimension() const override;
    
    virtual int NStateVariables() const override;  
    
    virtual int MatId();
    
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, std::map<int, TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    virtual void ContributeBCInterface(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    // Contribution for each state variable
    virtual void ApplyUxD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplyUyD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplySigmaN    (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplyQnD       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplyPN        (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplySin       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplySout      (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplySin       (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    virtual void ApplySout      (TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec) override;
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec) override;
    
    
    virtual int VariableIndex(const std::string &name) override;
    
    virtual int NSolutionVariables(int var) override;
    
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
    
    virtual void Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * Left, TPZCompEl * Right) override
    {
        TPZMaterial::Solution(data,dataleftvec,datarightvec,var,Solout,Left,Right);
    }
    
    
    // Here is needed to spefify which models are used in this bi-phasic model see Advanced-Petroleum-Reservoir-Simulation M. Rafiqul Islam
    
    /** @brief Characteristic length - m */
    REAL fLref; 
    
    /** @brief Permeability reference - m2 */
    REAL fKref;
    
    /** @brief Pressure reference - Pa */   
    REAL fQref;    
    
    /** @brief Pressure reference - Pa */   
    REAL fPref;
    
    /** @brief density reference - kg/m3 */
    REAL fRhoref;
    
    /** @brief viscosity reference - Pa s */
    REAL fEtaref;
    
    
    /** @brief K map */ 
    TPZStack< TPZFMatrix<REAL> > fKabsoluteMap; 

    /** @brief plane stress condition */    
    int fPlaneStress;
    
    /** @brief Use or not K map */      
    bool fYorN;
    
    /** @brief Simulation time step */
    REAL fDeltaT;
    
    /** @brief Simulation current time */
    REAL fTime;
    
    /** @brief Parameter representing temporal scheme for transport equation */
    REAL fTheta;
    
    /** @brief Parameter representing temporal scheme for conservation equation */
    REAL fGamma;    
    
    /** @brief Defines simulation bc numeric balance. */
    void SetXiBalance(REAL xi){ this->fxi = xi;}
        
    /** @brief Defines simulation time step. */
    void SetTimeStep(REAL timestep){ this->fDeltaT = timestep;}
    
    /** @brief Defines simulation time step. */
    void SetTime(REAL time){ this->fTime = time;}
    
    /** @brief Defines stemporal scheme. */
    void SetTScheme(REAL timegamma, REAL timetheta){ this->fTheta = timetheta; this->fGamma = timegamma;}
    
    void SetLastState(){ gState = ELastState;}
    
    void SetCurrentState(){ gState = ECurrentState;}
    
    /** @brief Defines simulation use a K map */
    void SetYorN(bool dummybool){ this->fYorN = dummybool;}     
    
    /**
     * @brief Lame First Parameter.
     * \f$ lamelamda \f$
     */
    STATE LameLambda();

    /**
     * @brief Undrained Lame First Parameter.
     * \f$ lamelamdaU \f$
     */
    STATE LameLambdaU();
    
    /**
     * @brief Lame Second Parameter.
     * \f$ lamemu \f$
     */
    STATE LameMu();

    /**
     * @brief Biot parameter Parameter.
     * \f$ lamemu \f$
     */
    STATE BiotAlpha();

    /**
     * @brief //Se o 1/M coeficiente poroelastico de armazenamento a volume constante.
     * \f$ Se \f$
     */
    STATE Se();
    
    
    /**
     * @name Getting Data
     * @{
     */
    
    void SetKMap(TPZStack< TPZFMatrix<REAL> > KabsoluteMap)
    {
        fKabsoluteMap = KabsoluteMap;
    }   
    
    /** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
    void CapillaryPressure(REAL So, REAL &pc, REAL &DpcDSo);
    
    /**
     * @brief Oil relative permeability.
     * \f$ Kro = Kro( Sw ) \f$
     */
    void Kro(REAL Sw, REAL &Kro, REAL &dKroDSw);
    
    /**
     * @brief Water relative permeability.
     * \f$ Krw = Krw( Sw ) \f$
     */
    void Krw(REAL Sw, REAL &Krw, REAL &dKrwSo);
    
    /**
     * @brief Water saturation maximum value of the fractional flow product function.
     * \f$ S* = Sw* \f$
     */
    void SWaterstar(REAL &Swstar, REAL &Po, REAL &Sw); 
    
    /**
     * @brief Fractional product function, water decrease direction (dSw/dt < 0).
     * \f$ fw*  \f$
     */
    void fWaterstar(REAL &fWstar, REAL Pw, REAL Sw, REAL &dfWstarDPw, REAL &dfWstarDSw);

    /**
     * @brief Fractional product function, oil decrease direction (dSw/dt > 0).
     * \f$ fo*  \f$
     */
    void fOilstar(REAL &fOstar, REAL Pw, REAL Sw, REAL &dfOstarDPw, REAL &dfOstarDSw);    

    /**
     * @brief Fractional product upwind function, Gdotn > 0 means water decrease (dSw/dt < 0), Gdotn < 0 means water increase (dSw/dt > 0).
     * \f$ f*  \f$
     */
    void fstar(REAL &fStar, REAL Pw, REAL Sw, REAL Gdotn, REAL &dfstarDPw, REAL &dfstarDSw);      
    
    /** 
     * @brief  Rock porosity. \f$ Phi = Phi( p ) \f$
     * @param po Refrence pressure
     */ 
    void Porosity(REAL po, REAL &poros, REAL &dPorosDp);
    
    /** 
     * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
     * @param po Refrence pressure
     */
    void RhoOil(REAL po, REAL &RhoOil, REAL &dRhoOilDpo);
    
    /** 
     * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
     * @param pw Refrence pressure
     */
    void RhoWater(REAL pw, REAL &RhoWater, REAL &dRhoWaterDpo);
    
    /** 
     * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
     * @param po Refrence pressure
     */
    void OilViscosity(REAL po, REAL &OilViscosity, REAL &dOilViscosityDpo);
    
    /** 
     * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
     * @param po Refrence pressure
     */ 
    void WaterViscosity(REAL po, REAL &WaterViscosity, REAL &dWaterViscosityDpo);
    
    /**
     * @brief Oil mobility.
     * \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$
     */
    void OilLabmda(REAL &OilLabmda, REAL Po, REAL Sw, REAL &dOilLabmdaDPo, REAL &dOilLabmdaDSw);
    
    /**
     * @brief Water mobility.
     * \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$
     */
    void WaterLabmda(REAL &WaterLabmda, REAL Pw, REAL Sw, REAL &dWaterLabmdaDPw, REAL &dWaterLabmdaDSw);
    
    /**
     * @brief Bulk mobility.
     * \f$ \lambda = \lambda( pw , Sw ) \f$
     */
    void Labmda(REAL &Labmda, REAL Pw, REAL Sw, REAL &dLabmdaDPw, REAL &dLabmdaDSw);
    
    /**
     * @brief Fractional oil flux.
     * \f$ f_{Oil} = f_{Oil}( po , Sw ) \f$
     */
    void fOil(REAL &fOil, REAL Pw, REAL Sw, REAL &dfOilDPw, REAL &dfOilDSw);
    
    /**
     * @brief Fractional water flux.
     * \f$ f_{Water} = f_{Water}( pw , Sw ) \f$
     */
    void fWater(REAL &fWater, REAL Pw, REAL Sw, REAL &dfWaterDPw, REAL &dfWaterDSw);
    

    // Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////
    
    /** @brief Capilar pressure. \f$ pc = pc( Sw ) \f$ */
    void CapillaryPressure(BFadREAL So, BFadREAL &pc);  
    
    /**
     * @brief Oil relative permeability.
     * \f$ Kro = Kro( Sw ) \f$
     */
    void Kro(BFadREAL Sw, BFadREAL &Kro);
    
    /**
     * @brief Water relative permeability.
     * \f$ Krw = Krw( Sw ) \f$
     */
    void Krw(BFadREAL Sw, BFadREAL &Krw);
    
    /** 
     * @brief Rock porosity. \f$ Phi = Phi( p ) \f$
     * @param po Refrence pressure
     */ 
    void Porosity(BFadREAL po, BFadREAL &poros);    
    
    /** 
     * @brief \f$ Oil density RhoOil = RhoOil( po ) \f$
     * @param po Refrence pressure
     */
    void RhoOil(BFadREAL po, BFadREAL &RhoOil); 
    
    /** 
     * @brief \f$ Water density RhoWater = RhoWater( pw ) \f$
     * @param pw Refrence pressure
     */
    void RhoWater(BFadREAL pw, BFadREAL &RhoWater);     
    
    /** 
     * @brief Oil viscosity. \f$ OilViscosity = ViscOleo( po ) \f$
     * @param po Refrence pressure
     */
    void OilViscosity(BFadREAL po, BFadREAL &OilViscosity); 
    
    /** 
     * @brief Water viscosity. \f$ WaterViscosity = WaterViscosity( pw ) \f$
     * @param po Refrence pressure
     */ 
    void WaterViscosity(BFadREAL po, BFadREAL &WaterViscosity); 
    
    /**
     * @brief Oil mobility.
     * \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$
     */
    void OilLabmda(BFadREAL OilLabmda, BFadREAL Po, BFadREAL &Sw);  
    
    /**
     * @brief Water mobility.
     * \f$ \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$
     */
    void WaterLabmda(BFadREAL WaterLabmda, BFadREAL Pw, BFadREAL &Sw);  
    
    
    /**
     * @brief Bulk mobility.
     * \f$ \lambda = \lambda( pw , Sw ) \f$
     */
    void Labmda(BFadREAL Labmda, BFadREAL Pw, BFadREAL &Sw);    
    
    /**
     * @brief Fractional oil flux.
     * \f$ f_{Oil} = f_{Oil}( po , Sw ) \f$
     */
    void fOil(BFadREAL fOil, BFadREAL Pw, BFadREAL &Sw);    
    
    
    /**
     * @brief Fractional water flux.
     * \f$ f_{Water} = f_{Water}( pw , Sw ) \f$
     */
    void fWater(BFadREAL fWater, BFadREAL Pw, BFadREAL &Sw);        
    

    // Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////  
    

    /** @brief Set characteristic length - m */
    void SetLreference(REAL &Lref){ fLref = Lref;}
    
    /** @brief Set permeability reference - m2 */
    void SetKreference(REAL &Kref){ fKref = Kref;}
    
    /** @brief Set pressure reference - Pa */
    void SetPreference(REAL &Pref){ fPref = Pref;}
    
    /** @brief Set flux reference - Pa */
    void SetQreference(REAL &Qref){ fQref = Qref;}     
    
    /** @brief Set density reference - kg/m3 */
    void SetRhoSCreference(REAL &Densityref){ fRhoref = Densityref;}
    
    /** @brief Set viscosity reference - Pa s */
    void SetEtaSCreference(REAL &Viscosityref){ fEtaref = Viscosityref;}
    
    /** @brief Oil density on standard conditions - kg/m3 */
    REAL RhoOilSC();
    
    /** @brief Water density on standard conditions - kg/m3 */
    REAL RhoWaterSC();
    
    /** @brief Gravity */
    TPZFMatrix<REAL> Gravity();
    
    /** @brief Absolute permeability. */
    void K(TPZFMatrix<REAL> &Kab);
    
    /** @brief Absolute permeability inverse. */
    TPZFMatrix<REAL>  Kinv(TPZFMatrix<REAL> &Kab);  
    
    /** @brief Absolute permeability. */
    void LoadKMap(std::string MaptoRead);   
    
    /** @} */
    public:
virtual int ClassId() const override;

};

#endif
