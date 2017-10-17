/*
 *  TPZAxiSymmetricDarcyFlow.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/11/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "pzdiscgal.h"
#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"

#include "tpzautopointer.h"
#include "SimulationData.h"
#include "ReservoirData.h"
#include "PetroPhysicData.h"
#include "Phase.h"


#ifndef TPZDARCYFLOW
#define TPZDARCYFLOW



class TPZAxiSymmetricDarcyFlow : public TPZDiscontinuousGalerkin {
    
private:
    
    REAL fepsilon;
    REAL fSalpha_max;
    
    TPZAutoPointer<SimulationData> fSimulationData;
    TPZAutoPointer<ReservoirData> fReservoirdata;
    TPZAutoPointer<PetroPhysicData> fPetrophysicdata;
    TPZAutoPointer<Phase> fluid_alpha;
    TPZAutoPointer<Phase> fluid_beta;
    TPZAutoPointer<Phase> fluid_gamma;
    
    // State variables used for weighted fluid blackoil formulation
    int fnstate_vars;
    
    TPZManVector<REAL,3> fBulkVelocity;
    REAL fAveragePressure;
    REAL fWaterSaturation;
    REAL fOilSaturation; //Alpha Saturation
    
    // Auxiliar variables require for the total bulk velocity -  weighted pressure black oil formulation
    
    
     // Velocity of the phase water, oil and gas
    TPZManVector<REAL,3> fWaterVelocity;
    TPZManVector<REAL,3> fOilVelocity;
    TPZManVector<REAL,3> fGasVelocity;
    
    // Pressure of the phase water, oil and gas
    TPZManVector<REAL,4> fWaterPressure;
    TPZManVector<REAL,4> fOilPressure;
    TPZManVector<REAL,4> fGasPressure;
    
    // Density of the phase water, oil and gas
    TPZManVector<REAL,4> fWaterDensity;
    TPZManVector<REAL,4> fOilDensity;
    TPZManVector<REAL,4> fGasDensity;

    // Fractional flow of the phase water, oil and gas
    TPZManVector<REAL,4> fFWater;
    TPZManVector<REAL,4> fFOil;
    TPZManVector<REAL,4> fFGas;
    
    // Mobility of the phase water, oil and gas
    TPZManVector<REAL,4> fWaterMobility;
    TPZManVector<REAL,4> fOilMobility;
    TPZManVector<REAL,4> fGasMobility;
    
    // Total mobility
    TPZManVector<REAL,4> fTotalMobility;
    
    // Total density (As function of the saturation)
    TPZManVector<REAL,4> fTotalDensity;
    
public:
    
    
    /**
     * Empty Constructor
     */
    TPZAxiSymmetricDarcyFlow();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZAxiSymmetricDarcyFlow(int matid);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZAxiSymmetricDarcyFlow(const TPZAxiSymmetricDarcyFlow &mat);
    
    /**
     * Destructor
     */
    ~TPZAxiSymmetricDarcyFlow();
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec);
    
    /** returns the name of the material */
    std::string Name() {
        return "TPZAxiSymmetricDarcyFlow";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const {return 2;}
    
    /** returns the number of state variables associated with the material */
    int NStateVariables() {return 1;} // for hdiv are 3
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout);
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name);
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var);
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){DebugStop();}
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){DebugStop();}
    
    
    // Contribute Methods being used
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    //   Contribute for Darcy system
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterfaceDarcy(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    
    //   Contribute for Transport of alpha phase
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    void ContributeAlpha(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeInterfaceAlpha(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef);
    
    
    /**
     * Unique identifier for serialization purposes
     */
    virtual int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    virtual void Write(TPZStream &buf, int withclassid) const;
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context);
    
    
    /**
     * Set the simulation data
     */
    void SetSimulationData(TPZAutoPointer<SimulationData> simulationdata){fSimulationData = simulationdata;}
    
    /**
     * Get the simulation data
     */
    TPZAutoPointer<SimulationData> GetSimulationData() {return fSimulationData;}
    
    /**
     * Set the reservoir data
     */
    void SetReservoirData(TPZAutoPointer<ReservoirData> ResData){fReservoirdata = ResData;}
    
    /**
     * Get the reservoir data
     */
    TPZAutoPointer<ReservoirData> GetReservoirData() {return fReservoirdata;}
    
    /**
     * Set the Petrophysics data
     */
    void SetPetroPhysicsData(TPZAutoPointer<PetroPhysicData> Petrophysicdata){fPetrophysicdata = Petrophysicdata;}
    
    /**
     * Get the Petrophysics data
     */
    TPZAutoPointer<PetroPhysicData> SetPetroPhysicsData() {return fPetrophysicdata;}
    
    /**
     * Set fluid alpha data
     */
    void SetFluidAlpha(TPZAutoPointer<Phase> Fluidmodeldata){fluid_alpha = Fluidmodeldata;}
    
    /**
     * Get fluid alpha data
     */
    TPZAutoPointer<Phase> GetFluidAlpha() {return fluid_alpha;}
    
    /**
     * Set fluid beta data
     */
    void SetFluidBeta(TPZAutoPointer<Phase> Fluidmodeldata){fluid_beta = Fluidmodeldata;}
    
    /**
     * Get fluid beta data
     */
    TPZAutoPointer<Phase> GetFluidBeta() {return fluid_beta;}
    
    /**
     * Set fluid gamma data
     */
    void SetFluidGamma(TPZAutoPointer<Phase> Fluidmodeldata){fluid_gamma = Fluidmodeldata;}
    
    /**
     * Get fluid gamma data
     */
    TPZAutoPointer<Phase> GetFluidGamma() {return fluid_gamma;}
    
    /**
     * Set number of state vars.
     */
    void SetNvars(int nvars){fnstate_vars = nvars;}
    
    /**
     * Get number of state vars.
     */
    int GetNvars() {return fnstate_vars;}
    
    
    /**
     * Compute Lambda coefficient
     */
    void ComputeLambda(TPZManVector<REAL> &lambda,TPZManVector<REAL> state_vars);
    
    /**
     * Compute Rho coefficient
     */
    void ComputeRho(TPZManVector<REAL> &rho,TPZManVector<REAL> state_vars);
    
    /**
     * Compute Rhof coefficient
     */
    void ComputeRhof(TPZManVector<REAL> &rhof,TPZManVector<REAL> state_vars);
    
    // Axuliar methods computed based on the current u, p, Sw and So values.
    
    /**
     * Compute the properties for computed state variables
     */
    void ComputeProperties(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL> > & props);
    
    /**
     * Compute the gravitational segregational fluxes
     */
    void GravitationalSegregation(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft,TPZVec<TPZMaterialData> &datavecright, TPZVec<TPZManVector<REAL> > & GravityFluxes, TPZManVector<REAL> & fstar);
    
    /**
     * Compute the capillary segregational fluxes
     */
    void CapillarySegregation(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft,TPZVec<TPZMaterialData> &datavecright, TPZVec<TPZManVector<REAL> > & GravitiFluxes, TPZManVector<REAL> & fstar);

    /**
     * Compute the gravitational segregational fluxes
     */
    void GravitationalSegregation(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL,3> >  & qg);

    /**
     * Compute the caplillary segregational fluxes
     */
    void CapillarySegregation(TPZVec<TPZMaterialData> &datavec, TPZVec<TPZManVector<REAL,3> > & qc, TPZManVector<REAL,3> & grads);
    
    /**
     * Compute the linearized version of the bubble function f
     */
    void f(REAL P, REAL Salpha, TPZManVector<REAL> & f_value);

    /**
     * Compute the linearized version of the bubble function f
     */
    void fLinear(REAL P, REAL Salpha, TPZManVector<REAL> & f_value);
    
    /**
     * Compute the linearized version of the expelling water
     */
    void fExp(REAL P, REAL Salpha, TPZManVector<REAL> & ExpL);
    
    /**
     * Compute the linearized version of the receive water
     */
    void fRec(REAL P, REAL Salpha, TPZManVector<REAL> & RecL);
    
    /**
     * Compute the linearized version of the expelling water
     */
    void fExpLinear(REAL P, REAL Salpha, TPZManVector<REAL> & ExpL);

    /**
     * Compute the linearized version of the receive water
     */
    void fRecLinear(REAL P, REAL Salpha, TPZManVector<REAL> & RecL);
    
    /**
     * Compute the properties for a given state variables
     */
    void ComputeProperties(TPZManVector<REAL> &state_vars, TPZVec<TPZManVector<REAL> > & props);
    
    void UpdateStateVariables(TPZManVector<REAL,3> u, REAL P, REAL Sw, REAL So);
    
    void PhasePressures();
    
    void PhaseDensities();
    
    void PhaseMobilities();
    
    void PhaseFractionalFlows();
    
    void TotalMobility();
    
    void TotalDensity();
    
    // System Properties
    
    void Rho_alpha(TPZVec<REAL> P_alpha);
    
    void Rho_beta(TPZVec<REAL> P_beta);
    
    void Rho_gamma(TPZVec<REAL> P_gamma);
    
    
};

#endif