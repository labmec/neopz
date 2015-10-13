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
#include "ReducedPVT.h"


#ifndef TPZDARCYFLOW
#define TPZDARCYFLOW



class TPZAxiSymmetricDarcyFlow : public TPZDiscontinuousGalerkin {
    
private:
    
    REAL fepsilon;
    
    TPZAutoPointer<SimulationData> fSimulationData;
    TPZAutoPointer<ReservoirData> fReservoirdata;
    TPZAutoPointer<PetroPhysicData> fPetrophysicdata;
    TPZAutoPointer<ReducedPVT> fluid_alpha;
    TPZAutoPointer<ReducedPVT> fluid_beta;
    TPZAutoPointer<ReducedPVT> fluid_gamma;
    
    // State variables used for weighted fluid blackoil formulation
    int fnstate_vars;
    
    TPZManVector<REAL,3> fBulkVelocity;
    REAL fAveragePressure;
    REAL fWaterSaturation;
    REAL fOilSaturation;
    
    // Auxiliar variables require for the total bulk velocity -  weighted pressure black oil formulation
    
    TPZManVector<REAL,3> fWaterVelocity;
    TPZManVector<REAL,3> fOilVelocity;
    TPZManVector<REAL,3> fGasVelocity;
    
    TPZManVector<REAL,4> fWaterPressure;
    TPZManVector<REAL,4> fOilPressure;
    TPZManVector<REAL,4> fGasPressure;
    
    TPZManVector<REAL,4> fWaterDensity;
    TPZManVector<REAL,4> fOilDensity;
    TPZManVector<REAL,4> fGasDensity;
    
    TPZManVector<REAL,4> fFWater;
    TPZManVector<REAL,4> fFOil;
    TPZManVector<REAL,4> fFGas;
    
    TPZManVector<REAL,4> fWaterMobility;
    TPZManVector<REAL,4> fOilMobility;
    TPZManVector<REAL,4> fGasMobility;
    
    TPZManVector<REAL,4> fTotalMobility;
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
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * It computes a contribution to the load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef);
    
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
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    void ContributeBCDarcy(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * Unique identifier for serialization purposes
     */
    int ClassId() const;
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid);
    
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
     * Set the simulation data
     */
    void SetReservoirData(TPZAutoPointer<ReservoirData> ResData){fReservoirdata = ResData;}
    
    /**
     * Get the simulation data
     */
    TPZAutoPointer<ReservoirData> GetReservoirData() {return fReservoirdata;}
    
    /**
     * Set the simulation data
     */
    void SetPetroPhysicsData(TPZAutoPointer<PetroPhysicData> Petrophysicdata){fPetrophysicdata = Petrophysicdata;}
    
    /**
     * Get the simulation data
     */
    TPZAutoPointer<PetroPhysicData> SetPetroPhysicsData() {return fPetrophysicdata;}
    
    /**
     * Set fluid alpha data
     */
    void SetFluidAlpha(TPZAutoPointer<ReducedPVT> Fluidmodeldata){fluid_alpha = Fluidmodeldata;}
    
    /**
     * Get fluid alpha data
     */
    TPZAutoPointer<ReducedPVT> GetFluidAlpha() {return fluid_alpha;}
    
    /**
     * Set fluid beta data
     */
    void SetFluidBeta(TPZAutoPointer<ReducedPVT> Fluidmodeldata){fluid_beta = Fluidmodeldata;}
    
    /**
     * Get fluid beta data
     */
    TPZAutoPointer<ReducedPVT> GetFluidBeta() {return fluid_beta;}
    
    /**
     * Set fluid gamma data
     */
    void SetFluidGamma(TPZAutoPointer<ReducedPVT> Fluidmodeldata){fluid_gamma = Fluidmodeldata;}
    
    /**
     * Get fluid gamma data
     */
    TPZAutoPointer<ReducedPVT> GetFluidGamma() {return fluid_gamma;}
    
    /**
     * Set state vars.
     */
    void SetNvars(int nvars){fnstate_vars = nvars;}
    
    /**
     * Get state vars.
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