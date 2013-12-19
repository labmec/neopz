#ifndef TPZPLASTICITYSIMULATION
#define TPZPLASTICITYSIMULATION

#include "TPZSandlerDimaggio.h"

class TPZPlasticitySimulation
{
protected:
    /// variable indicating if the maximum z stress is given
    bool fZStressKnown;
    
    /// variable indicating if the maximum r stress is given
    bool fRStressKnown;

    /// index of the initial pore closure stress
    int fPoreClosureIndex;
    
//    /// stress at which pores are closed
//    TPZManVector<STATE,2> fPoreStressRZ;

    /// simulated strain corresponding to porestress
    TPZManVector<STATE,2> fPoreStrainRZ;

    /// stress as read from the laboratory test
    TPZFMatrix<STATE> fStressRZInput;
    
    /// strain as read from the laboratory test
    TPZFMatrix<STATE> fStrainRZInput;
    
    /// stress as computed by the simulation
    TPZFMatrix<STATE> fStressRZSimulated;
    
    /// strain as computed by the simulation
    TPZFMatrix<STATE> fStrainRZSimulated;
    
    /// Sandler DiMaggio 
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> fSandler;
    
    /// Number of steps between both states
    int fNumSteps;
    
public:
    
    /// Default constructor
    TPZPlasticitySimulation();
    
    /// either the stress is determined or the deformation
    void SetZStressKnown(bool zstressknown)
    {
        fZStressKnown = zstressknown;
    }
    
    /// either the stress is determined or the deformation
    void SetRStressKnown(bool rstressknown)
    {
        fRStressKnown = rstressknown;
    }
    
    void SetSimulationInitialStep(int istep)
    {
        if (istep <0 || istep >= fStressRZInput.Rows()) {
            DebugStop();
        }
        fPoreClosureIndex = istep;
//        fPoreStressRZ[0] = fStressRZInput(istep,0);
//        fPoreStressRZ[1] = fStressRZInput(istep,1);
    }
    
    
    /// read the input strain and stress from the laboratory file
    void ReadInputStrainStress(const std::string &filename);
    
    /// read the input strain and stress from vectors
    void ReadInputStrainStress(const TPZVec<REAL> &sigax, const TPZVec<REAL> &epsax, const TPZVec<REAL> &sigr, const TPZVec<REAL> &epsr);
    
    /// set the SandlerDimaggio object
    void SetSandlerDimaggio(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj)
    {
        fSandler = obj;
    }
    
    /// compute the stress strain curve
    void PerformSimulation();
    
    /// get the simulated data
    void GetSimulatedStrainStress(TPZVec<REAL> &sigax, TPZVec<REAL> &epsax, TPZVec<REAL> &sigr, TPZVec<REAL> &epsr);
    
protected:
    
    /// Get the stress to the pore stress
    void ApplyInitialStress();
    
    /// Evoluate the stress and strain to step ist
    // Be aware that this method modifies the data fSandler
    void EvoluateToStep(TPZVec<STATE> &strainRZ, TPZVec<STATE> &stressRZ);
    
    
};

#endif
