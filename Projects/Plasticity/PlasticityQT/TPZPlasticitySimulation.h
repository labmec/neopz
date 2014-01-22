
#ifndef TPZPLASTICITYSIMULATION
#define TPZPLASTICITYSIMULATION

#include "TPZSandlerDimaggio.h"
#include <math.h>

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
		
		/// Compute the position of the ellips which contains the point
		REAL FindL(REAL I1, REAL sqrtJ2, STATE Lini);
    
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

    void SetElasticTransition(int istep) {
        //////////////////////////////////////////
        REAL sigrc = fStressRZInput(fPoreClosureIndex,0);
				REAL sigaxc = fStressRZInput(fPoreClosureIndex,1);
				REAL epsrc = fStrainRZInput(fPoreClosureIndex,0);
				REAL epsaxc = fStrainRZInput(fPoreClosureIndex,1);
				REAL I1C = 2.*sigrc+sigaxc;
        REAL sigr = fStressRZInput(istep,0);
				REAL sigax = fStressRZInput(istep,1);
				REAL epsr = fStrainRZInput(istep,0);
				REAL epsax = fStrainRZInput(istep,1);
				
				REAL DSigR = sigr-sigrc;
				REAL DEpsR = epsr - epsrc;
				REAL DSigAx = sigax-sigaxc;
				REAL DEpsAx = epsax-epsaxc;
				
				REAL elast = ((DSigAx - DSigR)*(DSigAx + 2*DSigR))/(-2*DEpsR*DSigR + DEpsAx*(DSigAx + DSigR));
				REAL nu = (-(DEpsR*DSigAx) + DEpsAx*DSigR)/(-2*DEpsR*DSigR + DEpsAx*(DSigAx + DSigR));
				
//				fSandler.fER.SetUp(elast,nu);
				
				REAL I1 = 2.*sigr+sigax;
				static const REAL sqrt3 = sqrt(3.);
				REAL sqrtJ2 = fabs(sigr-sigax)/sqrt3;
				REAL L = FindL(I1,sqrtJ2,0.);
				
				for(int i=fPoreClosureIndex; i<istep; i++)
				{
					REAL sigr = fStressRZInput(i,0);
					REAL sigax = fStressRZInput(i,1);
					REAL epsr = fStrainRZInput(i,0);
					REAL epsax = fStrainRZInput(i,1);
					REAL I1 = 2.*sigr+sigax;
					REAL sqrtJ2 = fabs(sigr-sigax)/sqrt3;
					REAL Lint = FindL(I1,sqrtJ2,L);
					if(Lint < L) 
					{
						L = Lint;						
					}
				}
				TPZElasticResponse ER;
				ER.SetUp(elast,nu);
				for(int i=fPoreClosureIndex; i<istep; i++)
				{
					REAL epsr = fStrainRZInput(i,0)-fStrainRZInput(fPoreClosureIndex,0);
					REAL epsax = fStrainRZInput(i,1)-fStrainRZInput(fPoreClosureIndex,1);
					TPZTensor<STATE> eps,sig;
					eps.XX() = epsr;
					eps.YY() = epsr;
					eps.ZZ() = epsax;
					ER.Compute<STATE>(eps,sig);
					sig.XX() += fStressRZInput(fPoreClosureIndex,0);
					sig.YY() += fStressRZInput(fPoreClosureIndex,0);
					sig.ZZ() += fStressRZInput(fPoreClosureIndex,1);
					REAL I1 = sig.I1();
					REAL sqrtJ2 = sqrt(sig.J2());
					REAL Lint = FindL(I1,sqrtJ2,L);
					if(Lint < L) 
					{
						L = Lint;						
					}

				}
				
				/// atualizar o objeto fsandler
				TPZPlasticState<STATE> locstate = fSandler.GetState();
				locstate.fAlpha = L;
				fSandler.SetState(locstate);
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
