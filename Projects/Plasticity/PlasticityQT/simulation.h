#ifndef SIMULATION_H
#define SIMULATION_H

#include <QVector>

#include "pzreal.h"
#include "TPZSandlerDimaggio.h"

class TPZPlasticityTest
{
    /// variable indicating if the maximum z stress is given
    bool fZStressKnown;

    /// variable indicating if the maximum r stress is given
    bool fRStressKnown;

    /// stress at which pores are closed
    TPZManVector<STATE,2> fPoreStressRZ;

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
    TPZPlasticityTest();

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

    /// set the stress at which pores are closed
    void SetPoreClosingStress(TPZVec<STATE> &porestress)
    {
        fPoreStressRZ = porestress;
    }

    /// read the input strain and stress from the laboratory file
    void ReadInputStrainStress(const std::string &filename);

    /// set the SandlerDimaggio object
    void SetSandlerDimaggio(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj)
    {
        fSandler = obj;
    }

    /// compute the stress strain curve
    void PerformSimulation();

    void SetUpSimulation ( REAL young, REAL poisson, REAL A, REAL B, REAL C, REAL D, REAL R, REAL W ) {
        this->_young = young;
        this->_poisson = poisson;
        this->_A = A;
        this->_B = B;
        this->_C = C;
        this->_D = D;
        this->_R = R;
        this->_W = W;
    }

    void RunSimulation () {
        // x y solution storage
        x_UCS.clear();
        y_UCS.clear();
        x_Triaxial.clear();
        y_Triaxial.clear();


        //simulation code here ...






    }

    double * get_UCS_Xvalues() {
        return x_UCS.data();
    }

    double * get_UCS_Yvalues() {
        return y_UCS.data();
    }

    double * get_Triaxial_Xvalues() {
        return x_Triaxial.data();
    }

    double * get_Triaxial_Yvalues() {
        return y_Triaxial.data();
    }

protected:

    /// Get the stress to the pore stress
    void ApplyInitialStress();

    /// Evoluate the stress and strain to step ist
    // Be aware that this method modifies the data fSandler
    void EvoluateToStep(TPZVec<STATE> &strainRZ, TPZVec<STATE> &stressRZ);

private:
    int _plasticityTest, _nofsteps, _unloadstep;
    REAL _young, _poisson, _E, _A, _B, _C, _D, _R, _W, _inttol, _epsx, _deltasigmaXX ;
    double _percent;

    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *sandler;

    QVector<double> x_UCS, y_UCS, x_Triaxial, y_Triaxial;


};
























//#include <TPZTensor.h>
//#include <TPZSandlerDimaggio.h>
//#include <TPZPlasticState.h>
//#include "TPZPlasticityTest.h"

//#include <QVector>

//class Simulation
//{
//private:
//    int _plasticityTest, _nofsteps, _unloadstep;
//    REAL _poisson, _E, _A, _B, _C, _D, _R, _W, _inttol, _epsx, _deltasigmaXX ;
//    double _percent;

//    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *sandler;

//    QVector<double> x, y;

//public:
//    Simulation();

//    void setupTriaxial();
//    void setupUCS();

//    double * getXvalues();
//    double * getYvalues();
//    int getXsize();
//};

#endif // SIMULATION_H
