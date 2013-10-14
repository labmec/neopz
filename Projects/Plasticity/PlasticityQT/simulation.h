#ifndef SIMULATION_H
#define SIMULATION_H

#include <QVector>
#include <QDebug>

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

    /// load the input strain and stress from GUI interface
    void LoadInputStrainStress (QVector<double> *sigmaAxialTotal, QVector<double> *sigmaLateral,
                                QVector<double> *defAxial, QVector<double> *defLateral) ;

    /// set the SandlerDimaggio object
    void SetSandlerDimaggio(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &obj)
    {
        fSandler = obj;
    }

    /// compute the stress strain curve
    void PerformSimulation();

    void SetUpSimulation ( REAL poisson, REAL young, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W ) {
        this->_young = young;
        this->_poisson = poisson;
        this->_A = A;
        this->_B = B;
        this->_C = C;
        this->_D = D;
        this->_R = R;
        this->_W = W;

//        qDebug() << "VALORES DA TELA 2.0!!!" <<poisson <<young <<A <<B <<C <<R <<D <<W;

        this->fSandler.SetUp(_poisson, _young, _A, _B, _C, _R, _D, _W);
    }

    void PrintResults () {
        std::cout << "Results: fStressRZSimulated= " << fStressRZSimulated << " fStrainRZSimulated= " << fStrainRZSimulated << std::endl;
    }

    const QVector<double> get_Stress_X () {
        for (int i=0; i<fStressRZSimulated.Rows(); i++) {
            Stress_X.insert(i, fStressRZSimulated.GetVal(i,0) * -1 );
        }
        return Stress_X;
    }
    const QVector<double> get_Stress_Y () {
        for (int i=0; i<fStressRZSimulated.Rows(); i++) {
            Stress_Y.insert(i, fStressRZSimulated.GetVal(i,1) * -1 );
        }
        return Stress_Y;
    }
    const QVector<double> get_Strain_X () {
        for (int i=0; i<fStrainRZSimulated.Rows(); i++) {
            Strain_X.insert(i, fStrainRZSimulated.GetVal(i,0) * -100 );
        }
        return Strain_X;
    }
    const QVector<double> get_Strain_Y () {
        for (int i=0; i<fStrainRZSimulated.Rows(); i++) {
            Strain_Y.insert(i, fStrainRZSimulated.GetVal(i,1) * -100 );
        }
        return Strain_Y;
    }


protected:

    /// Get the stress to the pore stress
    void ApplyInitialStress();

    /// Evoluate the stress and strain to step ist
    // Be aware that this method modifies the data fSandler
    void EvoluateToStep(TPZVec<STATE> &strainRZ, TPZVec<STATE> &stressRZ);

private:
    int _plasticityTest, _nofsteps, _unloadstep;
    REAL _young, _poisson, _A, _B, _C, _D, _R, _W, _inttol, _epsx, _deltasigmaXX ;
    double _percent;
    QVector<double> Stress_X, Stress_Y, Strain_X, Strain_Y;


};

#endif // SIMULATION_H
