#ifndef SIMULATION_H
#define SIMULATION_H

#include <QVector>
#include <QDebug>

#include "pzreal.h"
#include "TPZSandlerDimaggio.h"

#include "../TPZPlasticityTest.h"

class TPBrPlasticityTest : public TPZPlasticityTest
{

protected:
    int _plasticityTest, _nofsteps, _unloadstep;
    REAL _young, _poisson, _A, _B, _C, _D, _R, _W, _inttol, _epsx, _deltasigmaXX ;
    double _percent;
    QVector<double> Stress_0, Stress_1, Strain_0, Strain_1;

public:
    /// load the input strain and stress from GUI interface
    void LoadInputStrainStress ( const QVector<double> &sigmaAxialTotal, const QVector<double> &sigmaLateral,
                                 const QVector<double> &defAxial, const QVector<double> &defLateral );


    void SetUpSimulation ( REAL poisson, REAL young, REAL A, REAL B, REAL C, REAL R, REAL D, REAL W ) {
        this->_young = young;
        this->_poisson = poisson;
        this->_A = A;
        this->_B = B;
        this->_C = C;
        this->_D = D;
        this->_R = R;
        this->_W = W;
        this->fSandler.SetUp(_poisson, _young, _A, _B, _C, _R, _D, _W);
    }

    void PrintResults () {
        std::cout << "Results: fStressRZSimulated= " << fStressRZSimulated << " fStrainRZSimulated= " << fStrainRZSimulated << std::endl;
    }

    const QVector<double> get_Stress_0 () {
        for (int i=0; i<fStressRZSimulated.Rows(); i++) {
            Stress_0.insert(i, fStressRZSimulated.GetVal(i,0) * -1 );
        }
        return Stress_0;
    }
    const QVector<double> get_Stress_1 () {
        for (int i=0; i<fStressRZSimulated.Rows(); i++) {
            Stress_1.insert(i, fStressRZSimulated.GetVal(i,1) * -1 );
        }
        return Stress_1;
    }
    const QVector<double> get_Strain_0 () {
        for (int i=0; i<fStrainRZSimulated.Rows(); i++) {
            Strain_0.insert(i, fStrainRZSimulated.GetVal(i,0) * -100 );
        }
        return Strain_0;
    }
    const QVector<double> get_Strain_1 () {
        for (int i=0; i<fStrainRZSimulated.Rows(); i++) {
            Strain_1.insert(i, fStrainRZSimulated.GetVal(i,1) * -100 );
        }
        return Strain_1;
    }
};

#endif // SIMULATION_H
