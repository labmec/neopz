#ifndef SIMULATION_H
#define SIMULATION_H

#include "pzreal.h"

#include <TPZTensor.h>
#include <TPZSandlerDimaggio.h>
#include <TPZPlasticState.h>
#include "TPZPlasticityTest.h"

#include <QVector>

class Simulation
{
private:
    int _plasticityTest, _nofsteps, _unloadstep;
    REAL _poisson, _E, _A, _B, _C, _D, _R, _W, _inttol, _epsx, _deltasigmaXX ;
    double _percent;

    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *sandler;

    QVector<double> x, y;

public:
    Simulation();

    void setupTriaxial();
    void setupUCS();

    double * getXvalues();
    double * getYvalues();
    int getXsize();
};

#endif // SIMULATION_H
