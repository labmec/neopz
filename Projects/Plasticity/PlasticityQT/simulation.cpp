#include <iostream>
#include <cstdlib>


#include "simulation.h"

void TPBrPlasticityTest::LoadInputStrainStress ( const QVector<double> &sigmaAxialTotal, const QVector<double> &sigmaLateral,
                                                const QVector<double> &defAxial, const QVector<double> &defLateral) {
    int numlines = 0;
    numlines = sigmaLateral.size();

    fStressRZInput.Resize(numlines, 2);
    fStrainRZInput.Resize(numlines, 2);

    for (int i=0; i< numlines; i++) {
        fStressRZInput(i,1) = - sigmaAxialTotal.value(i);  //-sig_ax_t;
        fStressRZInput(i,0) = - sigmaLateral.value(i);     //-sig_r;
        fStrainRZInput(i,1) = - defAxial.value(i) / 100.;         //-eps_ax/100.;
        fStrainRZInput(i,0) = - defLateral.value(i) / 100.;       //-eps_r/100.;
    }

    qDebug() << "Fim do LoadInputStrainStress..." ;
}
