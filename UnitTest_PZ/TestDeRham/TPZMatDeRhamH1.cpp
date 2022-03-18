#include "TPZMatDeRhamH1.h"
#include "TPZMaterialDataT.h"
#include "pzaxestools.h"

void
TPZMatDeRhamH1::Contribute(const TPZMaterialDataT<STATE> &data, REAL weight,
                           TPZFMatrix<STATE> &ek,
                           TPZFMatrix<STATE> &ef)
{
    const TPZFMatrix<REAL> &phi = data.phi;
    const TPZFMatrix<REAL> &gradPhiaxes = data.dphix;
    TPZFNMatrix<3, REAL> gradPhi(3, phi.Rows(), 0.);
    TPZAxesTools<REAL>::Axes2XYZ(gradPhiaxes, gradPhi, data.axes);
    const int nFuncs = phi.Rows();
    for (int i = 0; i < nFuncs; i++) {
      for (int j = 0; j < nFuncs; j++) {
        for (int x = 0; x < 3; x++)
          ek(i, j) += gradPhi(x, i) * gradPhi(x, j) * weight;
      }
    }
}