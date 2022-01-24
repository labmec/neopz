
//
//  TPZSBFemVolumeL2.cpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#include "TPZSBFemVolumeL2.h"
#include "pzgeoelside.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "TPZMatSingleSpace.h"
#include "TPZMatErrorSingleSpace.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "pzcmesh.h"
#include "TPZGeoLinear.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

// #ifdef LOG4CXX
// static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
// #endif

TPZSBFemVolumeL2::TPZSBFemVolumeL2(TPZCompMesh & mesh, TPZGeoEl * gel) : TPZSBFemVolume(mesh, gel)
{
    fElementVec1D.Resize(3);
}

void TPZSBFemVolumeL2::LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd)
{
    fCoeficients = coef;
    fCoeficientsD = coefd;
}

void TPZSBFemVolumeL2::SetElementGroupIndex(int64_t index)
{
    fElementGroupIndex = index;
    TPZCompEl *celgr = Mesh()->Element(index);
    fElementGroup = celgr;
}

/**
 * @brief Computes solution and its derivatives in the local coordinate qsi.
 * @param qsi master element coordinate
 * @param sol finite element solution
 * @param dsol solution derivatives
 * @param axes axes associated with the derivative of the solution
 */
void TPZSBFemVolumeL2::ReallyComputeSolution(TPZMaterialDataT<STATE> & data)
{
    TPZVec<REAL> &qsi = data.xParametric;
    TPZSolVec<STATE> &sol = data.sol;
    TPZGradSolVec<STATE> &dsol = data.dsol;
    TPZFMatrix<REAL> &axes = data.axes;

    TPZCompMesh *cmesh = this->Mesh();
    data.sol.Resize(fCoeficients.Cols());
    data.dsol.Resize(fCoeficients.Cols());

    TPZGeoEl *Ref2D = this->Reference();
    int matid = Ref2D->MaterialId();
    TPZMaterial *mat2d = cmesh->FindMaterial(matid);

    int dim = Ref2D->Dimension();
    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
    if (sbfemparam < 0.) {
        std::cout << "sbfemparam " << sbfemparam << std::endl;
        sbfemparam = 0.;
    }
    if (IsZero(sbfemparam)) {
        for (int i = 0; i < dim - 1; i++) {
            qsi[i] = 0.;
        }
        if (dim == 2) {
            sbfemparam = 1.e-6;
            qsi[dim - 1] = 1. - 2.e-6;
        } else {
            sbfemparam = 1.e-4;
            qsi[dim - 1] = 1. - 2.e-4;
        }
    }

    auto CSkeleton = dynamic_cast<TPZInterpolationSpace *> (fElementVec1D[2]);

    TPZMaterialDataT<STATE> data1d, data2d;
    // compute the lower dimensional shape functions
    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);

    TPZGeoEl *Ref1D = CSkeleton->Reference();

    CSkeleton->InitMaterialData(data1d);
    CSkeleton->ComputeRequiredData(data1d, qsilow);

    Ref1D->Jacobian(qsilow, data1d.jacobian, data1d.axes, data1d.detjac, data1d.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    axes = data2d.axes;

    int nshape = data1d.fPhi.Rows();
    int nstate = mat2d->NStateVariables();
#ifdef PZDEBUG
    if (fPhi.Cols() != fCoeficients.Rows()) {
        DebugStop();
    }
#endif

    for (int s = 0; s < sol.size(); s++) {
        TPZManVector<std::complex<double>, 10> uh_xi(fPhi.Rows(), 0.), Duh_xi(fPhi.Rows(), 0.);
        int nphixi = fPhi.Rows();
        int numeig = fPhi.Cols();
        for (int c = 0; c < numeig; c++) {
            std::complex<double> xiexp;
            std::complex<double> xiexpm1;
            if (IsZero(fEigenvalues[c] + 0.5 * (dim - 2))) {
                xiexp = 1;
                xiexpm1 = 0;
            } else if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = sbfemparam;
                xiexpm1 = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] - 0.5 * (dim - 2));
                xiexpm1 = pow(sbfemparam, -fEigenvalues[c] - 1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixi; i++) {
                uh_xi[i] += fCoeficients(c, s) * xiexp * fPhi(i, c);
                Duh_xi[i] += -fCoeficients(c, s)*(fEigenvalues[c] + 0.5 * (dim - 2)) * xiexpm1 * fPhi(i, c);
            }
        }
        
#ifdef LOG4CXX2
        if (s == 0 && logger->isDebugEnabled()) {
            std::stringstream sout;
            sout << "uh_xi " << uh_xi << std::endl;
            sout << "Duh_xi " << Duh_xi << std::endl;
            data1d.fPhi.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        data.sol[s].Resize(nstate);
        data.sol[s].Fill(0.);
        TPZFNMatrix<9, STATE> dsollow(dim - 1, nstate, 0.), dsolxieta(dim, nstate, 0.);
        TPZManVector<STATE, 3> dsolxi(nstate, 0.);
        for (int ishape = 0; ishape < nshape; ishape++) {
            for (int istate = 0; istate < nstate; istate++) {
                sol[s][istate] += data1d.fPhi(ishape) * uh_xi[ishape * nstate + istate].real();
                dsolxi[istate] += data1d.fPhi(ishape) * Duh_xi[ishape * nstate + istate].real();
                for (int d = 0; d < dim - 1; d++) {
                    dsollow(d, istate) += data1d.fDPhi(d, ishape) * uh_xi[ishape * nstate + istate].real();
                }
            }
        }
        for (int istate = 0; istate < nstate; istate++) {
            for (int d = 0; d < dim - 1; d++) {
                dsolxieta(d, istate) = dsollow(d, istate);
            }
            dsolxieta(dim - 1, istate) = -dsolxi[istate] / 2.;
        }
        data.dsol[s].Resize(dim, nstate);
        data.dsol[s].Zero();
        for (int istate = 0; istate < nstate; istate++) {
            for (int d1 = 0; d1 < dim; d1++) {
                for (int d2 = 0; d2 < dim; d2++) {
                    data.dsol[s](d1, istate) += data2d.jacinv(d2, d1) * dsolxieta(d2, istate);
                }
            }
        }
    }
}

#include "pzaxestools.h"

TPZCompEl * CreateSBFemPressureCompEl(TPZCompMesh &mesh, TPZGeoEl *gel)
{
    return new TPZSBFemVolumeL2(mesh, gel);
}
