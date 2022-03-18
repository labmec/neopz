
//
//  TPZSBFemVolumeHdiv.cpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#include "TPZSBFemVolumeHdiv.h"
#include "pzgeoelside.h"
#include "TPZSBFemElementGroup.h"
#include "pzintel.h"
#include "TPZMatCombinedSpaces.h"
#include "TPZMatErrorCombinedSpaces.h"
#include "pzelmat.h"
#include "TPZBndCondT.h"
#include "pzcmesh.h"
#include "TPZGeoLinear.h"
#include "pzmultiphysicscompel.h"
#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

// #ifdef LOG4CXX
// static LoggerPtr logger(Logger::getLogger("pz.mesh.sbfemvolume"));
// #endif

TPZSBFemVolumeHdiv::TPZSBFemVolumeHdiv(TPZCompMesh & mesh, TPZGeoEl * gel) : TPZInterpolationSpace(mesh, gel), fElementGroupIndex(-1), fSkeleton(-1)
{
    fElementVec1D.Resize(3);
}

TPZCompEl * TPZSBFemVolumeHdiv::Clone(TPZCompMesh &mesh) const
{
    // till I remember how this works
    DebugStop();
    return 0;
}

void TPZSBFemVolumeHdiv::AddElement1D(TPZCompEl * cel, int localindex)
{
    fElementVec1D[localindex] = cel;
    auto ncon = fConnectIndexes.size();
    auto nconcel = cel->NConnects();
    
    fConnectIndexes.Resize(ncon+nconcel);
    for (auto i = 0; i < nconcel; i++)
    {
        fConnectIndexes[i+ncon] = cel->ConnectIndex(i);
    }
}

void TPZSBFemVolumeHdiv::SetPhiFlux(TPZFMatrix<std::complex<double>> &phiflux, TPZManVector<std::complex<double>> &eigval)
{        
    fEigenvalues = eigval;
    int nrow = fLocalIndicesFluxInt.size();
    fPhiFluxInt.Resize(nrow, phiflux.Cols());
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < phiflux.Cols(); j++)
        {
            fPhiFluxInt(i, j) = phiflux(fLocalIndicesFluxInt[i], j);
        }
    }    

    nrow = fLocalIndicesFlux.size();
    fPhiFlux.Resize(nrow, phiflux.Cols());
    for (int i = 0; i < nrow; i++)
    {
        for (int j = 0; j < phiflux.Cols(); j++)
        {
            fPhiFlux(i, j) = phiflux(fLocalIndicesFlux[i], j);
        }
    }
}

void TPZSBFemVolumeHdiv::SetLocalIndicesFlux(TPZManVector<int64_t> &localindicesfluxint, TPZManVector<int64_t> &localindicesflux)
{
    fLocalIndicesFluxInt = localindicesfluxint;
    fLocalIndicesFlux = localindicesflux;
}

void TPZSBFemVolumeHdiv::LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd)
{
    fCoeficients = coef;
    fCoeficientsD = coefd;
}

void TPZSBFemVolumeHdiv::ReallyComputeSolution(TPZMaterialDataT<STATE> & data)
{
    TPZCompMesh *cmesh = this->Mesh();
    data.sol.Resize(fCoeficients.Cols());
    data.dsol.Resize(fCoeficients.Cols());

    TPZManVector<REAL> qsi = data.xParametric;

    TPZGeoEl *Ref2D = this->Reference();
    TPZMaterial *mat2d = cmesh->FindMaterial(Ref2D->MaterialId());
    int dim = Ref2D->Dimension();

    REAL sbfemparam = (1. - qsi[dim - 1]) / 2.;
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
    
    auto CFlux = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> *> (fElementVec1D[1]);
#ifdef PZDEBUG
    if (!CFlux) DebugStop();
#endif
    
    TPZMaterialDataT<STATE> data2d;

    TPZManVector<REAL, 3> qsilow(qsi);
    qsilow.Resize(dim - 1);

    TPZGeoEl *Ref1D = CFlux->Reference();

    Ref1D->Jacobian(qsilow, data.jacobian, data.axes, data.detjac, data.jacinv);
    Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
    auto axes = data2d.axes;

    CFlux->ComputeRequiredData(data, qsilow);

    // TPZVec<TPZTransform<> > trvec;
    // CFlux->AffineTransform(trvec);

    // TPZManVector<REAL> locpt(CFlux->Reference()->Dimension());
    // trvec[CFlux->Index()].Apply(qsilow, locpt);
    // ComputeRequiredData(data,locpt);

    int nshape = fPhiFlux.Rows() + fPhiFluxInt.Rows();
    int nstate = mat2d->NStateVariables();
    
    data.divsol.Resize(data.sol.size());
    
    for (int s = 0; s < data.sol.size(); s++)
    {
        TPZManVector<std::complex<double>, 10> uh_xi(nshape, 0.), Duh_xi(nshape, 0.);
        int nphixiint = fPhiFluxInt.Rows();
        int nphixiext = fPhiFlux.Rows();
        int numeig = fEigenvalues.size();
        std::complex<double> xiexp;
        for (int c = 0; c < numeig; c++) {
            if (IsZero(fEigenvalues[c] + 1. + 0.5 * (dim - 2))) {
                xiexp = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c] -1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixiint; i++) {
                uh_xi[i] += fPhiFluxInt(i, c) * xiexp * fCoeficientsD(c, s);
            }
            for (int i = 0; i < nphixiext; i++) {
                uh_xi[i+nphixiint] +=  fPhiFlux(i, c) * xiexp * fCoeficientsD(c, s);    
            }
        }
        for (int c = numeig; c < 2*numeig; c++) {
            if (IsZero(fEigenvalues[c-numeig] + 1. + 0.5 * (dim - 2))) { 
                xiexp = 1;
            } else {
                xiexp = pow(sbfemparam, -fEigenvalues[c-numeig] -1. - 0.5 * (dim - 2));
            }
            for (int i = 0; i < nphixiint; i++) {
                uh_xi[i] += fPhiFluxInt(i, c) * xiexp * fCoeficients(c-numeig, s);
            }
            for (int i = 0; i < nphixiext; i++) {
                uh_xi[i+nphixiint] += fPhiFlux(i, c) * xiexp * fCoeficients(c-numeig, s);
            }
        }

        data.sol[s].Resize(3);
        data.sol[s].Fill(0.);

        auto multiphysicsel = dynamic_cast<TPZMultiphysicsElement * >(CFlux);
        auto collapsedel = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear> * >(fElementVec1D[1]);

        int ishape = 0;
        for (int ic = 0; ic < collapsedel->NConnects()-1; ic++)
        {
            auto nshape = collapsedel->NConnectShapeF(ic, collapsedel->GetPreferredOrder());
            for (int is = 0; is < nshape; is++)
            {
                for (int istate = 0; istate < nstate; istate++)
                {
                    for (int d = 0; d < dim; d++)
                    {
                        int id = (ishape + is)* nstate + istate;
                        data.sol[s][d] += data.fPhi(id) * data.fDeformedDirections(d,ic) * uh_xi[id].real();
                    }
                }
            }
            ishape += nshape;
        }
        {
            data.dsol[s].Redim(dim*nstate, dim);
            data.divsol[s].Resize(nstate);
            data.divsol[s].Fill(0.);
        }
    }
}

int TPZSBFemVolumeHdiv::NShapeF() const
{
    int nc = fElementVec1D[1]->NConnects();
    int nshape = 0;
    for (int ic=0; ic<nc; ic++) {
        TPZConnect &c = Connect(ic);
        nshape += c.NShape();
    }
    return nshape;
}

int TPZSBFemVolumeHdiv::NConnectShapeF(int icon, int order) const
{
    DebugStop();
    return 0;
}

TPZIntPoints & TPZSBFemVolumeHdiv::GetIntegrationRule()
{
    if(!fIntRule) InitializeIntegrationRule();
    return *fIntRule;
}

const TPZIntPoints & TPZSBFemVolumeHdiv::GetIntegrationRule() const
{
    if(!fIntRule) DebugStop();
    return *fIntRule;
}

void TPZSBFemVolumeHdiv::SetPreferredOrder ( int order )
{
    fPreferredOrder = order;
}

void TPZSBFemVolumeHdiv::BuildCornerConnectList(std::set<int64_t> &connectindexes) const
{
    fElementVec1D[1]->BuildCornerConnectList(connectindexes);
}

int TPZSBFemVolumeHdiv::Dimension() const
{
    return Reference()->Dimension();
}

int64_t TPZSBFemVolumeHdiv::ConnectIndex(int i) const
{
    if (i > fConnectIndexes.size())
    {
        DebugStop();
    }        
    return fConnectIndexes[i];
}

int TPZSBFemVolumeHdiv::NConnects() const
{
    return fConnectIndexes.size();   
}

TPZCompEl * CreateSBFemFluxCompEl(TPZCompMesh &mesh, TPZGeoEl *gel)
{
    return new TPZSBFemVolumeHdiv(mesh, gel);    
}
