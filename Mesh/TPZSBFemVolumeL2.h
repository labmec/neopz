//
//  TPZSBFemVolumeL2.hpp
//  PZ
//
//  Created by Karolinne Coelho on 25/01/2021.
//
//

#pragma once

#include <stdio.h>
#include "pzcompel.h"
#include "pzelmat.h"
#include "TPZSBFemVolume.h"
#include "TPZCompElHDivSBFem.h"
#include "TPZMultiphysicsCompMesh.h"
#include "TPZSBFemMultiphysicsElGroup.h"
#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzelchdiv.h"

#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

using namespace std;
using namespace pzshape;

class TPZSBFemVolumeL2 : public TPZSBFemVolume
{
    // TPZSBFemVolume has the fSkeleton data

	/** @brief List of pointers to computational elements */
    // Order of the elements: fDifPressure, fInternal, fSkeleton
	TPZManVector<TPZCompEl* ,3> fElementVec1D;
    
    /// vector of local indices of the internal pressure
    TPZManVector<int64_t> fLocalIndicesInt;
    
    /// Section of the phi vector associated the internal pressure
    TPZFNMatrix<30,std::complex<double>> fPhiInt;
    
    /// Multiplier coeficients associated with the internal pressure
    TPZFNMatrix<30, std::complex<double>> fCoeficientsInt;
    
    /// vector of local indices of the internal pressure
    TPZManVector<int64_t> fLocalIndicesDifPr;
    
    /// Section of the phi vector associated the internal pressure
    TPZFNMatrix<30,std::complex<double>> fPhiDifPr;
    
    /// Multiplier coeficients associated with the internal pressure
    TPZFNMatrix<30, std::complex<double>> fCoeficients;
    TPZFNMatrix<30, std::complex<double>> fCoeficientsD;

	/// Indexes of the connects of all pressure elements
	TPZVec<int64_t> fConnectIndexes;

public:
    
    TPZSBFemVolumeL2(TPZCompMesh & mesh, TPZGeoEl * gel);
    
    virtual ~TPZSBFemVolumeL2()
    {
        // Reference()->ResetReference();
    }

    void SetLocalIndices(TPZManVector<int64_t> &localindices)
    {
        fLocalIndices = localindices;
    }


    void InitMaterialData(TPZMaterialData & data) override
    {
        TPZSBFemVolume::InitMaterialData(data);
    }

    virtual void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override
    {
        // auto CSkeleton = dynamic_cast<TPZInterpolationSpace *> (fElementVec1D[2]);

        // // compute the lower dimensional shape functions
        // TPZManVector<REAL, 3> qsilow(qsi);
        // qsilow.Resize(this->Reference()->Dimension() - 1);
        // CSkeleton->ComputeRequiredData(data, qsilow);

        // auto Ref2D = this->Reference();
        // TPZMaterialDataT<STATE> data2d;
        // Ref2D->Jacobian(qsi, data2d.jacobian, data2d.axes, data2d.detjac, data2d.jacinv);
        // data.axes = data2d.axes;

        // data.xParametric = qsi;

        data.xParametric = qsi;
        // if (data.fNeedsSol)
        // {
            ReallyComputeSolution(data);
        // }
    }

    void LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd);

    virtual void ReallyComputeSolution(TPZMaterialDataT<STATE> & data) override;

    void AddElement1D(TPZCompEl * cel, int localindex)
    {
        fElementVec1D[localindex] = cel;
        auto ncon = fConnectIndexes.size();
        auto nconcel =cel->NConnects();
        fConnectIndexes.Resize(ncon+nconcel);
        for (auto i = 0; i < nconcel; i++)
        {
            fConnectIndexes[i+ncon] = cel->ConnectIndex(i);
        }
    }

    void SetElementGroupIndex(int64_t index);

    void SetPhiEigVal(TPZFMatrix<std::complex<double>> &phi, TPZManVector<std::complex<double>> &eigval)
    {
        fEigenvalues = eigval;
        
        int nrow = fLocalIndices.size();
        fPhi.Resize(nrow, phi.Cols());
        for (int i = 0; i < nrow; i++) {
            for (int j = 0; j < phi.Cols(); j++) {
                fPhi(i, j) = phi(fLocalIndices[i], j);
            }
        }
    };

    /** @brief Returns the number of nodes of the element */
    virtual int NConnects() const override
    {
        return fConnectIndexes.size();   
    }

    virtual int64_t ConnectIndex(int i) const override
    {
        if (i > fConnectIndexes.size())
        {
            DebugStop();
        }        
        return fConnectIndexes[i];
    }

    virtual int Dimension() const override
    {
        TPZGeoEl *reference = Reference();
        return reference->Dimension();
    }

    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override
    {
        DebugStop();
    }

    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        fElementVec1D[2]->BuildCornerConnectList(connectindexes);
    }

    virtual int NSideConnects(int iside) const override
    {
        DebugStop();
        return 0;
    }

    virtual int SideConnectLocId(int icon,int is) const override
    {
        DebugStop();
        return 0;
    }

    virtual int NShapeF() const override
    {
        int nc = fElementVec1D[2]->NConnects();
        int nshape = 0;
        for (int ic=0; ic<nc; ic++) {
            TPZConnect &c = Connect(ic);
            nshape += c.NShape();
        }
        return nshape;
    }

    virtual TPZCompEl *Element(int elindex)
    {
        return fElementVec1D[elindex];
    }

    virtual TPZManVector<TPZCompEl *,3> & ElementVec()
    {
        return fElementVec1D;
    }

    virtual void SetConnectIndexes(TPZVec<int64_t> &indexes)
    {
        fConnectIndexes = indexes;
    }
};

TPZCompEl * CreateSBFemPressureCompEl(TPZCompMesh &mesh, TPZGeoEl *gel);
