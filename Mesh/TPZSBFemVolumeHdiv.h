//
//  TPZSBFemVolumeHdiv.hpp
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

#include "pzshapelinear.h"
#include "pzshapetriang.h"
#include "pzshapequad.h"
#include "pzelchdiv.h"

#include "pzgraphelq2dd.h"
#include "pzgraphelq3dd.h"
#include "tpzgraphelprismmapped.h"

using namespace std;
using namespace pzshape;

class TPZSBFemVolumeHdiv : public TPZInterpolationSpace
{
    /// index of element group
    int64_t fElementGroupIndex = -1;

    int fSkeleton = -1;

    TPZCompEl *fElementGroup = 0;

	/** @brief List of pointers to computational elements */
    // Order of the elements: fLeftFlux, fRightFlux, fSkeleton
	TPZManVector<TPZCompEl* ,3> fElementVec1D;
    	
	/** @brief Indexes of the connects of the element */
	TPZVec<int64_t> fConnectIndexes;
    
    /// pointer to the integration rule
    TPZIntPoints *fIntRule = 0;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndicesFluxInt;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndicesFlux;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<double>> fPhiFluxInt;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<double>> fPhiFlux;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<std::complex<double>> fEigenvalues;
    
    /// Multiplier coeficients associated with the solution
    TPZFNMatrix<30, std::complex<double>> fCoeficients;
    TPZFNMatrix<30, std::complex<double>> fCoeficientsD;

public:
    
    TPZSBFemVolumeHdiv(TPZCompMesh & mesh, TPZGeoEl * gel);
    
    virtual ~TPZSBFemVolumeHdiv()
    {
    }

    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override;

    void AddElement1D(TPZCompEl * cel, int localindex);

    void SetPhiFlux(TPZFMatrix<std::complex<double>> &phiflux, TPZManVector<std::complex<double>> &eigval);

    void SetLocalIndicesFlux(TPZManVector<int64_t> &localindicesfluxint, TPZManVector<int64_t> &localindicesflux);

    void LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd);

    int NShapeF() const override;

    virtual int NConnectShapeF(int icon, int order) const override;

    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual const TPZIntPoints &GetIntegrationRule() const override;
    
    /** @brief Returns a reference to an integration rule suitable for integrating the interior of the element */
    virtual TPZIntPoints &GetIntegrationRule() override;

    virtual void SetPreferredOrder ( int order ) override;

    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override;

    virtual int Dimension() const override;

    virtual int64_t ConnectIndex(int i) const override;

    virtual int NConnects() const override;

    virtual void ReallyComputeSolution(TPZMaterialDataT<STATE> & data) override;

    void InitMaterialData(TPZMaterialData & data) override
    {
        auto datat = dynamic_cast<TPZMaterialDataT<STATE>*>(&data);
        InitMaterialDataT(*datat);
    }

    void InitMaterialDataT(TPZMaterialDataT<STATE> &data)
    {
        auto collapsedcompel = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear>*>(fElementVec1D[1]);
        collapsedcompel->InitMaterialData(data);
    }

    virtual void ComputeRequiredData(TPZMaterialDataT<STATE> &data, TPZVec<REAL> &qsi) override
    {
//         TPZManVector<REAL, 3> qsilow(qsi);
//         qsilow.Resize(this->Reference()->Dimension() - 1);

//         auto collapsedcompel = dynamic_cast<TPZCompElHDivSBFem<pzshape::TPZShapeLinear>*>(fElementVec1D[1]);
// #ifdef PZDEBUG  
//         if (!collapsedcompel) DebugStop();
// #endif
        
//         collapsedcompel->ComputeRequiredData(data, qsilow);

        data.xParametric = qsi;
//         if (data.fNeedsSol)
//         {
            ReallyComputeSolution(data);
        // }
    };

    void Shape(TPZVec<REAL> &qsi,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphidxi) override
    {
        DebugStop();
    };

    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const override
    {
        DebugStop();
        return 0;
    };

    virtual void SetConnectIndex(int inode, int64_t index) override
    {
        DebugStop();
    };

    virtual int NSideConnects(int iside) const override
    {
        DebugStop();
        return 0;
    };

    virtual int SideConnectLocId(int icon,int is) const override
    {
        DebugStop();
        return 0;
    };

	virtual void PRefine ( int order ) override
    {
        DebugStop();
    };
};

TPZCompEl * CreateSBFemFluxCompEl(TPZCompMesh &mesh, TPZGeoEl *gel);
