//
//  TPZSBFemVolumeMultiphysics.hpp
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

template <class TGeometry>
class TPZSBFemVolumeMultiphysics : public TPZMultiphysicsCompEl<TGeometry>
{
    /// index of element group
    int64_t fElementGroupIndex = -1;

    int fSkeleton = -1;

    TPZCompEl *fElementGroup = 0;

	TPZManVector<TPZCompElSide ,5> fElementVec;

	/** @brief List of pointers to computational elements */
    // Order of the elements: fDifPressure, fInterface, fLeftFlux, fRightFlux, fAverPressure, fInterface, fSkeleton
	TPZManVector<TPZCompEl* ,7> fElementVec1D;
    	
	/** @brief Indexes of the connects of the element */
	TPZVec<int64_t> fConnectIndexes;
    
    /// pointer to the integration rule
    TPZIntPoints *fIntRule = 0;
    
    /// vector of local indices of multipliers in the group
    TPZManVector<int64_t> fLocalIndices;
    
    /// Section of the phi vector associated with this volume element
    TPZFNMatrix<30,std::complex<double>> fPhi;
    
    /// Eigenvlues associated with the internal shape functions
    TPZManVector<std::complex<double>> fEigenvalues;

public:
    
    TPZSBFemVolumeMultiphysics(TPZCompMesh & mesh, TPZGeoEl * gel);
    
    virtual ~TPZSBFemVolumeMultiphysics()
    {
    }

    void AddElement1D(TPZCompEl * cel, int localindex);

    virtual void AddElement(TPZCompEl *cel, int64_t meshindex) override;

    virtual TPZManVector<TPZCompElSide,5> &ElementVec() override { return fElementVec; }

    virtual TPZCompEl *Element(int64_t elindex) override;

    /**@brief Returns referred element of this*/
	virtual TPZCompEl *ReferredElement(int64_t mesh) override;

    void SetSkeleton(int64_t skeleton);
    
    int64_t SkeletonIndex();
    
    void LoadCoef(TPZFMatrix<std::complex<double>> &coef, TPZFMatrix<std::complex<double>> &coefd);

    void AffineTransform(TPZVec<TPZTransform<> > &trVec) const override;

    void SetElementGroupIndex(int64_t index);

    /** @brief Returns the number of nodes of the element */
    virtual int NConnects() const override;

    virtual int64_t ConnectIndex(int i) const override;

    virtual int Dimension() const override;

    virtual void CreateGraphicalElement(TPZGraphMesh &graphmesh, int dimension) override;

    void EvaluateError(TPZVec<REAL> &errors,bool store_error) override;

    virtual TPZIntPoints & GetIntegrationRule() const override;

    virtual TPZIntPoints & GetIntegrationRule() override;

    void InitializeIntegrationRule() override;

    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override;

    virtual int NShapeF() const;

    virtual TPZCompEl *Element(int elindex);

    TPZManVector<TPZCompEl * ,7> & SBFemElementVec();

    void SetPhiEigVal(TPZFMatrix<std::complex<double>> &phi, TPZManVector<std::complex<double>> &eigval);

    void SetPhiFlux(TPZFMatrix<std::complex<double>> &phiflux, TPZManVector<std::complex<double>> &eigval);

    //@{
	/**
	 * @brief Post processing method which computes the solution for the var post processed variable.
	 * @param qsi coordinate of the point in master element space where the solution will be evaluated
	 * @param var variable which will be computed
	 * @param sol (output) solution computed at the given point
	 * @see TPZMaterial::VariableIndex
	 * @see TPZMaterial::NSolutionVariables
	 * @see TPZMaterial::Solution
	 */
	/** The var index is obtained by calling the TPZMaterial::VariableIndex method with a post processing name */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;

	void InitMaterialData(TPZVec<TPZMaterialDataT<STATE>> &dataVec, TPZVec<int64_t> *indices = 0) override;

    void ComputeRequiredData(TPZVec<TPZMaterialDataT<STATE>> &dataVec,TPZVec<REAL> &qsi);

    void SetLocalIndices(TPZManVector<int64_t> &localindices, TPZManVector<int64_t> &localindicesint, TPZManVector<int64_t> &localindicesflux);

    virtual void SetConnectIndexes(TPZVec<int64_t> &indexes) override
    {
        fConnectIndexes = indexes;
    }

    /** @brief Method for creating a copy of the element */
    virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        DebugStop();
        return 0;
    }

    virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
                                    std::map<int64_t,int64_t> & gl2lcConMap,
                                    std::map<int64_t,int64_t> & gl2lcElMap) const override
    {
        DebugStop();
        return 0;
    }

    virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override
    {
        DebugStop();
    }

    virtual void SetConnectIndex(int inode, int64_t index) override
    {
        DebugStop();
    }
};

TPZCompEl * CreateSBFemMultiphysicsLinearEl(TPZGeoEl *gel, TPZCompMesh &mesh);

TPZCompEl * CreateSBFemMultiphysicsQuadEl(TPZGeoEl *gel, TPZCompMesh &mesh);

TPZCompEl * CreateSBFemMultiphysicsCubeEl(TPZGeoEl *gel, TPZCompMesh &mesh);

TPZCompEl * CreateSBFemMultiphysicsPrismaEl(TPZGeoEl *gel, TPZCompMesh &mesh);
