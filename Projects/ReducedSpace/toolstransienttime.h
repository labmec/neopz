//
//  toolstransienttime.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_toolstransienttime_h
#define PZ_toolstransienttime_h

#include "pzcmesh.h"
#include "pzanalysis.h"
#include "tpzcompmeshreferred.h"
#include "pznlfluidstructure2d.h"


class TPZElastoPlasticAnalysis;


class ToolsTransient
{
    public:
    
    ToolsTransient();
    ToolsTransient(int pOrder);
    
    ~ToolsTransient();
    
    //---------------------------------------------------------------
    
    void Run();
    void RunPlasticity();
    
    void InitializeUncoupledMeshesAttributes();
    TPZCompMesh * ElastCMeshReferenceProcessed();
    void Mesh2D();
    TPZCompMesh * CMeshElastic();
    void SetSigmaNStripeNum(TPZCompMesh * cmeshref, int actStripe);
    void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh);
    TPZCompMeshReferred * CMeshReduced(TPZCompMesh * cmeshref);
    TPZCompMesh * CMeshPressure();
    void CMeshMultiphysics();
    
    bool SolveSistTransient(TPZAnalysis *an, bool initialElasticKickIsNeeded);
    
    void TransferSolutions(TPZCompMesh * lastMPhysicsCMesh, TPZCompMesh * lastElastReferredCMesh);
    void TransferElasticSolution(TPZCompMesh * cmeshFrom);
    REAL IntegrateSolution(TPZCompMesh * cmesh, int variable);//0 = meshvec[0] ; 1 = meshvec[1]
    void TransferLeakoff(TPZCompMesh * oldMphysicsCMesh);
    
    //---------------------------------------------------------------
    
    void MassMatrix(TPZFMatrix<REAL> & Un);
    
    void StiffMatrixLoadVec(TPZAnalysis *an,
                            TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);
    
    void PostprocessPressure();
    void PostProcessAcumVolW();
    void PostProcessVolLeakoff();
    
    REAL ComputeKIPlaneStrain();
    
    void CheckConv();
    
    void SolveInitialElastoPlasticity(TPZElastoPlasticAnalysis &analysis, TPZCompMesh *Cmesh);
    
    TPZCompMesh * CMeshElastoPlastic(TPZGeoMesh *gmesh, REAL SigmaN);    
    
    //---------------------------------------------------------------
    
    int fpOrder;
    bool fMustStop;
    
    TPZNLFluidStructure2d * fCouplingMaterial1;
    TPZNLFluidStructure2d * fCouplingMaterial2;
    TPZGeoMesh * fgmesh;
    
    /** fmeshvec[0] = Malha computacional elastica do tipo referred */
    /** fmeshvec[1] = Malha computacional de fluxo 1D */
    TPZVec<TPZCompMesh *> fmeshvec;
    TPZCompMesh * fmphysics;
};



template<class TVar>
class TElastSolFunction : public TPZFunction<TVar>
{
    public:
    
    /**
     * Class constructor
     */
    TElastSolFunction();
    
    TElastSolFunction(TPZCompMesh * cmesh);
    
    /**
     * Class destructor
     */
    ~TElastSolFunction();
    
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df);
    
    /** Returns number of functions.
     */
    virtual int NFunctions() const;
    
    /** Polynomial order of this function. In case of non-polynomial
     * function it can be a reasonable approximation order.
     */
    virtual int PolynomialOrder() const;
    
    TPZCompMesh * fcmesh;
    int64_t fIniElIndex;
    
};

template<class TVar>
class TLeakoffFunction : public TPZFunction<TVar>
{
public:
    
    /**
     * Class constructor
     */
    TLeakoffFunction();
    
    TLeakoffFunction(TPZCompMesh * cmesh);
    
    /**
     * Class destructor
     */
    ~TLeakoffFunction();
    
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df);
    
    /** Returns number of functions.
     */
    virtual int NFunctions() const;
    
    /** Polynomial order of this function. In case of non-polynomial
     * function it can be a reasonable approximation order.
     */
    virtual int PolynomialOrder() const;
    
    TPZCompMesh * fcmesh;
    int64_t fIniElIndex;
};

#endif
