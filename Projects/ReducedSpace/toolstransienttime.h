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
    
        TPZGeoMesh * Mesh2D(REAL lmax);
        TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh);
        void SetSigmaNStripeNum(TPZCompMesh * cmesh, int actStripe);
        void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh);
        TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
        TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh);
        TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial);
        
        bool SolveSistTransient(REAL & deltaT, REAL & actTime, REAL maxTime, TPZAnalysis *an,
                                TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int & step,
                                REAL Jradius, std::stringstream & outP, std::stringstream & outW, std::stringstream & outJ,
                                std::string & outputfile);
    
    void TransferElasticSolution(TPZCompMeshReferred * cmeshFrom, TPZCompMeshReferred * cmeshTo);
    REAL IntegrateSolution(TPZCompMesh * cmesh, int variable);
    std::map<int,REAL> TransferLeakoff(TPZCompMesh * oldMphysicsCMesh, TPZCompMesh * newFluidCMesh, std::stringstream & outVl);
    
    //---------------------------------------------------------------
    
    void MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics, TPZFMatrix<REAL> & Un);
    
    void StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis *an,
                                   TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);

    void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    void SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial, std::stringstream & outP);

    void FillPositionPressure(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::map<REAL,REAL> & pos_pressure);
    
    void PosProcessMult(TPZAnalysis *an, std::string plotfile);
    
    void PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum);
    
    REAL ComputeKIPlaneStrain(TPZCompMesh * elastMesh, REAL young, REAL poisson,
                              REAL radius, std::stringstream & outFile, int cent = -1,
                              REAL TimeValue = -1, bool firstCall = true);
    
    void CheckConv(TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial,
                          TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    void SolveInitialElastoPlasticity(TPZElastoPlasticAnalysis &analysis, TPZCompMesh *Cmesh);
    
    TPZCompMesh * CMeshElastoPlastic(TPZGeoMesh *gmesh, REAL SigmaN);    
    
    //---------------------------------------------------------------
    
    int fpOrder;
    bool fMustStop;
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
    virtual int NFunctions();
    
    /** Polynomial order of this function. In case of non-polynomial
     * function it can be a reasonable approximation order.
     */
    virtual int PolynomialOrder();
    
    TPZCompMesh * fcmesh;
    int fIniElIndex;
    
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
    virtual int NFunctions();
    
    /** Polynomial order of this function. In case of non-polynomial
     * function it can be a reasonable approximation order.
     */
    virtual int PolynomialOrder();
    
    TPZCompMesh * fcmesh;
    int fIniElIndex;
    std::map<int,REAL> fleakoffMap;
};

#endif
