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





class ToolsTransient {
    
    public:
    
    ToolsTransient();
    ToolsTransient(int pOrder,
                   REAL Lx, REAL Ly, REAL Lf, REAL Hf, REAL E, REAL Poiss, REAL Fx, REAL Fy, REAL Visc, REAL SigN,
                   REAL Qinj, REAL Ttot, REAL Nsteps, REAL Cl, REAL Pe, REAL SigmaConf, REAL Pref, REAL vsp, REAL KIc);
    
    ~ToolsTransient();
    
    //---------------------------------------------------------------
    
    void Run();
    
        TPZGeoMesh * Mesh2D(REAL lmax);
        TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh);
        void SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh);
        TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh);
        TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh);
        TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial);
        
        bool SolveSistTransient(REAL & deltaT, REAL & actTime, REAL maxTime, TPZFMatrix<REAL> & InitialSolution, TPZAnalysis *an,
                                TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, int & step,
                                REAL Jradius, std::stringstream & outP, std::stringstream & outW, std::stringstream & outJ,
                                std::string & outputfile);
    
    //---------------------------------------------------------------
    
    void MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics, TPZFMatrix<REAL> & Un);
    
    void StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis *an,
                                   TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);

    void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    void SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial, std::stringstream & outP);

    void FillPositionPressure(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::map<REAL,REAL> & pos_pressure);
    
    void PosProcessMult(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis *an, std::string plotfile);
    
    void PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum);
    
    REAL ComputeKIPlaneStrain(TPZCompMesh * elastMesh, REAL young, REAL poisson,
                              REAL radius, std::stringstream & outFile, int cent = -1,
                              REAL TimeValue = -1, bool firstCall = true);
    
    void CheckConv(TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial,
                          TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    //---------------------------------------------------------------
    
    int fpOrder;
    
    InputDataStruct * fInputData;
};

#endif
