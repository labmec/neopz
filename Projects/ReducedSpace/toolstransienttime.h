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
    
    ~ToolsTransient();
    
    //---------------------------------------------------------------
    
    static TPZGeoMesh * Mesh2D(REAL lf, REAL ldom, REAL hdom, REAL lmax);
    
    static TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh, int pOrder, REAL E, REAL poisson, REAL sigN);
    static TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder, REAL E, REAL poisson);
    static TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh, int pOrder, REAL Qinj);
    static TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial,
                                              REAL ED, REAL nu, REAL fx, REAL fy, REAL Hf, REAL Lf, REAL visc,
                                              REAL Qinj, REAL Cl, REAL Pe, REAL SigmaConf, REAL sigN, REAL Pref, REAL vsp, REAL KIc, REAL Lx);
    
    static void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    static void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);
    
    //---------------------------------------------------------------
    
    static void SolveSistTransient(REAL deltaT, REAL maxTime, TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an,
                                   TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    static void MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics, TPZFMatrix<REAL> & Un);
    
    static void StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis *an,
                                   TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec);
    
    //static TPZCompMesh *CMeshProjectionL2(TPZGeoMesh *gmesh, int dim, int matId, int pOrder, TPZVec<STATE> &solini);
    
    //retorna pressao na metade do comprimento da fratura
    static void SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial, std::stringstream & outP);
    static void FillPositionPressure(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::map<REAL,REAL> & pos_pressure);
    
    //Retona a solucao inicial referente a uma malha computacional
    //static TPZFMatrix<REAL> InitialSolution(TPZGeoMesh *gmesh, TPZCompMesh * cmesh, int matId, int porder, REAL valsol);
    
    static void PosProcessMult(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis *an, std::string plotfile);
    
    //static TPZFMatrix<REAL> SetSolution(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder, int matId, REAL valIni);
    
    static void PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum);
    
    static REAL ComputeKIPlaneStrain(TPZCompMesh * elastMesh, REAL young, REAL poisson,
                                     REAL radius, std::stringstream & outFile, int cent = -1, REAL TimeValue = -1, bool firstCall = true);
    
    static bool PropagateOneStep(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZNLFluidStructure2d * &mymaterial);
    
    static void CheckConv(TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial,
                          TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
};

#endif
