//
//  toolstransienttime.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef PZ_toolstransienttime_h
#define PZ_toolstransienttime_h

#include <iostream>

#include "pzelastpressure.h"
#include "pznlfluidstructure2d.h"
#include "pzcmesh.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "tpzautopointer.h"
#include "pzanalysis.h"
#include "pzfmatrix.h"
#include "pzl2projection.h"

class ToolsTransient {
    
    public:
    
    ToolsTransient();
    
    ~ToolsTransient();
    
    static void SolveSistTransient(REAL deltaT,REAL maxTime, TPZFMatrix<REAL> InitialSolution, TPZAnalysis *an, TPZNLFluidStructure2d * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    static TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZNLFluidStructure2d *mymaterial, TPZCompMesh *mphysics);
    
    static void StiffMatrixLoadVec(TPZNLFluidStructure2d *mymaterial, TPZCompMesh* mphysics, TPZAnalysis *an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);
    
    static TPZCompMesh *CMeshProjectionL2(TPZGeoMesh *gmesh, int dim, int matId, int pOrder, TPZVec<STATE> &solini);
    
    //retorna pressao na metade do comprimento da fratura
    static REAL SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, std::stringstream & outP);
    
    //Retona a solucao inicial referente a uma malha computacional
    static TPZFMatrix<REAL> InitialSolution(TPZGeoMesh *gmesh, TPZCompMesh * cmesh, int matId, int porder, REAL valsol);
    
    static void PosProcessMult(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis *an, std::string plotfile);
    
    static TPZFMatrix<REAL> SetSolution(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder, int matId, REAL valIni);
    
    static void PlotWIntegral(TPZCompMesh *cmesh, std::stringstream & outW, int solNum);
};

#endif
