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
    
    static void SolveSistTransient(REAL deltaT,REAL maxTime, TPZElastPressure * &mymaterial, TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
    static TPZAutoPointer <TPZMatrix<REAL> > MassMatrix(TPZElastPressure *mymaterial, TPZCompMesh *mphysics);
    
    static void StiffMatrixLoadVec(TPZElastPressure *mymaterial, TPZCompMesh* mphysics, TPZAnalysis &an, TPZFMatrix<REAL> &matK1, TPZFMatrix<REAL> &fvec);
    
    static TPZCompMesh *CMeshProjectionL2(TPZGeoMesh *gmesh, int pOrder, TPZVec<STATE> &solini);
    
    static void SaidaMathPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);
    
};

#endif
