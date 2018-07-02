//
//  ToolsOptimal.h
//  PZ
//
//  Created by Denise de Siqueira on 6/18/15.
//
//

#ifndef __PZ__ToolsOptimal__
#define __PZ__ToolsOptimal__

#include <iostream>

#endif /* defined(__PZ__ToolsOptimal__) */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "pzmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "TPZMaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "TPZRefPattern.h"


//Criacao de Malhas Geometricas

// Malha Triangular

TPZGeoMesh *OptimalGeoMesh(bool ftriang, REAL Lx, REAL Ly);
//Malha computacional convencional
TPZCompMesh *OptimalCompMesh(TPZGeoMesh *gmesh,int porder);
// Variavel de estado factivel
void StateVar(const TPZVec<REAL> &loc, TPZVec<STATE> &u, TPZFMatrix<STATE> &du);
void OptForcing(const TPZVec<REAL> &pt, TPZVec<STATE> &res);
void SolveLUOpt ( TPZAnalysis &an );
