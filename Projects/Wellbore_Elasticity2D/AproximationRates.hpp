//
//  AproximationRates.hpp
//  PZ
//
//  Created by Nathalia Batalha on 9/8/17.
//
//

#ifndef AproximationRates_hpp
#define AproximationRates_hpp

#include <stdio.h>

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "pzvec.h"
#include "pzadmchunk.h"
#include "pzcmesh.h"
#include "pzvec_extras.h"
#include "pzdebug.h"
#include "pzcheckgeom.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"
#include "pzgeoelbc.h"

#include "pzintel.h"
#include "pzcompel.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzmaterial.h"
#include "pzbndcond.h"
#include "pzelasmat.h"
#include "pzplaca.h"
#include "pzmat2dlin.h"
#include "pzmathyperelastic.h"
#include "pzmattest3d.h"
#include "pzmatplaca2.h"

#include "pzfunction.h"
#include "TPZVTKGeoMesh.h"

#include <time.h>
#include <stdio.h>
#include <fstream>

#include "pzlog.h"
#include "pzgmesh.h"
#include "TPZMatElasticity2D.h"
#include <iostream>
#include <fstream>
#include <string>
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"

#include "pzstepsolver.h"
#include "math.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "TPZReadGIDGrid.h"
#include "TPZParFrontStructMatrix.h"

#include <cmath>
#include <set>

#include "pzelast3d.h"
#include <pzgengrid.h>

#include "pzrandomfield.h"
#include "GeometricMesh.hpp"
#include "ComputationalMesh.hpp"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
#endif

/** @brief Compute approximation rates for Inclined wellbore analytic solution */
int ApproximationRates();

/** @brief Compute L2 error for a given computational mesh */
void ComputeErrorL2(TPZCompMesh * cmesh, REAL &error);

/** @brief compute uniform refinement */
void UniformRefinement(TPZGeoMesh *gmesh, int nh);

/** @brief compute inner product between tensors */
REAL Inner_Product(TPZFNMatrix<4,REAL> S, TPZFNMatrix<4,REAL> T);

#endif /* AproximationRates_hpp */
