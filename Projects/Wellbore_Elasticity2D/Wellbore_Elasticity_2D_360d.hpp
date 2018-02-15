//
//  Wellbore_Elasticity_2D_360d.hpp
//  PZ
//
//  Created by Nathalia on 6/7/16.
//
//

#ifndef Wellbore_Elasticity_2D_360d_hpp
#define Wellbore_Elasticity_2D_360d_hpp

#include <stdio.h>

#endif /* Wellbore_Elasticity_2D_360d_hpp */

#ifdef HAVE_CONFIG_H
#include <pz_config.h>
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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.elasticity"));
//static LoggerPtr loggeradap(Logger::getLogger("pz.adaptivity"));
//static LoggerPtr loggerconv(Logger::getLogger("pz.adaptivity.conv"));
//static LoggerPtr loggerpoint(Logger::getLogger("pz.adaptivity.points"));
#endif

#include <time.h>
#include <stdio.h>
#include <fstream>

using namespace std;

class X{
    
    X(){
        
        DebugStop();
        
    }
    
    ~X(){
        
        DebugStop();        
        
    }
    
    // Cria Malha Geometrica
    TPZGeoMesh *CircularGeoMesh(REAL rwb, REAL re, int ncirc, int nrad, REAL DrDcirc);
    
    // Cria Malha COmputacional
    TPZCompMesh *CircularCMesh(TPZGeoMesh *gmesh, int pOrder);
    
    // Executa o codigo
    int Problem2D();
    
};

