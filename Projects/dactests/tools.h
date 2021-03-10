//
//  tools.h
//  PZ
//
//  Created by Douglas Castro on 1/29/15.
//
//

#ifndef __PZ__tools__
#define __PZ__tools__

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

//#include "pzhdivfull.h"
#include "pzelchdiv.h"

#include "pzgeopyramid.h"

#include "pznumeric.h"

#include "TPZExtendGridDimension.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


#include "pzstack.h"

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>



class tools{
    
    
public:
    
    tools();
    ~tools();
    
    static void RotateGeomesh(TPZGeoMesh *gmesh, REAL CounterClockwiseAngle, int &Axis);
    
    static void RotateNode(TPZVec<REAL> &iCoords, REAL CounterClockwiseAngle, int &Axis);
    
    static void UniformRefinement(TPZGeoMesh* gmesh, int nDiv);
    
//    //criar elementos esqueleto
//    static void AddWrap(TPZMultiphysicsElement *mfcel, int matskeleton, TPZStack< TPZStack< TPZMultiphysicsElement *, 7> > &ListGroupEl);
    
    static void PrintLS(TPZAnalysis *an);
    
    static void SolveSyst(TPZAnalysis &an, TPZCompMesh *fCmesh, REAL &assemble_time, REAL &solving_time);
    
    static void PosProcessMultphysics(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics, TPZAnalysis &an, std::string plotfile, int dim);
    
    /* pos processamento para H1 */
    static void PosProcess(TPZAnalysis &an, std::string plotfile, int dim);
    
    static void PrintDebugMapForMathematica(std::string filenameHdiv, std::string filenameL2, std::map<REAL, REAL> &fDebugMapL2, std::map<REAL, REAL> &fDebugMapHdiv);
    
    static bool MyDoubleComparer(REAL a, REAL b);
    
};

#endif /* defined(__PZ__tools__) */
