#ifndef TPZSurfaceH
#define TPZSurfaceH

/*
 *  TPZCurve.h
 *  PZ
 *
 *  Created by Omar Duran Triana on 03/27/15.
 *  Copyright 2015 __Labmec__. All rights reserved.
 *
 */


#include "tpzautopointer.h"
#include "pzfmatrix.h"
#include <math.h>

#include "TPZVTKGeoMesh.h"

#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "tpzquadraticquad.h"
#include "tpzquadratictrig.h"

#include "TPZQuadSphere.h"
#include "TPZTriangleSphere.h"

#include "TPZQuadTorus.h"
#include "TPZTriangleTorus.h"


#include "tpzgeoblend.h"

#include "tpzgeoelrefpattern.h"
#include "pzgmesh.h"


class TPZSurface {
    
private:
    
    int fdimension;
    STATE fradius;
    STATE Pi;
    bool fIsclosed;
    
    TPZGeoMesh * fgeometricmesh;
    
public:
    
    TPZSurface();
    
    ~TPZSurface();
    
    void SetRadius(STATE r) { fradius = r; }
    STATE GetRadius() {  return fradius; }
    
    TPZGeoMesh * GetGeometry() {  return fgeometricmesh; }
    
    void MakeCube();
    void MakeRhombohedron();    
    void MakeSphereFromQuadrilateral();
    void MakeSphereFromTriangle();
    
    
    void PrintMe();
    void RefineMe(int i);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    TPZManVector<STATE,3> ParametricCircle(STATE t);
    TPZManVector<STATE,3> ParametricSphere(STATE phi,STATE theta);
    
};


#endif