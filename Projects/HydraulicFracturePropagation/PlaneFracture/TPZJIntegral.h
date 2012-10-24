//
//  TPZJIntegral.h
//  PZ
//
//  Created by Labmec on 25/04/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#ifndef PZ_TPZJIntegral_h
#define PZ_TPZJIntegral_h

#include <iostream>
#include "TPZGeoElement.h"
#include "pzcmesh.h"
#include "pzelasmat.h"
#include "pzelast3d.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "TPZPlaneFracture.h"

enum pathType {ELinearPath, EExternalArcPath, EInternalArcPath};

const REAL Pi = 3.1415926535897932384626433832795;


struct LinearPath
{
public:
    
    LinearPath();//It is not to be used
    LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim);
    ~LinearPath();
    
    void X(REAL t, TPZVec<REAL> & xt);
    void dXdt(REAL t, TPZVec<REAL> & dxdt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    TPZVec<REAL> Func(REAL t);
    
    TPZVec<REAL> BoundaryFunc(TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
protected:
    
    /** Origin of coupled arc */
    TPZVec<REAL> fOrigin;
    
    /** This direction defines the coupled arc plane.
     *  (this direction is orthogonal to coupled arc plane and defines
     *   the right hand convention for the coupled arc direction)
     */
    TPZVec<REAL> fNormalDirection;
    
    /** Radius of coupled arc */
    REAL fradius;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
    
    /**
     *  The Func() method need to call ComputeXInverse method to get solutions on compMesh.
     *  With the Id of initial element provided, near the searched point, the search algorithm could be optimized.
     */
    int fInitialElementId;
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** For 2D problems (plane strain or plane stress), fMeshDim=2 */
    /** For 3D problems, fMeshDim=3 */
    int fMeshDim;
    
    /** pressure applied inside fracture */
    REAL fcrackPressure;
};


struct ArcPath
{
public:
    
    ArcPath();//It is not to be used
    ArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim);
    ~ArcPath();
    
    void X(REAL t, TPZVec<REAL> & xt);
    void dXdt(REAL t, TPZVec<REAL> & dxdt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    TPZVec<REAL> Func(REAL t);
    
    TPZVec<REAL> BoundaryFunc(TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
protected:
    
    /** Origin of arc */
    TPZVec<REAL> fOrigin;
    
    /** This direction defines the arc plane.
     *  (this direction is orthogonal to arc plane and defines
     *   the right hand convention for the arc direction)
     */
    TPZVec<REAL> fNormalDirection;
    
    /** Radius of arc */
    REAL fradius;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
    
    /**
     *  The Func() method need to call ComputeXInverse method to get solutions on compMesh.
     *  With the Id of initial element provided, near the searched point, the search algorithm could be optimized.
     */
    int fInitialElementId;
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** For 2D problems (plane strain or plane stress), fMeshDim=2 */
    /** For 3D problems, fMeshDim=3 */
    int fMeshDim;
};

/**
 *  ITS ALWAYS GOOD TO REMEMBER:
 *          THE CLASS PATH CONSIDERS THAT THE NORMAL DIRECTION IS IN X,Z PLANE (JUST LIKE FRACTURE PLANE) AND
 *          THE ORIENTATION OF ARC and LINEAR stretch is:
 *              ARC : RIGHT HAND DIRECTION WITH RESPECT TO GIVEN NORMAL AXE (axe that defines the (orthogonal) arc plane).
 *              LINEAR: FROM THE END OF ARC (supposed to be inside crack plane) TO ORIGIN.
 * SO, ITS ALWAYS GOOD DEFINE NORMAL AXE TANGENT TO THE CRACK BOUNDARY, FOLLOWING RIGHT HAND FROM OUTSIDE CRACK TO INSIDE CRACK
 */
class Path
{
    
public:

    Path();

    /**
     * Given unidimensional element reffers to the cracktip element that will be used
     * to compute J-integral around it.
     * Obs.: normal direction must be in xz plane and the arcs (internal and external) will be in (y>0).
     */
    Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim);
    ~Path();
    
    LinearPath * GetLinearPath()
    {
        return fLinearPath;
    }
    
    ArcPath * GetArcPath()
    {
        return fArcPath;
    }
    
    int MeshDim()
    {
        return fMeshDim;
    }
    
    
protected:
    
    LinearPath * fLinearPath;
    ArcPath * fArcPath;
    int fMeshDim;
};



class JIntegral
{
public:
    
    JIntegral();
    ~JIntegral();
    
    void PushBackPath(Path *pathElem);
    
    TPZVec<REAL> IntegratePath(int p);
    
private:
    
    TPZVec<Path*> fPathVec;
};



#endif
