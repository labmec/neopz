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

enum pathType {ELinearPath, EExternalArcPath, EInternalArcPath};

const REAL Pi = 3.1415926535897932384626433832795;

/**
 *  ITS ALWAYS GOOD TO REMEMBER:
 *          THE CLASS PATH CONSIDERS THAT THE NORMAL DIRECTION IS IN X,Z PLANE (JUST LIKE FRACTURE PLANE) AND
 *          THE ORIENTATION OF THIS NORMAL DETERMINE THE ARC DIRECTION THAT RULES THE X_arc AND dXdt_arc METHODS, AS IT FOLLOWS:
 *              EXTERNAL ARC : RIGHT HAND DIRECTION WITH RESPECT TO GIVEN NORMAL AXE.
 *              INTERNAL ARC : LEFT HAND DIRECTION WITH RESPECT TO GIVEN NORMAL AXE.
 * SO, ITS ALWAYS GOOD DEFINE NORMAL AXE TANGENT TO THE CRACK BOUNDARY, FOLLOWING RIGHT HAND FOR THE EXTERNAL ARC BEGINNING AT OUTSIDE CRACK.
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
    Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim);
    Path(const Path &copy)
    {
        this->operator =(copy);
    }
    ~Path();
    
    void operator=(const Path &cpPath)
    {
        fOrigin = cpPath.fOrigin;
        fNormalDirection = cpPath.fNormalDirection;
        fr_int = cpPath.fr_int;
        fr_ext = cpPath.fr_ext;
        fDETdxdt = cpPath.fDETdxdt;
        fInitial2DElementId = cpPath.fInitial2DElementId;
        fcmesh = cpPath.fcmesh;
        fMeshDim = cpPath.fMeshDim;
    }
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    virtual TPZVec<REAL> Func(REAL t);
    
    /////////////
    
    virtual void X(REAL t, TPZVec<REAL> & xt);
    virtual void dXdt(REAL t, TPZVec<REAL> & dxdt);
    virtual void normalVec(REAL t, TPZVec<REAL> & n);    
    
    TPZVec<REAL> & Origin()
    {
        return fOrigin;
    }
    
    TPZVec<REAL> & NormalDirection()
    {
        return fNormalDirection;
    }
    
    REAL R_int()
    {
        return fr_int;
    }
    
    REAL R_ext()
    {
        return fr_ext;
    }
    
    REAL DETdxdt()
    {
        return fDETdxdt;
    }
    
    int Initial2DElementId()
    {
        return fInitial2DElementId;
    }
    
    TPZAutoPointer<TPZCompMesh> Cmesh()
    {
        return fcmesh;
    }
    
    int MeshDim()
    {
        return fMeshDim;
    }
    
protected:
    
    /** initial node of given unidimensional element (this element define an axis) */
    TPZVec<REAL> fOrigin;
    
    /** final node of given unidimensional element (this element define an axis) */
    TPZVec<REAL> fNormalDirection;
    
    /** radius of internal and external arcs */
    REAL fr_int;
    REAL fr_ext;
    
    /** determinant of [jacobian matrix](3x1) */
    REAL fDETdxdt;
    
    /** 
     *  The Func() method need to call ComputeXInverse method to get solutions on compMesh.
     *  With the Id of initial element provided, near the searched point, the search algorithm could be optimized.
     *  This element must be an 2D element of PlaneFracture Mesh (see TPZPlaneFracture)
     */
    int fInitial2DElementId;
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** For 2D problems (plane strain or plane stress), fMeshDim=2 */
    /** For 3D problems, fMeshDim=2 */
    int fMeshDim;
};


class linearPath : public Path
{
public:
    
    linearPath();//It is not to be used
    linearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim);
    linearPath(Path * thePath);
    
    virtual void X(REAL t, TPZVec<REAL> & xt);
    virtual void dXdt(REAL t, TPZVec<REAL> & dxdt);
    virtual void normalVec(REAL t, TPZVec<REAL> & n);
};


class externalArcPath : public Path
{
public:
    
    externalArcPath();//It is not to be used
    externalArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim);
    externalArcPath(Path * thePath);

    virtual void X(REAL t, TPZVec<REAL> & xt);
    virtual void dXdt(REAL t, TPZVec<REAL> & dxdt);
    virtual void normalVec(REAL t, TPZVec<REAL> & n);
};


class internalArcPath : public Path
{
public:
    
    internalArcPath();//It is not to be used
    internalArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim);
    internalArcPath(Path * thePath);

    virtual void X(REAL t, TPZVec<REAL> & xt);
    virtual void dXdt(REAL t, TPZVec<REAL> & dxdt);
    virtual void normalVec(REAL t, TPZVec<REAL> & n);
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
