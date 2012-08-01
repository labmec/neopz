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

// *** ALL THIS INFORMATIONS IS WITH RESPECT TO THE AXIS STABLISHED BY THE UNIDIMENTIONAL ELEMENT ***
//by consequence:
//  - clockwise:
//      - external arc direction is ruled by lefthand
//      - internal arc direction is ruled by righthand
//      - linear is in fracture plane from internal arc to external arc
//
//  - counterclockwise:
//      - external arc direction is ruled by righthand
//      - internal arc direction is ruled by lefthand
//      - linear is in fracture plane from external arc to internal arc
enum pathDirection {EClockwise, ECounterclockwise};


//This internal variable could be changed that the code remains consintent
const int __defaultDirection = ECounterclockwise;

/**
 *  ITS ALWAYS GOOD TO REMEMBER:
 *          THE CLASS PATH CONSIDERS THAT THE 1D ELEMENT IS IN X,Z PLANE (JUST LIKE FRACTURE PLANE).
 *          THE ORIENTATION OF THIS ELEMENT DETERMINE THE ARC DIRECTION, 
 *          USED IN X_arc AND dXdt_arc METHODS AS A BOOLEAN PARAMETER.
 *              IF TRUE : RIGHT HAND DIRECTION WITH RESPECT TO AXES DEFINED BY 1D ELEMENT.
 *              IF FALSE : LEFT HAND DIRECTION WITH RESPECT TO AXES DEFINED BY 1D ELEMENT.
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
    Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, double r_int, double r_ext, int meshDim);
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
        fInitial2DElementId = cpPath.fInitial2DElementId;
        fcmesh = cpPath.fcmesh;
        fMeshDim = cpPath.fMeshDim;
    }
    
    TPZVec<REAL> operator()(double t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    virtual TPZVec<REAL> Func(double t);

    /////////////
    
    virtual void X(double t, TPZVec<REAL> & xt);
    virtual void dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt);
    virtual void normalVec(double t, TPZVec<REAL> & n);
    
        /**
         * For a given value of parametric variable t[-1,+1], this method computes the X(t)
         * that belongs to the line in plane fracture that connect the external and internal arcs.
         */
        void X_line(double t, TPZVec<REAL> & xt);
        void dXdt_line(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt);
        
        /**
         * For a given value of parametric variable t[-1,+1], this method computes the X(t)
         * that belongs to the arc internal or external (decided by pathT).
         */
        void X_arc(double t, TPZVec<REAL> & xt, pathType pathT);
        void dXdt_arc(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt, pathType pathT);
    
    
    /** initial node of given unidimensional element (this element define an axis) */
    TPZVec<REAL> fOrigin;
    
    /** final node of given unidimensional element (this element define an axis) */
    TPZVec<REAL> fNormalDirection;
    
    /** radius of internal and external arcs */
    double fr_int;
    double fr_ext;
    
    /** The Func() method need to ComputeXInverse to get solutions on compMesh.
     *  With the Id of initial element provided, near the searched point, the search algorithm could be optimized.
     *  This element must be an 2D element of PlaneFracture Mesh (see TPZPlaneFracture)*/
    int fInitial2DElementId;
    
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    int fMeshDim;
};


class linearPath : public Path
{
public:
    linearPath(Path * thePath)
    {
        fOrigin = thePath->fOrigin;
        fNormalDirection = thePath->fNormalDirection;
        fr_int = thePath->fr_int;
        fr_ext = thePath->fr_ext;
        fInitial2DElementId = thePath->fInitial2DElementId;
        fcmesh = thePath->fcmesh;
        fMeshDim = thePath->fMeshDim;
    }
    virtual void X(double t, TPZVec<REAL> & xt);
    virtual void dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt);
    virtual void normalVec(double t, TPZVec<REAL> & n);
};


class externalArcPath : public Path
{
public:
    externalArcPath(Path * thePath)
    {
        fOrigin = thePath->fOrigin;
        fNormalDirection = thePath->fNormalDirection;
        fr_int = thePath->fr_int;
        fr_ext = thePath->fr_ext;
        fInitial2DElementId = thePath->fInitial2DElementId;
        fcmesh = thePath->fcmesh;
        fMeshDim = thePath->fMeshDim;
    }
    virtual void X(double t, TPZVec<REAL> & xt);
    virtual void dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt);
    virtual void normalVec(double t, TPZVec<REAL> & n);
};


class internalArcPath : public Path
{
public:
    internalArcPath(Path * thePath)
    {
        fOrigin = thePath->fOrigin;
        fNormalDirection = thePath->fNormalDirection;
        fr_int = thePath->fr_int;
        fr_ext = thePath->fr_ext;
        fInitial2DElementId = thePath->fInitial2DElementId;
        fcmesh = thePath->fcmesh;
        fMeshDim = thePath->fMeshDim;
    }
    virtual void X(double t, TPZVec<REAL> & xt);
    virtual void dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt);
    virtual void normalVec(double t, TPZVec<REAL> & n);
};


class JIntegral
{
public:
    JIntegral();
    ~JIntegral();
    
    void PushBackPath(Path *pathElem);
    TPZVec<REAL> IntegratePath(int p);
    
    TPZVec<Path*> fPathVec;
};


#endif
