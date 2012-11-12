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


const REAL Pi = 3.1415926535897932384626433832795;


class LinearPath
{
public:
    
    LinearPath();//It is not to be used
    LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure);
    LinearPath(LinearPath * cp);
    ~LinearPath();
    
    void X(REAL t, TPZVec<REAL> & xt);
    //void dXdt(REAL t, TPZVec<REAL> & dxdt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    REAL Radius();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    TPZVec<REAL> Func(REAL t);
    
    virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
protected:
    
    /** Initial point of Line */
    TPZVec<REAL> fInitialPoint;
    
    /** Origin of coupled arc */
    TPZVec<REAL> fFinalPoint;
    
    /** This direction defines the coupled arc plane.
     *  (this direction is orthogonal to coupled arc plane and defines
     *   the right hand convention for the coupled arc direction)
     */
    TPZVec<REAL> fNormalDirection;
    
    /** Radius of coupled arc */
    REAL fradius;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** pressure applied inside fracture */
    REAL fcrackPressure;
    
    /** map that holds t and respective elId and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIdqsi;
};


class ArcPath
{
public:
    
    ArcPath();//It is not to be used
    ArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
    ArcPath(ArcPath * cp);
    ~ArcPath();
    
    void X(REAL t, TPZVec<REAL> & xt);
    //void dXdt(REAL t, TPZVec<REAL> & dxdt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    REAL Radius();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    TPZVec<REAL> Func(REAL t);
    
    virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
    void SetRadius(REAL radius);
    
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
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** map that holds t and respective elId and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIdqsi;
};



class AreaPath
{
public:
    
    AreaPath();//It is not to be used
    AreaPath(LinearPath * givenLinearPath);
    ~AreaPath();
    
    REAL DETdxdt();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    TPZVec<REAL> Func(REAL t);
    
protected:
    
    struct LinearPath2 : public LinearPath
    {
        public:
        LinearPath2();
        LinearPath2(LinearPath * cp);
        ~LinearPath2();
        
        virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
        
        struct ArcPath2 : public ArcPath
        {
            public:
            ArcPath2();
            ArcPath2(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
            ~ArcPath2();
            
            virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
            
            TPZVec<REAL> ComputeFiniteDifference(REAL t, TPZVec<REAL> xt, REAL halfDelta = 0.1);
            
            TPZVec<REAL> FunctionAux(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & direction);
        };
        
        ArcPath2 * fArcPath;
    };
    
    LinearPath2 * fLinearPath;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
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
    Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure);
    ~Path();
    
    LinearPath * GetLinearPath()
    {
        return fLinearPath;
    }
    
    ArcPath * GetArcPath()
    {
        return fArcPath;
    }
    
    AreaPath * GetAreaPath()
    {
        return fAreaPath;
    }
    
    
protected:
    
    LinearPath * fLinearPath;
    ArcPath * fArcPath;
    AreaPath * fAreaPath;
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
