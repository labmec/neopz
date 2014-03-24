/*
 *  TPZJIntegral.h
 *
 *  Created by Cesar Lucci on 25/04/12.
 *  Copyright 2010 LabMeC. All rights reserved.
 */

#ifndef PZ_TPZJIntegral_h
#define PZ_TPZJIntegral_h

#include <iostream>
#include "TPZGeoElement.h"
#include "pzcmesh.h"
#include "pzelasmat.h"
#include "pzelast3d.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"


/**
 * @brief Plane of the fracture.
 * @author Cesar Lucci (Caju)
 * @since 09/08/2010
 */


class LinearPath3D
{
public:
    
    LinearPath3D();//It is not to be used
    LinearPath3D(TPZCompMesh * cmeshElastic,
                 TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius);
    LinearPath3D(LinearPath3D * cp);
    ~LinearPath3D();
    
    void X(REAL t, TPZVec<REAL> & xt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    REAL Radius();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    virtual TPZVec<REAL> Func(REAL t);
    
    virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
    virtual REAL ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZVec<STATE> & Sigma_n);
    
    TPZCompMesh * CMeshElastic()
    {
        return this->fcmeshElastic;
    }
    
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
    
    /** CMesh that constains elastic data */
    TPZCompMesh * fcmeshElastic;
    
    /** map that holds t and respective elIndex from ElasticMesh and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIndexqsi_Elastic;
};


class ArcPath3D
{
public:
    
    ArcPath3D();//It is not to be used
    ArcPath3D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
    ArcPath3D(ArcPath3D * cp);
    ~ArcPath3D();
    
    void X(REAL t, TPZVec<REAL> & xt);
    void normalVec(REAL t, TPZVec<REAL> & n);
    
    REAL DETdxdt();
    
    TPZVec<REAL> & Origin()
    {
        return this->fOrigin;
    }
    TPZVec<REAL> & NormalPlane()
    {
        return this->fNormalDirection;
    }
    REAL Radius();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    virtual TPZVec<REAL> Func(REAL t);
    
    virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
    
    virtual REAL ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZFMatrix<STATE> & Sigma, TPZFMatrix<STATE> & strain);
    
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
    
    /** CMesh that constains elastic data */
    TPZCompMesh * fcmeshElastic;
    
    /** map that holds t and respective elIndex from ElasticMesh and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIndexqsi_Elastic;
};


class AreaPath3D
{
public:
    
    AreaPath3D();//It is not to be used
    AreaPath3D(LinearPath3D * givenLinearPath3D);
    ~AreaPath3D();
    
    REAL DETdxdt();
    
    TPZVec<REAL> operator()(REAL t)
    {
        TPZVec<REAL> Vval = Func(t);
        return Vval;
    }
    
    virtual TPZVec<REAL> Func(REAL t);
    
protected:
    
    struct LinearPath3D_2 : public LinearPath3D
    {
        public:
        LinearPath3D_2();
        LinearPath3D_2(LinearPath3D * cp);
        ~LinearPath3D_2();
        
        virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
        
        struct ArcPath3D_2 : public ArcPath3D
        {
            public:
            ArcPath3D_2();
            ArcPath3D_2(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
            ~ArcPath3D_2();
            
            virtual TPZVec<REAL> Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt);
            
            TPZVec<REAL> ComputeFiniteDifference(REAL t, TPZVec<REAL> xt, REAL halfDelta = 0.1);
            
            TPZVec<REAL> FunctionAux(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & direction);
        };
        
        ArcPath3D_2 * fArcPath3D;
    };
    
    LinearPath3D_2 * fLinearPath3D;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
};


/**
 *  ITS ALWAYS GOOD TO REMEMBER:
 *          THE CLASS Path3D CONSIDERS THAT THE NORMAL DIRECTION IS IN X,Z PLANE (JUST LIKE FRACTURE PLANE) AND
 *          THE ORIENTATION OF ARC and LINEAR stretch is:
 *              ARC : RIGHT HAND DIRECTION WITH RESPECT TO GIVEN NORMAL AXE (axe that defines the (orthogonal) arc plane).
 *              LINEAR: FROM THE END OF ARC (supposed to be inside crack plane) TO ORIGIN.
 * SO, ITS ALWAYS GOOD DEFINE NORMAL AXE TANGENT TO THE CRACK BOUNDARY, FOLLOWING RIGHT HAND FROM OUTSIDE CRACK TO INSIDE CRACK
 */
class Path3D
{
public:

    Path3D();

    /**
     * Given unidimensional element reffers to the cracktip element that will be used
     * to compute J-integral around it.
     * Obs.: normal direction must be in xz plane and the arcs (internal and external) will be in (y>0).
     */
    Path3D(TPZCompMesh * cmeshElastic,
           TPZVec<REAL> &Origin, REAL &KIc, int &myLayer, int &myStripe,
           TPZVec<REAL> &normalDirection, REAL radius);
    
    ~Path3D();
    
    void ComputeJIntegral();
    
    REAL OriginZcoord()
    {
        return this->fOriginZcoord;
    }
    
    TPZVec<REAL> & Origin()
    {
        return this->fArcPath3D->Origin();
    }
    
    void SetKI(REAL KI)
    {
        this->fKI = KI;
    }
    
    REAL KI()
    {
        return this->fKI;
    }
    
    REAL KIc()
    {
        return this->fKIc;
    }
    
    int MyLayer()
    {
        return this->fmyLayer;
    }
    
    int MyStripe()
    {
        return this->fmyStripe;
    }
    
    TPZVec<REAL> & JDirection()
    {
        return this->fJDirection;
    }
    
    TPZVec<REAL> & NormalPlane()
    {
        return this->fArcPath3D->NormalPlane();
    }
    
    REAL Jintegral()
    {
        return this->fJintegral;
    }
    
protected:
    
    LinearPath3D * fLinearPath3D;
    ArcPath3D * fArcPath3D;
    AreaPath3D * fAreaPath3D;
    
    TPZVec<REAL> fNormalDirection;
    REAL fOriginZcoord;
    
    REAL fKI;
    REAL fKIc;
    int fmyLayer;
    int fmyStripe;
    TPZVec<REAL> fJDirection;
    REAL fJintegral;
};


class JIntegral3D
{
public:
    
    JIntegral3D();
    ~JIntegral3D();
    
    void Reset();
    
    int NPaths();
    
    virtual void PushBackPath3D(Path3D *Path3DElem);
    
    virtual void IntegratePath3D();
    
    Path3D * Path(int p)
    {
        return this->fPath3DVec[p];
    }
    
    virtual void IntegratePath3D(int p);
    
private:
    
    TPZVec<Path3D*> fPath3DVec;
};


#endif
