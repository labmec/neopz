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
    LinearPath3D(TPZVec<REAL> &FinalPoint,
                 TPZVec<REAL> &normalDirection,
                 REAL radius,
                 TPZCompMesh * cmeshElastic,
                 TPZCompMesh * cmeshFluid);
    
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
    
    virtual REAL ComputeNetPressure(REAL t, TPZVec<REAL> & xt, REAL prestress);
    
    bool ThereIsNegativeNetPressure();//Eh Utilizado pelo kernel para verificar se a solucao de pressao eh valida.
    
    REAL Jradius()
    {
        return this->fradius;
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
    TPZVec<REAL> fPlaneNormalDirection;
    
    /** Radius of coupled arc */
    REAL fradius;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
    
    /** CMesh that constains elastic data */
    TPZCompMesh * fcmeshElastic;
    
    /** CMesh that constains fluid data */
    TPZCompMesh * fcmeshFluid;
    
    /** map that holds t and respective elIndex from ElasticMesh and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIndexqsi_Elastic;
    
    /** map that holds t and respective elIndex from FluidMesh and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIndexqsi_Fluid;
};


class ArcPath3D
{
public:
    
    ArcPath3D();//It is not to be used
    ArcPath3D(TPZVec<REAL> &Origin,
              TPZVec<REAL> &normalDirection,
              REAL radius,
              TPZCompMesh * cmeshElastic);
    
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
        return this->fPlaneNormalDirection;
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
    
protected:
    
    /** Origin of arc */
    TPZVec<REAL> fOrigin;
    
    /** This direction defines the arc plane.
     *  (this direction is orthogonal to arc plane and defines
     *   the right hand convention for the arc direction)
     */
    TPZVec<REAL> fPlaneNormalDirection;
    
    /** Radius of arc */
    REAL fradius;
    
    /** Determinant of dXdt(3x1) */
    REAL fDETdxdt;
    
    /** CMesh that constains elastic data */
    TPZCompMesh * fcmeshElastic;
    
    /** map that holds t and respective elIndex from ElasticMesh and qsi for ComputeXInverse optimization */
    std::map< REAL , std::pair< int , TPZVec<REAL> > > f_t_elIndexqsi_Elastic;
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
    Path3D(TPZVec<REAL> &Origin, REAL &KIc,
           TPZVec<REAL> &normalDirection, REAL radius,
           TPZCompMesh * cmeshElastic,
           TPZCompMesh * cmeshFluid);
    
    ~Path3D();
    
    void ComputeJIntegral();
    
    bool ThereIsNegativeNetPressure();
    
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
    
    TPZVec<REAL> fPlaneNormalDirection;
    REAL fOriginZcoord;
    
    REAL fKI;
    REAL fKIc;
    
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
    
    bool ThereIsNegativeNetPressure();
    
    Path3D * Path(int p)
    {
        return this->fPath3DVec[p];
    }
    
    virtual void IntegratePath3D(int p);
    
private:
    
    TPZVec<Path3D*> fPath3DVec;
};


#endif
