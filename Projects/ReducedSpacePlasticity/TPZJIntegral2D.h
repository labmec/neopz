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


class LinearPath2D// : public LinearPath3D
{
public:
    
    LinearPath2D();//It is not to be used
    LinearPath2D(TPZCompMesh * cmeshElastic,
                 TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius);
    LinearPath2D(LinearPath2D * cp);
    ~LinearPath2D();
    
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
    virtual REAL ComputePressure(REAL t, TPZVec<REAL> & xt);
    
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



class ArcPath2D// : public ArcPath3D
{
public:
    
    ArcPath2D();//It is not to be used
    ArcPath2D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
    ArcPath2D(ArcPath2D * cp);
    ~ArcPath2D();
    
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


class Path2D
{
public:
    
    Path2D();
    
    /**
     * Given unidimensional element reffers to the cracktip element that will be used
     * to compute J-integral around it.
     * Obs.: normal direction must be in xz plane and the arcs (internal and external) will be in (y>0).
     */
    Path2D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius);
    ~Path2D();
    
    void ComputeJIntegral();
    
    REAL Jintegral()
    {
        return this->fJintegral;
    }
    
protected:
    
    LinearPath2D * fLinearPath2D;
    ArcPath2D * fArcPath2D;
    REAL fJintegral;
};




class JIntegral2D
{
public:
    
    JIntegral2D();
    ~JIntegral2D();
    
    virtual void SetPath2D(Path2D *Path2DElem);
    
    virtual void IntegratePath2D();
    
    Path2D * Path()
    {
        return this->fPath2D;
    }
    
private:
    
    Path2D * fPath2D;
};



#endif
