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
    LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim);
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
     *  This element must be an 2D element of PlaneFracture Mesh (see TPZPlaneFracture)
     */
    int fInitial2DElementId;
    
    /** CMesh that constains data */
    TPZAutoPointer<TPZCompMesh> fcmesh;
    
    /** For 2D problems (plane strain or plane stress), fMeshDim=2 */
    /** For 3D problems, fMeshDim=3 */
    int fMeshDim;
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
     *  This element must be an 2D element of PlaneFracture Mesh (see TPZPlaneFracture)
     */
    int fInitial2DElementId;
    
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
    Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim);
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




inline
TPZVec<REAL> BoundaryFunc(TPZVec<REAL> & xt, TPZVec<REAL> & nt, REAL DETdxdt, int meshDim, TPZAutoPointer<TPZCompMesh> cmesh, int & initial2DElementId)
{
    TPZVec<REAL> qsi(0);
    
    TPZGeoEl * geoEl = NULL;
    if(meshDim == 2)
    {
        qsi.Resize(2, 0.);
        int axe0 = 0;//axe X
        int axe1 = 1;//axe Y
        int axeNormal = 2;//axe Z
        int elFoundId = TPZPlaneFracture::PointElementOnPlaneMesh(cmesh->Reference(), initial2DElementId, xt, qsi, axe0, axe1, axeNormal, false);
        
        geoEl = cmesh->Reference()->ElementVec()[elFoundId];
    }
    
    else if(meshDim == 3)
    {
        qsi.Resize(3, 0.);
        geoEl = TPZPlaneFracture::PointElementOnFullMesh(xt, qsi, initial2DElementId, cmesh->Reference());
    }
    else
    {
        std::cout << "Mesh dimension must be 2 or 3! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    
    TPZCompEl * compEl = geoEl->Reference();
    
#ifdef DEBUG
    if(!compEl)
    {
        std::cout << "Null compEl!\nSee " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
#endif
    
    TPZInterpolationSpace * intpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
    TPZMaterialData data;
    intpEl->InitMaterialData(data);
    
    intpEl->ComputeShape(qsi, data);
    intpEl->ComputeSolution(qsi, data);
    
    TPZFMatrix<REAL> Sigma(meshDim,meshDim), strain(meshDim,meshDim), GradUtxy(meshDim,meshDim);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    if(meshDim == 2)
    {
        TPZFMatrix<REAL> GradUtax(meshDim,meshDim);
        GradUtax = data.dsol[0];
        GradUtxy(0,0) = GradUtax(0,0)*data.axes(0,0) + GradUtax(1,0)*data.axes(1,0);
        GradUtxy(1,0) = GradUtax(0,0)*data.axes(0,1) + GradUtax(1,0)*data.axes(1,1);
        GradUtxy(0,1) = GradUtax(0,1)*data.axes(0,0) + GradUtax(1,1)*data.axes(1,0);
        GradUtxy(1,1) = GradUtax(0,1)*data.axes(0,1) + GradUtax(1,1)*data.axes(1,1);
        
        TPZElasticityMaterial * elast2D = dynamic_cast<TPZElasticityMaterial *>(compEl->Material());
        
#ifdef DEBUG
        if(!elast2D)
        {
            std::cout << "This material might be TPZElasticityMaterial type!\nSee " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
#endif
        
        TPZVec<REAL> Solout(3);
        int var;
        
        var = 10;//Stress Tensor
        elast2D->Solution(data, var, Solout);
        Sigma(0,0) = Solout[0];
        Sigma(1,1) = Solout[1];
        Sigma(0,1) = Solout[2];
        Sigma(1,0) = Solout[2];
        
        var = 11;//Strain Tensor
        elast2D->Solution(data, var, Solout);
        strain(0,0) = Solout[0];
        strain(1,1) = Solout[1];
        strain(0,1) = Solout[2];
        strain(1,0) = Solout[2];
    }
    else if(meshDim == 3)
    {
        GradUtxy = data.dsol[0];
        
        TPZElasticity3D * elast3D = dynamic_cast<TPZElasticity3D *>(compEl->Material());
        
#ifdef DEBUG
        if(!elast3D)
        {
            std::cout << "This material might be TPZElastMat3D type!\nSee " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
#endif
        
        elast3D->ComputeStressTensor(Sigma, data);
        elast3D->ComputeStrainTensor(strain, GradUtxy);
    }
    
    TPZFMatrix<REAL> GradUt_Sigma(meshDim,meshDim,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < meshDim; r++)
    {
        for(int c = 0; c < meshDim; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<REAL> W_I(meshDim,meshDim,0.);
    for(int d = 0; d < meshDim; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<REAL> W_I_minus_GradUt_Sigma(meshDim,meshDim,0.);
    W_I_minus_GradUt_Sigma = W_I - GradUt_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUt_Sigma__n(meshDim,0.);
    for(int r = 0; r < meshDim; r++)
    {
        for(int c = 0; c < meshDim; c++)
        {
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*nt[c]) * DETdxdt;
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}


#endif
