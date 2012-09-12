//
//  TPZJIntegral.cpp
//  PZ
//
//  Created by Labmec on 25/04/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "TPZJIntegral.h"

#include "adapt.h"
#include "TPZPlaneFracture.h"
#include "pzelasmat.h"
#include "pzelast3d.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"
#include "pzaxestools.h"

#include "TPZVTKGeoMesh.h"

const int dim3D = 3;


//--------------------------------------------------------class JPath


Path::Path()
{
    fOrigin.Resize(0);
    fNormalDirection.Resize(0);
    
    fr_int = 0.;
    fr_ext = 0.;
    fDETdxdt = 0.;
    
    fInitial2DElementId = 0;
    
    fcmesh = NULL;
    
    fMeshDim = 0;
}


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    
    #ifdef DEBUG
    if(fabs(normalDirection[1]) > 1.E-8)
    {
        std::cout << "\nThe normal direction of J integral arc must be in XZ plane!!!\n";
        DebugStop();
    }
    if(r_ext <= r_int)
    {
        std::cout << "\nExternal radius must be higher than internal radius!!!\n";
        DebugStop();
    }
    #endif
    
    fr_int = r_int;
    fr_ext = r_ext;
    fDETdxdt = 0.;
    
    fInitial2DElementId = 0;
    
    fcmesh = cmesh;
    
    fMeshDim = meshDim;
}


Path::~Path()
{
    fOrigin.Resize(0);
    fNormalDirection.Resize(0);
    
    fr_int = 0.;
    fr_ext = 0.;
    fDETdxdt = 0.;
    
    fInitial2DElementId = -1;
    fcmesh = NULL;
    fMeshDim = -1;
}


TPZVec<REAL> Path::Func(REAL t)
{    
    TPZVec<REAL> xt(dim3D), nt(dim3D);
    X(t,xt);
    normalVec(t, nt);

    TPZVec<REAL> qsi;
    
    TPZGeoEl * geoEl = NULL;
    if(fMeshDim == 2)
    {
        qsi.Resize(2, 0.);
        int axe0 = 0;//axe X
        int axe1 = 1;//axe Y
        int axeNormal = 2;//axe Z
        int elFoundId = TPZPlaneFracture::PointElementOnPlaneMesh(fcmesh->Reference(), fInitial2DElementId, xt, qsi, axe0, axe1, axeNormal, false);

        geoEl = fcmesh->Reference()->ElementVec()[elFoundId];
    }
    else if(fMeshDim == 3)
    {
        qsi.Resize(3, 0.);
        geoEl = TPZPlaneFracture::PointElementOnFullMesh(xt, qsi, fInitial2DElementId, fcmesh->Reference());
    }
    else
    {
        std::cout << "Mesh dimension must be 2 or 3! See " << __PRETTY_FUNCTION__ << " !!!\n";
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
    
    TPZFMatrix<REAL> Sigma(fMeshDim,fMeshDim), strain(fMeshDim,fMeshDim), GradUtxy(fMeshDim,fMeshDim);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    if(fMeshDim == 2)
    {
        TPZFMatrix<REAL> GradUtax(fMeshDim,fMeshDim);
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
    else if(fMeshDim == 3)
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
    
    TPZFMatrix<REAL> GradUt_Sigma(fMeshDim,fMeshDim,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < fMeshDim; r++)
    {
        for(int c = 0; c < fMeshDim; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<REAL> W_I(fMeshDim,fMeshDim,0.);
    for(int d = 0; d < fMeshDim; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<REAL> W_I_minus_GradUt_Sigma(fMeshDim,fMeshDim,0.);
    W_I_minus_GradUt_Sigma = W_I - GradUt_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUt_Sigma__n(fMeshDim,0.);
    for(int r = 0; r < fMeshDim; r++)
    {
        for(int c = 0; c < fMeshDim; c++)
        {
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*nt[c]) * fDETdxdt;
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}


void Path::X(REAL t, TPZVec<REAL> & xt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::normalVec(REAL t, TPZVec<REAL> & n)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


//--------------------------------------------------------class linearPath


linearPath::linearPath() : Path()
{
    DebugStop();//Nao eh para usar construtor vazio
}


linearPath::linearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim) :
                                                                                        Path(cmesh, Origin, normalDirection, r_int, r_ext, meshDim)
{
    fDETdxdt = (fr_ext - fr_int)/2.;
}


linearPath::linearPath(Path * thePath)
{
    fOrigin = thePath->Origin();
    fNormalDirection = thePath->NormalDirection();
    fr_int = thePath->R_int();
    fr_ext = thePath->R_ext();
    fInitial2DElementId = thePath->Initial2DElementId();
    fcmesh = thePath->Cmesh();
    fMeshDim = thePath->MeshDim();
    
    fDETdxdt = (fr_ext - fr_int)/2.;
}

void linearPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    TPZVec<REAL> initialPoint(dim3D);
    initialPoint[0] = (fOrigin[0] - fr_ext*sin((Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    initialPoint[1] = fr_ext*cos((Pi)/2.);
    initialPoint[2] = (fOrigin[2] + fr_ext*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi)/2.));
    
    TPZVec<REAL> finalPoint(dim3D);
    finalPoint[0] = (fOrigin[0] + fr_int*sin((-Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    finalPoint[1] = fr_int*cos((-Pi)/2.);
    finalPoint[2] = (fOrigin[2] - fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((-Pi)/2.));
    
    for(int c = 0; c < dim3D; c++)
    {
        xt[c] = (1.-t)/2.*initialPoint[c] + (1.+t)/2.*finalPoint[c];
    }
}


void linearPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    dxdt.Resize(dim3D,0.);
    
    TPZVec<REAL> x0(dim3D);
    TPZVec<REAL> x1(dim3D);
    
    TPZVec<REAL> initialPoint(dim3D);
    initialPoint[0] = (fOrigin[0] - fr_ext*sin((Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    initialPoint[1] = fr_ext*cos((Pi)/2.);
    initialPoint[2] = (fOrigin[2] + fr_ext*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi)/2.));
    
    TPZVec<REAL> finalPoint(dim3D);
    finalPoint[0] = (fOrigin[0] + fr_int*sin((-Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    finalPoint[1] = fr_int*cos((-Pi)/2.);
    finalPoint[2] = (fOrigin[2] - fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((-Pi)/2.));
    
    for(int c = 0; c < dim3D; c++)
    {
        dxdt[c] = (x1[c]-x0[c])/2.;
    }
}


void linearPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] = 0.;
    n[1] = -1.;
    n[2] = 0.;
}


//--------------------------------------------------------class externalArcPath


externalArcPath::externalArcPath() : Path()
{
    DebugStop();//Nao eh para usar construtor vazio
}


externalArcPath::externalArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim) :
                                                                                                Path(cmesh, Origin, normalDirection, r_int, r_ext, meshDim)
{
    fDETdxdt = Pi*fr_ext/2.;
}


externalArcPath::externalArcPath(Path * thePath)
{
    fOrigin = thePath->Origin();
    fNormalDirection = thePath->NormalDirection();
    fr_int = thePath->R_int();
    fr_ext = thePath->R_ext();
    fInitial2DElementId = thePath->Initial2DElementId();
    fcmesh = thePath->Cmesh();
    fMeshDim = thePath->MeshDim();
    
    fDETdxdt = Pi*fr_ext/2.;
}

void externalArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    xt[0] = (fOrigin[0] - fr_ext*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fr_ext*cos((Pi*t)/2.);
    xt[2] = (fOrigin[2] + fr_ext*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
}


void externalArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    dxdt.Resize(dim3D,0.);
    
    dxdt[0] = -(Pi*fr_ext*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
    dxdt[1] = -(Pi*fr_ext*sin((Pi*t)/2.))/2.;
    dxdt[2] = +(Pi*fr_ext*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[0],fNormalDirection[2])))/2.;
}


void externalArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(dim3D);
    X(t, xt);
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = (xt[i] - fOrigin[i])/fr_ext;
    }
}


//--------------------------------------------------------class internalArcPath


internalArcPath::internalArcPath() : Path()
{
    DebugStop();//Nao eh para usar construtor vazio
}


internalArcPath::internalArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim) :
                                                                                                    Path(cmesh, Origin, normalDirection, r_int, r_ext, meshDim)
{
    fDETdxdt = Pi*fr_int/2.;
}


internalArcPath::internalArcPath(Path * thePath)
{
    fOrigin = thePath->Origin();
    fNormalDirection = thePath->NormalDirection();
    fr_int = thePath->R_int();
    fr_ext = thePath->R_ext();
    fInitial2DElementId = thePath->Initial2DElementId();
    fcmesh = thePath->Cmesh();
    fMeshDim = thePath->MeshDim();
    
    fDETdxdt = Pi*fr_int/2.;
}

void internalArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    xt[0] = (fOrigin[0] + fr_int*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fr_int*cos((Pi*t)/2.);
    xt[2] = (fOrigin[2] - fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
}


void internalArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    dxdt.Resize(dim3D,0.);
    
    dxdt[0] = +(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
    dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
    dxdt[2] = -(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
}


void internalArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(dim3D);
    X(t, xt);
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = (fOrigin[i] - xt[i])/fr_int;
    }
}


//--------------------------------------------------------class JIntegral


JIntegral::JIntegral()
{
    fPathVec.Resize(0);
}

JIntegral::~JIntegral()
{
    fPathVec.Resize(0);
}

void JIntegral::PushBackPath(Path * pathElem)
{
    int oldSize = fPathVec.NElements();
    fPathVec.Resize(oldSize+1);
    fPathVec[oldSize] = pathElem;
}

TPZVec<REAL> JIntegral::IntegratePath(int p)
{
    Path * jpathElem = fPathVec[p];

    REAL precisionIntegralRule = 1.E-15;
    Adapt intRule(precisionIntegralRule);
    
    externalArcPath * _ExtArcPath = new externalArcPath(jpathElem);
    //internalArcPath * _IntArcPath = new internalArcPath(jpathElem);
    //linearPath * _LinPath = new linearPath(jpathElem);
    
    ////////////////
    ////////////////
    
    int meshDim = 3;
    TPZVec<REAL> integrExtArc  = intRule.Vintegrate(*_ExtArcPath,meshDim,-1.,+1.);
    //TPZVec<REAL> integrIntArc  = intRule.Vintegrate(*_IntArcPath,meshDim,-1.,+1.);
    //TPZVec<REAL> lin  = intRule.Vintegrate(*_LinPath,meshDim,-1.,+1.);
    
    TPZVec<REAL> answ(meshDim);
    if(meshDim == 2)
    {
        answ[0] = 2.*integrExtArc[0];
        answ[1] = 0.;
    }
    else if(meshDim == 3)
    {
        //Pela simetria do problema em relacao ao plano xz, deve-se somar a este vetor seu espelho em relacao ao plano xz.
        answ[0] = 2.*integrExtArc[0];
        answ[1] = 0.;
        answ[2] = 2.*integrExtArc[2];
    }
    
    return answ;
}





