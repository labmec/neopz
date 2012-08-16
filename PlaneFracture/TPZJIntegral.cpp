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

const REAL Pi = 3.1415926535897932384626433832795;

const int dim3D = 3;


//--------------------------------------------------------class JPath


Path::Path()
{
    fOrigin.Resize(0);
    fNormalDirection.Resize(0);
    
    fr_int = 0.;
    fr_ext = 0.;
    
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
    #endif
    
    fr_int = r_int;
    fr_ext = r_ext;
    
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
}


TPZVec<REAL> Path::Func(REAL t)
{    
    TPZVec<REAL> xt(dim3D), dxdt(dim3D), nt(dim3D);
    REAL DETdxdt;
    this->X(t,xt);
    this->dXdt(t, dxdt, DETdxdt);
    this->normalVec(t, nt);

    TPZVec<REAL> qsi;
    
    TPZGeoEl * geoEl = NULL;
    if(fMeshDim == 2)
    {
        qsi.Resize(2, 0.);
        int axe0 = 0;//axe X
        int axe1 = 1;//axe Y
        int axeNormal = 2;//axe Z
        int elFoundId = TPZPlaneFracture::PointElementOnPlaneMesh(this->fcmesh->Reference(), fInitial2DElementId, xt, qsi, axe0, axe1, axeNormal, false);

        geoEl = this->fcmesh->Reference()->ElementVec()[elFoundId];
    }
    else if(fMeshDim == 3)
    {
        qsi.Resize(3, 0.);
        geoEl = TPZPlaneFracture::PointElementOnFullMesh(xt, qsi, fInitial2DElementId, this->fcmesh->Reference());
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
    
    REAL W = 0.;
    TPZFMatrix<REAL> Sigma(fMeshDim,fMeshDim), strain(fMeshDim,fMeshDim), GradU(fMeshDim,fMeshDim);
    GradU = data.dsol[0];
    
    if(fMeshDim == 2)
    {
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
        
        REAL poisson = elast2D->Nu();
        REAL young = elast2D->E();
        
        W = (((poisson - 1.)*strain(0,0)*strain(0,0) + (2.*poisson - 1.)*(strain(0,1)*strain(0,1) + strain(1,0)*strain(1,0)) - 4.*poisson*strain(0,0)*strain(1,1) + (poisson - 1.)*strain(1,1)*strain(1,1))*young)/(2.*(poisson - 1. + 2.*poisson*poisson));
    }
    else if(fMeshDim == 3)
    {
        TPZElasticity3D * elast3D = dynamic_cast<TPZElasticity3D *>(compEl->Material());
        
        #ifdef DEBUG
        if(!elast3D)
        {
            std::cout << "This material might be TPZElastMat3D type!\nSee " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
        #endif
        
        elast3D->ComputeStressTensor(Sigma, data);
        elast3D->ComputeStrainTensor(strain, GradU);
        REAL poisson = elast3D->GetPoisson();
        REAL young = elast3D->GetE();
        
        W = (((poisson - 1.)*strain(0,0)*strain(0,0) + (2.*poisson - 1.)*strain(0,1)*strain(0,1) - strain(0,2)*strain(0,2) - strain(1,0)*strain(1,0) - strain(1,1)*strain(1,1) - strain(1,2)*strain(1,2) - strain(2,0)*strain(2,0) - strain(2,1)*strain(2,1) + poisson*(2.*strain(0,2)*strain(0,2) + 2.*strain(1,0)*strain(1,0) + strain(1,1)*strain(1,1) + 2.*(strain(1,2)*strain(1,2) + strain(2,0)*strain(2,0) + strain(2,1)*strain(2,1))) - 4.*poisson*strain(1,1)*strain(2,2) + (poisson - 1.)*strain(2,2)*strain(2,2) - 4.*poisson*strain(0,0)*(strain(1,1) + strain(2,2)))*young)/(2.*(poisson - 1. + 2.*poisson*poisson));
    }
    TPZFMatrix<REAL> W_I(fMeshDim,fMeshDim,0.);
    for(int d = 0; d < fMeshDim; d++)
    {
        W_I(d,d) = W;
    }
    GradU.Transpose();
    TPZFMatrix<REAL> GradUtranspose_Sigma(fMeshDim,fMeshDim,0.);
    GradU.Multiply(Sigma, GradUtranspose_Sigma);
    
    TPZFMatrix<REAL> W_I_minus_GradUtranspose_Sigma(fMeshDim,fMeshDim,0.);
    W_I_minus_GradUtranspose_Sigma = W_I - GradUtranspose_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUtranspose_Sigma__n(fMeshDim,0.);
    for(int r = 0; r < fMeshDim; r++)
    {
        for(int c = 0; c < fMeshDim; c++)
        {
            W_I_minus_GradUtranspose_Sigma__n[r] += (W_I_minus_GradUtranspose_Sigma(r,c)*nt[c]) * DETdxdt;
        }
    }
    
    return W_I_minus_GradUtranspose_Sigma__n;
}


void Path::X(REAL t, TPZVec<REAL> & xt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::dXdt(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::normalVec(REAL t, TPZVec<REAL> & n)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}

//------------

void Path::X_line(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    TPZVec<REAL> initialPoint(dim3D);
    TPZVec<REAL> finalPoint(dim3D);
    
    if(__defaultDirection == EClockwise)
    {
        this->X_arc(+1., initialPoint, EInternalArcPath);
        this->X_arc(-1., finalPoint, EExternalArcPath);
    }
    else if(__defaultDirection == ECounterclockwise)
    {
        this->X_arc(+1., initialPoint, EExternalArcPath);
        this->X_arc(-1., finalPoint, EInternalArcPath);
    }
    
    for(int c = 0; c < dim3D; c++)
    {
        xt[c] = (1.-t)/2.*initialPoint[c] + (1.+t)/2.*finalPoint[c];
    }
}


void Path::dXdt_line(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    dxdt.Resize(dim3D,0.);
    
    TPZVec<REAL> x0(dim3D);
    TPZVec<REAL> x1(dim3D);
    
    
    if(__defaultDirection == EClockwise)
    {
        this->X_arc(+1., x0, EInternalArcPath);
        this->X_arc(-1., x1, EExternalArcPath);
    }
    else if(__defaultDirection == ECounterclockwise)
    {
        this->X_arc(+1., x0, EExternalArcPath);
        this->X_arc(-1., x1, EInternalArcPath);
    }
    
    DETdxdt = 0.;
    for(int c = 0; c < dim3D; c++)
    {
        dxdt[c] = (x1[c]-x0[c])/2.;
        DETdxdt += dxdt[c]*dxdt[c];
    }
    DETdxdt = sqrt(DETdxdt);
}


void Path::X_arc(REAL t, TPZVec<REAL> & xt, pathType pathT)
{
    xt.Resize(dim3D,0.);
    
    if(__defaultDirection == EClockwise)
    {
        if(pathT == EInternalArcPath)
        {
            xt[0] = (fOrigin[0] - fr_int*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
            xt[1] = fr_int*cos((Pi*t)/2.);
            xt[2] = (fOrigin[2] + fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
        }
        else if(pathT == EExternalArcPath)
        {                   
            xt[0] = (fOrigin[0] + fr_ext*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
            xt[1] = fr_ext*cos((Pi*t)/2.);
            xt[2] = (fOrigin[2] - fr_ext*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
        }
    }
    else if(__defaultDirection == ECounterclockwise)
    {
        if(pathT == EExternalArcPath)
        {
            xt[0] = (fOrigin[0] - fr_ext*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
            xt[1] = fr_ext*cos((Pi*t)/2.);
            xt[2] = (fOrigin[2] + fr_ext*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
        }
        else if(pathT == EInternalArcPath)
        {                   
            xt[0] = (fOrigin[0] + fr_int*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
            xt[1] = fr_int*cos((Pi*t)/2.);
            xt[2] = (fOrigin[2] - fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
        }
    }
}


void Path::dXdt_arc(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt, pathType pathT)
{
    dxdt.Resize(dim3D,0.);

    if(__defaultDirection == EClockwise)
    {
        if(pathT == EInternalArcPath)
        {
            dxdt[0] = -(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
            dxdt[2] = +(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            DETdxdt = Pi*fr_int/2.;
        }
        else if(pathT == EExternalArcPath)
        {                   
            dxdt[0] = +(Pi*fr_ext*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            dxdt[1] = -(Pi*fr_ext*sin((Pi*t)/2.))/2.;
            dxdt[2] = -(Pi*fr_ext*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            DETdxdt = Pi*fr_ext/2.;
        }
    }
    if(__defaultDirection == ECounterclockwise)
    {
        if(pathT == EExternalArcPath)
        {
            dxdt[0] = -(Pi*fr_ext*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            dxdt[1] = -(Pi*fr_ext*sin((Pi*t)/2.))/2.;
            dxdt[2] = +(Pi*fr_ext*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[0],fNormalDirection[2])))/2.;
            DETdxdt = Pi*fr_ext/2.;
        }
        else if(pathT == EInternalArcPath)
        {                   
            dxdt[0] = +(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
            dxdt[2] = -(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
            DETdxdt = Pi*fr_int/2.;
        }
    }
}


//--------------------------------------------------------class linearPath


void linearPath::X(REAL t, TPZVec<REAL> & xt)
{
    X_line(t, xt);
}


void linearPath::dXdt(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    dXdt_line(t, dxdt, DETdxdt);
}


void linearPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] = 0.;
    n[1] = -1.;
    n[2] = 0.;
}


//--------------------------------------------------------class externalArcPath


void externalArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    pathType pathT = EExternalArcPath;
    X_arc(t, xt, pathT);
}


void externalArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    pathType pathT = EExternalArcPath;
    dXdt_arc(t,dxdt, DETdxdt, pathT);
}


void externalArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> x(dim3D);
    X(t, x);
    
    REAL normaN = 0.;
    for(int i = 0; i < dim3D; i++)
    {
        normaN += (x[i] - fOrigin[i]) * (x[i] - fOrigin[i]);
    }
    normaN = sqrt(normaN);
    
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = +1. * fabs(x[i] - fOrigin[i])/normaN;
    }
}


//--------------------------------------------------------class internalArcPath


void internalArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    pathType pathT = EInternalArcPath;
    X_arc(t, xt, pathT);
}


void internalArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    pathType pathT = EInternalArcPath;
    dXdt_arc(t,dxdt, DETdxdt, pathT);
}


void internalArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> x(dim3D);
    X(t, x);
    
    REAL normaN = 0.;
    for(int i = 0; i < dim3D; i++)
    {
        normaN += (x[i] - fOrigin[i]) * (x[i] - fOrigin[i]);
    }
    normaN = sqrt(normaN);
    
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = -1. * fabs(x[i] - fOrigin[i])/normaN;
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
    
    linearPath * _LinearPath = new linearPath(jpathElem);
    externalArcPath * _ExtArcPath = new externalArcPath(jpathElem);
    internalArcPath * _IntArcPath = new internalArcPath(jpathElem);
    
    //////////////////////
//    TPZVec<REAL> funcAnsw(2,0.);
//    for(int i = -100; i <= 100; i++)
//    {
//        double tt = i/100.;
//        
//        funcAnsw = _LinearPath->Func(tt);
//        funcAnsw = _ExtArcPath->Func(tt);
//        funcAnsw = _IntArcPath->Func(tt);
//    }
    //////////////////////
    
    int meshDim = 2;
    TPZVec<REAL> integrLinPath = intRule.Vintegrate(*_LinearPath,meshDim,-1.,+1.);
    TPZVec<REAL> integrExtArc  = intRule.Vintegrate(*_ExtArcPath,meshDim,-1.,+1.);
    TPZVec<REAL> integrIntArc  = intRule.Vintegrate(*_IntArcPath,meshDim,-1.,+1.);

    TPZVec<REAL> answ(meshDim);
    if(meshDim == 2)
    {
        answ[0] = 2.*(integrLinPath[0] + integrExtArc[0] + integrIntArc[0]);
        answ[1] = 0.;
    }
    else if(meshDim == 3)
    {
        //Pela simetria do problema em relacao ao plano xz, deve-se somar a este vetor seu espelho em relacao ao plano xz.
        answ[0] = 2.*(integrLinPath[0] + integrExtArc[0] + integrIntArc[0]);
        answ[1] = 0.;
        answ[2] = 2.*(integrLinPath[2] + integrExtArc[2] + integrIntArc[2]);
    }
    
    return answ;
}





