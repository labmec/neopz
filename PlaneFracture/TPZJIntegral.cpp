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
#include "pzelast3d.h"
#include "pzcompel.h"
#include "pzinterpolationspace.h"

const double Pi = 3.1415926535897932384626433832795;

const int dim3D = 3;


//--------------------------------------------------------class JPath


Path::Path()
{
    fInitialNode.Resize(0);
    fFinalNode.Resize(0);
    
    fr_int = 0.;
    fr_ext = 0.;
    
    fInitial2DElementId = 0;
    
    fcmesh = NULL;
}


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZGeoEl * el1D, double r_int, double r_ext)
{
#ifdef DEBUG
    if(!el1D || el1D->Dimension() != 1)
    {
        std::cout << "Null cracktip element or cracktip element given is not 1D!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " method!\n";
        DebugStop();
    }
    if(fabs(el1D->NodePtr(0)->Coord(1) > 1.E-9) || fabs(el1D->NodePtr(1)->Coord(1) > 1.E-9))
    {
        std::cout << "Given cracktip element is not on x,z plane (y coordinate = 0)!!\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " method!\n";
        DebugStop();
    }
#endif
    
    fInitialNode.Resize(dim3D);
    fInitialNode[0] = el1D->NodePtr(0)->Coord(0);
    fInitialNode[1] = el1D->NodePtr(0)->Coord(1);
    fInitialNode[2] = el1D->NodePtr(0)->Coord(2);
    
    fFinalNode.Resize(dim3D);
    fFinalNode[0] = el1D->NodePtr(1)->Coord(0);
    fFinalNode[1] = el1D->NodePtr(1)->Coord(1);
    fFinalNode[2] = el1D->NodePtr(1)->Coord(2);
    
    fr_int = r_int;
    fr_ext = r_ext;
    
    fcmesh = cmesh;
}


Path::~Path()
{
    fInitialNode.Resize(0);
    fFinalNode.Resize(0);
    
    fr_int = 0.;
    fr_ext = 0.;
}

TPZVec<REAL> Path::Func(double t)
{
    TPZVec<REAL> xt(dim3D), dxdt(dim3D), nt(dim3D);
    REAL DETdxdt;
    this->X(t,xt);
    this->dXdt(t, dxdt, DETdxdt);
    this->normalVec(t, nt);
    
    TPZGeoEl * geoEl = TPZPlaneFracture::PointElementOnFullMesh(xt, fInitial2DElementId, this->fcmesh->Reference());
    TPZVec<REAL> qsi(dim3D);
    bool isInsideDomain = geoEl->ComputeXInverse(xt, qsi);
    
    #ifdef DEBUG
    if(!isInsideDomain)
    {
        std::cout << "ComputeXInverse doen't work on " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    #endif
    
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

    TPZElasticity3D * elast3D = dynamic_cast<TPZElasticity3D *>(compEl->Material());
    
    #ifdef DEBUG
    if(!elast3D)
    {
        std::cout << "This material might be TPZElastMat3D type!\nSee " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
    }
    #endif
    
    TPZFMatrix<REAL> sigma, gradU = data.dsol[0];
    elast3D->ComputeStressTensor(sigma, data);
    
    
    double lambda = elast3D->GetLambda();
    double mu = elast3D->GetMU();
    
    double W = mu*(sigma(0,0)*sigma(0,0) + sigma(0,1)*sigma(0,1) + sigma(0,2)*sigma(0,2) + sigma(1,0)*sigma(1,0) 
                   + sigma(1,1)*sigma(1,1) + sigma(1,2)*sigma(1,2) + sigma(2,0)*sigma(2,0) + sigma(2,1)*sigma(2,1) + sigma(2,2)*sigma(2,2)) + 
    (lambda*(sigma(0,0)*sigma(0,0) + sigma(1,1)*sigma(1,1) + 4.*sigma(1,1)*sigma(2,2) + sigma(2,2)*sigma(2,2) + 4.*sigma(0,0)*(sigma(1,1) + sigma(2,2))))/2.;
    
    TPZVec<REAL> Vansw(0);
    
    return Vansw;
}


void Path::Print(std::ostream& out)
{
    out << "\n===========================================================\n";
    out << "Axis initial coordinates: ";
    out << "( " << fInitialNode[0] << " , " << fInitialNode[1] << " , " << fInitialNode[2] << " )\n";
    out << "Axis final coordinates: ";
    out << "( " << fFinalNode[0] << " , " << fFinalNode[1] << " , " << fFinalNode[2] << " )\n";
    out << "Internal radius = " << fr_int << std::endl;
    out << "External radius = " << fr_ext << std::endl;
}


void Path::X(double t, TPZVec<REAL> & xt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}


void Path::normalVec(double t, TPZVec<REAL> & n)
{
    std::cout << "The code should not enter here, but in derivated class!\n";
    DebugStop();
}

//------------

void Path::X_line(double t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D);
    
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


void Path::dXdt_line(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    dxdt.Resize(dim3D);
    
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


void Path::X_arc(double t, TPZVec<REAL> & xt, pathType pathT)
{
    xt.Resize(dim3D);
    
    double n0x = fInitialNode[0];
    double n0z = fInitialNode[2];
    
    double n1x = fFinalNode[0];
    double n1z = fFinalNode[2];
    
    if(__defaultDirection == EClockwise)
    {
        if(pathT == EInternalArcPath)
        {
            xt[0] = (n0x + n1x - 2.*fr_int*sin((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            xt[1] = fr_int*cos((Pi*t)/2.);
            xt[2] = (n0z + n1z + 2.*fr_int*cos(atan2(n1z-n0z,n1x-n0x))*sin((Pi*t)/2.))/2.;
        }
        else if(pathT == EExternalArcPath)
        {                   
            xt[0] = (n0x + n1x + 2.*fr_ext*sin((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            xt[1] = fr_ext*cos((Pi*t)/2.);
            xt[2] = (n0z + n1z - 2.*fr_ext*cos(atan2(n1z-n0z,n1x-n0x))*sin((Pi*t)/2.))/2.;
        }
    }
    else if(__defaultDirection == ECounterclockwise)
    {
        if(pathT == EExternalArcPath)
        {
            xt[0] = (n0x + n1x - 2.*fr_ext*sin((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            xt[1] = fr_ext*cos((Pi*t)/2.);
            xt[2] = (n0z + n1z + 2.*fr_ext*cos(atan2(n1z-n0z,n1x-n0x))*sin((Pi*t)/2.))/2.;
        }
        else if(pathT == EInternalArcPath)
        {                   
            xt[0] = (n0x + n1x + 2.*fr_int*sin((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            xt[1] = fr_int*cos((Pi*t)/2.);
            xt[2] = (n0z + n1z - 2.*fr_int*cos(atan2(n1z-n0z,n1x-n0x))*sin((Pi*t)/2.))/2.;
        }
    }
}


void Path::dXdt_arc(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt, pathType pathT)
{
    dxdt.Resize(dim3D);
    
    double n0x = fInitialNode[0];
    double n0z = fInitialNode[2];
    
    double n1x = fFinalNode[0];
    double n1z = fFinalNode[2];
    
    if(__defaultDirection == EClockwise)
    {
        if(pathT == EInternalArcPath)
        {
            dxdt[0] = -(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
            dxdt[2] = +(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(n1z-n0z,n1x-n0x)))/2.;
            DETdxdt = Pi*fr_int/2.;
        }
        else if(pathT == EExternalArcPath)
        {                   
            dxdt[0] = +(Pi*fr_ext*cos((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            dxdt[1] = -(Pi*fr_ext*sin((Pi*t)/2.))/2.;
            dxdt[2] = -(Pi*fr_ext*cos((Pi*t)/2.)*cos(atan2(n1z-n0z,n1x-n0x)))/2.;
            DETdxdt = Pi*fr_ext/2.;
        }
    }
    if(__defaultDirection == ECounterclockwise)
    {
        if(pathT == EExternalArcPath)
        {
            dxdt[0] = -(Pi*fr_ext*cos((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            dxdt[1] = -(Pi*fr_ext*sin((Pi*t)/2.))/2.;
            dxdt[2] = +(Pi*fr_ext*cos((Pi*t)/2.)*cos(atan2(n1x-n0x,n1z-n0z)))/2.;
            DETdxdt = Pi*fr_ext/2.;
        }
        else if(pathT == EInternalArcPath)
        {                   
            dxdt[0] = +(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(n1z-n0z,n1x-n0x)))/2.;
            dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
            dxdt[2] = -(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(n1z-n0z,n1x-n0x)))/2.;
            DETdxdt = Pi*fr_int/2.;
        }
    }
}


//--------------------------------------------------------class linearPath


void linearPath::X(double t, TPZVec<REAL> & xt)
{
    X_line(t, xt);
}


void linearPath::dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    dXdt_line(t, dxdt, DETdxdt);
}


void linearPath::normalVec(double t, TPZVec<REAL> & n)
{
    n[0] = 0.;
    n[1] = -1.;
    n[2] = 0.;
}


//--------------------------------------------------------class externalArcPath


void externalArcPath::X(double t, TPZVec<REAL> & xt)
{
    pathType pathT = EExternalArcPath;
    X_arc(t, xt, pathT);
}


void externalArcPath::dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    pathType pathT = EExternalArcPath;
    dXdt_arc(t,dxdt, DETdxdt, pathT);
}


void externalArcPath::normalVec(double t, TPZVec<REAL> & n)
{
    TPZVec<REAL> centerPoint(dim3D);
    
    for(int c = 0; c < dim3D; c++)
    {
        centerPoint[c] = (fInitialNode[c] + fFinalNode[c])/2.;
    }
    
    TPZVec<REAL> x(dim3D);
    X(t, x);
    
    double normaN = 0.;
    for(int i = 0; i < dim3D; i++)
    {
        normaN += (x[i] - centerPoint[i]) * (x[i] - centerPoint[i]);
    }
    normaN = sqrt(normaN);
    
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = +1. * fabs(x[i] - centerPoint[i])/normaN;
    }
}


//--------------------------------------------------------class internalArcPath


void internalArcPath::X(double t, TPZVec<REAL> & xt)
{
    pathType pathT = EInternalArcPath;
    X_arc(t, xt, pathT);
}


void internalArcPath::dXdt(double t, TPZVec<REAL> & dxdt, REAL & DETdxdt)
{
    pathType pathT = EInternalArcPath;
    dXdt_arc(t,dxdt, DETdxdt, pathT);
}


void internalArcPath::normalVec(double t, TPZVec<REAL> & n)
{
    TPZVec<REAL> centerPoint(dim3D);
    
    for(int c = 0; c < dim3D; c++)
    {
        centerPoint[c] = (fInitialNode[c] + fFinalNode[c])/2.;
    }
    
    TPZVec<REAL> x(dim3D);
    X(t, x);
    
    double normaN = 0.;
    for(int i = 0; i < dim3D; i++)
    {
        normaN += (x[i] - centerPoint[i]) * (x[i] - centerPoint[i]);
    }
    normaN = sqrt(normaN);
    
    for(int i = 0; i < dim3D; i++)
    {
        n[i] = -1. * fabs(x[i] - centerPoint[i])/normaN;
    }
}


//--------------------------------------------------------class JPath

JPath::JPath()
{
    //nothing to do
}

JPath::~JPath()
{
    //nothing to do
}

JPath::JPath(linearPath linearArc, externalArcPath extArc, internalArcPath intArc)
{
    this->fLinearPath = linearArc;
    this->fExtArcPath = extArc;
    this->fIntArcPath = intArc;    
}


//--------------------------------------------------------class JIntegral


JIntegral::JIntegral()
{
    fJPathVec.Resize(0);
}

JIntegral::~JIntegral()
{
    fJPathVec.Resize(0);
}

void JIntegral::PushBackJPath(JPath jpathElem)
{
    int oldSize = fJPathVec.NElements();
    fJPathVec.Resize(oldSize+1);
    fJPathVec[oldSize] = jpathElem;
}

TPZVec<REAL> JIntegral::IntegrateJPath(int p)
{
    JPath jpathElem = fJPathVec[p];

    double precisionIntegralRule = 1.E-100;
    Adapt intRule(precisionIntegralRule);
    
    TPZVec<REAL> integrLinPath = intRule.Vintegrate(jpathElem.fLinearPath,dim3D,Path::leftLimit(),Path::rightLimit());
    TPZVec<REAL> integrExtArc  = intRule.Vintegrate(jpathElem.fExtArcPath,dim3D,Path::leftLimit(),Path::rightLimit());
    TPZVec<REAL> integrIntArc  = intRule.Vintegrate(jpathElem.fIntArcPath,dim3D,Path::leftLimit(),Path::rightLimit());
    
    //multiplicado por 2 pois estou aproveitando a simetria do domínio em relacao ao plano xz.
    TPZVec<REAL> answ(dim3D);
    for(int d = 0; d < dim3D; d++)
    {
        answ[d] = 2.*(integrLinPath[d] + integrExtArc[d] + integrIntArc[d]);
    }
    
    return answ;
}





