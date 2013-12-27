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
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"
#include "TPZTimer.h"

const REAL gIntegrPrecision = 1.e-4;
const REAL gScaleFactor = 1.e5;


//--------------------------------------------------------class LinearPath3D


LinearPath3D::LinearPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


LinearPath3D::LinearPath3D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
                           TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius)
{    
    fFinalPoint = FinalPoint;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmeshElastic = cmeshElastic;
    fcmeshFluid = cmeshFluid;
    //Se crackPressure acabar sendo positivo, teremos ponta da fratura fechando (imporei KI = 0.)!!!
    
    fInitialPoint.Resize(3, 0.);
    fInitialPoint[0] = (fFinalPoint[0] - fradius*sin((M_PI)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    fInitialPoint[1] = fradius*cos((M_PI)/2.); 
    fInitialPoint[2] = (fFinalPoint[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((M_PI)/2.));
    
    f_t_elIdqsi_Elastic.clear();
    f_t_elIdqsi_Fluid.clear();
}

LinearPath3D::LinearPath3D(LinearPath3D * cp)
{
    fInitialPoint = cp->fInitialPoint;
    fFinalPoint = cp->fFinalPoint;
    fNormalDirection = cp->fNormalDirection;
    fradius = cp->fradius;
    
    fDETdxdt = cp->fDETdxdt;
    
    fcmeshElastic = cp->fcmeshElastic;
    fcmeshFluid = cp->fcmeshFluid;
    
    f_t_elIdqsi_Elastic.clear();
    f_t_elIdqsi_Fluid.clear();
}

LinearPath3D::~LinearPath3D()
{
    fcmeshElastic = NULL;
}


void LinearPath3D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    for(int c = 0; c < 3; c++)
    {
        xt[c] = (1.-t)/2.*fInitialPoint[c] + (1.+t)/2.*fFinalPoint[c];
    }
}

void LinearPath3D::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] =  0.;
    n[1] = -1.;
    n[2] =  0.;
}

REAL LinearPath3D::DETdxdt()
{
    return fDETdxdt;
}

REAL LinearPath3D::Radius()
{
    return fradius;
}


#ifdef print_Jx_Linear
std::map<REAL,REAL> functionLinJx;
#endif

TPZVec<REAL> LinearPath3D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> linContribution(3,0.);
    linContribution = Function(t, xt, nt);
    
    linContribution[0] = linContribution[0] * gScaleFactor;
    linContribution[1] = 0.;
    linContribution[2] = linContribution[2] * gScaleFactor;
    
    #ifdef print_Jx_Linear
    functionLinJx[xt[0]] = linContribution[0];
    #endif
    
    return linContribution;
}


TPZVec<REAL> LinearPath3D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> GradUtxy(3,3);
    TPZVec<STATE> Sigma_n(3,0.);
    ComputeElasticData(t, xt, GradUtxy, Sigma_n);
    
    TPZVec<REAL> minusGradUt_Sigma__n(3,0.);
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            minusGradUt_Sigma__n[r] += -(GradUtxy(r,c)*Sigma_n[c]);
        }
    }
    
    return minusGradUt_Sigma__n;
}

void LinearPath3D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZVec<STATE> & Sigma_n)
{
    fcmeshElastic->LoadReferences();
    
    TPZVec<REAL> qsi(3,0.);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Elastic.lower_bound(t);
    if(it != f_t_elIdqsi_Elastic.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
//    else
//    {
//        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
//    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementId, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Elastic[t] = std::make_pair(geoEl->Id(), qsi);
    
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
    
    //remember: Dont need to multiply with axes once this
    //          last is IDENTITY in PZ (because is 3D element).
    GradUtxy = data.dsol[0];
    Sigma_n[1] = ComputePressure(t, xt);
}


REAL LinearPath3D::ComputePressure(REAL t, TPZVec<REAL> & xt)
{
    fcmeshFluid->LoadReferences();
    
    TPZVec<REAL> qsi(2,0.);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Fluid.lower_bound(t);
    if(it != f_t_elIdqsi_Fluid.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Fluid.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
//    else
//    {
//        qsi.Resize(fcmeshFluid->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
//    }
    TPZGeoEl * geoEl = fcmeshFluid->Reference()->FindElement(xt, qsi, InitialElementId, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Fluid[t] = std::make_pair(geoEl->Id(), qsi);
    
    TPZCompEl * cel = geoEl->Reference();
    if(!cel)
    {
        DebugStop();
    }
    TPZInterpolationSpace * sp = dynamic_cast <TPZInterpolationSpace*>(cel);
    if(!sp)
    {
        DebugStop();
    }
    
    TPZMaterialData data;
    sp->InitMaterialData(data);
    
    sp->ComputeShape(qsi, data);
    sp->ComputeSolution(qsi, data);
    
    REAL press = data.sol[0][0];
    
    return press;
}

//-------------------------class LinearPath2D

LinearPath2D::LinearPath2D() : LinearPath3D()
{
    
}

LinearPath2D::LinearPath2D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
                           TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius) :
              LinearPath3D(cmeshElastic,cmeshFluid,FinalPoint,normalDirection,radius)
{
    
}

LinearPath2D::LinearPath2D(LinearPath2D * cp) : LinearPath3D(cp)
{
    
}

LinearPath2D::~LinearPath2D()
{
    fcmeshElastic = NULL;
}
    
TPZVec<REAL> LinearPath2D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> linContribution(2,0.);
    linContribution = Function(t, xt, nt);
    
    linContribution[0] = linContribution[0] * gScaleFactor;
    linContribution[1] = 0.;
    
#ifdef print_Jx_Linear
    functionLinJx[xt[0]] = linContribution[0];
#endif
    
    return linContribution;
}

TPZVec<REAL> LinearPath2D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> GradUtxy(2,2);
    TPZVec<STATE> Sigma_n(2,0.);
    ComputeElasticData(t, xt, GradUtxy, Sigma_n);
    
    TPZVec<REAL> minusGradUt_Sigma__n(2,0.);
    for(int r = 0; r < 2; r++)
    {
        for(int c = 0; c < 2; c++)
        {
            minusGradUt_Sigma__n[r] += -(GradUtxy(r,c)*Sigma_n[c]);
        }
    }
    
    return minusGradUt_Sigma__n;
}

void LinearPath2D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZVec<STATE> & Sigma_n)
{
    fcmeshElastic->LoadReferences();
    
    TPZVec<REAL> qsi(0);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Elastic.lower_bound(t);
    if(it != f_t_elIdqsi_Elastic.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementId, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Elastic[t] = std::make_pair(geoEl->Id(), qsi);
    
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
    
    TPZFMatrix<STATE> GradUtax(2,2);
    GradUtax = data.dsol[0];
    GradUtxy(0,0) = GradUtax(0,0)*data.axes(0,0) + GradUtax(1,0)*data.axes(1,0);
    GradUtxy(1,0) = GradUtax(0,0)*data.axes(0,1) + GradUtax(1,0)*data.axes(1,1);
    GradUtxy(0,1) = GradUtax(0,1)*data.axes(0,0) + GradUtax(1,1)*data.axes(1,0);
    GradUtxy(1,1) = GradUtax(0,1)*data.axes(0,1) + GradUtax(1,1)*data.axes(1,1);
    
    Sigma_n[1] = ComputePressure(t, xt);
}


REAL LinearPath2D::ComputePressure(REAL t, TPZVec<REAL> & xt)
{
    fcmeshFluid->LoadReferences();
    
    TPZVec<REAL> qsi(0);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Fluid.lower_bound(t);
    if(it != f_t_elIdqsi_Fluid.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Fluid.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshFluid->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshFluid->Reference()->FindElement(xt, qsi, InitialElementId, 1);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Fluid[t] = std::make_pair(geoEl->Id(), qsi);
    
    TPZCompEl * cel = geoEl->Reference();
    if(!cel)
    {
        DebugStop();
    }
    TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*>(cel);
    if(!sp)
    {
        DebugStop();
    }

    TPZMaterialData data;
    sp->InitMaterialData(data);
        
    sp->ComputeShape(qsi, data);
    sp->ComputeSolution(qsi, data);
    
    REAL press = data.sol[0][0];
    
    return press;
}


//--------------------------------------------------------class ArcPath3D


ArcPath3D::ArcPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


ArcPath3D::ArcPath3D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = M_PI*fradius/2.;
    
    fcmeshElastic = cmeshElastic;
    
    f_t_elIdqsi_Elastic.clear();
}


ArcPath3D::ArcPath3D(ArcPath3D * cp)
{
    fOrigin = cp->fOrigin;
    fNormalDirection = cp->fNormalDirection;
    fradius = cp->fradius;
    
    fDETdxdt = cp->fDETdxdt;
    
    fcmeshElastic = cp->fcmeshElastic;
    
    f_t_elIdqsi_Elastic.clear();
}


ArcPath3D::~ArcPath3D()
{
    fcmeshElastic = NULL;
}


void ArcPath3D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    xt[0] = (fOrigin[0] - fradius*sin((M_PI*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fradius*cos((M_PI*t)/2.);
    xt[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((M_PI*t)/2.));
}


//void ArcPath3D::dXdt(REAL t, TPZVec<REAL> & dxdt)
//{
//    dxdt.Resize(3,0.);
//
//    dxdt[0] = -(M_PI*fradius*cos((M_PI*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
//    dxdt[1] = -(M_PI*fradius*sin((M_PI*t)/2.))/2.;
//    dxdt[2] = +(M_PI*fradius*cos((M_PI*t)/2.)*cos(atan2(fNormalDirection[0],fNormalDirection[2])))/2.;
//}


void ArcPath3D::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(3);
    X(t, xt);
    for(int i = 0; i < 3; i++)
    {
        n[i] = (xt[i] - fOrigin[i])/fradius;
    }
}


REAL ArcPath3D::DETdxdt()
{
    return fDETdxdt;
}


TPZVec<REAL> ArcPath3D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(3,0.);
    arcContribution = Function(t, xt, nt);
    
    arcContribution[0] = arcContribution[0] * gScaleFactor;
    arcContribution[1] = 0.;
    arcContribution[2] = arcContribution[2] * gScaleFactor;
    
    return arcContribution;
}


TPZVec<REAL> ArcPath3D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> Sigma(3,3), strain(3,3), GradUtxy(3,3);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    
    ComputeElasticData(t, xt, GradUtxy, Sigma, strain);
    
    TPZFMatrix<STATE> GradUt_Sigma(3,3,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<STATE> W_I(3,3,0.);
    for(int d = 0; d < 3; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<STATE> W_I_minus_GradUt_Sigma(3,3,0.);
    W_I_minus_GradUt_Sigma = W_I - GradUt_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUt_Sigma__n(3,0.);
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*nt[c]);
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}

void ArcPath3D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZFMatrix<STATE> & Sigma, TPZFMatrix<STATE> & strain)
{
    TPZVec<REAL> qsi(0);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Elastic.lower_bound(t);
    if(it != f_t_elIdqsi_Elastic.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementId, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Elastic[t] = std::make_pair(geoEl->Id(), qsi);
    
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

void ArcPath3D::SetRadius(REAL radius)
{
    fradius = radius;
    fDETdxdt = M_PI*radius/2.;
    
    f_t_elIdqsi_Elastic.clear();
}


REAL ArcPath3D::Radius()
{
    return fradius;
}

//-------------------------class ArcPath2D

ArcPath2D::ArcPath2D() : ArcPath3D()
{
    
}

ArcPath2D::ArcPath2D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius) :
           ArcPath3D(cmeshElastic,Origin,normalDirection,radius)
{
    
}

ArcPath2D::ArcPath2D(ArcPath2D * cp) : ArcPath3D(cp)
{
    
}

ArcPath2D::~ArcPath2D()
{
    fcmeshElastic = NULL;
}

TPZVec<REAL> ArcPath2D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(2,0.);
    arcContribution = Function(t, xt, nt);
    
    arcContribution[0] = arcContribution[0] * gScaleFactor;
    arcContribution[1] = 0.;
    
    return arcContribution;
}

TPZVec<REAL> ArcPath2D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> Sigma(2,2), strain(2,2), GradUtxy(2,2);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    
    ComputeElasticData(t, xt, GradUtxy, Sigma, strain);

    TPZFMatrix<STATE> GradUt_Sigma(2,2,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < 2; r++)
    {
        for(int c = 0; c < 2; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<STATE> W_I(2,2,0.);
    for(int d = 0; d < 2; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<STATE> W_I_minus_GradUt_Sigma(2,2,0.);
    W_I_minus_GradUt_Sigma = W_I - GradUt_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUt_Sigma__n(2,0.);
    for(int r = 0; r < 2; r++)
    {
        for(int c = 0; c < 2; c++)
        {
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*nt[c]);
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}

void ArcPath2D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZFMatrix<STATE> & Sigma, TPZFMatrix<STATE> & strain)
{
    TPZVec<REAL> qsi(0);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Elastic.lower_bound(t);
    if(it != f_t_elIdqsi_Elastic.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementId, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Elastic[t] = std::make_pair(geoEl->Id(), qsi);
    
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
    
    TPZFMatrix<STATE> GradUtax(2,2);
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
    
    TPZVec<STATE> Solout(3);
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


//--------------------------------------------------------class AreaPath3D


AreaPath3D::LinearPath3D_2::LinearPath3D_2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath3D::LinearPath3D_2::LinearPath3D_2(LinearPath3D * cp) : LinearPath3D(cp)
{
    fArcPath3D = new ArcPath3D_2(this->fcmeshElastic, this->fFinalPoint, this->fNormalDirection, this->fradius);
    
#ifdef DEBUG
    if(!fArcPath3D)
    {
        DebugStop();
    }
#endif
}

AreaPath3D::LinearPath3D_2::~LinearPath3D_2()
{
    fcmeshElastic = NULL;
}

TPZVec<REAL> AreaPath3D::LinearPath3D_2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> arcIntegral(3,0.);
    REAL arcRadius = 0.;
    for(int c = 0; c < 3; c++)
    {
        arcRadius += (fFinalPoint[c] - xt[c])*(fFinalPoint[c] - xt[c]);
    }
    arcRadius = sqrt(arcRadius);
    fArcPath3D->SetRadius(arcRadius);
    
    Adapt intRule(1.e-1);
    arcIntegral = intRule.Vintegrate(*(fArcPath3D),3,-1.,+1.);
    
    for(int i = 0; i < 3; i++)
    {
        arcIntegral[i] = arcIntegral[i] * fArcPath3D->DETdxdt();
    }
    
    return arcIntegral;
}

//------------

AreaPath3D::LinearPath3D_2::ArcPath3D_2::ArcPath3D_2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath3D::LinearPath3D_2::ArcPath3D_2::ArcPath3D_2(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius) :
ArcPath3D(cmeshElastic,Origin,normalDirection,radius)
{
    
}

AreaPath3D::LinearPath3D_2::ArcPath3D_2::~ArcPath3D_2()
{
    fcmeshElastic = NULL;
}

TPZVec<REAL> AreaPath3D::LinearPath3D_2::ArcPath3D_2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> answ(3,0.);
    answ = ComputeFiniteDifference(t, xt, 0.001);
    
    return answ;
}

TPZVec<REAL> AreaPath3D::LinearPath3D_2::ArcPath3D_2::ComputeFiniteDifference(REAL t, TPZVec<REAL> xt, REAL halfDelta)
{
    TPZVec<REAL> direction = this->fNormalDirection;
    REAL norm = 0;
    for(int i = 0; i < 3; i++)
    {
        norm += direction[i]*direction[i];
    }
    norm = sqrt(norm);
    
#ifdef DEBUG
    if(norm < 1.8-12)
    {
        DebugStop();
    }
#endif
    
    TPZVec<REAL> upPos(3,0.), downPos(3,0.);
    for (int i = 0; i < 3; i++)
    {
        upPos[i]    = xt[i] + halfDelta*direction[i]/norm;
        downPos[i]  = xt[i] - halfDelta*direction[i]/norm;
    }
    
    TPZVec<REAL> upSol(3,0.), downSol(3,0.);
    upSol = FunctionAux(t,upPos,direction);
    downSol = FunctionAux(t,downPos,direction);
    
    TPZVec<REAL> answ(3,0.);
    for(int i = 0; i < 3; i++)
    {
        answ[i] = (upSol[i] - downSol[i])/(2.*halfDelta);
    }
    
    return answ;
}

TPZVec<REAL> AreaPath3D::LinearPath3D_2::ArcPath3D_2::FunctionAux(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & direction)
{
    TPZVec<REAL> qsi(0);
    
    long InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi_Elastic.lower_bound(t);
    if(it != f_t_elIdqsi_Elastic.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementId, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIdqsi_Elastic[t] = std::make_pair(geoEl->Id(), qsi);
    
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
    
    TPZFMatrix<STATE> Sigma(3,3), strain(3,3), GradUtxy(3,3);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    
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
    
    TPZFMatrix<STATE> GradUt_Sigma(3,3,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<STATE> W_I(3,3,0.);
    for(int d = 0; d < 3; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<STATE> W_I_minus_GradUt_Sigma(3,3,0.);
    W_I_minus_GradUt_Sigma = W_I - GradUt_Sigma;
    
    TPZVec<REAL> W_I_minus_GradUt_Sigma__n(3,0.);
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*direction[c]);
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}

//------------

AreaPath3D::AreaPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


AreaPath3D::AreaPath3D(LinearPath3D * givenLinearPath3D)
{
    fLinearPath3D = new LinearPath3D_2(givenLinearPath3D);
    
    fDETdxdt = fLinearPath3D->DETdxdt();//O ArcPath3D_2 jah incluirah seu Detdxdt(), soh faltarah do LinearPath3D_2::Detdxdt()
    
    #ifdef DEBUG
    if(!fLinearPath3D)
    {
        DebugStop();
    }
    #endif
}


AreaPath3D::~AreaPath3D()
{
    fLinearPath3D = NULL;
}


REAL AreaPath3D::DETdxdt()
{
    return fDETdxdt;
}


TPZVec<REAL> AreaPath3D::Func(REAL t)
{
    TPZVec<REAL> areaContribution(3,0.);
    areaContribution = fLinearPath3D->Func(t);
    
    return areaContribution;
}


//--------------------------------------------------------class JPath3D


Path3D::Path3D()
{
    fLinearPath3D = NULL;
    fArcPath3D = NULL;
    fAreaPath3D = NULL;
    
    fOriginZcoord = 0.;
    
    fJDirection.Resize(3,0.);
    fJintegral = 0.;
}


Path3D::Path3D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
               TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fLinearPath3D = new LinearPath3D(cmeshElastic,cmeshFluid,Origin,normalDirection,radius);
    fArcPath3D = new ArcPath3D(cmeshElastic,Origin,normalDirection,radius);
    fAreaPath3D = new AreaPath3D(fLinearPath3D);
    
    fOriginZcoord = Origin[2];
    
    fJDirection.Resize(3,0.);
    fJintegral = 0.;
    
#ifdef DEBUG
    if(fabs(normalDirection[1]) > 1.E-8)
    {
        std::cout << "\nThe normal direction of J integral arc must be in XZ plane!!!\n";
        DebugStop();
    }
    if(radius < 1.e-3)
    {
        std::cout << "\nRadius too small!!!\n";
        DebugStop();
    }
#endif
}


Path3D::~Path3D()
{
    fLinearPath3D = NULL;
    fArcPath3D = NULL;
    fAreaPath3D = NULL;
    
    fOriginZcoord = 0.;
    
    fJDirection.Resize(0);
    fJintegral = 0.;
}

void Path3D::ComputeJIntegral()
{
    Adapt intRule(gIntegrPrecision);
    
    //------------------ LINE
    TPZTimer linInt("LinearIntegration");
    linInt.start();
    
    TPZVec<REAL> linJintegral(3,0.);
    linJintegral = intRule.Vintegrate(*(fLinearPath3D),3,-1.,+1.);
    
    //Simetry in xz plane
    linJintegral[0] = 2. * linJintegral[0] * fLinearPath3D->DETdxdt() / gScaleFactor;
    linJintegral[1] = 0.;
    linJintegral[2] = 2. * linJintegral[2] * fLinearPath3D->DETdxdt() / gScaleFactor;
    
    linInt.stop();
    
    //------------------ ARC
    TPZTimer arcInt("ArcIntegration");
    arcInt.start();
    
    TPZVec<REAL> arcJintegral(3,0.);
    arcJintegral = intRule.Vintegrate(*(fArcPath3D),3,-1.,+1.);
    
    //Simetry in xz plane
    arcJintegral[0] = 2. * arcJintegral[0] * fArcPath3D->DETdxdt() / gScaleFactor;
    arcJintegral[1] = 0.;
    arcJintegral[2] = 2. * arcJintegral[2] * fArcPath3D->DETdxdt() / gScaleFactor;
    
    arcInt.stop();
    
    //------------------ AREA
    TPZTimer areaInt("AreaIntegration");
    areaInt.start();
    
    TPZVec<REAL> areaJIntegral(3,0.);
    intRule.SetPrecision(1.e-2);
    //areaJIntegral = intRule.Vintegrate(*(fAreaPath3D),3,-1.,+1.);
    std::cout << "\nAQUICAJU : nao estah integrando na area!!!\n";
    
    //Simetry in xz plane
    areaJIntegral[0] = 2. * areaJIntegral[0] * fAreaPath3D->DETdxdt() / (gScaleFactor * gScaleFactor);
    areaJIntegral[1] = 0.;
    areaJIntegral[2] = 2. * areaJIntegral[2] * fAreaPath3D->DETdxdt() / (gScaleFactor * gScaleFactor);
    
    areaInt.stop();
    
    //------------------ COMBINING
    REAL Jx = linJintegral[0] + arcJintegral[0] + areaJIntegral[0];
    REAL Jy = linJintegral[1] + arcJintegral[1] + areaJIntegral[1];
    REAL Jz = linJintegral[2] + arcJintegral[2] + areaJIntegral[2];
    
    fJintegral = sqrt(Jx*Jx + Jy*Jy + Jz*Jz);
    
    //Normalizing
    fJDirection[0] = Jx/fJintegral;
    fJDirection[1] = Jy/fJintegral;
    fJDirection[2] = Jz/fJintegral;
    
    //std::cout << "J = { " << fJDirection[0] << " , " << fJDirection[1] << " , " << fJDirection[2] << " } --> " << "(" << fJintegral << ") : ";
}


Path2D::Path2D()
{
    fLinearPath2D = NULL;
    fArcPath2D = NULL;
    fJintegral = 0.;
}


Path2D::Path2D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
               TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fLinearPath2D = new LinearPath2D(cmeshElastic,cmeshFluid,Origin,normalDirection,radius);
    fArcPath2D = new ArcPath2D(cmeshElastic,Origin,normalDirection,radius);
    fJintegral = 0.;
    
#ifdef DEBUG
    if(fabs(normalDirection[1]) > 1.E-8)
    {
        std::cout << "\nThe normal direction of J integral arc must be in XZ plane!!!\n";
        DebugStop();
    }
    if(radius < 1.e-3)
    {
        std::cout << "\nRadius too small!!!\n";
        DebugStop();
    }
#endif
}


Path2D::~Path2D()
{
    fLinearPath2D = NULL;
    fArcPath2D = NULL;
    fJintegral = 0.;
}

void Path2D::ComputeJIntegral()
{
    Adapt intRule(gIntegrPrecision);
    
    //------------------ LINE
    TPZTimer linInt("LinearIntegration");
    linInt.start();
    
    TPZVec<REAL> linJintegral(2,0.);
    linJintegral = intRule.Vintegrate(*(fLinearPath2D),2,-1.,+1.);
    
    //Simetry in xz plane
    linJintegral[0] = 2. * linJintegral[0] * fLinearPath2D->DETdxdt() / gScaleFactor;
    linJintegral[1] = 0.;
    
    linInt.stop();
    
    //------------------ ARC
    TPZTimer arcInt("ArcIntegration");
    arcInt.start();
    
    TPZVec<REAL> arcJintegral(2,0.);
    arcJintegral = intRule.Vintegrate(*(fArcPath2D),2,-1.,+1.);
    
    //Simetry in xz plane
    arcJintegral[0] = 2. * arcJintegral[0] * fArcPath2D->DETdxdt() / gScaleFactor;
    arcJintegral[1] = 0.;
    
    arcInt.stop();

    //------------------ COMBINIG
    fJintegral = linJintegral[0] + arcJintegral[0];
    
    //    std::cout << "DeltaT integracao linha = " << linInt.seconds() << " s" << std::endl;
    //    std::cout << "DeltaT integracao arco = " << arcInt.seconds() << " s" << std::endl;
    //    std::cout << "J = " << answ[0] << "\n";
}


//--------------------------------------------------------class JIntegral


JIntegral3D::JIntegral3D()
{
    fPath3DVec.Resize(0);
}

JIntegral3D::~JIntegral3D()
{
    Reset();
}

void JIntegral3D::Reset()
{
    for(int j = 0; j < NPaths(); j++)
    {
        fPath3DVec[j] = NULL;
        delete fPath3DVec[j];
    }
    fPath3DVec.Resize(0);
}

int JIntegral3D::NPaths()
{
    return fPath3DVec.NElements();
}

void JIntegral3D::PushBackPath3D(Path3D * Path3DElem)
{
    int oldSize = fPath3DVec.NElements();
    fPath3DVec.Resize(oldSize+1);
    fPath3DVec[oldSize] = Path3DElem;
}

void JIntegral3D::IntegratePath3D()
{
    std::cout << "Computing J-integral...\n";
    for(int p = 0; p < NPaths(); p++)
    {
        IntegratePath3D(p);
        std::cout << p+1 << " of " << NPaths() << " computed!\n";
    }
}

void JIntegral3D::IntegratePath3D(int p)
{
    fPath3DVec[p]->ComputeJIntegral();
}

JIntegral2D::JIntegral2D()
{
    fPath2D = NULL;
}

JIntegral2D::~JIntegral2D()
{
    fPath2D = NULL;
}

void JIntegral2D::SetPath2D(Path2D * Path2DElem)
{
    fPath2D = Path2DElem;
}

void JIntegral2D::IntegratePath2D()
{
    fPath2D->ComputeJIntegral();
}



