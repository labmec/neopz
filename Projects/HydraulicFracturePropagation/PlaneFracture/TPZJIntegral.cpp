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
#include "TPZPlaneFractureMesh.h"


const REAL gIntegrPrecision = 1.e-4;


//--------------------------------------------------------class LinearPath3D


LinearPath3D::LinearPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


LinearPath3D::LinearPath3D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
                           TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius)
{    
    this->fFinalPoint = FinalPoint;
    this->fNormalDirection = normalDirection;
    this->fradius = radius;
    
    this->fDETdxdt = this->fradius/2.;
    
    this->fcmeshElastic = cmeshElastic;
    this->fcmeshFluid = cmeshFluid;
    
    this->fInitialPoint.Resize(3, 0.);
    this->fInitialPoint[0] =
                (this->fFinalPoint[0] - this->fradius*sin((M_PI)/2.)*sin(atan2(this->fNormalDirection[2],this->fNormalDirection[0])));
    
    this->fInitialPoint[1] =
                this->fradius*cos((M_PI)/2.);
    
    this->fInitialPoint[2] =
                (this->fFinalPoint[2] + this->fradius*cos(atan2(this->fNormalDirection[2],this->fNormalDirection[0]))*sin((M_PI)/2.));
    
    this->f_t_elIndexqsi_Elastic.clear();
    this->f_t_elIndexqsi_Fluid.clear();
}

LinearPath3D::LinearPath3D(LinearPath3D * cp)
{
    this->fInitialPoint = cp->fInitialPoint;
    this->fFinalPoint = cp->fFinalPoint;
    this->fNormalDirection = cp->fNormalDirection;
    this->fradius = cp->fradius;
    
    this->fDETdxdt = cp->fDETdxdt;
    
    this->fcmeshElastic = cp->fcmeshElastic;
    this->fcmeshFluid = cp->fcmeshFluid;
    
    this->f_t_elIndexqsi_Elastic.clear();
    this->f_t_elIndexqsi_Fluid.clear();
}

LinearPath3D::~LinearPath3D()
{
    this->fcmeshElastic = NULL;
    this->fcmeshFluid = NULL;
}

void LinearPath3D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    for(int c = 0; c < 3; c++)
    {
        xt[c] = (1.-t)/2.*this->fInitialPoint[c] + (1.+t)/2.*this->fFinalPoint[c];
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
    return this->fDETdxdt;
}

REAL LinearPath3D::Radius()
{
    return this->fradius;
}

TPZVec<REAL> LinearPath3D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    normalVec(t, nt);
    
    TPZVec<REAL> linContribution(3,0.);
    linContribution = Function(t, xt, nt);
    
    linContribution[0] = linContribution[0];
    linContribution[1] = 0.;
    linContribution[2] = linContribution[2];
    
    return linContribution;
}


TPZVec<REAL> LinearPath3D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> GradUtxy(3,3,0.);
    TPZVec<STATE> Sigma_n(3,0.);
    REAL young = ComputeElasticData(t, xt, GradUtxy, Sigma_n);
    
    TPZVec<REAL> minusGradUt_Sigma__n(3,0.);
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            minusGradUt_Sigma__n[r] += -young*(GradUtxy(r,c)*Sigma_n[c]);
        }
    }
    
    return minusGradUt_Sigma__n;
}

REAL LinearPath3D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZVec<STATE> & Sigma_n)
{
    if(this->fcmeshElastic->Reference()->Reference() != this->fcmeshElastic)
    {
        this->fcmeshElastic->LoadReferences();
    }
    
    TPZVec<REAL> qsi(3,0.);
    
    long InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = this->f_t_elIndexqsi_Elastic.lower_bound(t);
    if(it != this->f_t_elIndexqsi_Elastic.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(this->f_t_elIndexqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }

    TPZGeoEl * geoEl = this->fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    this->f_t_elIndexqsi_Elastic[t] = std::make_pair(geoEl->Index(), qsi);
    
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
    
    {
        TPZGeoEl * geoEl2D = this->fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 2);
        int insideMatId = geoEl2D->MaterialId();
        
        if(globMaterialIdGen.IsInsideFractMat(insideMatId) == false)
        {
            bool found = false;
            for(int s = geoEl2D->NNodes(); s < geoEl2D->NSides(); s++)
            {
                if(found)
                {
                    break;
                }
                
                TPZGeoElSide elSide(geoEl2D,s);
                TPZGeoElSide neighSide(elSide.Neighbour());
                while(neighSide != elSide)
                {
                    if(globMaterialIdGen.IsInsideFractMat(neighSide.Element()->MaterialId()))
                    {
                        found = true;
                        insideMatId = neighSide.Element()->MaterialId();
                        break;
                    }
                    neighSide = neighSide.Neighbour();
                }
            }
        }
        if(globMaterialIdGen.IsInsideFractMat(insideMatId) == false)
        {
            std::cout << "\n\nNao achou vizinho dentro da fratura!!!";
            DebugStop();
        }
        
        int layer = globMaterialIdGen.WhatLayer(insideMatId);

/** Nesta abordagem, a pressao aplicada na parede da fratura eh a pressao de fluido : PIOR */
//        REAL prestress = -globLayerStruct.GetLayer(layer).fSigYY;
//        Sigma_n[1] = this->ComputeNetPressure(t,xt,prestress);
//        Se for utilizar isso, vá no método LinearPath3D::ComputeElasticData e descomente a linha "this->fcmeshElastic->LoadReferences();"

/** Nesta abordagem, a pressao aplicada na parede da fratura eh a pressao media : MELHOR */
        int stripe = globMaterialIdGen.WhatStripe(insideMatId);
        Sigma_n[1] = globElastReducedSolution.GetNetPressure(layer, stripe);
    }
    
    TPZElasticity3D * mat3d = dynamic_cast<TPZElasticity3D*>(compEl->Material());
    REAL young = mat3d->GetE();
    
    return young;
}


REAL LinearPath3D::ComputeNetPressure(REAL t, TPZVec<REAL> & xt, REAL prestress)
{
    if(this->fcmeshFluid->Reference()->Reference() != this->fcmeshFluid)
    {
        this->fcmeshFluid->LoadReferences();
    }
    
    TPZVec<REAL> qsi(2,0.);

    long InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = this->f_t_elIndexqsi_Fluid.lower_bound(t);
    if(it != this->f_t_elIndexqsi_Fluid.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(this->f_t_elIndexqsi_Fluid.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(this->fcmeshFluid->Reference()->ElementVec()[InitialElementIndex]->Dimension(),0.);
    }

    TPZGeoEl * geoEl = this->fcmeshFluid->Reference()->FindElement(xt, qsi, InitialElementIndex, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    this->f_t_elIndexqsi_Fluid[t] = std::make_pair(geoEl->Index(), qsi);
    
    TPZInterpolationSpace * intpEl = dynamic_cast<TPZInterpolationSpace *>(geoEl->Reference());
    TPZMaterialData data;
    intpEl->InitMaterialData(data);
    
    intpEl->ComputeShape(qsi, data);
    intpEl->ComputeSolution(qsi, data);
    
    REAL pressure = data.sol[0][0];
    
    return (pressure - prestress);
}

//--------------------------------------------------------class ArcPath3D


ArcPath3D::ArcPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


ArcPath3D::ArcPath3D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    this->fOrigin = Origin;
    this->fNormalDirection = normalDirection;
    this->fradius = radius;
    
    this->fDETdxdt = M_PI*this->fradius/2.;
    
    this->fcmeshElastic = cmeshElastic;
    
    this->f_t_elIndexqsi_Elastic.clear();
}


ArcPath3D::ArcPath3D(ArcPath3D * cp)
{
    this->fOrigin = cp->fOrigin;
    this->fNormalDirection = cp->fNormalDirection;
    this->fradius = cp->fradius;
    
    this->fDETdxdt = cp->fDETdxdt;
    
    this->fcmeshElastic = cp->fcmeshElastic;
    
    this->f_t_elIndexqsi_Elastic.clear();
}


ArcPath3D::~ArcPath3D()
{
    this->fcmeshElastic = NULL;
}


void ArcPath3D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    xt[0] = (this->fOrigin[0] - this->fradius*sin((M_PI*t)/2.)*sin(atan2(this->fNormalDirection[2],this->fNormalDirection[0])));
    xt[1] = this->fradius*cos((M_PI*t)/2.);
    xt[2] = (this->fOrigin[2] + this->fradius*cos(atan2(this->fNormalDirection[2],this->fNormalDirection[0]))*sin((M_PI*t)/2.));
}


void ArcPath3D::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(3);
    X(t, xt);
    for(int i = 0; i < 3; i++)
    {
        n[i] = (xt[i] - this->fOrigin[i])/this->fradius;
    }
}


REAL ArcPath3D::DETdxdt()
{
    return this->fDETdxdt;
}


TPZVec<REAL> ArcPath3D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(3,0.);
    arcContribution = Function(t, xt, nt);
    
    arcContribution[0] = arcContribution[0];
    arcContribution[1] = 0.;
    arcContribution[2] = arcContribution[2];
    
    return arcContribution;
}


TPZVec<REAL> ArcPath3D::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZFMatrix<STATE> Sigma(3,3), strain(3,3), GradUtxy(3,3);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    
    REAL young = ComputeElasticData(t, xt, GradUtxy, Sigma, strain);
    
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
            W_I_minus_GradUt_Sigma__n[r] += young * (W_I_minus_GradUt_Sigma(r,c)*nt[c]);
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}

REAL ArcPath3D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZFMatrix<STATE> & Sigma, TPZFMatrix<STATE> & strain)
{
    if(this->fcmeshElastic->Reference()->Reference() != this->fcmeshElastic)
    {
        this->fcmeshElastic->LoadReferences();
    }
    
    TPZVec<REAL> qsi(0);
    
    long InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = this->f_t_elIndexqsi_Elastic.lower_bound(t);
    if(it != this->f_t_elIndexqsi_Elastic.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(this->f_t_elIndexqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(this->fcmeshElastic->Reference()->ElementVec()[InitialElementIndex]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = this->fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    this->f_t_elIndexqsi_Elastic[t] = std::make_pair(geoEl->Index(), qsi);
    
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
    
    REAL young = elast3D->GetE();
    
    return young;
}

void ArcPath3D::SetRadius(REAL radius)
{
    this->fradius = radius;
    this->fDETdxdt = M_PI*radius/2.;
    
    this->f_t_elIndexqsi_Elastic.clear();
}


REAL ArcPath3D::Radius()
{
    return this->fradius;
}

//--------------------------------------------------------class AreaPath3D


AreaPath3D::LinearPath3D_2::LinearPath3D_2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath3D::LinearPath3D_2::LinearPath3D_2(LinearPath3D * cp) : LinearPath3D(cp)
{
    this->fArcPath3D = new ArcPath3D_2(this->fcmeshElastic, this->fFinalPoint, this->fNormalDirection, this->fradius);
    
#ifdef DEBUG
    if(!this->fArcPath3D)
    {
        DebugStop();
    }
#endif
}

AreaPath3D::LinearPath3D_2::~LinearPath3D_2()
{
    this->fcmeshElastic = NULL;
}

TPZVec<REAL> AreaPath3D::LinearPath3D_2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> arcIntegral(3,0.);
    REAL arcRadius = 0.;
    for(int c = 0; c < 3; c++)
    {
        arcRadius += (this->fFinalPoint[c] - xt[c])*(this->fFinalPoint[c] - xt[c]);
    }
    arcRadius = sqrt(arcRadius);
    this->fArcPath3D->SetRadius(arcRadius);
    
    Adapt intRule(1.e-1);
    arcIntegral = intRule.Vintegrate(*(this->fArcPath3D),3,-1.,+1.);
    
    for(int i = 0; i < 3; i++)
    {
        arcIntegral[i] = arcIntegral[i] * this->fArcPath3D->DETdxdt();
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
    this->fcmeshElastic = NULL;
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
    if(this->fcmeshElastic->Reference()->Reference() != this->fcmeshElastic)
    {
        this->fcmeshElastic->LoadReferences();
    }
    
    TPZVec<REAL> qsi(0);
    
    long InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = this->f_t_elIndexqsi_Elastic.lower_bound(t);
    if(it != this->f_t_elIndexqsi_Elastic.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(this->f_t_elIndexqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(this->fcmeshElastic->Reference()->ElementVec()[InitialElementIndex]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = this->fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 3);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    this->f_t_elIndexqsi_Elastic[t] = std::make_pair(geoEl->Index(), qsi);
    
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
    this->fLinearPath3D = new LinearPath3D_2(givenLinearPath3D);
    
    this->fDETdxdt = this->fLinearPath3D->DETdxdt();//O ArcPath3D_2 jah incluirah seu Detdxdt(), soh faltarah do LinearPath3D_2::Detdxdt()
    
    #ifdef DEBUG
    if(!this->fLinearPath3D)
    {
        DebugStop();
    }
    #endif
}


AreaPath3D::~AreaPath3D()
{
    this->fLinearPath3D = NULL;
}


REAL AreaPath3D::DETdxdt()
{
    return this->fDETdxdt;
}


TPZVec<REAL> AreaPath3D::Func(REAL t)
{
    TPZVec<REAL> areaContribution(3,0.);
    areaContribution = this->fLinearPath3D->Func(t);
    
    return areaContribution;
}


//--------------------------------------------------------class JPath3D


Path3D::Path3D()
{
    this->fLinearPath3D = NULL;
    this->fArcPath3D = NULL;
    this->fAreaPath3D = NULL;
    
    this->fNormalDirection.Resize(0,0.);
    this->fOriginZcoord = 0.;
    this->fKI = 0.;
    this->fKIc = 0.;
    this->fmyLayer = -1;
    this->fmyStripe = -1;
    
    this->fJDirection.Resize(3,0.);
    this->fJintegral = 0.;
}


Path3D::Path3D(TPZCompMesh * cmeshElastic, TPZCompMesh * cmeshFluid,
               TPZVec<REAL> &Origin, REAL &KIc, int &myLayer, int &myStripe,
               TPZVec<REAL> &normalDirection, REAL radius)
{
    this->fLinearPath3D = new LinearPath3D(cmeshElastic,cmeshFluid,Origin,normalDirection,radius);
    this->fArcPath3D = new ArcPath3D(cmeshElastic,Origin,normalDirection,radius);
    this->fAreaPath3D = new AreaPath3D(this->fLinearPath3D);
    
    this->fNormalDirection = normalDirection;
    this->fOriginZcoord = Origin[2];
    
    this->fKI = 0.;
    this->fKIc = KIc;
    this->fmyLayer = myLayer;
    this->fmyStripe = myStripe;
    this->fJDirection.Resize(3,0.);
    this->fJintegral = 0.;
    
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
    this->fLinearPath3D = NULL;
    this->fArcPath3D = NULL;
    this->fAreaPath3D = NULL;
    
    this->fOriginZcoord = 0.;
    
    this->fJDirection.Resize(0);
    this->fJintegral = 0.;
}

void Path3D::ComputeJIntegral()
{
    this->fJintegral = 0.;
    this->fKI = 0.;
    this->fJDirection.Fill(0.);
    
    Adapt intRule(gIntegrPrecision);
    
    //------------------ LINE
    TPZVec<REAL> linJintegral(3,0.);
    linJintegral = intRule.Vintegrate(*(this->fLinearPath3D),3,-1.,+1.);
    
    //Simetry in xz plane
    linJintegral[0] = 2. * linJintegral[0] * this->fLinearPath3D->DETdxdt();
    linJintegral[1] = 0.;
    linJintegral[2] = 2. * linJintegral[2] * this->fLinearPath3D->DETdxdt();
    
    //------------------ ARC
    TPZVec<REAL> arcJintegral(3,0.);
    arcJintegral = intRule.Vintegrate(*(this->fArcPath3D),3,-1.,+1.);
    
    //Simetry in xz plane
    arcJintegral[0] = 2. * arcJintegral[0] * this->fArcPath3D->DETdxdt();
    arcJintegral[1] = 0.;
    arcJintegral[2] = 2. * arcJintegral[2] * this->fArcPath3D->DETdxdt();
    
    //------------------ AREA
    TPZVec<REAL> areaJIntegral(3,0.);
    intRule.SetPrecision(1.e-2);
    //areaJIntegral = intRule.Vintegrate(*(this->fAreaPath3D),3,-1.,+0.9);
    //std::cout << "Nao estah integrando na area!!!\n";
    
    //Simetry in xz plane
    areaJIntegral[0] = 2. * areaJIntegral[0] * this->fAreaPath3D->DETdxdt();
    areaJIntegral[1] = 0.;
    areaJIntegral[2] = 2. * areaJIntegral[2] * this->fAreaPath3D->DETdxdt();
    
    //------------------ COMBINING
    REAL Jx = linJintegral[0] + arcJintegral[0] + areaJIntegral[0];
    REAL Jy = linJintegral[1] + arcJintegral[1] + areaJIntegral[1];
    REAL Jz = linJintegral[2] + arcJintegral[2] + areaJIntegral[2];
    
    //------------------- PROJETANDO O J VECTOR PARA A NORMAL EXTERNA
    REAL ex = this->fNormalDirection[0];
    REAL ez = this->fNormalDirection[2];
    ex = ex/sqrt(ex*ex + ez*ez);
    ez = ez/sqrt(ex*ex + ez*ez);
    
    REAL Jxproj = Jx - ex*(ex*Jx + ez*Jz);
    REAL Jyproj = Jy;
    REAL Jzproj = Jz - ez*(ex*Jx + ez*Jz);
    
    //------------------- SETANDO VALORES
    REAL Jprojnorm = sqrt(Jxproj*Jxproj + Jyproj*Jyproj + Jzproj*Jzproj);
    this->fJDirection[0] = Jxproj/Jprojnorm;
    this->fJDirection[1] = Jyproj/Jprojnorm;
    this->fJDirection[2] = Jzproj/Jprojnorm;
    this->fJintegral = Jprojnorm;
    
    this->fKI = sqrt(this->fJintegral);//<< in fact, this->fJintegral is (young * J-integral)
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
    std::cout << "\n>>>> Computing J-integral (Tot = " << NPaths() << ")\n";
    
    for(int p = 0; p < NPaths(); p++)
    {
        IntegratePath3D(p);
    }
}

void JIntegral3D::IntegratePath3D(int p)
{
    fPath3DVec[p]->ComputeJIntegral();
    std::cout << p+1 << " of " << NPaths() << " ok!\n";
    
    std::cout << "Jvec" << p << " = { "
              << fPath3DVec[p]->JDirection()[0] << " , "
              << fPath3DVec[p]->JDirection()[2] << " };\n"
              << "normJvec" << p << " = " << fPath3DVec[p]->Jintegral() << ";\n"
              << "KI = " << fPath3DVec[p]->KI() << "\n"
              << "KI/KIc = " << fPath3DVec[p]->KI()/fPath3DVec[p]->KIc() << "\n\n";
}


