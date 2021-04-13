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


LinearPath3D::LinearPath3D(TPZVec<REAL> &FinalPoint,
                           TPZVec<REAL> &normalDirection,
                           REAL radius,
                           TPZCompMesh * cmeshElastic,
                           TPZCompMesh * cmeshFluid)
{    
    this->fFinalPoint = FinalPoint;
    this->fPlaneNormalDirection = normalDirection;
    this->fradius = radius;
    
    this->fDETdxdt = this->fradius/2.;
    
    this->fcmeshElastic = cmeshElastic;
    this->fcmeshFluid = cmeshFluid;
    
    this->fInitialPoint.Resize(3, 0.);
    this->fInitialPoint[0] =
                (this->fFinalPoint[0] - this->fradius*sin((M_PI)/2.)*sin(atan2(this->fPlaneNormalDirection[2],this->fPlaneNormalDirection[0])));
    
    this->fInitialPoint[1] =
                this->fradius*cos((M_PI)/2.);
    
    this->fInitialPoint[2] =
                (this->fFinalPoint[2] + this->fradius*cos(atan2(this->fPlaneNormalDirection[2],this->fPlaneNormalDirection[0]))*sin((M_PI)/2.));
    
    this->f_t_elIndexqsi_Elastic.clear();
    this->f_t_elIndexqsi_Fluid.clear();
}

LinearPath3D::LinearPath3D(LinearPath3D * cp)
{
    this->fInitialPoint = cp->fInitialPoint;
    this->fFinalPoint = cp->fFinalPoint;
    this->fPlaneNormalDirection = cp->fPlaneNormalDirection;
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
    
    int64_t InitialElementIndex = 0;
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

    TPZGeoEl * geoEl = this->fcmeshElastic->Reference()->FindElementCaju(xt, qsi, InitialElementIndex, 3);
    
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
    
#ifdef PZDEBUG
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
        TPZGeoEl * geoEl2D = this->fcmeshElastic->Reference()->FindElementCaju(xt, qsi, InitialElementIndex, 2);
        int insideMatId = geoEl2D->MaterialId();
        
        if(globMaterialIdGen.IsInsideFractMat(insideMatId) == false)
        {//Quando t = +1, pode achar o elemento 2D de fora da fratura, vizinho do de dentro da fratura
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

        REAL prestress = -globLayerStruct.GetLayer(layer).fSigYY;
        
        bool useConstantP = true;
        if(useConstantP)
        {
            const TPZFMatrix<STATE> &cmeshElasticSol = this->fcmeshElastic->Solution();
            REAL totalPress = cmeshElasticSol.Get(1,0);//pressao constante = solucao elastica
            Sigma_n[1] = totalPress - prestress;
        }
        else//calcula utilizando pressao de fluidos
        {
            Sigma_n[1] = MAX(this->ComputeNetPressure(t, xt, prestress),0.);
        }
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
    
    int64_t InitialElementIndex = 0;
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

    TPZGeoEl * geoEl = this->fcmeshFluid->Reference()->FindElementCaju(xt, qsi, InitialElementIndex, 2);
    
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

bool LinearPath3D::ThereIsNegativeNetPressure()
{
    if(this->fcmeshFluid->Reference()->Reference() != this->fcmeshFluid)
    {
        this->fcmeshFluid->LoadReferences();
    }
    
    REAL t = 0.95;
    TPZVec<REAL> xt(3,0.);
    this->X(t,xt);
    
    REAL prestress = globLayerStruct.GetLowerPreStress();
    {
        REAL z = xt[2];
        int lay = globLayerStruct.WhatLayer(z);
        REAL SigYY = -globLayerStruct.GetLayer(lay).fSigYY;
        
        if(fabs(SigYY - prestress) > 1.E-3)
        {
            return false;
        }
    }
    REAL netPress = this->ComputeNetPressure(t,xt,prestress);
    
    if(netPress < 0.)
    {//net pressure is negative
        return true;
    }
    //else...
    return false;//net pressure is not negative
}

//--------------------------------------------------------class ArcPath3D

ArcPath3D::ArcPath3D()
{
    DebugStop();//Nao eh para usar construtor vazio
}


ArcPath3D::ArcPath3D(TPZVec<REAL> &Origin,
                     TPZVec<REAL> &normalDirection,
                     REAL radius,
                     TPZCompMesh * cmeshElastic)
{
    this->fOrigin = Origin;
    this->fPlaneNormalDirection = normalDirection;
    this->fradius = radius;
    
    this->fDETdxdt = M_PI*this->fradius/2.;
    
    this->fcmeshElastic = cmeshElastic;
    
    this->f_t_elIndexqsi_Elastic.clear();
}


ArcPath3D::ArcPath3D(ArcPath3D * cp)
{
    this->fOrigin = cp->fOrigin;
    this->fPlaneNormalDirection = cp->fPlaneNormalDirection;
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
    
    xt[0] = (this->fOrigin[0] - this->fradius*sin((M_PI*t)/2.)*sin(atan2(this->fPlaneNormalDirection[2],this->fPlaneNormalDirection[0])));
    xt[1] = this->fradius*cos((M_PI*t)/2.);
    xt[2] = (this->fOrigin[2] + this->fradius*cos(atan2(this->fPlaneNormalDirection[2],this->fPlaneNormalDirection[0]))*sin((M_PI*t)/2.));
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
    
    int64_t InitialElementIndex = 0;
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
    TPZGeoEl * geoEl = this->fcmeshElastic->Reference()->FindElementCaju(xt, qsi, InitialElementIndex, 3);
    
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
    
#ifdef PZDEBUG
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
    
#ifdef PZDEBUG
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

REAL ArcPath3D::Radius()
{
    return this->fradius;
}

//--------------------------------------------------------class JPath3D


Path3D::Path3D()
{
    this->fLinearPath3D = NULL;
    this->fArcPath3D = NULL;
    
    this->fPlaneNormalDirection.Resize(0,0.);
    this->fOriginZcoord = 0.;
    
    this->fKI = 0.;
    this->fKIc = 0.;
    
    this->fJDirection.Resize(3,0.);
    this->fJintegral = 0.;
}


Path3D::Path3D(TPZVec<REAL> &Origin, REAL &KIc,
               TPZVec<REAL> &normalDirection, REAL radius,
               TPZCompMesh * cmeshElastic,
               TPZCompMesh * cmeshFluid)
{
    this->fLinearPath3D = new LinearPath3D(Origin,normalDirection,radius,cmeshElastic,cmeshFluid);
    this->fArcPath3D = new ArcPath3D(Origin,normalDirection,radius,cmeshElastic);
    
    this->fPlaneNormalDirection = normalDirection;
    this->fOriginZcoord = Origin[2];
    
    this->fKI = 0.;
    this->fKIc = KIc;

    this->fJDirection.Resize(3,0.);
    this->fJintegral = 0.;
    
#ifdef PZDEBUG
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
    
    //------------------ COMBINING
    REAL Jx = linJintegral[0] + arcJintegral[0];
    REAL Jy = linJintegral[1] + arcJintegral[1];
    REAL Jz = linJintegral[2] + arcJintegral[2];
    
    //------------------- PROJETANDO O J VECTOR PARA A NORMAL EXTERNA
    REAL ex_ = this->fPlaneNormalDirection[0];
    REAL ez_ = this->fPlaneNormalDirection[2];
    REAL ex = ex_/sqrt(ex_*ex_ + ez_*ez_);
    REAL ez = ez_/sqrt(ex_*ex_ + ez_*ez_);
    
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

bool Path3D::ThereIsNegativeNetPressure()
{
    return this->fLinearPath3D->ThereIsNegativeNetPressure();
}

//--------------------------------------------------------class JIntegral


JIntegral3D::JIntegral3D()
{
    this->fPath3DVec.Resize(0);
}

JIntegral3D::~JIntegral3D()
{
    this->Reset();
}

void JIntegral3D::Reset()
{
    this->fPath3DVec.Resize(0);
}

int JIntegral3D::NPaths()
{
    return this->fPath3DVec.NElements();
}

void JIntegral3D::PushBackPath3D(Path3D * Path3DElem)
{
    int oldSize = this->fPath3DVec.NElements();
    this->fPath3DVec.Resize(oldSize+1);
    this->fPath3DVec[oldSize] = Path3DElem;
}

void JIntegral3D::IntegratePath3D()
{
    std::cout << "\n>>>> Computing J-integral (Tot = " << NPaths() << ")\n";
    
    for(int p = 0; p < NPaths(); p++)
    {
        this->IntegratePath3D(p);
    }
}

bool JIntegral3D::ThereIsNegativeNetPressure()
{
    for(int p = 0; p < NPaths(); p++)
    {
        bool thereIsNegNetPress = this->fPath3DVec[p]->ThereIsNegativeNetPressure();
        if(thereIsNegNetPress)
        {
            return true;
        }
    }
    
    return false;
}

void JIntegral3D::IntegratePath3D(int p)
{
    this->fPath3DVec[p]->ComputeJIntegral();
    std::cout << p+1 << " of " << NPaths() << " ok!\n";
    
    std::cout << "Jvec" << p << " = { "
              << this->fPath3DVec[p]->JDirection()[0] << " , "
              << this->fPath3DVec[p]->JDirection()[2] << " };\n"
              << "normJvec" << p << " = " << this->fPath3DVec[p]->Jintegral() << ";\n"
              << "KI = " << this->fPath3DVec[p]->KI() << "\n"
              << "KI/KIc = " << this->fPath3DVec[p]->KI()/this->fPath3DVec[p]->KIc() << "\n\n";
}


