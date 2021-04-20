 //
//  TPZJIntegral.cpp
//  PZ
//
//  Created by Labmec on 25/04/12.
//  Copyright (c) 2012 __MyCompanyName__. All rights reserved.
//

#include <iostream>

#include "TPZJIntegral2D.h"

#include "adapt.h"
#include "TPZVTKGeoMesh.h"
#include "pzintel.h"

const REAL gIntegrPrecision = 1.e-4;

//-------------------------class LinearPath2D

LinearPath2D::LinearPath2D()
{
    
}

LinearPath2D::LinearPath2D(TPZCompMesh * cmeshElastic,
                           TPZVec<REAL> &FinalPoint,
                           TPZVec<REAL> &normalDirection,
                           REAL radius)

{
    fFinalPoint = FinalPoint;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmeshElastic = cmeshElastic;
    //Se crackPressure acabar sendo positivo, teremos ponta da fratura fechando (imporei KI = 0.)!!!
    
    fInitialPoint.Resize(3, 0.);
    fInitialPoint[0] = (fFinalPoint[0] - fradius*sin((M_PI)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    fInitialPoint[1] = fradius*cos((M_PI)/2.);
    fInitialPoint[2] = (fFinalPoint[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((M_PI)/2.));
    
    f_t_elIndexqsi_Elastic.clear();
}

LinearPath2D::LinearPath2D(LinearPath2D * cp)
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
    
    linContribution[0] = linContribution[0];
    linContribution[1] = 0.;
    
    return linContribution;
}

void LinearPath2D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    for(int c = 0; c < 3; c++)
    {
        xt[c] = (1.-t)/2.*fInitialPoint[c] + (1.+t)/2.*fFinalPoint[c];
    }
}

void LinearPath2D::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] =  0.;
    n[1] = -1.;
    n[2] =  0.;
}

REAL LinearPath2D::DETdxdt()
{
    return this->fDETdxdt;
}

REAL LinearPath2D::Radius()
{
    return this->fradius;
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

REAL LinearPath2D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZVec<STATE> & Sigma_n)
{
    fcmeshElastic->LoadReferences();
    
    TPZVec<REAL> qsi(0);
    
    int64_t InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIndexqsi_Elastic.lower_bound(t);
    if(it != f_t_elIndexqsi_Elastic.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIndexqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementIndex]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIndexqsi_Elastic[t] = std::make_pair(geoEl->Index(), qsi);
    
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
    
    TPZFMatrix<STATE> GradUtax(2,2);
    GradUtax = data.dsol[0];
    GradUtxy(0,0) = GradUtax(0,0)*data.axes(0,0) + GradUtax(1,0)*data.axes(1,0);
    GradUtxy(1,0) = GradUtax(0,0)*data.axes(0,1) + GradUtax(1,0)*data.axes(1,1);
    GradUtxy(0,1) = GradUtax(0,1)*data.axes(0,0) + GradUtax(1,1)*data.axes(1,0);
    GradUtxy(1,1) = GradUtax(0,1)*data.axes(0,1) + GradUtax(1,1)*data.axes(1,1);
    
    Sigma_n[1] = ComputePressure(t, xt);
    
    return 0.;
}


REAL LinearPath2D::ComputePressure(REAL t, TPZVec<REAL> & xt)
{
    TPZFMatrix<STATE> &sol = this->fcmeshElastic->Solution();
    REAL press = sol(1,0);
    
    return press;
}

//-------------------------class ArcPath2D

ArcPath2D::ArcPath2D()
{
    
}

ArcPath2D::ArcPath2D(TPZCompMesh * cmeshElastic, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = M_PI*fradius/2.;
    
    fcmeshElastic = cmeshElastic;
    
    f_t_elIndexqsi_Elastic.clear();
}

ArcPath2D::ArcPath2D(ArcPath2D * cp)
{
    
}

ArcPath2D::~ArcPath2D()
{
    fcmeshElastic = NULL;
}

void ArcPath2D::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    xt[0] = (fOrigin[0] - fradius*sin((M_PI*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fradius*cos((M_PI*t)/2.);
    xt[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((M_PI*t)/2.));
}


void ArcPath2D::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(3);
    X(t, xt);
    for(int i = 0; i < 3; i++)
    {
        n[i] = (xt[i] - fOrigin[i])/fradius;
    }
}


REAL ArcPath2D::DETdxdt()
{
    return this->fDETdxdt;
}


TPZVec<REAL> ArcPath2D::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(2,0.);
    arcContribution = Function(t, xt, nt);
    
    arcContribution[0] = arcContribution[0];
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

REAL ArcPath2D::ComputeElasticData(REAL t, TPZVec<REAL> & xt, TPZFMatrix<STATE> & GradUtxy, TPZFMatrix<STATE> & Sigma, TPZFMatrix<STATE> & strain)
{
    TPZVec<REAL> qsi(0);
    
    int64_t InitialElementIndex = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIndexqsi_Elastic.lower_bound(t);
    if(it != f_t_elIndexqsi_Elastic.end())
    {
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIndexqsi_Elastic.size() > 0)
    {
        it--;
        InitialElementIndex = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmeshElastic->Reference()->ElementVec()[InitialElementIndex]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmeshElastic->Reference()->FindElement(xt, qsi, InitialElementIndex, 2);
    
    if(!geoEl)
    {
        std::cout.precision(15);
        std::cout << "\n\ngeoEl not found!\n";
        std::cout << "xt={ " << xt[0] << " , " << xt[1] << " , " << xt[2] << "};\n";
        std::cout << "See " << __PRETTY_FUNCTION__ << " !!!\n\n";
        DebugStop();
    }
    
    f_t_elIndexqsi_Elastic[t] = std::make_pair(geoEl->Index(), qsi);
    
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
    
    TPZFMatrix<STATE> GradUtax(2,2);
    GradUtax = data.dsol[0];
    GradUtxy(0,0) = GradUtax(0,0)*data.axes(0,0) + GradUtax(1,0)*data.axes(1,0);
    GradUtxy(1,0) = GradUtax(0,0)*data.axes(0,1) + GradUtax(1,0)*data.axes(1,1);
    GradUtxy(0,1) = GradUtax(0,1)*data.axes(0,0) + GradUtax(1,1)*data.axes(1,0);
    GradUtxy(1,1) = GradUtax(0,1)*data.axes(0,1) + GradUtax(1,1)*data.axes(1,1);
    
    TPZElasticityMaterial * elast2D = dynamic_cast<TPZElasticityMaterial *>(compEl->Material());
    
#ifdef PZDEBUG
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
    
    return 0.;
}


Path2D::Path2D()
{
    fLinearPath2D = NULL;
    fArcPath2D = NULL;
    fJintegral = 0.;
}


Path2D::Path2D(TPZCompMesh * cmeshElastic,
               TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fLinearPath2D = new LinearPath2D(cmeshElastic,Origin,normalDirection,radius);
    fArcPath2D = new ArcPath2D(cmeshElastic,Origin,normalDirection,radius);
    fJintegral = 0.;
    
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
    TPZVec<REAL> linJintegral(2,0.);
    linJintegral = intRule.Vintegrate(*(fLinearPath2D),2,-1.,+1.);
    
    //Simetry in xz plane
    linJintegral[0] = 2. * linJintegral[0] * fLinearPath2D->DETdxdt();
    linJintegral[1] = 0.;
    
    //------------------ ARC
    TPZVec<REAL> arcJintegral(2,0.);
    arcJintegral = intRule.Vintegrate(*(fArcPath2D),2,-1.,+1.);
    
    //Simetry in xz plane
    arcJintegral[0] = 2. * arcJintegral[0] * fArcPath2D->DETdxdt();
    arcJintegral[1] = 0.;

    //------------------ COMBINIG
    fJintegral = linJintegral[0] + arcJintegral[0];//<< in fact, fJintegral is (young * J-integral)
    
    //    std::cout << "DeltaT integracao linha = " << linInt.seconds() << " s" << std::endl;
    //    std::cout << "DeltaT integracao arco = " << arcInt.seconds() << " s" << std::endl;
    //    std::cout << "J = " << answ[0] << "\n";
}


//--------------------------------------------------------class JIntegral


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



