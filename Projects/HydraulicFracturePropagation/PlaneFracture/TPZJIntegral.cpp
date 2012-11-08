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


REAL precisionIntegralRule = 1.E-7;

//--------------------------------------------------------class LinearPath


LinearPath::LinearPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


LinearPath::LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmesh = cmesh;
    fMeshDim = meshDim;
    fcrackPressure = fabs(pressure);
    
    f_t_elIdqsi.clear();
}


LinearPath::~LinearPath()
{
    fcmesh = NULL;
}


void LinearPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    TPZVec<REAL> initialPoint(3);
    initialPoint[0] = (fOrigin[0] - fradius*sin((Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    initialPoint[1] = fradius*cos((Pi)/2.);
    initialPoint[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi)/2.));
    
    for(int c = 0; c < 3; c++)
    {
        xt[c] = (1.-t)/2.*initialPoint[c] + (1.+t)/2.*fOrigin[c];
    }
}


//void LinearPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
//{
//    dxdt.Resize(3,0.);
//    
//    TPZVec<REAL> xarc(3);
//    t = 1;
//    xarc[0] = (fOrigin[0] - fradius*sin(Pi/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
//    xarc[1] = fradius*cos((Pi*t)/2.);
//    xarc[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
//    
//    for(int c = 0; c < 3; c++)
//    {
//        dxdt[c] = (fOrigin[c] - xarc[c])/2.;
//    }
//}


void LinearPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] =  0.;
    n[1] = -1.;
    n[2] =  0.;
}


REAL LinearPath::DETdxdt()
{
    return fDETdxdt;
}

//std::map<REAL,REAL> functionLinJx;
TPZVec<REAL> LinearPath::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);

    TPZVec<REAL> linContribution(fMeshDim,0.);
    linContribution = Function(t, xt, nt);
    
    if(fMeshDim == 2)
    {
        linContribution[0] = 2.*linContribution[0];
        linContribution[1] = 0.;
    }
    else if(fMeshDim == 3)
    {
        //Simetry in xz plane
        linContribution[0] = 2.*linContribution[0];
        linContribution[1] = 0.;
        linContribution[2] = 2.*linContribution[2];
    }
    
    //functionLinJx[xt[0]] = linContribution[0];
    
    return linContribution;
}


TPZVec<REAL> LinearPath::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> qsi(0);
    
    int InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi.lower_bound(t);
    if(it != f_t_elIdqsi.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmesh->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, InitialElementId, fMeshDim);

    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    
    f_t_elIdqsi[t] = std::make_pair(geoEl->Id(), qsi);
    
    TPZVec<REAL> minusGradUt_Sigma__n(fMeshDim,0.);
    
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
    
    TPZFMatrix<REAL> GradUtxy(fMeshDim,fMeshDim);
    GradUtxy.Zero();
    if(fMeshDim == 2)
    {
        TPZFMatrix<REAL> GradUtax(fMeshDim,fMeshDim);
        GradUtax = data.dsol[0];
        GradUtxy(0,0) = GradUtax(0,0)*data.axes(0,0) + GradUtax(1,0)*data.axes(1,0);
        GradUtxy(1,0) = GradUtax(0,0)*data.axes(0,1) + GradUtax(1,0)*data.axes(1,1);
        GradUtxy(0,1) = GradUtax(0,1)*data.axes(0,0) + GradUtax(1,1)*data.axes(1,0);
        GradUtxy(1,1) = GradUtax(0,1)*data.axes(0,1) + GradUtax(1,1)*data.axes(1,1);
    }
    else if(fMeshDim == 3)
    {
        GradUtxy = data.dsol[0];
    }
    
    TPZVec<REAL> Sigma_n(fMeshDim,0.);
    Sigma_n[1] = fcrackPressure;
    for(int r = 0; r < fMeshDim; r++)
    {
        for(int c = 0; c < fMeshDim; c++)
        {
            minusGradUt_Sigma__n[r] += -(GradUtxy(r,c)*Sigma_n[c]) * DETdxdt();
        }
    }
    
    return minusGradUt_Sigma__n;
}


REAL LinearPath::Origin(int c)
{
    REAL originCoordC = fOrigin[c];
    
    return originCoordC;
}


//--------------------------------------------------------class ArcPath


ArcPath::ArcPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


ArcPath::ArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = Pi*fradius/2.;
    
    fcmesh = cmesh;
    fMeshDim = meshDim;
    
    f_t_elIdqsi.clear();
}


ArcPath::ArcPath(ArcPath * cp)
{
    fOrigin = cp->fOrigin;
    fNormalDirection = cp->fNormalDirection;
    fradius = cp->fradius;
    
    fDETdxdt = cp->fDETdxdt;
    
    fcmesh = cp->fcmesh;
    fMeshDim = cp->fMeshDim;
    
    f_t_elIdqsi.clear();
}


ArcPath::~ArcPath()
{
    fcmesh = NULL;
}


void ArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    xt[0] = (fOrigin[0] - fradius*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fradius*cos((Pi*t)/2.);
    xt[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
}


//void ArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
//{
//    dxdt.Resize(3,0.);
//    
//    dxdt[0] = -(Pi*fradius*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
//    dxdt[1] = -(Pi*fradius*sin((Pi*t)/2.))/2.;
//    dxdt[2] = +(Pi*fradius*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[0],fNormalDirection[2])))/2.;
//}


void ArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(3);
    X(t, xt);
    for(int i = 0; i < 3; i++)
    {
        n[i] = (xt[i] - fOrigin[i])/fradius;
    }
}


REAL ArcPath::DETdxdt()
{
    return fDETdxdt;
}


TPZVec<REAL> ArcPath::Func(REAL t)
{
    TPZVec<REAL> xt(3), nt(3);
    
    X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(fMeshDim);
    arcContribution = Function(t, xt, nt);
    
    //Simetry in xz plane
    if(fMeshDim == 2)
    {
        arcContribution[0] = 2.*arcContribution[0];
        arcContribution[1] = 0.;
    }
    else if(fMeshDim == 3)
    {
        arcContribution[0] = 2.*arcContribution[0];
        arcContribution[1] = 0.;
        arcContribution[2] = 2.*arcContribution[2];
    }
    
    return arcContribution;
}


TPZVec<REAL> ArcPath::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> qsi(0);
    
    int InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi.lower_bound(t);
    if(it != f_t_elIdqsi.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmesh->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, InitialElementId, fMeshDim);

    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }

    f_t_elIdqsi[t] = std::make_pair(geoEl->Id(), qsi);
    
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
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*nt[c]) * DETdxdt();
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}


void ArcPath::SetRadius(REAL radius)
{
    fradius = radius;
    fDETdxdt = Pi*radius/2.;
    
    f_t_elIdqsi.clear();
}


REAL ArcPath::Radius()
{
    return fradius;
}


//--------------------------------------------------------class AreaPath


AreaPath::ArcPath2::ArcPath2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath::ArcPath2::ArcPath2(ArcPath * cp) : ArcPath(cp)
{
    
}

AreaPath::ArcPath2::~ArcPath2()
{
    fcmesh = NULL;
}

TPZVec<REAL> AreaPath::ArcPath2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> answ(3,0.);
    answ = ComputeFiniteDifference(t, xt, 0.005);
    
    return answ;
}

TPZVec<REAL> AreaPath::ArcPath2::ComputeFiniteDifference(REAL t, TPZVec<REAL> xt, REAL halfDelta)
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

TPZVec<REAL> AreaPath::ArcPath2::FunctionAux(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & direction)
{
    TPZVec<REAL> qsi(0);
    
    int InitialElementId = 0;
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi.lower_bound(t);
    if(it != f_t_elIdqsi.end())
    {
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else if(f_t_elIdqsi.size() > 0)
    {
        it--;
        InitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmesh->Reference()->ElementVec()[InitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, InitialElementId, 3);
    
    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    
    f_t_elIdqsi[t] = std::make_pair(geoEl->Id(), qsi);
    
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
    
    TPZFMatrix<REAL> Sigma(3,3), strain(3,3), GradUtxy(3,3);
    Sigma.Zero();
    strain.Zero();
    GradUtxy.Zero();
    
    if(fMeshDim == 2)
    {
        DebugStop();//Integral de area em simulacao 2D?
    }

    //else
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
            W_I_minus_GradUt_Sigma__n[r] += (W_I_minus_GradUt_Sigma(r,c)*direction[c]) * DETdxdt();
        }
    }
    
    return W_I_minus_GradUt_Sigma__n;
}


AreaPath::AreaPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


AreaPath::AreaPath(LinearPath * givenLinearPath, ArcPath * givenArcPath)
{
    fLinearPath = givenLinearPath;
    fArcPath = new ArcPath2(givenArcPath);
    
    #ifdef DEBUG
    if(!fLinearPath || !fArcPath)
    {
        DebugStop();//deu errado
    }
    #endif
}


AreaPath::~AreaPath()
{
    fLinearPath = NULL;
    fArcPath = NULL;
}


TPZVec<REAL> AreaPath::Integrate(int linNDiv)
{
    TPZVec<REAL> arcIntegr(3,0.), areaContribution(3,0.);
    
    REAL radius = fArcPath->Radius();
    REAL faixa = radius/linNDiv;
    REAL actRadius = 0.;
    
    for(int pos = 0; pos < linNDiv; pos++)
    {
        REAL t = -1. + (1. + 2.*pos)/linNDiv;
        TPZVec<REAL> xt(3,0.);
        fLinearPath->X(t, xt);
        
        actRadius = 0.;
        for(int c = 0; c < 3; c++)
        {
            actRadius += (xt[c] - fLinearPath->Origin(c))*(xt[c] - fLinearPath->Origin(c));
        }
        actRadius = sqrt(actRadius);
        
        fArcPath->SetRadius(actRadius);
        
        Adapt intRule(precisionIntegralRule);
        arcIntegr = intRule.Vintegrate(*fArcPath, 3, -1., 1.);
        
        areaContribution[0] += arcIntegr[0]*faixa;
        areaContribution[1]  = 0.;
        areaContribution[2] += arcIntegr[2]*faixa;
    }
    
    //CUIDADO!!!    A integral do arco jah contempla a simetria (jah multiplica x*2,y*0,z*2),
    //              nao precisando portanto fazer aqui novamente!!!
    return areaContribution;
}



//--------------------------------------------------------class JPath


Path::Path()
{
    fLinearPath = NULL;
    fArcPath = NULL;
    fMeshDim = 0;
}


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim)
{
    fLinearPath = new LinearPath(cmesh,Origin,normalDirection,radius,pressure,meshDim);
    fArcPath = new ArcPath(cmesh,Origin,normalDirection,radius,meshDim);
    fAreaPath = new AreaPath(fLinearPath,fArcPath);
    
    fMeshDim = meshDim;
    
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


Path::~Path()
{
    fLinearPath = NULL;
    fArcPath = NULL;
    
    fMeshDim = 0;
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

    
    Adapt intRule(precisionIntegralRule);
    
    int meshDim = jpathElem->MeshDim();
    
    TPZVec<REAL> linJintegral(meshDim,0.);
    linJintegral = intRule.Vintegrate(*(jpathElem->GetLinearPath()),meshDim,-1.,+1.);
    
    //4debug Jlin_x
    {
//        std::map<REAL,REAL>::iterator ii;
//        std::ofstream plotLinJx("PlotLinJx.txt");
//        plotLinJx.precision(10);
//        plotLinJx << "Jx = {";
//        for(ii = functionLinJx.begin(); ii != functionLinJx.end(); ii++)
//        {
//            plotLinJx << "{" << ii->first << "," << ii->second << "},";
//        }
//        plotLinJx << "};\n";
//        plotLinJx << "Jxplot = ListPlot[Jx,AxisOrigin->{0,0}]\n";
//        plotLinJx.close();
    }
    
    TPZVec<REAL> arcJintegral(meshDim,0.);
    arcJintegral = intRule.Vintegrate(*(jpathElem->GetArcPath()),meshDim,-1.,+1.);
    //
    TPZVec<REAL> areaJIntegral(meshDim,0.);
    //areaJIntegral = jpathElem->GetAreaPath()->Integrate(5);
    
    TPZVec<REAL> answ(meshDim);
    for(int i = 0; i < meshDim; i++)
    {
        answ[i] = linJintegral[i] + arcJintegral[i] + areaJIntegral[i];
    }
    
    if(meshDim==2)
    {
        std::cout << "J = " << answ[0] << "\n";
    }
    else
    {
        std::cout << "J = " << answ[0] << "\n";
    }
    
    return answ;
}

/*
 Franc2D : somente pressão de 5 no interior da fratura resulta em J = 0.0027
 */





