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
#include "TPZTimer.h"

const REAL gIntegrPrecision = 1.e-5;
const REAL gScaleFactor = 1.e5;


//--------------------------------------------------------class LinearPath


LinearPath::LinearPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


LinearPath::LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &FinalPoint, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure)
{    
    fFinalPoint = FinalPoint;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmesh = cmesh;
    fcrackPressure = fabs(pressure);
    
    fInitialPoint.Resize(3, 0.);
    fInitialPoint[0] = (fFinalPoint[0] - fradius*sin((Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    fInitialPoint[1] = fradius*cos((Pi)/2.);
    fInitialPoint[2] = (fFinalPoint[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi)/2.));
    
    f_t_elIdqsi.clear();
}

LinearPath::LinearPath(LinearPath * cp)
{
    fInitialPoint = cp->fInitialPoint;
    fFinalPoint = cp->fFinalPoint;
    fNormalDirection = cp->fNormalDirection;
    fradius = cp->fradius;
    
    fDETdxdt = cp->fDETdxdt;
    
    fcmesh = cp->fcmesh;
    fcrackPressure = cp->fcrackPressure;
    
    f_t_elIdqsi.clear();
}

LinearPath::~LinearPath()
{
    fcmesh = NULL;
}


void LinearPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(3,0.);
    
    for(int c = 0; c < 3; c++)
    {
        xt[c] = (1.-t)/2.*fInitialPoint[c] + (1.+t)/2.*fFinalPoint[c];
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

REAL LinearPath::Radius()
{
    return fradius;
}


#ifdef print_Jx_Linear
std::map<REAL,REAL> functionLinJx;
#endif

TPZVec<REAL> LinearPath::Func(REAL t)
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
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, InitialElementId, 3);
    
    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    
    f_t_elIdqsi[t] = std::make_pair(geoEl->Id(), qsi);
    
    TPZVec<REAL> minusGradUt_Sigma__n(3,0.);
    
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
    
    TPZFMatrix<REAL> GradUtxy(3,3);
    GradUtxy.Zero();
    
    GradUtxy = data.dsol[0];
    
    TPZVec<REAL> Sigma_n(3,0.);
    Sigma_n[1] = fcrackPressure;
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            minusGradUt_Sigma__n[r] += -(GradUtxy(r,c)*Sigma_n[c]);
        }
    }
    
    return minusGradUt_Sigma__n;
}


//--------------------------------------------------------class ArcPath


ArcPath::ArcPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


ArcPath::ArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = Pi*fradius/2.;
    
    fcmesh = cmesh;
    
    f_t_elIdqsi.clear();
}


ArcPath::ArcPath(ArcPath * cp)
{
    fOrigin = cp->fOrigin;
    fNormalDirection = cp->fNormalDirection;
    fradius = cp->fradius;
    
    fDETdxdt = cp->fDETdxdt;
    
    fcmesh = cp->fcmesh;
    
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
    
    TPZVec<REAL> arcContribution(3);
    arcContribution = Function(t, xt, nt);
    
    arcContribution[0] = arcContribution[0] * gScaleFactor;
    arcContribution[1] = 0.;
    arcContribution[2] = arcContribution[2] * gScaleFactor;
    
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
    
    TPZFMatrix<REAL> GradUt_Sigma(3,3,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<REAL> W_I(3,3,0.);
    for(int d = 0; d < 3; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<REAL> W_I_minus_GradUt_Sigma(3,3,0.);
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

AreaPath::LinearPath2::LinearPath2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath::LinearPath2::LinearPath2(LinearPath * cp) : LinearPath(cp)
{
    fArcPath = new ArcPath2(this->fcmesh, this->fFinalPoint, this->fNormalDirection, this->fradius);
    
#ifdef DEBUG
    if(!fArcPath)
    {
        DebugStop();
    }
#endif
}

AreaPath::LinearPath2::~LinearPath2()
{
    fcmesh = NULL;
}

TPZVec<REAL> AreaPath::LinearPath2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> arcIntegral(3,0.);
    REAL arcRadius = 0.;
    for(int c = 0; c < 3; c++)
    {
        arcRadius += (fFinalPoint[c] - xt[c])*(fFinalPoint[c] - xt[c]);
    }
    arcRadius = sqrt(arcRadius);
    fArcPath->SetRadius(arcRadius);
    
    Adapt intRule(1.e-1);
    arcIntegral = intRule.Vintegrate(*(fArcPath),3,-1.,+1.);
    
    for(int i = 0; i < 3; i++)
    {
        arcIntegral[i] = arcIntegral[i] * fArcPath->DETdxdt();
    }
    
    return arcIntegral;
}

//------------

AreaPath::LinearPath2::ArcPath2::ArcPath2()
{
    DebugStop();//Nao eh para usar construtor vazio
}

AreaPath::LinearPath2::ArcPath2::ArcPath2(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius) :
ArcPath(cmesh,Origin,normalDirection,radius)
{
    
}

AreaPath::LinearPath2::ArcPath2::~ArcPath2()
{
    fcmesh = NULL;
}

TPZVec<REAL> AreaPath::LinearPath2::ArcPath2::Function(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> answ(3,0.);
    answ = ComputeFiniteDifference(t, xt, 0.001);
    
    return answ;
}

TPZVec<REAL> AreaPath::LinearPath2::ArcPath2::ComputeFiniteDifference(REAL t, TPZVec<REAL> xt, REAL halfDelta)
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

TPZVec<REAL> AreaPath::LinearPath2::ArcPath2::FunctionAux(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & direction)
{
//    double x = xt[0];
//    double y = xt[1];
//    double z = xt[2];
//    TPZVec<REAL> answTeste(3,0.);
//    answTeste[0] = -sin(2.*x*y*z);
//    return answTeste; //Para Z = -5 e r = 0.6, a integral final (multiplicada por 2) = 0.0669331
    
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
    
    TPZFMatrix<REAL> GradUt_Sigma(3,3,0.);
    GradUtxy.Multiply(Sigma, GradUt_Sigma);
    
    REAL W = 0.;
    for(int r = 0; r < 3; r++)
    {
        for(int c = 0; c < 3; c++)
        {
            W += 0.5*Sigma(r,c)*strain(r,c);
        }
    }
    
    TPZFMatrix<REAL> W_I(3,3,0.);
    for(int d = 0; d < 3; d++)
    {
        W_I(d,d) = W;
    }
    
    TPZFMatrix<REAL> W_I_minus_GradUt_Sigma(3,3,0.);
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

AreaPath::AreaPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


AreaPath::AreaPath(LinearPath * givenLinearPath)
{
    fLinearPath = new LinearPath2(givenLinearPath);
    
    fDETdxdt = fLinearPath->DETdxdt();//O ArcPath2 jah incluirah seu Detdxdt(), soh faltarah do LinearPath2::Detdxdt()
    
    #ifdef DEBUG
    if(!fLinearPath)
    {
        DebugStop();
    }
    #endif
}


AreaPath::~AreaPath()
{
    fLinearPath = NULL;
}


REAL AreaPath::DETdxdt()
{
    return fDETdxdt;
}


TPZVec<REAL> AreaPath::Func(REAL t)
{
    TPZVec<REAL> areaContribution(3,0.);
    areaContribution = fLinearPath->Func(t);
    
    return areaContribution;
}


//--------------------------------------------------------class JPath


Path::Path()
{
    fLinearPath = NULL;
    fArcPath = NULL;
}


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure)
{
    fLinearPath = new LinearPath(cmesh,Origin,normalDirection,radius,pressure);
    fArcPath = new ArcPath(cmesh,Origin,normalDirection,radius);
    fAreaPath = new AreaPath(fLinearPath);
    
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
    
    
    Adapt intRule(gIntegrPrecision);
    
    TPZTimer linInt("LinearIntegration"); linInt.start();
    TPZVec<REAL> linJintegral(3,0.);
    linJintegral = intRule.Vintegrate(*(jpathElem->GetLinearPath()),3,-1.,+1.);

    //Simetry in xz plane
    linJintegral[0] = 2. * linJintegral[0] * jpathElem->GetLinearPath()->DETdxdt() / gScaleFactor;
    linJintegral[1] = 0.;
    linJintegral[2] = 2. * linJintegral[2] * jpathElem->GetLinearPath()->DETdxdt() / gScaleFactor;

    linInt.stop();
    
    //4debug Jlin_x
    {
        #ifdef print_Jx_Linear
        std::map<REAL,REAL>::iterator ii;
        std::ofstream plotLinJx("PlotLinJx.txt");
        plotLinJx.precision(10);
        plotLinJx << "Jx = {";
        for(ii = functionLinJx.begin(); ii != functionLinJx.end(); ii++)
        {
            plotLinJx << "{" << ii->first << "," << ii->second << "},";
        }
        plotLinJx << "};\n";
        plotLinJx << "Jxplot = ListPlot[Jx,AxisOrigin->{0,0}]\n";
        plotLinJx.close();
        functionLinJx.clear();
        #endif
    }
    
    TPZTimer arcInt("ArcIntegration"); arcInt.start();
    TPZVec<REAL> arcJintegral(3,0.);
    arcJintegral = intRule.Vintegrate(*(jpathElem->GetArcPath()),3,-1.,+1.);

    //Simetry in xz plane
    arcJintegral[0] = 2. * arcJintegral[0] * jpathElem->GetArcPath()->DETdxdt() / gScaleFactor;
    arcJintegral[1] = 0.;
    arcJintegral[2] = 2. * arcJintegral[2] * jpathElem->GetArcPath()->DETdxdt() / gScaleFactor;

    arcInt.stop();
    //
    TPZTimer areaInt("AreaIntegration"); areaInt.start();
    TPZVec<REAL> areaJIntegral(3,0.);
    intRule.SetPrecision(1.e-3);
    areaJIntegral = intRule.Vintegrate(*(jpathElem->GetAreaPath()),3,-1.,+1.);

    //Simetry in xz plane
    areaJIntegral[0] = 2. * areaJIntegral[0] * jpathElem->GetAreaPath()->DETdxdt() / (gScaleFactor * gScaleFactor);
    areaJIntegral[1] = 0.;
    areaJIntegral[2] = 2. * areaJIntegral[2] * jpathElem->GetAreaPath()->DETdxdt() / (gScaleFactor * gScaleFactor);

    areaInt.stop();
    
    TPZVec<REAL> answ(3);
    for(int i = 0; i < 3; i++)
    {
        answ[i] = linJintegral[i] + arcJintegral[i] + areaJIntegral[i];
    }
    
    std::cout << "DeltaT integracao linha = " << linInt.seconds() << " s" << std::endl;
    std::cout << "DeltaT integracao arco = " << arcInt.seconds() << " s" << std::endl;
    std::cout << "DeltaT integracao area = " << areaInt.seconds() << " s" << std::endl;
    std::cout << "Jx AREA = " << areaJIntegral[0] << "\n";
    std::cout << "J = " << answ[0] << "\n";
    
    return answ;
}

/*
 Franc2D : somente pressão de 5 no interior da fratura resulta em J = 0.0027
 */





