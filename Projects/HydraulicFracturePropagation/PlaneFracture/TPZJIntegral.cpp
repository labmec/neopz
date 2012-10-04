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

const int dim3D = 3;



//--------------------------------------------------------class LinearPath


LinearPath::LinearPath()
{
    DebugStop();//Nao eh para usar construtor vazio
}


LinearPath::LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmesh = cmesh;
    fMeshDim = meshDim;
    fInitial2DElementId = 0;
}


void LinearPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    TPZVec<REAL> initialPoint(dim3D);
    initialPoint[0] = (fOrigin[0] - fradius*sin((Pi)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    initialPoint[1] = fradius*cos((Pi)/2.);
    initialPoint[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi)/2.));
    
    for(int c = 0; c < dim3D; c++)
    {
        xt[c] = (1.-t)/2.*initialPoint[c] + (1.+t)/2.*fOrigin[c];
    }
}


void LinearPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    dxdt.Resize(dim3D,0.);
    
    TPZVec<REAL> xarc(3);
    t = 1;
    xarc[0] = (fOrigin[0] - fradius*sin(Pi/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xarc[1] = fradius*cos((Pi*t)/2.);
    xarc[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
    
    for(int c = 0; c < dim3D; c++)
    {
        dxdt[c] = (fOrigin[c] - xarc[c])/2.;
    }
}


void LinearPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    n[0] = 0.;
    n[1] = -1.;
    n[2] = 0.;
}


REAL LinearPath::DETdxdt()
{
    return fDETdxdt;
}


TPZVec<REAL> LinearPath::Func(REAL t)
{
    TPZVec<REAL> xt(dim3D), nt(dim3D);
    
    this->X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> linContribution(fMeshDim);
    linContribution = BoundaryFunc(xt, nt, DETdxdt(), fMeshDim, fcmesh, fInitial2DElementId);
    
    return linContribution;
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
    fInitial2DElementId = 0;
}


void ArcPath::X(REAL t, TPZVec<REAL> & xt)
{
    xt.Resize(dim3D,0.);
    
    xt[0] = (fOrigin[0] - fradius*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
    xt[1] = fradius*cos((Pi*t)/2.);
    xt[2] = (fOrigin[2] + fradius*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
}


void ArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
{
    dxdt.Resize(dim3D,0.);
    
    dxdt[0] = -(Pi*fradius*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
    dxdt[1] = -(Pi*fradius*sin((Pi*t)/2.))/2.;
    dxdt[2] = +(Pi*fradius*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[0],fNormalDirection[2])))/2.;
}


void ArcPath::normalVec(REAL t, TPZVec<REAL> & n)
{
    TPZVec<REAL> xt(dim3D);
    X(t, xt);
    for(int i = 0; i < dim3D; i++)
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
    TPZVec<REAL> xt(dim3D), nt(dim3D);
    
    this->X(t,xt);
    this->normalVec(t, nt);
    
    TPZVec<REAL> arcContribution(fMeshDim);
    arcContribution = BoundaryFunc(xt, nt, DETdxdt(), fMeshDim, fcmesh, fInitial2DElementId);
    
    return arcContribution;
}


//--------------------------------------------------------class internalArcPath


//internalArcPath::internalArcPath() : Path()
//{
//    DebugStop();//Nao eh para usar construtor vazio
//}
//
//
//internalArcPath::internalArcPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL r_int, REAL r_ext, int meshDim) :
//                                                                                                    Path(cmesh, Origin, normalDirection, r_int, r_ext, meshDim)
//{
//    fDETdxdt = Pi*fr_int/2.;
//}
//
//
//internalArcPath::internalArcPath(Path * thePath)
//{
//    fOrigin = thePath->Origin();
//    fNormalDirection = thePath->NormalDirection();
//    fr_int = thePath->R_int();
//    fradius = thePath->R_ext();
//    fInitial2DElementId = thePath->Initial2DElementId();
//    fcmesh = thePath->Cmesh();
//    fMeshDim = thePath->MeshDim();
//
//    fDETdxdt = Pi*fr_int/2.;
//}
//
//void internalArcPath::X(REAL t, TPZVec<REAL> & xt)
//{
//    xt.Resize(dim3D,0.);
//
//    xt[0] = (fOrigin[0] + fr_int*sin((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])));
//    xt[1] = fr_int*cos((Pi*t)/2.);
//    xt[2] = (fOrigin[2] - fr_int*cos(atan2(fNormalDirection[2],fNormalDirection[0]))*sin((Pi*t)/2.));
//}
//
//
//void internalArcPath::dXdt(REAL t, TPZVec<REAL> & dxdt)
//{
//    dxdt.Resize(dim3D,0.);
//
//    dxdt[0] = +(Pi*fr_int*cos((Pi*t)/2.)*sin(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
//    dxdt[1] = -(Pi*fr_int*sin((Pi*t)/2.))/2.;
//    dxdt[2] = -(Pi*fr_int*cos((Pi*t)/2.)*cos(atan2(fNormalDirection[2],fNormalDirection[0])))/2.;
//}
//
//
//void internalArcPath::normalVec(REAL t, TPZVec<REAL> & n)
//{
//    TPZVec<REAL> xt(dim3D);
//    X(t, xt);
//    for(int i = 0; i < dim3D; i++)
//    {
//        n[i] = (fOrigin[i] - xt[i])/fr_int;
//    }
//}


//--------------------------------------------------------class JPath


Path::Path()
{
    fLinearPath = NULL;
    fArcPath = NULL;
    fMeshDim = 0;
}


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, int meshDim)
{
    fLinearPath = new LinearPath(cmesh,Origin,normalDirection,radius,meshDim);
    fArcPath = new ArcPath(cmesh,Origin,normalDirection,radius,meshDim);
    
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

    REAL precisionIntegralRule = 1.E-15;
    Adapt intRule(precisionIntegralRule);
    
    int meshDim = jpathElem->MeshDim();
    TPZVec<REAL> linJintegral(meshDim,0.);
    //linJintegral = intRule.Vintegrate(*(jpathElem->GetLinearPath()),meshDim,-1.,+1.);
    //
    TPZVec<REAL> arcJintegral(meshDim,0.);
    arcJintegral = intRule.Vintegrate(*(jpathElem->GetArcPath()),meshDim,-1.,+1.);
    
    TPZVec<REAL> answ(meshDim);
    if(meshDim == 2)
    {
        answ[0] = 2.*(linJintegral[0] + arcJintegral[0]);
        answ[1] = 0.;
    }
    else if(meshDim == 3)
    {
        //Pela simetria do problema em relacao ao plano xz, deve-se somar a este vetor seu espelho em relacao ao plano xz.
        answ[0] = 2.*(linJintegral[0] + arcJintegral[0]);
        answ[1] = 0.;
        answ[2] = 2.*(linJintegral[2] + arcJintegral[2]);
    }
    
    return answ;
}





