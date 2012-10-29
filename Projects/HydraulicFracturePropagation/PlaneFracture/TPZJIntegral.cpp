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


LinearPath::LinearPath(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim)
{
    fOrigin = Origin;
    fNormalDirection = normalDirection;
    fradius = radius;
    
    fDETdxdt = fradius/2.;
    
    fcmesh = cmesh;
    fMeshDim = meshDim;
    fInitialElementId = 0;
    fcrackPressure = fabs(pressure);
    
    f_t_elIdqsi.clear();
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
    n[0] =  0.;
    n[1] = -1.;
    n[2] =  0.;
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
    linContribution = BoundaryFunc(t, xt, nt);
    
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
    
    return linContribution;
}


TPZVec<REAL> LinearPath::BoundaryFunc(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> qsi(0);
    
    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi.lower_bound(t);
    if(it != f_t_elIdqsi.end())
    {
        fInitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmesh->Reference()->ElementVec()[fInitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, fInitialElementId, fMeshDim);

    if(!geoEl)
    {
        std::cout << "geoEl not found! See " << __PRETTY_FUNCTION__ << " !!!\n";
        DebugStop();
    }
    
    f_t_elIdqsi[t] = std::make_pair(geoEl->Id(), qsi);
    
    TPZVec<REAL> minusGradUt_Sigma__n(fMeshDim,0.);
    if(geoEl->MaterialId() == __3DrockMat_quarterPoint)
    {
        //For a while will not compute solution in quarter point elements
        //return minusGradUt_Sigma__n;
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
    fInitialElementId = 0;
    
    f_t_elIdqsi.clear();
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
    arcContribution = BoundaryFunc(t, xt, nt);
    
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


TPZVec<REAL> ArcPath::BoundaryFunc(REAL t, TPZVec<REAL> & xt, TPZVec<REAL> & nt)
{
    TPZVec<REAL> qsi(0);

    std::map< REAL , std::pair< int , TPZVec<REAL> > >::iterator it = f_t_elIdqsi.lower_bound(t);
    if(it != f_t_elIdqsi.end())
    {
        fInitialElementId = it->second.first;
        qsi = it->second.second;
    }
    else
    {
        qsi.Resize(fcmesh->Reference()->ElementVec()[fInitialElementId]->Dimension(),0.);
    }
    TPZGeoEl * geoEl = fcmesh->Reference()->FindElement(xt, qsi, fInitialElementId, fMeshDim);

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
//    fInitialElementId = thePath->InitialElementId();
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


Path::Path(TPZAutoPointer<TPZCompMesh> cmesh, TPZVec<REAL> &Origin, TPZVec<REAL> &normalDirection, REAL radius, REAL pressure, int meshDim)
{
    fLinearPath = new LinearPath(cmesh,Origin,normalDirection,radius,pressure,meshDim);
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

    REAL precisionIntegralRule = 1.E-7;
    Adapt intRule(precisionIntegralRule);
    
    int meshDim = jpathElem->MeshDim();
    
//    {////////////4debug
//        TPZVec<REAL> JintegralTemp(meshDim,0.);
//        std::stringstream linear0;
//        std::stringstream linear1;
//        std::stringstream linear2;
//        std::stringstream arc0;
//        std::stringstream arc1;
//        std::stringstream arc2;
//        linear0 << "lin0 = {";
//        linear1 << "lin1 = {";
//        linear2 << "lin2 = {";
//        arc0 << "arc0 = {";
//        arc1 << "arc1 = {";
//        arc2 << "arc2 = {";
//
//        double nn = 50;
//        for(int i = int(-nn); i <= int(+nn); i++)
//        {
//            double tt = i/nn;
//            TPZVec<REAL> xt(3,0.);
//            
//            JintegralTemp = jpathElem->GetLinearPath()->Func(tt);
//            jpathElem->GetLinearPath()->X(tt, xt);
//            linear0 << "{" << xt[0] << "," << JintegralTemp[0] << "}";
//            linear1 << "{" << xt[0] << "," << JintegralTemp[1] << "}";
//            linear2 << "{" << xt[0] << "," << JintegralTemp[2] << "}";
//            
//            JintegralTemp = jpathElem->GetArcPath()->Func(tt);
//            jpathElem->GetArcPath()->X(tt, xt);
//            arc0 << "{" << xt[0] << "," << JintegralTemp[0] << "}";
//            arc1 << "{" << xt[0] << "," << JintegralTemp[1] << "}";
//            arc2 << "{" << xt[0] << "," << JintegralTemp[2] << "}";
//            
//            if(tt != 1.)
//            {
//                linear0 << ",";
//                linear1 << ",";
//                linear2 << ",";
//                arc0 << ",";
//                arc1 << ",";
//                arc2 << ",";
//            }
//        }
//        linear0 << "};\n";
//        linear1 << "};\n";
//        linear2 << "};\n";
//        arc0 << "};\n";
//        arc1 << "};\n";
//        arc2 << "};\n";
//        
//        std::ofstream saidaAll("funcs.txt");
//        saidaAll << linear0.str() << linear1.str() << linear2.str() << arc0.str() << arc1.str() << arc2.str() << "\n";
//        saidaAll << "l0=ListPlot[lin0,AxisOrigin->{0,0}]\n";
//        saidaAll << "l1=ListPlot[lin1,AxisOrigin->{0,0}]\n";
//        saidaAll << "l2=ListPlot[lin2,AxisOrigin->{0,0}]\n";
//        saidaAll << "a0=ListPlot[arc0,AxisOrigin->{0,0}]\n";
//        saidaAll << "a1=ListPlot[arc1,AxisOrigin->{0,0}]\n";
//        saidaAll << "a2=ListPlot[arc2,AxisOrigin->{0,0}]\n";
//    }
    
    TPZVec<REAL> linJintegral(meshDim,0.);
    linJintegral = intRule.Vintegrate(*(jpathElem->GetLinearPath()),meshDim,-1.,+1.);
    //
    TPZVec<REAL> arcJintegral(meshDim,0.);
    arcJintegral = intRule.Vintegrate(*(jpathElem->GetArcPath()),meshDim,-1.,+1.);
    
    TPZVec<REAL> answ(meshDim);
    for(int i = 0; i < meshDim; i++)
    {
        answ[i] = linJintegral[i] + arcJintegral[i];
    }
    
    return answ;
}

/*
 Franc2D : somente pressão de 5 no interior da fratura resulta em J = 0.0027
 */





