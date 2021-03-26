//
//  pzgeom_utility.h
//  AcoplamentoH1Hdiv
//
//  Created by Philippe Devloo on 31/08/19.
//

#ifndef pzgeom_utility_h
#define pzgeom_utility_h

#include <stdio.h>
class TPZGeoMesh;
class TPZGeoElSide;
#include "pzgeoel.h"
#include "pzgnode.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpzpyramid.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

#include "TPZGeoLinear.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzgeopoint.h"
#include "pzgeoprism.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"


/** @brief Verifies if pt (in parametric domain of the side) is within boundaries */
template<class Topology>
bool IsInSideParametricDomain(int side, const TPZVec<REAL> &pt, REAL tol);

//template<typename T,
//typename std::enable_if<std::is_same<T,Fad<REAL>>::value>::type* = nullptr>
//inline bool IsInSideParametricDomain(int side, const TPZVec<T> &pt, REAL tol){
//    TPZVec<REAL> qsiReal(pt.size(),-1);
//    for(int i = 0; i < qsiReal.size(); i++) qsiReal[i] = pt[i].val();
//    return IsInSideParametricDomain(side,qsiReal,tol);
//}


//template<class Topology, class T>
//void GetSideShapeFunction(int side, TPZVec<T> &qsiSide, TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi );


/**
 * @brief Projects point pt (in parametric coordinate system) in the element parametric domain.
 * @return Returns the side where the point was projected.
 * @note Observe that if the point is already in the parametric domain, the method will return
 * \f$ NSides() - 1 \f$
 */
template<class Topology>
int ProjectInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain){
    const int nsides = Topology::NSides;
    if(Topology::IsInParametricDomain(qsi,0.)){///it is already in the domain
        qsiInDomain = qsi;
        return nsides-1;
    }//if
    
    int winnerSide = -1;
    REAL winnerDistance = 1e12;
    TPZManVector<REAL,3> pt1(qsi.NElements()), pt2(qsi.NElements());
    for(int is = 0; is < nsides-1; is++){
        
        ///Go from NSides-1 to side is
        TPZTransform<> T1 = Topology::SideToSideTransform(nsides-1, is);
        pt1.Resize(Topology::SideDimension(is));
        T1.Apply(qsi,pt1);
        
        ///Check if the point is within side boundaries
        bool IsInSideDomain = IsInSideParametricDomain<Topology>(is,pt1,0.);
        if(!IsInSideDomain) continue;
        
        ///Come back from side is to \f$ NSides-1 \f$
        TPZTransform<> T2 = Topology::SideToSideTransform(is,nsides-1);
        T2.Apply(pt1,pt2);
        
        ///Compare original to mapped point
        REAL distance = 0.;
        for(int i = 0; i < qsi.NElements(); i++){
            REAL val = qsi[i]-pt2[i];
            distance += val*val;
        }//i
        distance = sqrt(distance);
        
        ///The closest side point to the original is the projected point
        if(distance < winnerDistance){
            winnerDistance = distance;
            winnerSide = is;
            qsiInDomain = pt2;
        }
    }//for is
    
    return winnerSide;
    
}//method

// project the point towards the center of the element
// find the intersecting side
template<class Topology>
int ProjectBissectionInParametricDomain(TPZVec<REAL> &qsi, TPZVec<REAL> &qsiInDomain){
    const int nsides = Topology::NSides;
    const int dim = Topology::Dimension;
    
    qsiInDomain.Resize(dim);
    
    if(Topology::IsInParametricDomain(qsi,0.))///it is already in the domain
    {
        qsiInDomain = qsi;
        return nsides-1;
    }
    
    REAL tol = 1.e-10;
    
    ///first, will be made a project to center direction
    TPZVec<REAL> OutPt(qsi), InnPt(dim);
    Topology::CenterPoint(nsides-1,InnPt);
    
    REAL dist = 0.;
    for(int c = 0; c < dim; c++)
    {
        dist += (InnPt[c] - OutPt[c])*(InnPt[c] - OutPt[c]);
    }
    dist = sqrt(dist);
    
    while(dist > tol)
    {
        for(int c = 0; c < dim; c++)
        {
            qsiInDomain[c] = (InnPt[c] + OutPt[c])/2.;
        }
        if(Topology::IsInParametricDomain(qsiInDomain,0.))
        {
            InnPt = qsiInDomain;
        }
        else
        {
            OutPt = qsiInDomain;
        }
        dist = 0.;
        for(int c = 0; c < dim; c++)
        {
            dist += (InnPt[c] - OutPt[c])*(InnPt[c] - OutPt[c]);
        }
        dist = sqrt(dist);
    }
    
    ///found in witch side the projection belongs
    int winnerSide = -1;
    TPZManVector<REAL,3> pt1(dim), pt2(dim);
    for(int is = 0; is < nsides-1; is++)
    {
        ///Go orthogonally from \f$ NSides-1 \f$ to side is
        TPZTransform<> T1 = Topology::SideToSideTransform(nsides-1, is);
        T1.Apply(qsiInDomain,pt1);
        
        ///Come back from side is to \f$ NSides-1 \f$
        TPZTransform<> T2 = Topology::SideToSideTransform(is,nsides-1);
        T2.Apply(pt1,pt2);
        
        ///Compare ptInDomain to transformed point
        dist = 0.;
        for(int c = 0; c < dim; c++)
        {
            dist += (qsiInDomain[c]-pt2[c]) * (qsiInDomain[c]-pt2[c]);
        }//i
        dist = sqrt(dist);
        
        ///Closest side
        if(dist < tol)
        {
            winnerSide = is;
            qsiInDomain = pt2;
            break;
        }
    }//for is
    
    return winnerSide;
    
}//method


template<class Topology>
inline void GetSideShapeFunction(int side, TPZVec<Fad<REAL> > &qsiSide, TPZFMatrix<Fad<REAL> > &phi,TPZFMatrix<Fad<REAL> > &dphi ){
    MElementType  sideType = Topology::Type(side);
#ifdef PZDEBUG
    std::stringstream sout;
    REAL tol = 1e-12;
    TPZManVector<REAL,3> qsi(qsiSide.size());
    for (int i=0; i<qsi.size(); i++) {
        qsi[i] = qsiSide[i].val();
    }
    if(!IsInSideParametricDomain<Topology>(side,qsi,tol)){
        sout<<"The method expects the coordinates in the side's parametric domain. Exiting..."<<std::endl;
        PZError << "\n" << sout.str() << "\n";
        DebugStop();
    }
#endif
    switch (sideType){
        case EPoint:
            return pzgeom::TPZGeoPoint::TShape(qsiSide,phi,dphi);
        case EOned:
            return pzgeom::TPZGeoLinear::TShape(qsiSide,phi,dphi);
        case ETriangle:
            return pzgeom::TPZGeoTriangle::TShape(qsiSide,phi,dphi);
        case EQuadrilateral:
            return pzgeom::TPZGeoQuad::TShape(qsiSide,phi,dphi);
        case ETetraedro:
            return pzgeom::TPZGeoTetrahedra::TShape(qsiSide,phi,dphi);
        case EPiramide:
            return pzgeom::TPZGeoPyramid::TShape(qsiSide,phi,dphi);
        case EPrisma:
            return pzgeom::TPZGeoPrism::TShape(qsiSide,phi,dphi);
        case ECube:
            return pzgeom::TPZGeoCube::TShape(qsiSide,phi,dphi);
        default:{
            std::stringstream sout;
            sout<<"Could not find associated shape function to the side. Details are as follows:"<<std::endl;
            sout<<"Element is of type "<<MElementType_Name(Topology::Type())<<std::endl;
            sout<<"Side\t"<<side<<" is of type\t"<< MElementType_Name(sideType)<<std::endl;
            PZError<<std::endl<<sout.str()<<std::endl;
            DebugStop();
        }
    }
}

template<class Topology>
inline void GetSideShapeFunction(int side, TPZVec<REAL> &qsiSide, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi ){
    MElementType  sideType = Topology::Type(side);
#ifdef PZDEBUG
    std::stringstream sout;
    REAL tol = pztopology::GetTolerance();
    if(!IsInSideParametricDomain<Topology>(side,qsiSide,tol)){
        sout<<"The method expects the coordinates in the side's parametric domain. Exiting..."<<std::endl;
        DebugStop();
    }
#endif
    switch (sideType){
        case EPoint:
            return pzgeom::TPZGeoPoint::TShape(qsiSide,phi,dphi);
        case EOned:
            return pzgeom::TPZGeoLinear::TShape(qsiSide,phi,dphi);
        case ETriangle:
            return pzgeom::TPZGeoTriangle::TShape(qsiSide,phi,dphi);
        case EQuadrilateral:
            return pzgeom::TPZGeoQuad::TShape(qsiSide,phi,dphi);
        case ETetraedro:
            return pzgeom::TPZGeoTetrahedra::TShape(qsiSide,phi,dphi);
        case EPiramide:
            return pzgeom::TPZGeoPyramid::TShape(qsiSide,phi,dphi);
        case EPrisma:
            return pzgeom::TPZGeoPrism::TShape(qsiSide,phi,dphi);
        case ECube:
            return pzgeom::TPZGeoCube::TShape(qsiSide,phi,dphi);
        default:{
            std::stringstream sout;
            sout<<"Could not find associated shape function to the side. Details are as follows:"<<std::endl;
            sout<<"Element is of type "<<MElementType_Name(Topology::Type())<<std::endl;
            sout<<"Side\t"<<side<<" is of type\t"<< MElementType_Name(sideType)<<std::endl;
            PZError<<std::endl<<sout.str()<<std::endl;
            DebugStop();
        }
    }
}

template<class Topology>
inline bool IsInSideParametricDomain(int side, const TPZVec<REAL> &pt, REAL tol){
    MElementType type = Topology::Type(side);
    switch(type){
        case EPoint:
            return pztopology::TPZPoint::IsInParametricDomain(pt,tol);
        case EOned:
            return pztopology::TPZLine::IsInParametricDomain(pt,tol);
        case ETriangle:
            return pztopology::TPZTriangle::IsInParametricDomain(pt,tol);
        case EQuadrilateral:
            return pztopology::TPZQuadrilateral::IsInParametricDomain(pt,tol);
            
        case ETetraedro:
            return pztopology::TPZTetrahedron::IsInParametricDomain(pt,tol);
            
        case EPiramide:
            return pztopology::TPZPyramid::IsInParametricDomain(pt,tol);
            
        case EPrisma:
            return pztopology::TPZPrism::IsInParametricDomain(pt,tol);
            
        case ECube:
            return pztopology::TPZCube::IsInParametricDomain(pt,tol);
            
        default:
            std::stringstream sout;
            sout << "Fatal error at " << __PRETTY_FUNCTION__ << " - Element type " << type << " not found";
            PZError << "\n" << sout.str() << "\n";
            DebugStop();
            return false;
            break;
    }//case
    
    return false;
}//method


#endif /* pzgeom_utility_h */
