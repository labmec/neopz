/**
 * @file
 * @brief Creates TPZNodeRep classes for several master elements. 
 */

#ifndef BORLAND
#include "pznoderep.h.h"
#endif

#ifdef BORLAND
#include "pznoderep.h"
#endif
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "pzgeoprism.h"

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpzpyramid.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

using namespace pztopology;
namespace pzgeom{

template<int N, class Topology>
bool TPZNodeRep<N,Topology>::IsLinearMapping() const
{
    return true;
}

template<int N, class Topology>
template<class T>
void TPZNodeRep<N,Topology>::GetSideShapeFunction(int side, TPZVec<T> &qsiSide, TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi ){
    MElementType  sideType = Topology::Type(side);
#ifdef PZDEBUG
    std::stringstream sout;
    REAL tol = 1e-12;
    if(!IsInSideParametricDomain(side,qsiSide,tol)){
        sout<<"The method expects the coordinates in the side's parametric domain. Exiting..."<<std::endl;
        DebugStop();
        PZError << "\n" << sout.str() << "\n";
        #ifdef LOG4CXX
        LOGPZ_FATAL(lognoderep,sout.str().c_str());
        #endif
    }
#endif
    switch (sideType){
        case EPoint:
            return TPZGeoPoint::TShape(qsiSide,phi,dphi);
        case EOned:
            return TPZGeoLinear::TShape(qsiSide,phi,dphi);
        case ETriangle:
            return TPZGeoTriangle::TShape(qsiSide,phi,dphi);
        case EQuadrilateral:
            return TPZGeoQuad::TShape(qsiSide,phi,dphi);
        case ETetraedro:
            return TPZGeoTetrahedra::TShape(qsiSide,phi,dphi);
        case EPiramide:
            return TPZGeoPyramid::TShape(qsiSide,phi,dphi);
        case EPrisma:
            return TPZGeoPrism::TShape(qsiSide,phi,dphi);
        case ECube:
            return TPZGeoCube::TShape(qsiSide,phi,dphi);
        default:{
            std::stringstream sout;
            sout<<"Could not find associated shape function to the side. Details are as follows:"<<std::endl;
            sout<<"Element is of type "<<MElementType_Name(Topology::Type())<<std::endl;
            sout<<"Side\t"<<side<<" is of type\t"<< MElementType_Name(sideType)<<std::endl;
            PZError<<std::endl<<sout.str()<<std::endl;
            #ifdef LOG4CXX
            LOGPZ_FATAL(lognoderep,sout.str().c_str());
            #endif
            DebugStop();
        }
    }
}

    template<int N, class Topology>
    template<class T>
    void TPZNodeRep<N,Topology>::GetPointInSideInfluence(const int &origSide, const int &subSide, const TPZVec<T> &qsiSide, T &correctionFactor){
        MElementType  sideType = Topology::Type(origSide);
#ifdef PZDEBUG
        std::stringstream sout;
        REAL tol = 1e-12;
        if(!IsInSideParametricDomain(origSide,qsiSide,tol)){
            sout<<"The method expects the coordinates in the side's parametric domain. Exiting..."<<std::endl;
            DebugStop();
            PZError << "\n" << sout.str() << "\n";
            #ifdef LOG4CXX
            LOGPZ_FATAL(lognoderep,sout.str().c_str());
            #endif
        }
#endif
        switch (sideType){
            case EPoint:
                return TPZGeoPoint::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case EOned:
                return TPZGeoLinear::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case ETriangle:
                return TPZGeoTriangle::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case EQuadrilateral:
                return TPZGeoQuad::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case ETetraedro:
                return TPZGeoTetrahedra::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case EPiramide:
                return TPZGeoPyramid::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case EPrisma:
                return TPZGeoPrism::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            case ECube:
                return TPZGeoCube::CalcSideInfluence(subSide,qsiSide,correctionFactor);
            default:{
//                std::stringstream sout;
//                sout<<"Could not find associated shape function to the side. Details are as follows:"<<std::endl;
//                sout<<"Element is of type "<<MElementType_Name(Topology::Type())<<std::endl;
//                sout<<"Side\t"<<side<<" is of type\t"<< MElementType_Name(sideType)<<std::endl;
//                PZError<<std::endl<<sout.str()<<std::endl;
//                #ifdef LOG4CXX
//                LOGPZ_FATAL(lognoderep,sout.str().c_str());
//                #endif
                DebugStop();
            }
        }
    }

template class TPZNodeRep<1,TPZPoint>;
template class TPZNodeRep<2,TPZLine>;
template class TPZNodeRep<3,TPZTriangle>;
template class TPZNodeRep<4,TPZQuadrilateral>;
template class TPZNodeRep<5,TPZPyramid>;
template class TPZNodeRep<4,TPZTetrahedron>;
template class TPZNodeRep<6,TPZPrism>;
template class TPZNodeRep<8,TPZCube>;

template void TPZNodeRep<1,TPZPoint>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<2,TPZLine>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<3,TPZTriangle>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<4,TPZQuadrilateral>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<5,TPZPyramid>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<4,TPZTetrahedron>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<6,TPZPrism>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);
template void TPZNodeRep<8,TPZCube>::GetSideShapeFunction<REAL>(int , TPZVec<REAL> &, TPZFMatrix<REAL> &,TPZFMatrix<REAL> &);

}
#ifdef _AUTODIFF
template<class T=REAL>
class Fad;
template void pzgeom::TPZNodeRep<1,TPZPoint>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<2,TPZLine>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<3,TPZTriangle>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<4,TPZQuadrilateral>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<5,TPZPyramid>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<4,TPZTetrahedron>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<6,TPZPrism>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
template void pzgeom::TPZNodeRep<8,TPZCube>::GetSideShapeFunction<Fad<REAL>>(int , TPZVec<Fad<REAL>> &, TPZFMatrix<Fad<REAL>> &,TPZFMatrix<Fad<REAL>> &);
#endif

