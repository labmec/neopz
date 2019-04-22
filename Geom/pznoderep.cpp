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
void TPZNodeRep<N,Topology>::GetSideShapeFunction(int side, TPZVec<REAL> &qsiSide, TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi ){
    MElementType  sideType = Topology::Type(side);
#ifdef PZDEBUG
    REAL tol = 1e-12;
    if(!IsInSideParametricDomain(side,qsiSide,tol)){
        PZError<<"The method expects the coordinates in the side's parametric domain. Exiting..."<<std::endl;
        DebugStop();
    }
#ifdef LOG4CXX
    LOGPZ_FATAL(lognoderep,sout.str().c_str());
#endif
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
        default:
            PZError<<"Could not find associated shape function to the side. Details are as follows:"<<std::endl;
            PZError<<"Element is of type "<<MElementType_Name(Topology::Type())<<std::endl;
            PZError<<"Side\t"<<side<<" is of type\t"<< MElementType_Name(sideType)<<std::endl;
        #ifdef LOG4CXX
            LOGPZ_FATAL(lognoderep,sout.str().c_str());
        #endif
            DebugStop();
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
}
