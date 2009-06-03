//
// C++ Interface: pznoderep
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "pznoderep.h"
#include "pzlog.h"
#include <sstream>

#include "tpzpoint.h"
#include "tpzline.h"
#include "tpztriangle.h"
#include "tpzquadrilateral.h"
#include "tpzpyramid.h"
#include "tpztetrahedron.h"
#include "tpzcube.h"
#include "tpzprism.h"

#ifdef LOG4CXX
static LoggerPtr loggernoderep(Logger::getLogger("pz.geom.noderep"));
#endif

//using namespace pztopology;

namespace pzgeom {

  /**
 * Constructor with node map
   */
template<int N, class Topology>
TPZNodeRep<N,Topology>::TPZNodeRep(const TPZNodeRep<N,Topology> &cp, std::map<int,int> & gl2lcNdMap)
{
  int i;
  for(i = 0; i < N; i++)
  {
    if (gl2lcNdMap.find(cp.fNodeIndexes[i]) == gl2lcNdMap.end())
    {
      std::stringstream sout;
      sout << "ERROR in - " << __PRETTY_FUNCTION__
          << " trying to clone a node " << i << " index " << cp.fNodeIndexes[i]
          << " wich is not mapped";
      LOGPZ_ERROR(loggernoderep,sout.str().c_str());
      fNodeIndexes[i] = -1;
      continue;
    }
    fNodeIndexes[i] = gl2lcNdMap [ cp.fNodeIndexes[i] ];
  }
}

template<int N, class Topology>
bool TPZNodeRep<N,Topology>::IsInSideParametricDomain(int side, TPZVec<REAL> &pt, REAL tol){
    MElementType type = Topology::Type(side);
    switch(type){
      case EPoint:
        return pztopology::TPZPoint::IsInParametricDomain(pt,tol);
        break;
      case EOned:
        return pztopology::TPZLine::IsInParametricDomain(pt,tol);
        break;
      case ETriangle:
        return pztopology::TPZTriangle::IsInParametricDomain(pt,tol);
        break;
      case EQuadrilateral:
        return pztopology::TPZQuadrilateral::IsInParametricDomain(pt,tol);
        break;
      case ETetraedro:
        return pztopology::TPZTetrahedron::IsInParametricDomain(pt,tol);
        break;
      case EPiramide:
        return pztopology::TPZPyramid::IsInParametricDomain(pt,tol);
        break;
      case EPrisma:
        return pztopology::TPZPrism::IsInParametricDomain(pt,tol);
        break;
      case ECube:
        return pztopology::TPZCube::IsInParametricDomain(pt,tol);
        break;
      default: 
        std::stringstream sout;
        sout << "Fatal error at " << __PRETTY_FUNCTION__ << " - Element type " << type << " not found";
        PZError << "\n" << sout.str() << "\n";
#ifdef LOG4CXX
        LOGPZ_FATAL(lognoderep,sout.str().c_str());
#endif
        DebugStop();
        return false;
        break;
    }///case
  
  return false;
}///method

}
