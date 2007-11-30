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

}
