//
// C++ Implementation: pznoderep
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
static LoggerPtr logger(Logger::getLogger("pz.geom.noderep"));
#endif

namespace pzgeom {

  /**
 * Constructor with node map
   */
template<int N>
TPZNodeRep<N>::TPZNodeRep(const TPZNodeRep<N> &cp,
                 std::map<int,int> & gl2lcNdMap)
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
      LOGPZ_ERROR(logger,sout.str().c_str());
      fNodeIndexes[i] = -1;
      continue;
    }
    fNodeIndexes[i] = gl2lcNdMap [ cp.fNodeIndexes[i] ];
  }
}

template class TPZNodeRep<1>;
template class TPZNodeRep<2>;
template class TPZNodeRep<3>;
template class TPZNodeRep<4>;
template class TPZNodeRep<5>;
template class TPZNodeRep<6>;
template class TPZNodeRep<8>;

};
