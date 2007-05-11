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

#ifndef PZNODEREPH
#define PZNODEREPH

#include "pzvec.h"

#include <map>
class TPZGeoMesh;

#include "pzlog.h"
#ifdef LOG4CXX
static LoggerPtr lognoderep(Logger::getLogger("pz.geom.tpznoderep"));
#endif

namespace pzgeom {

template<int N, class Topology>
class TPZNodeRep : public Topology

{
public:
  int fNodeIndexes[N];
  /**
   * Constructor with list of nodes
   */
  TPZNodeRep(TPZVec<int> &nodeindexes)
  {
    int nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
    memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int));
    int i;
    for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
  }

  /**
   * Empty constructor
   */
  TPZNodeRep()
  {
    int i;
    for(i=0; i<N; i++) fNodeIndexes[i]=-1;
  }

  /**
   * Constructor with node map
   */
  TPZNodeRep(const TPZNodeRep &cp,
                 std::map<int,int> & gl2lcNdMap);

  /**
   * Copy constructor
   */
  TPZNodeRep(const TPZNodeRep<N,Topology> &cp)
  {
    memcpy(fNodeIndexes,cp.fNodeIndexes,N*sizeof(int));
  }

  void Initialize(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh)
  {
    int nn = nodeindexes.NElements() < N ? nodeindexes.NElements() : N;
#ifndef NDEBUG
  if(nodeindexes.NElements() != N)
  {
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << " Nodeindexes have wrong size " << nodeindexes.NElements() << " but should be " << N;
#ifdef LOG4CXX
    LOGPZ_ERROR(lognoderep,sout.str().c_str());
#else
    std::cout << sout.str().c_str() << std::endl;
#endif
  }
#endif
    memcpy(fNodeIndexes,&nodeindexes[0],nn*sizeof(int));
    int i;
    for(i=nn; i<N; i++) fNodeIndexes[i]=-1;
  }

};

};

#endif
