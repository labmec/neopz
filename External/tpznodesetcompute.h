//
// C++ Interface: tpznodesetcompute
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZNODESETCOMPUTE_H
#define TPZNODESETCOMPUTE_H
#include "pzvec.h"
#include "pzstack.h"
#include <set>

/**
This class will compute the cardinality of a nodegraph, identifying nodes as vertices, lines, faces or volumes
It will also compress the nodegraph for nodes with identical connectivity graph

@author Philippe R. B. Devloo
*/
class TPZNodesetCompute{
public:
    TPZNodesetCompute();

    ~TPZNodesetCompute();
    /**
    Group the node graph as passed by the parameters
    @param nodegraph [in] node connectivity graph (large array)
    @param nodegraphindex [in] pointer in the nodegraph (small array)
    */
    void AnalyseGraph();
    
    void BuildNodeGraph(TPZVec<int> &blockgraph, TPZVec<int> &blockgraphindex);
    
    void BuildVertexGraph(TPZStack<int> &blockgraph, TPZVec<int> &blockgraphindex);
    
    void BuildElementGraph(TPZStack<int> &blockgraph, TPZStack<int> &blockgraphindex);

void BuildNodeSet(int node, std::set<int> &nodeset);  

TPZVec<int> &Nodegraph() {
return fNodegraph;
}

TPZVec<int> &Nodegraphindex() {
return fNodegraphindex;
}

private:
  TPZVec<int> fNodegraph;
  TPZVec<int> fNodegraphindex;
  int fMaxSeqNum;
  int fMaxLevel;
  TPZVec<int> fSeqNumber;
  TPZStack<int> fSeqCard;
  TPZVec<int> fLevel;

/**
 * This method will analyse the set inclusion of the current node, calling the method
 * recursively if another node need to be analysed first
 */
void AnalyseNode(int node, TPZVec< std::set<int> > &nodeset);  

/**
 * Build the set of nodes which are vertices
 */
void BuildVertexSet(int node, std::set<int> &nodeset);

/**
 * Look for elements formed by vertices, intersecting with the intersectvertices, one by one
 * If the intersection does not remove any of the intersectvertices, we found an element!
 */
 void AnalyseForElements(std::set<int> &vertices, std::set< std::set<int> > &elements);
/**
working a set of vertex nodes with nodes which have to be intersected (tested)
first parameter the set of nodes which need to form elements
second parameter nodes whose intersection need to be considered
third parameter vertex sets which have already been considered
*/

void SubstractLowerNodes(int node, std::set<int> &nodeset);
};

#endif
