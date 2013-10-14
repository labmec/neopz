/**
 * @file
 * @brief Contains the TPZNodesetCompute class which computes the cardinality of a nodegraph.
 */
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

template<class TVar>
class TPZBlock;

/**
 * @brief Computes the cardinality of a nodegraph, identifying nodes as vertices, lines, faces or volumes. \ref util "Utility"
 * @author Philippe R. B. Devloo
 */
 /**
 * It will also compress the nodegraph for nodes with identical connectivity graph
 */
class TPZNodesetCompute {
public:
    TPZNodesetCompute();
	
    ~TPZNodesetCompute();
    /**
	 * @brief Group the node graph as passed by the parameters
	 */
    void AnalyseGraph();
    
    /** @brief Build the graph which groups the equations of each node */
    void BuildNodeGraph(TPZVec<long> &blockgraph, TPZVec<long> &blockgraphindex);
    
    /** @brief build the graph which builds the equations linked to vertices */
    void BuildVertexGraph(TPZStack<long> &blockgraph, TPZVec<long> &blockgraphindex);
    
    /** @brief Build the  graph which groups the equations grouped by elements */
    void BuildElementGraph(TPZStack<long> &blockgraph, TPZStack<long> &blockgraphindex);
	
	void BuildNodeSet(long node, std::set<long> &nodeset);
	
	/** @brief Expand the graph acording to the block structure */
	static void ExpandGraph(TPZVec<long> &graph, TPZVec<long> &graphindex, TPZBlock<STATE> &block,
							TPZVec<long> &expgraph, TPZVec<long> &expgraphindex);
	/** @brief Color the graph into mutually independent blocks */
	static int ColorGraph(TPZVec<long> &graph, TPZVec<long> &graphindex, long neq,
						  TPZVec<int> &colors);
	
	/** @brief Returns the level of the nodes */
	TPZVec<int> &Levels()
	{
		return fLevel;
	}
	
	/** @brief Returns the maximum level */
	int MaxLevel()
	{
		return fMaxLevel;
	}
    
	void Print(std::ostream &file) const;
	
	static void Print(std::ostream &file, const TPZVec<long> &graphindex, const TPZVec<long> &graph);
	
	static void Print(std::ostream &file, const std::set<long> &nodeset, const char *text);
	
	TPZVec<long> &Nodegraph() {
		return fNodegraph;
	}
	
	TPZVec<long> &Nodegraphindex() {
		return fNodegraphindex;
	}
	
	TPZVec<int> &IsIncluded()
	{
		return fIsIncluded;
	}
	
private:
	/**
	 * @brief The node graph as passed on by the finite element mesh \n
	 * His node graph is organized by sequence numbers
	 */
	TPZManVector<long> fNodegraph;
	TPZVec<long> fNodegraphindex;
	/** @brief Counter for the condensed node graph */
	long fMaxSeqNum;
	int fMaxLevel;
	
	/** @brief Sequence number associated with each node after condensing */
	TPZVec<long> fSeqNumber;
	/** @brief Number of nodes associated with each sequence number */
	TPZStack<long> fSeqCard;
	/** @brief Inclusion relation ship between nodes */
	TPZVec<int> fLevel;
	/** @brief Vector indicating whether a node connectivity is included in another one */
	TPZVec<int> fIsIncluded;
	
	/**
	 * @brief This method will analyse the set inclusion of the current node, calling the method \n
	 * recursively if another node need to be analysed first
	 */
	void AnalyseNode(long node, TPZVec< std::set<long> > &nodeset);  
	
	/** @brief Look for elements formed by vertices, intersecting with the intersectvertices, one by one */
	/** If the intersection does not remove any of the intersectvertices, we found an element! */
	void AnalyseForElements(std::set<long> &vertices, std::set< std::set<long> > &elements);
	/**
	 * @brief working a set of vertex nodes with nodes which have to be intersected (tested)
	 * @param node node whose intersection need to be considered
	 * @param nodeset the set of nodes which need to form elements
	 */
	void SubstractLowerNodes(long node, std::set<long> &nodeset);
};

#endif
