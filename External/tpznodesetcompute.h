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
    void BuildNodeGraph(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex);
    
    /** @brief build the graph which builds the equations linked to vertices */
    void BuildVertexGraph(TPZStack<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex);
    
    /** @brief Build the  graph which groups the equations grouped by elements */
    void BuildElementGraph(TPZStack<int64_t> &blockgraph, TPZStack<int64_t> &blockgraphindex);
	
	void BuildNodeSet(int64_t node, std::set<int64_t> &nodeset);
	
	/** @brief Expand the graph acording to the block structure */
	static void ExpandGraph(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, TPZBlock &block,
							TPZVec<int64_t> &expgraph, TPZVec<int64_t> &expgraphindex);
	/** @brief Color the graph into mutually independent blocks */
	static int ColorGraph(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, int64_t neq,
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
	
	static void Print(std::ostream &file, const TPZVec<int64_t> &graphindex, const TPZVec<int64_t> &graph);
	
	static void Print(std::ostream &file, const std::set<int64_t> &nodeset, const char *text);
	
	TPZManVector<int64_t> &Nodegraph() {
		return fNodegraph;
	}
	
	TPZManVector<int64_t> &Nodegraphindex() {
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
	TPZManVector<int64_t> fNodegraph;
	TPZManVector<int64_t> fNodegraphindex;
	/** @brief Counter for the condensed node graph */
	int64_t fMaxSeqNum;
	int fMaxLevel;
	
	/** @brief Sequence number associated with each node after condensing */
	TPZVec<int64_t> fSeqNumber;
	/** @brief Number of nodes associated with each sequence number */
	TPZStack<int64_t> fSeqCard;
	/** @brief Inclusion relation ship between nodes */
	TPZVec<int> fLevel;
	/** @brief Vector indicating whether a node connectivity is included in another one */
	TPZVec<int> fIsIncluded;
	
	/**
	 * @brief This method will analyse the set inclusion of the current node, calling the method \n
	 * recursively if another node need to be analysed first
	 */
	void AnalyseNode(int64_t node, TPZVec< std::set<int64_t> > &nodeset);  
	
	/** @brief Look for elements formed by vertices, intersecting with the intersectvertices, one by one */
	/** If the intersection does not remove any of the intersectvertices, we found an element! */
	void AnalyseForElements(std::set<int64_t> &vertices, std::set< std::set<int64_t> > &elements);
	/**
	 * @brief working a set of vertex nodes with nodes which have to be intersected (tested)
	 * @param node node whose intersection need to be considered
	 * @param nodeset the set of nodes which need to form elements
	 */
	void SubtractLowerNodes(int64_t node, std::set<int64_t> &nodeset);
};

#endif
