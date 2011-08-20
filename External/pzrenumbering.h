/**
 * @file
 * @brief Contains the TPZRenumbering class which defines the behaviour to implementing node sequence numbering optimization.
 */

#ifndef TPZRENUMBERING_H
#define TPZRENUMBERING_H

#include "pzvec.h"
#include <set>
/** 
 * @brief This abstract class which defines the behaviour which derived classes need to implement \n
 * for implementing node sequence numbering optimization. \ref util "Utility"
 * @ingroup util
 */
class TPZRenumbering {
public:
	
	int fHDivPermute;
	
	TPZRenumbering() : fNElements(0), fNNodes(0)
	{
	}
	
	TPZRenumbering(int NElements, int NNodes);
	
	virtual ~TPZRenumbering()
	{
	}
	
	void SetElementsNodes(int NElements, int NNodes)
	{
		this->fNElements = NElements;
		this->fNNodes = NNodes;
	}
	/**
	 * @brief This method declares the element graph to the object
	 */
	/**
	 * The first vector contains the element node number \n
	 * The second vector contains the index where to find the first node number
	 * of each element \n
	 * The size of second vector is fNElements+1
	 */
	void SetElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex);
	
	/** @brief Sets the number of equations associated with each node */
	/** The derived class may or may not take this data into consideration */
	void SetNodeWeights(TPZVec<int> &weights)
	{
		fNodeWeights = weights;
	}
	
	/** @brief This will reset all datastructures the object may contain. */
	/** Node resequencing algorithms may require a possibly large amount of temporary data */
	virtual void ClearDataStructures();
	
	virtual void Resequence(TPZVec<int> &perm, TPZVec<int> &iperm)
	{
		/*if(fHDivPermute)
		 {
		 AdjustHDivPermutation(perm, iperm);
		 }*/
		std::cout << "Resequence not implemented\n";
		DebugStop();
	}
	
	
	/**
	 * @brief Will convert an element graph defined by elgraph and elgraphindex
	 * into a node graph defined by nodegraph and nodegraphindex
	 */
	void ConvertGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex, TPZVec<int> &nodegraph, TPZVec<int> &nodegraphindex);
	
	/** @brief Convert a traditional elgraph to an element to element graph */
	void ConvertToElementoToElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex, TPZVec<int> &eltotelgraph, TPZVec<int> &eltoelweight, TPZVec<int> &eltoelgraphindex);

	/** @brief Stores the graph of nodes to elements */
	void NodeToElGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex, TPZVec<int> &nodetoelgraph, TPZVec<int> &nodetoelgraphindex);
	
	/**
	 * @brief Will assign a color to the nodes in the graph such that no two connected nodes have the same color
	 * the return value indicates the number of colors in the graph
	 */
	int ColorNodes(TPZVec<int> &nodegraph, TPZVec<int> &nodegraphindex, TPZVec<int> &family, TPZVec<int> &colors);
	
	/** @brief Prints graph */
	void Print(TPZVec<int> &grapho, TPZVec<int> &graphoindex, const char *name = 0, std::ostream &out = std::cout);
	
	/**
	 * @brief Analyse the graph, find the corner nodes \n
	 * Number of elements which should be considered for determining corner nodes
	 */
	void CornerEqs(int mincorners, int nelconsider, std::set<int> &cornernodes);
	
protected:
	/** @brief Number of elements in the graph */
	int fNElements;
	
	/** @brief Number of nodes in the graph */
	int fNNodes;
	/** @brief Number of equations associated with each node */
	TPZVec<int> fNodeWeights;

	/** @brief Node number of each element*/
	TPZVec<int> fElementGraph;

	/** @brief Indicates for each element the index of the first entry with
	 * fElementGraph for that element
	 */

    /** Size of the vector fNElements+1 */
	TPZVec<int> fElementGraphIndex;
	
};

class TPZCompMesh;
/** @brief Makes resequence to renumbering */
void ResequenceByGeometry(TPZCompMesh *cmesh,const TPZVec<REAL> &normal);

#endif //TPZRENUMBERING_H











