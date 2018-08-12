/**
 * @file
 * @brief Contains the TPZRenumbering class which defines the behavior to implementing node sequence numbering optimization.
 */

#ifndef TPZRENUMBERING_H
#define TPZRENUMBERING_H

#include "pzvec.h"
#include <set>
#include "TPZSavable.h"

class TPZCompMesh;

/** 
 * @brief This abstract class which defines the behavior which derived classes need to implement \n
 * for implementing node sequence numbering optimization. \ref util "Utility"
 * @ingroup util
 */
class TPZRenumbering : public TPZSavable {
public:
	
	int fHDivPermute;
	
	TPZRenumbering() : fNElements(0), fNNodes(0)
	{
	}
	
	TPZRenumbering(int64_t NElements, int64_t NNodes);
	
	virtual ~TPZRenumbering()
	{
	}
        
        int ClassId() const;
        void Read(TPZStream& buf, void* context);
        void Write(TPZStream& buf, int withclassid) const;

	void SetElementsNodes(int64_t NElements, int64_t NNodes)
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
	void SetElementGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex);
	
	/** @brief Sets the number of equations associated with each node */
	/** The derived class may or may not take this data into consideration */
	void SetNodeWeights(TPZVec<int> &weights)
	{
		fNodeWeights = weights;
	}
	
	/** @brief This will reset all datastructures the object may contain. */
	/** Node resequencing algorithms may require a possibly large amount of temporary data */
	virtual void ClearDataStructures();
	
	virtual void Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &iperm)
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
	void ConvertGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, TPZVec<int64_t> &nodegraph, TPZVec<int64_t> &nodegraphindex);
	
	/** @brief Convert a traditional elgraph to an element to element graph */
	void ConvertToElementoToElementGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, TPZVec<int64_t> &eltotelgraph, TPZVec<int> &eltoelweight, TPZVec<int64_t> &eltoelgraphindex);

	/** @brief Stores the graph of nodes to elements */
	void NodeToElGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, TPZVec<int64_t> &nodetoelgraph, TPZVec<int64_t> &nodetoelgraphindex);
	
	/**
	 * @brief Will assign a color to the nodes in the graph such that no two connected nodes have the same color
	 * the return value indicates the number of colors in the graph
	 */
	int64_t ColorNodes(TPZVec<int64_t> &nodegraph, TPZVec<int64_t> &nodegraphindex, TPZVec<int> &family, TPZVec<int> &colors);
	
        /**
	 * @brief Assigns a color to the elements in the elementIndices list such that 
         * two elements that share a connect have different colors.
	 * The return value indicates the number of colors.
	 */
        static int64_t ColorElements(const TPZCompMesh *cmesh, const TPZVec<int64_t> &elementIndices, TPZVec<int64_t> &elementColors);
	
	/** @brief Prints graph */
	void Print(TPZVec<int64_t> &grapho, TPZVec<int64_t> &graphoindex, const char *name = 0, std::ostream &out = std::cout);
	
	/**
	 * @brief Analyzes the graph, finds the corner nodes \n
	 * Number of elements which should be considered for determining corner nodes
	 */
	void CornerEqs(unsigned int mincorners, int64_t nelconsider, std::set<int> &eligible, std::set<int> &cornernodes);
	
protected:
	/** @brief Number of elements in the graph */
	int64_t fNElements;
	
	/** @brief Number of nodes in the graph */
	int64_t fNNodes;
        
	/** @brief Number of equations associated with each node */
	TPZVec<int> fNodeWeights;

	/** @brief Node number of each element*/
	TPZVec<int64_t> fElementGraph;

	/** @brief Indicates for each element the index of the first entry with
	 * fElementGraph for that element
	 * The size of this vector is fNElements+1 
         */
	TPZVec<int64_t> fElementGraphIndex;
	
};

class TPZCompMesh;
/** @brief Makes resequence to renumbering */
void ResequenceByGeometry(TPZCompMesh *cmesh,const TPZVec<REAL> &normal);

#endif //TPZRENUMBERING_H











