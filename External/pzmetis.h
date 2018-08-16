/**
 * @file
 * @brief Contains TPZMetis class which implements the renumbering for elements of a mesh to minimize the band.
 */

#ifndef TPZMETIS_H
#define TPZMETIS_H

#include "pzstack.h"
#include "TPZRenumbering.h"
#include <iostream>

/**
 * @ingroup util
 * @brief Implements renumbering for elements of a mesh. \ref util "Utility"
 */
class TPZMetis : public TPZRenumbering {
public:
	/**
	 * @brief Perform the renumbering of elements. The aim of this operation is to minimize the
	 * band of the resulting stiffeness matrix.
	 */
	virtual void Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm);
	
    /** @brief Constructor. */
	/** 
	 * Instantiates an object which will compute the resequencing
	 * scheme of the metis package.
	 */
    TPZMetis(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes)
    {
        
    }
    
    TPZMetis();
	/** @brief Destructor */
	virtual ~TPZMetis() {}
	/** @brief Prints the current object data structure. */
	void Print(std::ostream &out);
	void Print(std::ostream &out,char * title);
	/**
	 * @brief Subdivides a Graph in nParts.
	 * @param nParts Number of subdomains the Original domain must be divided to.
	 * @param Domains A vector the subdomain index for each vertex
	 */
	/** The Adjacency list works according to the MeTiS specification as found on MeTiS manual. \n
	 * The graph information is stored using the CSR (Compressed Storage Format) \n
	 * The CSR for a graph with 'n' vertices and 'm' edges is represented using two vector XAdj \n
	 * and Adjacency. XAdj is of size 'n + 1' while Adjacency is of size '2 * m'. \n
	 * The graph structures is stored as follows: \n
	 * For the ith vertex, its Adjacency list is stored in the Adjacency vector positions \n
	 * startint in XAdj[i] until XAdj[i+1]-1.
	 */
	void Subdivide(int nParts, TPZVec < int >  & Domains);
private:
	
};

#endif //TPZMETIS_H
