/**
 * @file
 * @brief Contains the TPZBoostGraph class which implements an interface to the BGL for graph operations.
 */

#ifndef TPZBOOSTGRAPH_H
#define TPZBOOSTGRAPH_H
#include "TPZRenumbering.h"

#ifdef USING_BOOST

#include "pzvec.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>
#include <boost/graph/properties.hpp>

/**
 * @brief Defining the graph type \n
 * From Boost Sloan_Ordering example
 */
typedef boost::adjacency_list<boost::setS,boost::vecS,boost::undirectedS,
boost::property<boost::vertex_color_t,boost::default_color_type,
boost::property<boost::vertex_degree_t,int,
boost::property<boost::vertex_priority_t,double > > > > Graph;

typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::vertices_size_type size_type;

typedef std::pair<std::size_t, std::size_t> Pair;


/**
 * @brief Implements an interface to the BGL (Boost Graph Library) for graph operations. \ref util "Utility"
 * @author Gustavo Camargo Longhin
 * @version 0.1
 * @since October 2007 
 */
/**
 * BGL is implemented on top of Boost C++
 * It behaves as close as possible to STL implementations, aproaches edges and vertices access by means of iterators.
 * Still bein evaluated
 */
class TPZBoostGraph : public TPZRenumbering {
    
public:
	
	enum GraphType { KMC, KMCExpensive, Sloan };
	
	
	
	TPZBoostGraph(GraphType tp) : TPZRenumbering(), fGType(tp)
	{
	}
    
	/** @brief Simple constructor */
	TPZBoostGraph(int64_t NElements, int64_t NNodes);

	virtual ~TPZBoostGraph();

	void CompressedResequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm);
	
	/**
	 * @brief Perform the renumbering of elements. The aim of this operation is to minimize the
	 * band of the resulting stiffeness matrix.
	 */
	void ResequenceOld(TPZVec<int> &perm, TPZVec<int> &inverseperm);
	void Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm);
	void setGType(GraphType M) { fGType = M; }
	/**
	 * @brief This will reset all datastructures the object may contain. \n
	 * Node resequencing algorithms may require a possibly large amount of temporary data
	 */
	virtual void ClearDataStructures();
private:
	/** @brief Defines the list of Edges for the current Graph */
	TPZVec<Pair> m_Edges;
	/** @brief Defines the Graph itself. This structure will be submitted to the Boost kernel */
	Graph m_Graph;
	//Creating two iterators over the vertices
	//graph_traits<Graph>::vertex_iterator ui, ui_end;
	
	/** @brief Creating a property_map with the degrees of the degrees of each vertex */
	boost::property_map<Graph,boost::vertex_degree_t>::type m_Degrees;
	//    = get(vertex_degree, G);
	//   for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
	//     deg[*ui] = degree(*ui, G);
	
	/** @brief Creating a property_map for the indices of a vertex */
	boost::property_map<Graph, boost::vertex_index_t>::type m_Index_map;
	// = get(vertex_index, G);  
	TPZVec<int64_t> m_Connects;
	
    GraphType fGType;
};



#endif // USING_BOOST

#endif //TPZBOOSTGRAPH_H
