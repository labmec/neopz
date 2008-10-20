/* Generated by Together */

#ifndef TPZBOOSTGRAPH_H
#define TPZBOOSTGRAPH_H
#include "pzrenumbering.h"

#ifdef USING_BOOST

#include <vector>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/sloan_ordering.hpp>
#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/bandwidth.hpp>
#include <boost/graph/profile.hpp>
#include <boost/graph/wavefront.hpp>

  //using namespace boost;
  //using namespace std;

  /**
   * Defining the graph type 
   * From Boost Sloan_Ordering example
   */
  typedef boost::adjacency_list<
    boost::setS, 
    boost::vecS, 
    boost::undirectedS,
    boost::property<
    boost::vertex_color_t, 
    boost::default_color_type,
    boost::property<
    boost::vertex_degree_t,
    int,
    boost::property<
    boost::vertex_priority_t,
    double > > > > Graph;
  
  typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
  typedef boost::graph_traits<Graph>::vertices_size_type size_type;

  typedef std::pair<std::size_t, std::size_t> Pair;


/**
 * Implements an interface to the BGL (Boost Graph Library) for graph operations.
 * BGL is implemented on top of Boost C++
 * It behaves as close as possible to STL implementations, aproaches edges and vertices access by means of iterators.
 * Still bein evaluated
 * @author Gustavo Camargo Longhin
 * @version 0.1
 * @since October 2007 
 */
class TPZBoostGraph : public TPZRenumbering {
public:

  enum GraphType { KMC, KMCExpensive, Sloan };
  
  GraphType fGType;
  
  TPZBoostGraph() : TPZRenumbering(), fGType(KMC)
  {
  }
  /**
    * Simple constructor 
    */
  TPZBoostGraph(int NElements, int NNodes);
  

  virtual ~TPZBoostGraph();

  /**
   * Perform the renumbering of elements. The aim of this operation is to minimize the
   * band of the resulting stiffeness matrix.
   */
  void ResequenceOld(TPZVec<int> &perm, TPZVec<int> &inverseperm);
  void Resequence(TPZVec<int> &perm, TPZVec<int> &inverseperm);
  /**
   * This method declares the element graph to the object
   * The first vector contains the element node number
   * The second vector contains the number of nodes of each element
   */
  virtual void SetElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex);
  virtual void SetElementGraphOld(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex);

  /**
   * Sets the number of equations associated with each node
   * The derived class may or may not take this data into
   * consideration
   */
  virtual void SetNodeWeights(TPZVec<int> &weights);

  /**
   * This will reset all datastructures the object may contain.
   * Node resequencing algorithms may require a possibly large
   * amount of temporary data
   */
  virtual void ClearDataStructures();
private:
  /**
   * Defines the list of Edges for the current Graph
   */
  std::vector<Pair> m_Edges;
  /**
   * Defines the Graph itself.
   * This structure will be submitted to the Boost kernel
   */
  Graph m_Graph;
  //Creating two iterators over the vertices
  //graph_traits<Graph>::vertex_iterator ui, ui_end;

  //Creating a property_map with the degrees of the degrees of each vertex
  boost::property_map<Graph,boost::vertex_degree_t>::type m_Degrees;
//    = get(vertex_degree, G);
//   for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
//     deg[*ui] = degree(*ui, G);

  //Creating a property_map for the indices of a vertex
  boost::property_map<Graph, boost::vertex_index_t>::type m_Index_map;
  // = get(vertex_index, G);
  /**Node number of each element*/
  TPZVec<int> fElementGraph;
  

  /**Indicates for each element the index of the first entry with
     fElementGraph for that element
     Size of the vector fNElements+1*/
  TPZVec<int> fElementGraphIndex;
  
  std::vector<int> m_Connects;

};

#endif // USING_BOOST

#endif //TPZBOOSTGRAPH_H
/*    Graph m_Graph;
template <class Graph, class OutputIterator,
                  class ColorMap, class DegreeMap, 
                  class PriorityMap>
  OutputIterator sloan_ordering(Graph& g,
                 OutputIterator permutation, 
                 ColorMap color, 
                 DegreeMap degree, 
                 PriorityMap priority );
*/
