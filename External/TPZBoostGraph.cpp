/**
 * @file
 * @brief Contains the implementation of the TPZBoostGraph methods.
 */

#include "TPZBoostGraph.h"
#include <fstream>
#include "pzmanvector.h"

#ifdef USING_BOOST

#include <boost/graph/properties.hpp>
#include <boost/graph/compressed_sparse_row_graph.hpp>
#include "pzlog.h"
//#include <sys/time.h>
#include <stdio.h>

#ifdef PZ_LOG
static TPZLogger logger("pz.external.boostgraph");
#endif

using namespace boost;
using namespace std;

TPZBoostGraph::TPZBoostGraph(int64_t NElements, int64_t NNodes) : TPZRenumbering(NElements,NNodes),fGType(Sloan)//,fGType(KMCExpensive)
{
    m_Graph.clear();
}


TPZBoostGraph::~TPZBoostGraph()
{
    //Performe termination tasks
}


void TPZBoostGraph::ClearDataStructures()
{
    m_Edges.clear();
	TPZRenumbering::ClearDataStructures();
}

void TPZBoostGraph::CompressedResequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm)
{
    

    /* if the graph is empty, trivial */
    if (this->fNNodes == 0)
    {
        perm.resize(0);
        inverseperm.resize(0);
        return;
    }
    /* define type graph */
    typedef boost::compressed_sparse_row_graph<boost::directedS> BoostGraph;
    
    /* this code is a copy modified from the method ConvertGraph. Used here to create a  Compressed Sparse Row Graph Boost */
    int64_t nod,el;
    TPZVec<int64_t> nodtoelgraphindex;
    TPZVec<int64_t> nodtoelgraph;
    
    NodeToElGraph(fElementGraph,fElementGraphIndex,nodtoelgraph,nodtoelgraphindex);
    
    std::vector<std::pair<std::size_t, std::size_t> > edges;
    
    size_t maxsize = 100000;
    int resize_times = 1;
    edges.reserve(1000000);
    
    for(nod=0; nod<fNNodes; nod++)
    {
        int64_t firstel = nodtoelgraphindex[nod];
        int64_t lastel = nodtoelgraphindex[nod+1];
        std::set<int64_t> nodecon;
        for(el=firstel; el<lastel; el++)
        {
            int64_t gel = nodtoelgraph[el];
            int64_t firstelnode = fElementGraphIndex[gel];
            int64_t lastelnode = fElementGraphIndex[gel+1];
            nodecon.insert(&fElementGraph[firstelnode],&fElementGraph[(lastelnode-1)]+1);
        }
        nodecon.erase(nod);
        
        std::set<int64_t>::iterator it;
        for(it = nodecon.begin(); it!= nodecon.end(); it++)
        {
            edges.push_back(std::make_pair(nod, *it));
            edges.push_back(std::make_pair(*it, nod));
            if (edges.size() > (edges.size()+maxsize*resize_times)) {
                edges.reserve(edges.size()+maxsize*(++resize_times));
            }
        }
    }
    
    BoostGraph G(boost::edges_are_unsorted_multi_pass, edges.begin(), edges.end(), fNNodes);
    
    boost::property_map<BoostGraph, boost::vertex_index_t>::type boost_index_map;
    boost_index_map = get(boost::vertex_index, G);
    
    // Compute graph re-ordering
    std::vector<std::size_t> inv_perm(fNNodes);
    
    boost::cuthill_mckee_ordering(G, inv_perm.begin());
    
    perm.Resize(fNNodes);
    inverseperm.Resize(fNNodes);
    
    for (std::size_t i = 0; i < fNNodes; ++i)
    {
        perm[inv_perm[i]] = i;
        inverseperm[i]=inv_perm[i];
    }

}
void TPZBoostGraph::Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm)
{
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "TPZBoostGraph::Resequence started \n";
        Print(fElementGraph,fElementGraphIndex,"Element graph when entering Resequence",sout);
        LOGPZ_DEBUG(logger,sout.str());
    }
#endif
    if (this->fNNodes == 0) {
        perm.resize(0);
        inverseperm.resize(0);
        return;
    }
    Graph G;
    size_type i;
    size_type elgraphsize = fElementGraphIndex.NElements()-1;
    TPZManVector<int64_t> nconnects(fNNodes,0);
    
    for(i=0; i < elgraphsize; i++)
    {
        int64_t first, second;
        first = fElementGraphIndex[i];
        second = fElementGraphIndex[i+1];
        int64_t j,k;
        for(j=first; j< second; j++)
        {
            for(k=j+1; k<second; k++)
            {
                add_edge(fElementGraph[j],fElementGraph[k],G);
                nconnects[fElementGraph[j]]++;
            }
        }
    }
    
    for(i=0; i< (size_type)nconnects.size(); i++)
    {
        if(!nconnects[i])
        {
            add_edge(i,i,G);
        }
    }
    
    graph_traits<Graph>::vertex_iterator ui, ui_end;
    
    property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
    for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
        deg[*ui] = degree(*ui, G);
    
#ifdef PZ_LOG
    if (logger.isDebugEnabled())
    {
        std::stringstream sout;
        sout << "NNodes " << fNNodes << std::endl;
        sout << "NElements " << fNElements << std::endl;
        
        sout << "Original Bandwidth: " << bandwidth(G) << std::endl;
        sout << " Original profile: "
        << profile(G)
        << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    /*  std::cout << " Original max_wavefront: "
     << max_wavefront(G)
     << std::endl;
     std::cout << " Original aver_wavefront: "
     << aver_wavefront(G)
     << std::endl;
     std::cout << " Original rms_wavefront: "
     << rms_wavefront(G)
     << std::endl;*/
    //   std::cout << "Number of Vertices " << num_vertices(G) << std::endl;
    //   std::cout << "Number of Edges " << num_edges(G) << std::endl;
    int64_t nVertices = num_vertices(G);
    //  std::cout << "Number of Vertices " << nVertices << std::endl;
    TPZVec<Vertex> inv_perm(nVertices);
    TPZVec<size_type> l_perm(nVertices);
    for(size_type i = 0; i < (size_type)l_perm.size(); i++) l_perm[i]=-1;
    for(size_type i = 0; i < (size_type)l_perm.size(); i++) inv_perm[i]=0;
    perm.Resize(fNNodes,-1);
    
    switch(fGType)
    {
        case KMC:
            
        {
            Vertex s = vertex(0, G);
            //reverse cuthill_mckee_ordering
            cuthill_mckee_ordering(G, s, inv_perm.begin(), get(vertex_color, G),
                                   get(vertex_degree, G));
#ifdef PZ_LOG
            if(logger.isDebugEnabled())
            {
                LOGPZ_DEBUG(logger,"Reverse Cuthill-McKee ordering:")
            }
#endif
            
        }
            break;
            
        case KMCExpensive:
        {
            cuthill_mckee_ordering(G, inv_perm.begin(), get(vertex_color, G),
                                   get(vertex_degree, G));
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger, "Reverse Expensive Cuthill-McKee ordering:")
            }
#endif
            
        }
            break;
        case Sloan:
        {
#ifdef PZ_LOG
            if (logger.isDebugEnabled()) {
                LOGPZ_DEBUG(logger,  "Sloan ordering:")
            }
#endif
            
            //Setting the start node
            Vertex s = vertex(0, G);
            int ecc;   //defining a variable for the pseudoperipheral radius
            
            //Calculating the pseudoeperipheral node and radius
            Vertex e = pseudo_peripheral_pair(G, s, ecc, get(vertex_color, G), get(vertex_degree, G) );
#ifdef PZ_LOG
            if (logger.isDebugEnabled())
            {
                std::stringstream sout;
                sout << "Sloan Starting vertex: " << s << endl;
                sout << "Sloan Pseudoperipheral vertex: " << e << endl;
                sout << "Sloan Pseudoperipheral radius: " << ecc << endl << endl;
                LOGPZ_DEBUG(logger, sout.str())
            }
#endif
            //Sloan ordering
            sloan_ordering(G, s, e, inv_perm.begin(), get(vertex_color, G),
                           get(vertex_degree, G), get(vertex_priority, G));
        }
            break;
            
    }
    
    for (size_type c = 0; c != (size_type)inv_perm.size(); ++c)
    {
        l_perm[inv_perm[c]] = c;
    }
    
    /*
     std::cout << "l_perm = ";
     for (size_type i = 0; i < l_perm.size(); ++i)
     {
     std::cout << i << "/" << l_perm[i] << " ";
     }
     std::cout << std::endl;
     */
    //    property_map<Graph, vertex_index_t>::type
    //    index_map = get(vertex_index, G);
    
    /*    std::cout << "  bandwidth: "
     << bandwidth(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
     << std::endl;
     std::cout << "  profile: "
     << profile(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
     << std::endl;
     std::cout << "  max_wavefront: "
     << max_wavefront(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
     << std::endl;
     std::cout << "  aver_wavefront: "
     << aver_wavefront(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
     << std::endl;
     std::cout << "  rms_wavefront: "
     << rms_wavefront(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
     << std::endl;
     std::cout << "  inverse map bandwidth: "
     << bandwidth(G, make_iterator_property_map(&inv_perm[0], index_map, inv_perm[0]))
     << std::endl;
     std::cout << "  inverse map profile: "
     << profile(G, make_iterator_property_map(&inv_perm[0], index_map, inv_perm[0]))
     << std::endl;
     std::cout << "  inverse map max_wavefront: "
     << max_wavefront(G, make_iterator_property_map(&inv_perm[0], index_map, inv_perm[0]))
     << std::endl;
     std::cout << "  inverse map aver_wavefront: "
     << aver_wavefront(G, make_iterator_property_map(&inv_perm[0], index_map, inv_perm[0]))
     << std::endl;
     std::cout << "  inverse map rms_wavefront: "
     << rms_wavefront(G, make_iterator_property_map(&inv_perm[0], index_map, inv_perm[0]))
     << std::endl;*/
    perm.Resize(l_perm.size());
    inverseperm.Resize(inv_perm.size());
    for(i=0; i<(size_type)l_perm.size(); i++)
    {
        perm[i] = l_perm[i];
        inverseperm[i] = inv_perm[i];
    }
#ifdef PZ_LOG
    if (logger.isDebugEnabled()) {
        LOGPZ_DEBUG(logger, "TPZBoostGraph::Resequence finished")
    }
#endif
    
}
/*
 
 std::cout << "l_perm[index_map[inv_perm[c]]] " <<
 l_perm[index_map[inv_perm[c]]] << " index_map[inv_perm[c]] " <<
 index_map[inv_perm[c]] <<
 " Inv_Perm[c] " << inv_perm[c] << "\n";
 
 
 void TPZBoostGraph::Resequence(TPZVec<int> &permuta, TPZVec<int> &inverseperm)
 {
 //Creating a graph and adding the edges from above into it
 Graph G;//(10);
 for (int i = 0; i < m_Edges.size(); ++i)
 add_edge(m_Edges[i].first, m_Edges[i].second, G);
 
 //Creating two iterators over the vertices
 graph_traits<Graph>::vertex_iterator ui, ui_end;
 
 //Creating a property_map with the degrees of the degrees of each vertex
 property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
 for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
 deg[*ui] = degree(*ui, G);
 
 //Creating a property_map for the indices of a vertex
 property_map<Graph, vertex_index_t>::type index_map = get(vertex_index, G);
 
 std::cout << "original bandwidth: " << bandwidth(G) << std::endl;
 std::cout << "original profile: " << profile(G) << std::endl;
 std::cout << "original max_wavefront: " << max_wavefront(G) << std::endl;
 std::cout << "original aver_wavefront: " << aver_wavefront(G) << std::endl;
 std::cout << "original rms_wavefront: " << rms_wavefront(G) << std::endl;
 
 
 //Creating a vector of vertices
 std::vector<Vertex> sloan_order(num_vertices(G));
 //Creating a vector of size_type
 std::vector<size_type> perm(num_vertices(G));
 permuta.Resize(perm.size());
 {
 //sloan_ordering
 sloan_ordering(G, sloan_order.begin(),
 get(vertex_color, G),
 make_degree_map(G),
 get(vertex_priority, G) );
 
 cout << endl << "Sloan ordering without a start-vertex:" << endl;
 cout << "  ";
 for (std::vector<Vertex>::const_iterator i=sloan_order.begin();
 i != sloan_order.end(); ++i)
 {
 cout << index_map[*i] << " ";
 permuta[*i]=index_map[*i];
 }
 cout << endl;
 
 for (size_type c = 0; c != sloan_order.size(); ++c)
 perm[index_map[sloan_order[c]]] = c;
 std::cout << "  bandwidth: "
 << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  profile: "
 << profile(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  max_wavefront: "
 << max_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  aver_wavefront: "
 << aver_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  rms_wavefront: "
 << rms_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 }
 
 {
 
 //Setting the start node
 Vertex s = vertex(0, G);
 int ecc;   //defining a variable for the pseudoperipheral radius
 
 //Calculating the pseudoeperipheral node and radius
 Vertex e = pseudo_peripheral_pair(G, s, ecc, get(vertex_color, G), get(vertex_degree, G) );
 
 cout << endl;
 cout << "Starting vertex: " << s << endl;
 cout << "Pseudoperipheral vertex: " << e << endl;
 cout << "Pseudoperipheral radius: " << ecc << endl << endl;
 
 
 
 //Sloan ordering
 sloan_ordering(G, s, e, sloan_order.begin(), get(vertex_color, G),
 get(vertex_degree, G), get(vertex_priority, G));
 
 cout << "Sloan ordering starting at: " << s << endl;
 cout << "  ";
 
 for (std::vector<Vertex>::const_iterator i = sloan_order.begin();
 i != sloan_order.end(); ++i)
 {
 cout << index_map[*i] << " ";
 permuta[*i]=index_map[*i];
 }
 cout << endl;
 
 for (size_type c = 0; c != sloan_order.size(); ++c)
 perm[index_map[sloan_order[c]]] = c;
 std::cout << "  bandwidth: "
 << bandwidth(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  profile: "
 << profile(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  max_wavefront: "
 << max_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  aver_wavefront: "
 << aver_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 std::cout << "  rms_wavefront: "
 << rms_wavefront(G, make_iterator_property_map(&perm[0], index_map, perm[0]))
 << std::endl;
 }
 
 }
 */
/*
 {
 //perm.Resize(fNNodes+1);
 Graph G;
 for (int i = 0; i < m_Edges.size(); ++i)
 add_edge(m_Edges[i].first, m_Edges[i].second, G);
 
 graph_traits<Graph>::vertex_iterator ui, ui_end;
 
 property_map<Graph,vertex_degree_t>::type deg = get(vertex_degree, G);
 for (boost::tie(ui, ui_end) = vertices(G); ui != ui_end; ++ui)
 deg[*ui] = degree(*ui, G);
 
 property_map<Graph, vertex_index_t>::type
 index_map = get(vertex_index, G);
 
 std::cout << "original bandwidth: " << bandwidth(G) << std::endl;
 std::cout << "Number of Vertices " << num_vertices(G) << std::endl;
 std::cout << "Number of Edges " << num_edges(G) << std::endl;
 
 std::vector<Vertex> inv_perm(num_vertices(G));
 std::vector<size_type> l_perm(num_vertices(G));
 
 {
 //reverse cuthill_mckee_ordering
 //    cuthill_mckee_ordering(G, inv_perm.rbegin(), get(vertex_color, G),
 //      make_degree_map(G));
 cuthill_mckee_ordering(G, inv_perm.begin(), get(vertex_color, G),
 make_degree_map(G));
 cout << "Reverse Cuthill-McKee ordering:" << endl;
 cout << "  ";
 perm.Resize(l_perm.size(),0);
 int orignode = 0;
 for (std::vector<Vertex>::const_iterator i=inv_perm.begin();i != inv_perm.end();++i)
 {
 cout << *i << ": "  << index_map[*i] << " ";
 perm[orignode] = index_map[*i];
 orignode++;
 }
 cout << endl;
 
 for (size_type c = 0; c != inv_perm.size(); ++c)
 {
 l_perm[index_map[inv_perm[c]]] = c;
 std::cout << "l_perm[index_map[inv_perm[c]]] " <<
 l_perm[index_map[inv_perm[c]]] << " index_map[inv_perm[c]] " <<
 index_map[inv_perm[c]] <<
 " Inv_Perm[c] " << inv_perm[c] << "\n";
 }
 std::cout << "  bandwidth: "
 << bandwidth(G, make_iterator_property_map(&l_perm[0], index_map, l_perm[0]))
 << std::endl;
 }
 
 
 
 }
 */

/*
 std::vector<Vertex> inv_perm(num_vertices(m_Graph));
 std::vector<size_type> l_perm(num_vertices(m_Graph));
 {
 Vertex s = vertex(0, m_Graph);
 //reverse cuthill_mckee_ordering
 cuthill_mckee_ordering(m_Graph, s, inv_perm.rbegin(), get(vertex_color, m_Graph),
 get(vertex_degree, m_Graph));
 for (std::vector<Vertex>::const_iterator i=inv_perm.begin();
 i != inv_perm.end(); ++i)
 cout << m_Index_map[*i] << " ";
 cout << endl;
 
 for (size_type c = 0; c != inv_perm.size(); ++c)
 l_perm[m_Index_map[inv_perm[c]]] = c;
 std::cout << "  bandwidth: "
 << bandwidth(m_Graph, make_iterator_property_map(&l_perm[0], m_Index_map, l_perm[0]))
 << std::endl;
 }
 int i;
 int nVertex = l_perm.size();
 perm.Resize(nVertex);
 inverseperm.Resize(nVertex);
 for(i = 0; i < nVertex; i++)
 {
 perm[i]=l_perm[i];
 inverseperm[i]=inv_perm[i];
 }
 
 */
#endif
