/**
 * @file
 * @brief Contains the implementation of the TPZNodesetCompute methods. 
 */
//
// C++ Implementation: tpznodesetcompute
//
// Description:
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2004
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpznodesetcompute.h"
#include "pzstack.h"
#include "pzblock.h"
#include <set>
#include <map>
#include <algorithm>
#include <iterator>
#include <fstream>
#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.mesh.nodesetcompute");
#endif

//static std::ofstream out("nodeset.txt");

TPZNodesetCompute::TPZNodesetCompute()
{
  fMaxSeqNum = -1;
  fMaxLevel = -1;
}


TPZNodesetCompute::~TPZNodesetCompute()
{
}


    /**
    Group the node graph as passed by the parameters
    */
void TPZNodesetCompute::AnalyseGraph()
{
  fMaxSeqNum = -1;
  fMaxLevel = -1;
  int64_t nnodes = fNodegraphindex.NElements()-1;
  fSeqNumber.Resize(nnodes);
  this->fLevel.Resize(nnodes);
  fSeqNumber.Fill(-1);
	fIsIncluded.Resize(nnodes);
	fIsIncluded.Fill(0);
  fLevel.Fill(0);
  TPZVec<std::set<int64_t> > nodeset(nnodes);
  int64_t in;
  for(in=0; in<nnodes; in++) 
  {
	  if(!(in%1000))
	  {
		  std::cout << "*";
		  std::cout.flush();
	  }
	  AnalyseNode(in,nodeset);
  }
	std::cout << std::endl;
  
}

/**
 * This method will analyse the set inclusion of the current node, calling the method
 * recursively if another node need to be analysed first
 */
void TPZNodesetCompute::AnalyseNode(int64_t node, TPZVec< std::set<int64_t> > &nodeset)
{
  if(fSeqNumber[node] != -1) return;
  if(! nodeset[node].size()) 
  {
    nodeset[node].insert(&fNodegraph[fNodegraphindex[node]],&fNodegraph[fNodegraphindex[node+1]]);
    nodeset[node].insert(node);
  }
  int minlevel = 0;
  std::set<int64_t>::iterator it;
  TPZStack<int64_t> equalnodes;
  for(it = nodeset[node].begin(); it != nodeset[node].end(); it++)
  {
    int64_t othernode = *it;
    if(othernode == node) continue;
    // build the data structure of the connected node
    if(! nodeset[othernode].size())
    {
      nodeset[othernode].insert(&fNodegraph[fNodegraphindex[othernode]],&fNodegraph[fNodegraphindex[othernode+1]]);
      nodeset[othernode].insert(othernode);
    }
    // the other node is included in my connectivity
    bool inc = includes(nodeset[node].begin(),nodeset[node].end(),nodeset[othernode].begin(),nodeset[othernode].end());
    // the other node has a different connectivity
    bool diff = nodeset[node] != nodeset[othernode];
    if( inc && diff) 
    {
      // Analyse the connectivity graph of the other node
      // Its graph is smaller than mine
		fIsIncluded[othernode] = 1;
      AnalyseNode(othernode,nodeset);
      minlevel = minlevel < fLevel[othernode]+1 ? fLevel[othernode]+1 : minlevel;
#ifdef PZ_LOG
		if(fLevel[othernode] >= 0)
		{
			std::stringstream sout;
			sout << "The level of " << node << " is increased because of " << othernode << " ";
			sout << "Level of othernode " << fLevel[othernode] << " Seqnumber othernode " << fSeqNumber[othernode];
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
    }
    else if(!diff)
    {
      // The graphs are equal. These nodes should be grouped
		if(fIsIncluded[node]) 
		{
			fIsIncluded[othernode] = 1;
		}
      equalnodes.Push(othernode);
    }
    // the other node has been analysed (and is probably not equal...)
    if(inc && diff && fSeqNumber[othernode] != -1)
    {
		fIsIncluded[othernode] = 1;
      minlevel = minlevel < fLevel[othernode]+1 ? fLevel[othernode]+1 : minlevel;
    } else if(!diff && fSeqNumber[othernode] != -1)
    {
      // the level should be at least the level of the other node
      minlevel = minlevel < fLevel[othernode] ? fLevel[othernode] : minlevel;
		// should not happen because if the nodes are equal, they both have the same sequence number
		DebugStop();
    }
  }
  // assign a sequence number to the node
  fMaxSeqNum++;
  fSeqNumber[node] = fMaxSeqNum;
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << "Assigning Seq Number " << fMaxSeqNum << " and level " << minlevel << " to nodes " << node << " " << equalnodes;
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  fSeqCard.Push(1);
  // the maximum level of the set of nodes
  fMaxLevel = fMaxLevel < minlevel ? minlevel : fMaxLevel;
  // assign the level of the node
  fLevel[node] = minlevel;
	nodeset[node].erase(node);
  
  // memory clean up
  for(it = nodeset[node].begin(); it != nodeset[node].end(); it++)
  {
    int64_t othernode = *it;
	  int level = fLevel[othernode];
	  int64_t seq = fSeqNumber[othernode];
    if(seq != -1 && level <= minlevel) 
	{
		nodeset[othernode].clear();
		if(othernode == node)
		{
			LOGPZ_ERROR(logger," othernode equal to node!!")
			DebugStop();
			break;
		}
	}
  }
  nodeset[node].clear();
  // initialize the datastructure of the nodes which have the same connectivity
  int64_t neq = equalnodes.NElements();
  int64_t ieq;
  for(ieq=0; ieq<neq; ieq++)
  {
    fSeqNumber[equalnodes[ieq]] = fMaxSeqNum;
    fLevel[equalnodes[ieq]] = minlevel;
    // I think this is overkill : this nodeset should already be empty...
    nodeset[equalnodes[ieq]].clear();
    // the number of nodes associated with this sequence number
    fSeqCard[fMaxSeqNum]++;
  }
}

void TPZNodesetCompute::BuildNodeGraph(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex)
{
  int64_t nnodes = fSeqNumber.NElements();
  blockgraphindex.Resize(fSeqCard.NElements()+1);
  blockgraph.Resize(nnodes);
  int64_t seq = 0;
  blockgraphindex[0] = 0;
  blockgraphindex[1] = 0;
  for(seq=2; seq<=fSeqCard.NElements(); seq++)
  {
    blockgraphindex[seq] = blockgraphindex[seq-1]+fSeqCard[seq-2];
  }
  int64_t in;
  for(in = 0; in< nnodes; in++)
  {
    int64_t seqn = fSeqNumber[in];
    blockgraph[blockgraphindex[seqn+1]] = in;
    blockgraphindex[seqn+1]++;
  }
}

void TPZNodesetCompute::BuildVertexGraph(TPZStack<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex)
{
  std::map<int64_t,int64_t> vertices;
  int64_t nnodes = fSeqNumber.NElements();
  int64_t in;
  for(in=0; in<nnodes; in++)
  {
    if(fLevel[in] == this->fMaxLevel) vertices[fSeqNumber[in]] = in;
  }
  int64_t nvert = vertices.size();
  blockgraphindex.Resize(nvert+1);
  blockgraphindex[0] = 0;
  int iv = 1;
  std::map<int64_t,int64_t>::iterator it;
  for(it=vertices.begin(); it != vertices.end(); it++)
  {
    int node = (*it).second;
    blockgraph.Push(node);
    std::set<int64_t> vertexset;
    std::set<int> included;
    std::set<int> notincluded;
    // The set will contain the connectivity of the node
    BuildNodeSet(node,vertexset);
    included.insert((*it).first);
    std::set<int64_t>::iterator versetit;
    for(versetit = vertexset.begin(); versetit != vertexset.end(); versetit++)
    {
      int64_t linkednode = *versetit;
      if(linkednode == node) continue;
      int64_t seq = fSeqNumber[linkednode];
      
      // if the sequence number of the equation is already included in the nodeset of the vertex, put it in the blockgraph
      // this node has equal connectivity with other node
      if(included.count(seq))
      {
        blockgraph.Push(linkednode);
      }
      // if the seq number is recognized as not to be included stop analysing
      // this node has equal connectivity with other node
      else if(notincluded.count(seq))
      {
        continue;
      }
      // the equation hasn t been analysed yet 
      else 
      {
        std::set<int64_t> locset;
        BuildNodeSet(linkednode,locset);
        // if its nodeset is included in the vertexset
        if(includes(vertexset.begin(),vertexset.end(),locset.begin(),locset.end()))
        {
          included.insert(seq);
          blockgraph.Push(linkednode);
        }
        // else this sequence number should not be analysed anymore
        else
        {
          notincluded.insert(seq);
        }
      }
    }
    blockgraphindex[iv] = blockgraph.NElements();
    iv++;
  }
}

void TPZNodesetCompute::BuildNodeSet(int64_t node, std::set<int64_t> &nodeset)
{
  nodeset.clear();
  nodeset.insert(&fNodegraph[fNodegraphindex[node]],&fNodegraph[fNodegraphindex[node+1]]);
  nodeset.insert(node);
}

/**
 * Look for elements formed by vertices, intersecting with the intersectvertices, one by one
 * If the intersection does not remove any of the intersectvertices, we found an element!
 */
void TPZNodesetCompute::AnalyseForElements(std::set<int64_t> &vertices, std::set< std::set<int64_t> > &elements)
{

  if(!vertices.size()) return;
  std::set<int64_t> elem;
  std::set<int64_t>::iterator intit,diffit;
#ifdef PZ_LOG
    if(logger.isDebugEnabled())
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " Original set of nodes ";
		Print(sout,vertices,0);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  for(intit = vertices.begin(); intit != vertices.end(); intit++)
  {
    std::set<int64_t> locset,diffset,interset,unionset,loclocset;
    // locset = all nodes connected to a vertex
    BuildNodeSet(*intit,locset);
    // only diffset with nodes which have no connection with lower nodes interest us
    if(*locset.begin() < *vertices.begin()) continue;
    // diffset contains all vertices in the mesh, except for those in the influence zone of intit????
		set_difference(vertices.begin(),vertices.end(),locset.begin(),locset.end(),inserter(diffset,diffset.begin()));
    // the influence zone of the vertex includes other vertices
    if(diffset.size())
    {
#ifdef PZ_LOG
		{
			std::stringstream sout;
			sout << "Difference after taking the intersection with " << *intit;
			Print(sout,diffset," Difference set");
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
    // some unions need to be made before calling this method
      for(diffit=diffset.begin(); diffit!= diffset.end(); diffit++) 
      {
      
        loclocset.clear();
        // locset now will contain the influence zone of each vertex node in diffset
        BuildNodeSet(*diffit,loclocset);
        if(*loclocset.begin() < *vertices.begin()) continue;
        unionset.insert(loclocset.begin(),loclocset.end());
      }
      // locset now contains the union of all influence zones of vertices influenced by intit
      diffset.clear();
      // diffset will now contain only vertex nodes
      set_intersection(unionset.begin(),unionset.end(),vertices.begin(),vertices.end(),inserter(diffset,diffset.begin()));
#ifdef PZ_LOG
		{
			std::stringstream sout;
			Print(sout,diffset,"First set to be reanalised");
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
      set_intersection(vertices.begin(),vertices.end(),locset.begin(),locset.end(),inserter(interset,interset.begin()));
#ifdef PZ_LOG
		{
			std::stringstream sout;
			Print(sout,interset,"Second set to be reanalised");
			LOGPZ_DEBUG(logger,sout.str())
		}
#endif
      AnalyseForElements(diffset,elements);
      AnalyseForElements(interset,elements);
      return;
    }
  }
  
  intit = vertices.begin();
  BuildNodeSet(*intit,elem);
  intit++;
//  elem = vertices;
  // this code doesnt make sense!!!
  for(;intit != vertices.end(); intit++)
  {
    std::set<int64_t> locset,interset;
    BuildNodeSet(*intit,locset);
    set_intersection(elem.begin(),elem.end(),locset.begin(),locset.end(),inserter(interset,interset.begin()));
    elem = interset;
  }
  if(vertices != elem)
  {
#ifdef PZ_LOG
	  {
		  std::stringstream sout;
		  sout << "Discarding a vertex set as incomplete";
		  Print(sout,vertices,0);
		  LOGPZ_DEBUG(logger,sout.str())
	  }
#endif
  }
  else if(elem.size())
  {
#ifdef PZ_LOG
	  {
		  std::stringstream sout;
		  Print(sout,elem,"Inserted element");
		  LOGPZ_DEBUG(logger,sout.str())
	  }
#endif
    elements.insert(elem);
  }
}

void TPZNodesetCompute::BuildElementGraph(TPZStack<int64_t> &blockgraph, TPZStack<int64_t> &blockgraphindex)
{
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__ << " entering build element graph\n";
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  blockgraph.Resize(0);
  blockgraphindex.Resize(1);
  blockgraphindex[0] = 0;
  int64_t in;
  int64_t nnodes = this->fNodegraphindex.NElements()-1;
  std::set<int64_t> nodeset;
  for(in=0; in<nnodes; in++)
  {
    std::set< std::set<int64_t> > elements;
    BuildNodeSet(in,nodeset);
#ifdef PZ_LOG
	  {
		  std::stringstream sout;
		  sout << "Nodeset for " << in << ' ';
		  Print(sout,nodeset,"Nodeset");
		  LOGPZ_DEBUG(logger,sout.str())
	  }
#endif
    SubtractLowerNodes(in,nodeset);
#ifdef PZ_LOG
		  {
			  std::stringstream sout;
			  Print(sout,nodeset,"LowerNodes result");
			  LOGPZ_DEBUG(logger,sout.str())
		  }
#endif
    AnalyseForElements(nodeset,elements);
    std::set< std::set<int64_t> >::iterator itel;
    for(itel = elements.begin(); itel != elements.end(); itel++)
    {
      std::set<int64_t>::const_iterator it;
      for(it = (*itel).begin(); it!= (*itel).end(); it++)
      {
        blockgraph.Push(*it);
      }
      blockgraphindex.Push(blockgraph.NElements());
    }
  }
}

void TPZNodesetCompute::SubtractLowerNodes(int64_t node, std::set<int64_t> &nodeset)
{
  std::set<int64_t> lownode,lownodeset,unionset;
  std::set<int64_t>::iterator it;
#ifdef PZ_LOG
	{
		std::stringstream sout;
		sout << __PRETTY_FUNCTION__;
		Print(sout,nodeset," Incoming nodeset");
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  for(it=nodeset.begin(); it != nodeset.end() && *it < node; it++)
  {
    BuildNodeSet(*it,lownodeset);
    unionset.insert(lownodeset.begin(),lownodeset.end());
    lownodeset.clear();
  }
  set_difference(nodeset.begin(),nodeset.end(),unionset.begin(),unionset.end(),
    inserter(lownode,lownode.begin()));
#ifdef PZ_LOG
	{
		std::stringstream sout;
		Print(sout,lownode," What is left after substracting the influence of lower numbered nodes ");
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  unionset.clear();
  for(it=lownode.begin(); it!=lownode.end(); it++)
  {
    BuildNodeSet(*it,lownodeset);
    unionset.insert(lownodeset.begin(),lownodeset.end());
    lownodeset.clear();
  }
  lownode.clear();
  set_intersection(unionset.begin(),unionset.end(),nodeset.begin(),nodeset.end(),
    inserter(lownode,lownode.begin()));
#ifdef PZ_LOG
	{
		std::stringstream sout;
		Print(sout,lownode," Resulting lower nodeset");
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
  nodeset = lownode;
}

void TPZNodesetCompute::Print(std::ostream &file) const
{
  file << "TPZNodesetCompute\n";
  file << "Node graph\n";
  Print(file,fNodegraphindex,fNodegraph);
  file << "Node sequence number and level\n";
  int64_t nnode = fNodegraphindex.NElements()-1;
  int64_t in;
  for(in=0; in<nnode; in++) file << in << "/" << fSeqNumber[in] << "/" << fLevel[in] << " ";
  file << std::endl;
}

void TPZNodesetCompute::Print(std::ostream &file, const TPZVec<int64_t> &graphindex, const TPZVec<int64_t> &graph)
{
  int64_t nnode = graphindex.NElements()-1;
  int64_t in;
  for(in =0; in<nnode; in++)
  {
    file << "Node Number " << in;
    int64_t first = graphindex[in];
    int64_t last = graphindex[in+1];
    int64_t jn;
    for(jn = first; jn< last; jn++) file << " " << graph[jn];
    file << std::endl;
  }
}

void TPZNodesetCompute::Print(std::ostream &file, const std::set<int64_t> &nodeset, const char *text)
{
  if(text) file << text;
  std::set<int64_t>::const_iterator it;
  for(it=nodeset.begin(); it!=nodeset.end(); it++) file << *it << ' ';
  file << std::endl;
}

  /**
  * Expand the graph acording to the block structure
  */
void TPZNodesetCompute::ExpandGraph(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, TPZBlock &block,
    TPZVec<int64_t> &expgraph, TPZVec<int64_t> &expgraphindex)
{
  int64_t expgraphsize = 0;
  int64_t nbl = graph.NElements();
  expgraphindex.Resize(graphindex.NElements());
  int64_t ibl;
  for(ibl=0; ibl<nbl; ibl++)
  {
    expgraphsize += block.Size(graph[ibl]);
  }
  expgraph.Resize(expgraphsize);
  int64_t counter = 0;
  int64_t numblocks = graphindex.NElements()-1;
  expgraphindex[0] = counter;
  int64_t blcounter = 0;
  for(ibl=0; ibl < numblocks; ibl++)
  {
    int64_t first = graphindex[ibl];
    int64_t last = graphindex[ibl+1];
    int64_t ieq;
    for(ieq = first; ieq<last; ieq++)
    {
      int64_t blsize =  block.Size(graph[ieq]);
      int64_t pos = block.Position(graph[ieq]);
      int64_t b;
      for(b=0; b<blsize; b++)
      {
        expgraph[counter++] =  pos+b;
      }
    }
    if(expgraphindex[blcounter] != counter) expgraphindex[++blcounter] = counter;
  }
  expgraphindex.Resize(blcounter+1);
}

  /**
  * Color the graph into mutually independent blocks
  */
int TPZNodesetCompute::ColorGraph(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, int64_t neq,
    TPZVec<int> &colors)
{
  TPZVec<int> eqcolor(neq);
  int color = 0;
  bool hasuncolored = true;
  int64_t nblocks = graphindex.NElements()-1;
  colors.Resize(nblocks);
  colors.Fill(-1);
  while(hasuncolored)
  {
    hasuncolored = false;
    eqcolor.Fill(-1);
    int64_t ibl;
    for(ibl=0; ibl<nblocks; ibl++)
    {
      if(colors[ibl] != -1) continue;
      int64_t first = graphindex[ibl];
      int64_t last = graphindex[ibl+1];
      int64_t ieq;
      for(ieq=first; ieq<last; ieq++)    
      {
        if(eqcolor[graph[ieq]] == color) break;
      }
      if(ieq != last)
      {
        hasuncolored = true;
      }
      else
      {
        colors[ibl] = color;
        for(ieq=first; ieq<last; ieq++)    
        {
          eqcolor[graph[ieq]] = color;
        }
      }
    }
    color++;
  }
  return color;
}
