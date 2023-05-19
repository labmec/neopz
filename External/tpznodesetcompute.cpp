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
#include "TPZEquationFilter.h"
#include <set>
#include <map>
#include <algorithm>
#include <iterator>
#include <fstream>
#include "pzlog.h"
#include "TPZSimpleTimer.h"
#include "pzvec_extras.h"
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
  TPZSimpleTimer analyse("AnalyseGraph");
  fMaxSeqNum = -1;
  fMaxLevel = -1;
  int64_t nnodes = fNodegraphindex.NElements()-1;
  fSeqNumber.Resize(nnodes);
  this->fLevel.Resize(nnodes);
  fSeqNumber.Fill(-1);
	fIsIncluded.Resize(nnodes);
	fIsIncluded.Fill(0);
  fLevel.Fill(0);
  int64_t in;
  int count{0};//for displaying info
  const char* spinner = "-\\|/";
  const int spin_len = 4;
  for(in=0; in<nnodes; in++) 
  {
	  if(!(in%1000))
	  {
		  std::cout <<'\r'<<"TPZNodesetCompute::AnalyseGraph "<<spinner[count++%spin_len];
		  std::cout.flush();
	  }
	  AnalyseNode(in);
  }
  if(count > 0){
    std::cout << '\r'<< "TPZNodesetCompute::AnalyseGraph done!"<<std::endl;
  }
  
}

/**
 * This method will analyse the set inclusion of the current node, calling the method
 * recursively if another node need to be analysed first
 */
void TPZNodesetCompute::AnalyseNode(const int64_t node)
{
  if(fSeqNumber[node] != -1) return;

  int minlevel = 0;
  TPZStack<int64_t> equalnodes;

  const auto mysz =  fNodegraphindex[node+1] - fNodegraphindex[node];
  TPZManVector<int64_t,1000> mynodeset(mysz+1);
  {
    const auto first = fNodegraphindex[node];
    for(auto i = 0; i < mysz; i++){
      mynodeset[i] = fNodegraph[first+i];
    }
    mynodeset[mysz] = node;
    std::sort(mynodeset.begin(),mynodeset.end());
  }
  
  for(auto othernode : mynodeset)
  {
    if(othernode == node) continue;
    // build the data structure of the connected node
    const auto othersz =  fNodegraphindex[othernode+1] - fNodegraphindex[othernode];
    TPZManVector<int64_t,1000> othernodeset(othersz+1);
    {
      const auto first = fNodegraphindex[othernode];
      for(auto i = 0; i < othersz; i++){
        othernodeset[i] = fNodegraph[first+i];
      }
      othernodeset[othersz] = othernode;
      std::sort(othernodeset.begin(),othernodeset.end());
    }
    
    
    // connectivity of the other node is included in my connectivity
    const bool inc = std::includes(mynodeset.begin(),mynodeset.end(), othernodeset.begin(), othernodeset.end());
    // the other node has a different connectivity
    const bool diff = !std::equal(mynodeset.begin(),mynodeset.end(), othernodeset.begin(), othernodeset.end());
    if(inc && diff) 
    {
      // Analyse the connectivity graph of the other node
      // Its graph is smaller than mine
      fIsIncluded[othernode] = 1;
      AnalyseNode(othernode);
#ifdef PZ_LOG
      if(fLevel[othernode] >= minlevel)
      {
        if(logger.isDebugEnabled()){
          std::stringstream sout;
          sout << "The level of " << node << " is increased because of " << othernode << " ";
          sout << "Level of othernode " << fLevel[othernode] << " Seqnumber othernode " << fSeqNumber[othernode];
          LOGPZ_DEBUG(logger,sout.str());
        }
      }
#endif
      minlevel = minlevel < fLevel[othernode]+1 ? fLevel[othernode]+1 : minlevel;
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
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      sout << "Assigning Seq Number " << fMaxSeqNum << " and level " << minlevel << " to nodes " << node << " " << equalnodes;
      LOGPZ_DEBUG(logger,sout.str());
    }
	}
#endif
  fSeqCard.Push(1);
  // the maximum level of the set of nodes
  fMaxLevel = fMaxLevel < minlevel ? minlevel : fMaxLevel;
  // assign the level of the node
  fLevel[node] = minlevel;
	
  // initialize the datastructure of the nodes which have the same connectivity
  int64_t neq = equalnodes.NElements();
  int64_t ieq;
  for(ieq=0; ieq<neq; ieq++)
  {
    fSeqNumber[equalnodes[ieq]] = fMaxSeqNum;
    fLevel[equalnodes[ieq]] = minlevel;
    // the number of nodes associated with this sequence number
    fSeqCard[fMaxSeqNum]++;
  }
}

void TPZNodesetCompute::BuildNodeGraph(TPZVec<int64_t> &blockgraph, TPZVec<int64_t> &blockgraphindex)
{
  TPZSimpleTimer timer("TPZNodesetCompute::BuildNodeGraph");
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
  TPZSimpleTimer timer("TPZNodesetCompute::BuildVertexGraph");
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
  const auto lastpos = fNodegraphindex[node+1];
  nodeset.insert(&fNodegraph[fNodegraphindex[node]],&fNodegraph[0]+lastpos);
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
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      sout << __PRETTY_FUNCTION__ << " Original set of nodes ";
      Print(sout,vertices,0);
      LOGPZ_DEBUG(logger,sout.str());
    }
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
			if(logger.isDebugEnabled()){
        std::stringstream sout;
        sout << "Difference after taking the intersection with " << *intit;
        Print(sout,diffset," Difference set");
        LOGPZ_DEBUG(logger,sout.str());
      }
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
			if(logger.isDebugEnabled()){
        std::stringstream sout;
        Print(sout,diffset,"First set to be reanalised");
        LOGPZ_DEBUG(logger,sout.str());
      }
		}
#endif
      set_intersection(vertices.begin(),vertices.end(),locset.begin(),locset.end(),inserter(interset,interset.begin()));
#ifdef PZ_LOG
		{
			if(logger.isDebugEnabled()){
        std::stringstream sout;
        Print(sout,interset,"Second set to be reanalised");
        LOGPZ_DEBUG(logger,sout.str());
      }
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
		  if(logger.isDebugEnabled()){
        std::stringstream sout;
        sout << "Discarding a vertex set as incomplete";
        Print(sout,vertices,0);
        LOGPZ_DEBUG(logger,sout.str());
      }
	  }
#endif
  }
  else if(elem.size())
  {
#ifdef PZ_LOG
	  {
		  if(logger.isDebugEnabled()){
        std::stringstream sout;
        Print(sout,elem,"Inserted element");
        LOGPZ_DEBUG(logger,sout.str());
      }
	  }
#endif
    elements.insert(elem);
  }
}

void TPZNodesetCompute::BuildElementGraph(TPZStack<int64_t> &blockgraph, TPZStack<int64_t> &blockgraphindex)
{
  TPZSimpleTimer timer("TPZNodesetCompute::BuildElementGraph");
#ifdef PZ_LOG
	{
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      sout << __PRETTY_FUNCTION__ << " entering build element graph\n";
      LOGPZ_DEBUG(logger,sout.str());
    }
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
		  if(logger.isDebugEnabled()){
        std::stringstream sout;
        sout << "Nodeset for " << in << ' ';
        Print(sout,nodeset,"Nodeset");
        LOGPZ_DEBUG(logger,sout.str());
      }
	  }
#endif
    SubtractLowerNodes(in,nodeset);
#ifdef PZ_LOG
    {
      if(logger.isDebugEnabled()){
        std::stringstream sout;
        Print(sout,nodeset,"LowerNodes result");
        LOGPZ_DEBUG(logger,sout.str());
      }
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
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      sout << __PRETTY_FUNCTION__;
      Print(sout,nodeset," Incoming nodeset");
      LOGPZ_DEBUG(logger,sout.str());
    }
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
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      Print(sout,lownode," What is left after substracting the influence of lower numbered nodes ");
      LOGPZ_DEBUG(logger,sout.str());
    }
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
		if(logger.isDebugEnabled()){
      std::stringstream sout;
      Print(sout,lownode," Resulting lower nodeset");
      LOGPZ_DEBUG(logger,sout.str());
    }
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
    file << "Node Number " << in << ':' ;
    int64_t first = graphindex[in];
    int64_t last = graphindex[in+1];
    int64_t jn;
    for(jn = first; jn< last; jn++) file << " " << graph[jn];
    file << std::endl;
  }
}

void TPZNodesetCompute::Print(std::ostream &file,
                              const TPZVec<int64_t> &graphindex,
                              const TPZVec<int64_t> &graph,
                              const TPZVec<int> &color)
{
  int64_t nnode = graphindex.NElements()-1;
  int64_t in;
  for(in =0; in<nnode; in++)
  {
    int64_t first = graphindex[in];
    int64_t last = graphindex[in+1];
    int64_t jn;
    file << "Node Number " << in << "(size "<<last-first<<" color "<<color[in]<<" ) :" ;
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
void TPZNodesetCompute::ExpandGraph(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex,
                                    TPZBlock &block, TPZVec<int64_t> &expgraph,
                                    TPZVec<int64_t> &expgraphindex, TPZVec<int64_t> &removed_blocks)
{
  TPZSimpleTimer timer("TPZNodesetCompute::ExpandGraph");
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
  removed_blocks.Resize(0);
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
    if(expgraphindex[blcounter] != counter) {
      const auto first = expgraphindex[blcounter];
      std::sort(&expgraph[0]+first,&expgraph[0]+counter);
      expgraphindex[++blcounter] = counter;
    }else{
      removed_blocks.push_back(ibl);
    }
  }
  expgraphindex.Resize(blcounter+1);
  if (counter != expgraphsize){
    DebugStop();
  }
}

void TPZNodesetCompute::FilterGraph(const TPZEquationFilter &eqfilt, TPZVec<int64_t> &graph,
                                    TPZVec<int64_t> &graphindex, TPZVec<int64_t> &removed_blocks)
{

  /*
    we need to make a copy of graphindex
    since we update graphindex[ibl+1] when
    analysing block ibl
  */
  const auto ebgindex = graphindex;
  const auto nblocks = ebgindex.size()-1;
  
  removed_blocks.Resize(0);
  int64_t indexcount = 0, blcount = 0;
  for(auto ibl = 0; ibl < nblocks; ibl++){
        
    /*
      the prefix o will stand for old (before eq filter),
      and n for new (after eq filter).
    */
    const int ofirst = ebgindex[ibl];
    const int olast = ebgindex[ibl+1];
    const int oneqbl = olast - ofirst;
    TPZVec<int64_t> orig(oneqbl), filt(oneqbl);
    for(auto i = 0; i < oneqbl; i++){
      orig[i]=filt[i]=ebg[ofirst+i];
    }
    eqfilt.Filter(orig,filt);
    const int nneqbl = filt.size();
    //check for empty block
    if(nneqbl){
      newgraphindex[blcount] = indexcount;
      for(int ieq = 0; ieq < nneqbl; ieq++){
        newgraph[indexcount+ieq] = filt[ieq];
      }
      indexcount += nneqbl;
      newgraphindex[blcount+1] = indexcount;
      blcount++;
    }
  }
  auto lasteq = newgraphindex[blcount];
  newgraphindex.Resize(blcount+1);
  newgraph.Resize(lasteq);
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
    nodecolor.Fill(-1);
    int64_t ibl;
    //we iterate over the blocks of the input graph
    for(ibl=0; ibl<nblocks; ibl++)
    {
      if(colors[ibl] != -1) continue;
      int64_t first = graphindex[ibl];
      int64_t last = graphindex[ibl+1];
      int64_t inode;
      //we look over the nodes of the graph block
      for(inode=first; inode<last; inode++)    
      {
        auto node = graph[inode];
        if(nodecolor[node] == color) break;
        //now we check the connected nodes
        const int fneigh = fNodegraphindex[node];
        const int lneigh = fNodegraphindex[node+1];
        int in;
        for(in = fneigh; in < lneigh; in++){
          auto neigh = fNodegraph[in];
          if(nodecolor[neigh] == color) break;
        }
        if(in != lneigh){break;}
      }
      if(inode != last)
      {
        hasuncolored = true;
      }
      else
      {
        colors[ibl] = color;
        for(inode=first; inode<last; inode++)    
        {
          nodecolor[graph[inode]] = color;
        }
      }
    }
    color++;
  }
  return color;
}
