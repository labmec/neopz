/**
 * @file
 * @brief Contains the implementation of the TPZRenumbering methods. 
 */

#include "TPZRenumbering.h"
#include "pzvec.h"
#include "pzerror.h"
#include "pzstack.h"
#include <map>
#include <set>
#include <algorithm>
#include "pzlog.h"
#include <algorithm>

#include "TPZTimer.h"
#include "Hash/TPZHash.h"
#include "TPZStream.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "TPZThreadPool.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.renumbering"));
#endif

using namespace std;

void TPZRenumbering::NodeToElGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, TPZVec<int64_t> &nodtoelgraph, TPZVec<int64_t> &nodtoelgraphindex) {
	
	TPZVec<int64_t> nelcon(fNNodes+1,0);
  	int64_t nod,last = elgraphindex[fNElements];
  	for(nod = 0; nod<last; nod++) {
		nelcon[elgraph[nod]]++;
  	}
	nodtoelgraphindex = nelcon;
	
  	for(nod=fNNodes; nod>0; nod--) nodtoelgraphindex[nod] = nodtoelgraphindex[nod-1];
  	nodtoelgraphindex[0] = 0;
  	for(nod=1;nod<=fNNodes;nod++) nodtoelgraphindex[nod] += nodtoelgraphindex[nod-1];
	
	nodtoelgraph.Resize(nodtoelgraphindex[fNNodes]);
	nodtoelgraph.Fill (-1);
	
    TPZVec<int64_t> nodtoelgraphposition(nodtoelgraphindex);
    
	int64_t el;
  	for(el=0; el<fNElements; el++) {
		int64_t firstnode = elgraphindex[el];
		int64_t lastnode = elgraphindex[el+1];
		for(nod=firstnode;nod<lastnode;nod++) {
			int64_t gnod = elgraph[nod];
            nodtoelgraph[nodtoelgraphposition[gnod]] = el;
            nodtoelgraphposition[gnod]++;
#ifdef PZDEBUG
            if (nodtoelgraphposition[gnod] > nodtoelgraphindex[gnod+1]) {
                PZError << __PRETTY_FUNCTION__ << " wrong data structure\n";
            }
#endif
            /*
			int firstel= nodtoelgraphindex[gnod];
			int lastel = nodtoelgraphindex[gnod+1];
			while(firstel<lastel && nodtoelgraph[firstel] != -1) firstel++;
			if(firstel == lastel) {
				PZError << "TPZCompMesh::ComputeConnecttoElGraph wrong data structure\n";
				continue;
			} else {
				nodtoelgraph[firstel] = el;
			}
            */
		}
  	}
}

void TPZRenumbering::ConvertGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex, TPZVec<int64_t> &nodegraph, TPZVec<int64_t> &nodegraphindex) {
    TPZTimer convert("Converting graph ");
    convert.start();
    int64_t nod, el;
    TPZVec<int64_t> nodtoelgraphindex;
    TPZVec<int64_t> nodtoelgraph;

    NodeToElGraph(elgraph, elgraphindex, nodtoelgraph, nodtoelgraphindex);

    nodegraphindex.Resize(fNNodes + 1);
    nodegraphindex.Fill(0);

    int64_t nodegraphincrement = 10000;
    nodegraph.Resize(nodegraphincrement);
    int64_t nextfreeindex = 0;
    for (nod = 0; nod < fNNodes; nod++) {
        int64_t firstel = nodtoelgraphindex[nod];
        int64_t lastel = nodtoelgraphindex[nod + 1];
        std::set<int64_t> nodecon;
        for (el = firstel; el < lastel; el++) {
            int64_t gel = nodtoelgraph[el];
            int64_t firstelnode = elgraphindex[gel];
            int64_t lastelnode = elgraphindex[gel + 1];
            nodecon.insert(&elgraph[firstelnode], &elgraph[(lastelnode - 1)] + 1);
        }
        nodecon.erase(nod);
        while (nextfreeindex + nodecon.size() >= nodegraph.NElements()) nodegraph.Resize(nodegraph.NElements() + nodegraphincrement);
        std::set<int64_t>::iterator it;
        for (it = nodecon.begin(); it != nodecon.end(); it++) nodegraph[nextfreeindex++] = *it;
        nodegraphindex[nod + 1] = nextfreeindex;
    }
    convert.stop();
    //    std::cout << convert.processName().c_str()  << convert << std::endl;
}

TPZRenumbering::TPZRenumbering(int64_t NElements, int64_t NNodes){
	fNElements = NElements;
	fNNodes = NNodes;
}

int TPZRenumbering::ClassId() const {
    return Hash("TPZRenumbering");
}

void TPZRenumbering::Read(TPZStream& buf, void* context) {
    buf.Read(&fHDivPermute);
    buf.Read(&fNElements);
    buf.Read(&fNNodes);
    buf.Read(fNodeWeights);
    buf.Read(fElementGraph);
    buf.Read(fElementGraphIndex);
}

void TPZRenumbering::Write(TPZStream& buf, int withclassid) const {
    buf.Write(&fHDivPermute);
    buf.Write(&fNElements);
    buf.Write(&fNNodes);
    buf.Write(fNodeWeights);
    buf.Write(fElementGraph);
    buf.Write(fElementGraphIndex);
}

int64_t TPZRenumbering::ColorNodes(TPZVec<int64_t> &nodegraph, TPZVec<int64_t> &nodegraphindex, TPZVec<int> &family, TPZVec<int> &colors) {
    TPZStack<int> usedcolors;
    TPZStack<int64_t> ncolorsbyfamily;
    if (nodegraph.NElements() - 1 != family.NElements()) {
        cout << "TPZRenumbering::ColorNodes inconsistent input parameters" << std::endl;
        DebugStop();
    }
    int64_t nnodes = nodegraphindex.NElements() - 1;
    colors.Resize(nnodes);
    colors.Fill(-1);
    int current_family = 0;
    int64_t nodeshandled = 0;
    int64_t ncolors = 0;
    while (nodeshandled < nnodes) {
        int64_t nod;
        current_family = 0;
        usedcolors.Resize(0);
        for (nod = 0; nod < nnodes; nod++) {
            int64_t firstnod = nodegraphindex[nod];
            int64_t lastnod = nodegraphindex[nod + 1];
            usedcolors.Fill(-1);
            for (int64_t ind = firstnod; ind < lastnod; ind++) {
                int64_t nodcon = nodegraph[ind];
                if (family[nodcon] != current_family) continue;
                if (colors[nodcon] != -1) usedcolors[colors[nodcon]] = 1;
            }
            int64_t ic;
            for (ic = 0; ic < usedcolors.NElements(); ic++){
                if (usedcolors[ic] != 1){
                    break;
                }
            }
            if (ic == usedcolors.NElements()) {
                usedcolors.Push(1);
            }
            colors[nod] = ic;
            nodeshandled++;
        }
        ncolorsbyfamily.Push(usedcolors.NElements());
        ncolors += usedcolors.NElements();
        current_family++;
    }
    return ncolors;
}

int64_t TPZRenumbering::ColorElements(const TPZCompMesh *cmesh, const TPZVec<int64_t> &elementIndices, TPZVec<int64_t> &elementColors) {
    const int64_t n_connects = cmesh->NConnects();
    const int64_t nel = cmesh->NElements();
    const int64_t nel_to_be_colored = elementIndices.size();

    if (nel == 0 || nel_to_be_colored == 0) return 0;
    
    elementColors.Resize(nel_to_be_colored);
    elementColors.Fill(-1);
    int64_t n_colors = 0;
    int64_t next_color_initial_index = 0;
    std::function<void(int64_t,int64_t,TPZTaskGroup *)> color = [n_connects, nel_to_be_colored, &elementIndices, &elementColors, cmesh, &color, &n_colors, &next_color_initial_index](int64_t initial_index, int64_t currentColor, TPZTaskGroup *group) {
        // determines whether each connect is in an element which has this color
        TPZStack<int64_t> connects_this_color;
        bool created_new_color = false; // flag: in the case we find an uncolored element 
        // that cannot be colored with this color, we should use another color,
        // but just for the first element we find.
        for (int64_t iel = initial_index; iel < nel_to_be_colored; ++iel) {
            auto elindex = elementIndices[iel];
            // if this element hasn't been computed in a previous pass
            if (elementColors[iel] == -1) {
                auto cel = cmesh->Element(elindex);
                if (!cel) continue;
                TPZStack<int64_t> connectlist;
                cel->BuildConnectList(connectlist);
                bool skip_element = false;
                for (auto connect : connectlist) {
                    if (std::find(connects_this_color.begin(), connects_this_color.end(), connect) != connects_this_color.end()) {
                        skip_element = true;
                        break;
                    }
                }
                if (skip_element){
                    if (!created_new_color) {
                        if (group->Active() < TPZThreadPool::globalInstance().threadCount()){
                            TPZThreadPool::globalInstance().run(1, group, color, iel, n_colors++, group);
                            created_new_color = true;
                        } else {
                            if (next_color_initial_index == -1) {
                                next_color_initial_index = iel;
                            } else {
                                next_color_initial_index = std::min(next_color_initial_index, iel);
                            }
                        }
                    }
                    continue;
                }
                for (auto connect : connectlist) {
                    connects_this_color.push_back(connect);
                }
                elementColors[iel] = currentColor;
            }
        }
    };
    do {
        int64_t initial_index = next_color_initial_index;
        next_color_initial_index = -1;
        TPZTaskGroup group;
        TPZThreadPool::globalInstance().run(1, &group, color, initial_index, n_colors++, &group);
        group.Wait();
    } while (next_color_initial_index != -1);
#ifdef PZDEBUG
#define CHECK_COLORING
#ifdef CHECK_COLORING
    for (unsigned int color = 0; color < n_colors; ++color) {
        TPZManVector<bool> color_connect(n_connects, false);
        for (int64_t iel = 0; iel < nel_to_be_colored; ++iel) {
            auto elindex = elementIndices[iel];
            if (elementColors[iel] == color) {
                auto cel = cmesh->Element(elindex);
                if (!cel) DebugStop(); //colored element is not in the mesh
                TPZStack<int64_t> connectlist;
                cel->BuildConnectList(connectlist);
                for (auto connect : connectlist) {
                    if (color_connect[connect]) { //this connect is in another element of the same color
                        DebugStop();
                    }
                }
                for (auto connect : connectlist) {
                    color_connect[connect] = true;
                }
            } else if (elementColors[iel] == -1){ 
                auto cel = cmesh->Element(elindex);
                if (cel) DebugStop(); //uncolored element
            }
        }
    }
#endif
#endif
    return n_colors;
}

void TPZRenumbering::Print(TPZVec<int64_t> &graph, TPZVec<int64_t> &graphindex, const char *name, std::ostream& out) {
    out << "Graph: " << name << endl;
    for (int64_t i = 0; i < graphindex.NElements() - 1; i++) {
        out << "Graph item: " << i << "\t";
        for (int64_t j = graphindex[i]; j < graphindex[i + 1]; j++) {
            if (j >= graph.NElements()) {
                out << "wrong graphindex graph.NElements = " << graph.NElements() << " i = " << i << "graphindex[i] = " << graphindex[i] << " " << graphindex[i + 1] << endl;
                break;
            } else {
                out << graph[j] << "\t";
            }
        }
        out << endl;
    }
}

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzcompel.h"

static REAL XMin(TPZGeoEl *gel, const TPZVec<REAL> &normal) {
	
	int nnode = gel->NNodes();
	REAL distmin = 1.e10;
	if(!nnode) {
		cout << "Geometric element without nodes??? \n";
		gel->Print(cout);
		return distmin;
	}
	TPZGeoNode *gn = gel->NodePtr(0);
	distmin = gn->Coord(0)*normal[0]+ gn->Coord(1)*normal[1]+ gn->Coord(2)*normal[2];
	int in;
	for(in=0; in<nnode; in++) {
		gn=gel->NodePtr(in);
		REAL dist = gn->Coord(0)*normal[0]+ gn->Coord(1)*normal[1]+ gn->Coord(2)*normal[2];
		distmin = distmin < dist ? distmin : dist;
	}
	return distmin;
	
}
void ResequenceByGeometry(TPZCompMesh *cmesh, const TPZVec<REAL> &normal) {
	
	TPZCompEl *cel;
	multimap<REAL, TPZCompEl *> elmap;
	int64_t nelem = cmesh->ElementVec().NElements();
	int64_t iel;
	for(iel=0; iel<nelem; iel++) {
		cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel) continue;
		REAL dist = XMin(gel,normal);
		elmap.insert(pair<REAL,TPZCompEl *>(dist,cel));
	}
	TPZVec<int64_t> Permute(cmesh->NConnects(),-1);
	multimap<REAL, TPZCompEl *>::iterator it;
	it = elmap.begin();
	int64_t iseq=0;
	while(it != elmap.end()) {
		cel = it->second;
		int nc = cel->NConnects();
		int ic;
		for(ic=0; ic<nc; ic++) {
			int64_t seqnum = cel->Connect(ic).SequenceNumber();
			if(Permute[seqnum] == -1) Permute[seqnum] = iseq++;
		}
		it++;
	}
	cmesh->Permute(Permute);
	cmesh->CleanUpUnconnectedNodes();
	
}

/**
 * Convert a traditional elgraph to an element to element graph
 */
void TPZRenumbering::ConvertToElementoToElementGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex,
													 TPZVec<int64_t> &eltoelgraph, TPZVec<int> &eltoelweight, TPZVec<int64_t> &eltoelgraphindex)
{
	TPZVec<int64_t> nodegraph;
	TPZVec<int64_t> nodegraphindex;
	LOGPZ_DEBUG(logger, "before NodeToElGraph")
	NodeToElGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
	LOGPZ_DEBUG(logger, "after NodeToElGraph")
	int64_t nelements = elgraphindex.NElements()-1;
	eltoelgraphindex.Resize(nelements+1);
	eltoelgraphindex[0] = 0;
	eltoelgraph.Resize(1000);
	eltoelweight.Resize(1000);
	int64_t iel;
	for(iel=0; iel<nelements; iel++)
	{
		map<int64_t,int> elset;
		int64_t firstnodeindex = elgraphindex[iel];
		int64_t lastnodeindex = elgraphindex[iel+1];
		int64_t nodeindex;
		for(nodeindex = firstnodeindex; nodeindex< lastnodeindex; nodeindex++)
		{
			int64_t node = elgraph[nodeindex];
			int64_t firstelindex = nodegraphindex[node];
			int64_t lastelindex = nodegraphindex[node+1];
			int64_t elindex;
			for(elindex = firstelindex; elindex < lastelindex; elindex++)
			{
				int64_t element = nodegraph[elindex];
				if(element != iel)
				{
					int nweight = 0;
					if(node < fNodeWeights.NElements()) nweight = fNodeWeights[node];
					elset[element] += nweight;
				}
			}
		}
		int64_t eltoelsize = eltoelgraph.NElements();
		if(eltoelgraphindex[iel]+elset.size() >= eltoelsize)
		{
			eltoelgraph.Resize(eltoelsize+elset.size()+1000);
			eltoelweight.Resize(eltoelgraph.NElements());
		}
		int64_t count = eltoelgraphindex[iel];
		map<int64_t,int>::iterator it;
		for(it=elset.begin(); it != elset.end(); it++,count++)
		{
			eltoelgraph[count] = it->first;
			eltoelweight[count] = it->second;
		}
		eltoelgraphindex[iel+1] = count;
	}
	eltoelgraph.Resize(eltoelgraphindex[nelements]);
	eltoelweight.Resize(eltoelgraph.NElements());
}

void TPZRenumbering::SetElementGraph(TPZVec<int64_t> &elgraph, TPZVec<int64_t> &elgraphindex){
#ifdef SLOANDEBUG
	Print(elgraph, elgraphindex, "original element graph", cout);
#endif
	fElementGraph = elgraph;
	fElementGraphIndex = elgraphindex;
}

void TPZRenumbering::ClearDataStructures(){
	fNodeWeights.Resize(0);
	fElementGraph.Resize(0);
	fElementGraphIndex.Resize(0);
}

/**
 * Analyse the graph, find the corner nodes
 */
void TPZRenumbering::CornerEqs(unsigned int mincorners, int64_t nelconsider, std::set<int> &eligible, std::set<int> &cornernodes)
{
	
	TPZVec<int64_t> nodtoelgraphindex;
	TPZVec<int64_t> nodtoelgraph;
	int64_t sub = 0;
	
	NodeToElGraph(fElementGraph,fElementGraphIndex,nodtoelgraph,nodtoelgraphindex);
	
	int64_t nelem = fElementGraphIndex.NElements()-1;
    if (nelconsider > nelem) {
        DebugStop();
    }
	int64_t element;
	for (element=0; element < nelconsider; element++) 
	{
		int64_t firstindex = fElementGraphIndex[element];
		int64_t lastindex = fElementGraphIndex[element+1];
		if (firstindex == lastindex) {
			continue;
		}
		TPZStack<int64_t> corners;
        std::set<int64_t> elcornernodes;
		// a vector of sets of element connections for each node
        // first key : number elements connected to the node
        // second key : pair of node and set of element indices
		std::multimap<int64_t,std::pair<int64_t,std::set<int64_t> > > connectivities;
		typedef std::multimap<int64_t,std::pair<int64_t,std::set<int64_t> > > map_type;
		int64_t ind;
		int64_t maxelcon = 0;
		for (ind=firstindex; ind<lastindex; ind++) {
			int64_t node = fElementGraph[ind];
			int64_t firstelind = nodtoelgraphindex[node];
			int64_t lastelind = nodtoelgraphindex[node+1];
			std::set<int64_t> elcon;
			elcon.insert(&nodtoelgraph[firstelind],(&nodtoelgraph[lastelind-1])+1);
			maxelcon = maxelcon < elcon.size() ? elcon.size() : maxelcon;
			connectivities.insert(map_type::value_type(elcon.size(), std::pair<int64_t, std::set<int64_t> >(node,elcon)));
		}
		
		map_type::reverse_iterator it = connectivities.rbegin();		
		// put all nodes with maximum connectivities on the stack
		int64_t maxconnect = it->first;
		std::pair<map_type::const_iterator, map_type::const_iterator> p = connectivities.equal_range(it->first);
		map_type::const_iterator ithead;
		for (ithead = p.first; ithead != p.second; ithead++) 
		{
            int64_t seqnum = ithead->second.first;
            if (eligible.find(seqnum) != eligible.end()) {
                corners.Push(seqnum);
                cornernodes.insert(seqnum);
                elcornernodes.insert(seqnum);
            }
		}
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Element " << element << " First stage corner nodes " << corners;
			if (logger->isDebugEnabled())
			{
				LOGPZ_DEBUG(logger, sout.str())
			}
        }
#endif
		// look the included sets
		map_type::iterator itf = connectivities.begin();		
		for (itf=connectivities.begin(); itf != connectivities.end(); itf++)
		{
			map_type::reverse_iterator it = connectivities.rbegin();
			while (it->first > itf->first && it!=connectivities.rend()) 
			{
				std::set<int64_t>::iterator smallsetbeg, smallsetend, largesetbeg, largesetend;
				smallsetbeg = itf->second.second.begin();
				smallsetend = itf->second.second.end();
				largesetbeg = it->second.second.begin();
				largesetend = it->second.second.end();
				
				if (includes(largesetbeg,largesetend,smallsetbeg,smallsetend)) {
					break;
				}
				it++;
			}
			// the set is not included in any
			// we allready put the maxconnect nodes on the stack
			if (it->first == itf->first && it->first != maxconnect) {
                int64_t seqnum = itf->second.first;
                if (eligible.find(seqnum) != eligible.end()) 
                {
                    corners.Push(seqnum);
                    cornernodes.insert(seqnum);
                    elcornernodes.insert(seqnum);
                }
			}
		}
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Element " << element << " cornernodes.size " << elcornernodes.size() << " Second stage corner nodes " << corners;
			if (logger->isDebugEnabled())
			{
				LOGPZ_DEBUG(logger, sout.str())
			}
        }
#endif
		if (elcornernodes.size() < mincorners) {
			it = connectivities.rbegin();
			int64_t numelconnected = it->first-1;
			int64_t ncorners = corners.NElements();
			while (numelconnected > 1 && ncorners < mincorners) {
				std::pair<map_type::const_iterator, map_type::const_iterator> p = connectivities.equal_range(numelconnected);
				map_type::const_iterator ithead;
				for (ithead = p.first; ithead != p.second; ithead++) 
				{
                    int64_t seqnum = ithead->second.first;
                    if (eligible.find(seqnum) != eligible.end()) 
                    {
                        corners.Push(seqnum);
                        cornernodes.insert(seqnum);
                        elcornernodes.insert(seqnum);
                    }
					ncorners = elcornernodes.size();
				}
				numelconnected--;
			}
		}
		int64_t ieq;
		for (ieq=0; ieq<corners.NElements(); ieq++) {
			if (cornernodes.find(corners[ieq]) == cornernodes.end()) {
				DebugStop();
			}
		}
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Element " << element << " sub " << sub << " Final stage corner nodes " << corners;
			if (logger->isDebugEnabled())
			{
				LOGPZ_DEBUG(logger, sout.str())
			}
        }
#endif
		
		sub++;
	}
}

template class TPZRestoreClass<TPZRenumbering>;
