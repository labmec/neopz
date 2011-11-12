/**
 * @file
 * @brief Contains the implementation of the TPZRenumbering methods. 
 */

#include "pzrenumbering.h"
#include "pzvec.h"
#include "pzerror.h"
#include "pzstack.h"
#include <map>
#include <set>
#include <algorithm>
#include "pzlog.h"
#include <algorithm>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.renumbering"));
#endif

using namespace std;

void TPZRenumbering::NodeToElGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex, TPZVec<int> &nodtoelgraph, TPZVec<int> &nodtoelgraphindex){
	
	TPZVec<int> nelcon(fNNodes+1,0);
  	int nod,last = elgraphindex[fNElements];
  	for(nod = 0; nod<last; nod++) {
		nelcon[elgraph[nod]]++;
  	}
	nodtoelgraphindex = nelcon;
	
  	for(nod=fNNodes; nod>0; nod--) nodtoelgraphindex[nod] = nodtoelgraphindex[nod-1];
  	nodtoelgraphindex[0] = 0;
  	for(nod=1;nod<=fNNodes;nod++) nodtoelgraphindex[nod] += nodtoelgraphindex[nod-1];
	
	nodtoelgraph.Resize(nodtoelgraphindex[fNNodes]);
	nodtoelgraph.Fill (-1);
	
	int el;
  	for(el=0; el<fNElements; el++) {
		int firstnode = elgraphindex[el];
		int lastnode = elgraphindex[el+1];
		for(nod=firstnode;nod<lastnode;nod++) {
			int gnod = elgraph[nod];
			int firstel= nodtoelgraphindex[gnod];
			int lastel = nodtoelgraphindex[gnod+1];
			while(firstel<lastel && nodtoelgraph[firstel] != -1) firstel++;
			if(firstel == lastel) {
				PZError << "TPZCompMesh::ComputeConnecttoElGraph wrong data structure\n";
				continue;
			} else {
				nodtoelgraph[firstel] = el;
			}
		}
  	}
}

void TPZRenumbering::ConvertGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex, TPZVec<int> &nodegraph, TPZVec<int> &nodegraphindex){
	int nod,el;
	TPZVec<int> nodtoelgraphindex;
	TPZVec<int> nodtoelgraph;
	
	NodeToElGraph(elgraph,elgraphindex,nodtoelgraph,nodtoelgraphindex);
	
	nodegraphindex.Resize(fNNodes+1);
  	nodegraphindex.Fill(0);
	
	int nodegraphincrement = 10000;
  	nodegraph.Resize(nodegraphincrement);
  	int nextfreeindex = 0;
  	for(nod=0; nod<fNNodes; nod++) {
		int firstel = nodtoelgraphindex[nod];
		int lastel = nodtoelgraphindex[nod+1];
		std::set<int> nodecon;
		for(el=firstel; el<lastel; el++) {
			int gel = nodtoelgraph[el];
			int firstelnode = elgraphindex[gel];
			int lastelnode = elgraphindex[gel+1];
            nodecon.insert(&elgraph[firstelnode],&elgraph[(lastelnode-1)]+1);
		}
        nodecon.erase(nod);
        while(nextfreeindex+(int)nodecon.size() >= nodegraph.NElements()) nodegraph.Resize(nodegraph.NElements()+nodegraphincrement);
        std::set<int>::iterator it;
        for(it = nodecon.begin(); it!= nodecon.end(); it++) nodegraph[nextfreeindex++] = *it;
		nodegraphindex[nod+1] = nextfreeindex;
  	}
}

TPZRenumbering::TPZRenumbering(int NElements, int NNodes){
	fNElements = NElements;
	fNNodes = NNodes;
}

int TPZRenumbering::ColorNodes(TPZVec<int> &nodegraph, TPZVec<int> &nodegraphindex, TPZVec<int> &family, TPZVec<int> &colors) {
	
	TPZStack<int> usedcolors;
	TPZStack<int> ncolorsbyfamily;
	if(nodegraph.NElements()-1 != family.NElements()) {
		cout << "TPZRenumbering::ColorNodes inconsistent input parameters\n";
	}
	int nnodes = nodegraphindex.NElements()-1;
	colors.Resize(nnodes);
	colors.Fill(-1);
	int curfam = 0;
	int nodeshandled = 0;
	int ncolors = 0;
	while(nodeshandled < nnodes) {
		int nod;
		curfam = 0;
		usedcolors.Resize(0);
		for(nod = 0; nod < nnodes; nod++) {
			int firstnod = nodegraphindex[nod];
			int lastnod = nodegraphindex[nod+1];
			usedcolors.Fill(-1);
			int ind, nodcon;
			for(ind= firstnod; ind<lastnod; ind++) {
				nodcon = nodegraph[ind];
				if(family[nodcon] != curfam) continue;
				if(colors[nodcon] != -1) usedcolors[colors[nodcon]] = 1;
			}
			int ic;
			for(ic=0; ic<usedcolors.NElements(); ic++) if(usedcolors[ic] != 1) break;
			if(ic == usedcolors.NElements()) usedcolors.Push(1);
			colors[nod] = ic;
			nodeshandled++;
		}
		ncolorsbyfamily.Push(usedcolors.NElements());
		ncolors += usedcolors.NElements();
		curfam++;
	}
	return ncolors;
}

void TPZRenumbering::Print(TPZVec<int> &grapho, TPZVec<int> &graphoindex, const char *name, std::ostream& out){
	
	int i,j;
	out << "Grapho: " << name << endl;
	for (i=0;i<graphoindex.NElements()-1;i++){
		out << "Grapho item: " << i << "\t";
		for(j=graphoindex[i];j<graphoindex[i+1];j++){
			if(j >= grapho.NElements()) {
				out << "graphoindex errado grapho.NElements = " << grapho.NElements() << " i = " << i << "graphoindex[i] = " << graphoindex[i] << " " << graphoindex[i+1] << endl;
				break;
			} else {
				out << grapho[j] <<"\t";
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
	int nelem = cmesh->ElementVec().NElements();
	int iel;
	for(iel=0; iel<nelem; iel++) {
		cel = cmesh->ElementVec()[iel];
		if(!cel) continue;
		TPZGeoEl *gel = cel->Reference();
		if(!gel) continue;
		REAL dist = XMin(gel,normal);
		elmap.insert(pair<REAL,TPZCompEl *>(dist,cel));
	}
	TPZVec<int> Permute(cmesh->NConnects(),-1);
	multimap<REAL, TPZCompEl *>::iterator it;
	it = elmap.begin();
	int iseq=0;
	while(it != elmap.end()) {
		cel = it->second;
		int nc = cel->NConnects();
		int ic;
		for(ic=0; ic<nc; ic++) {
			int seqnum = cel->Connect(ic).SequenceNumber();
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
void TPZRenumbering::ConvertToElementoToElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex,
													 TPZVec<int> &eltoelgraph, TPZVec<int> &eltoelweight, TPZVec<int> &eltoelgraphindex)
{
	TPZVec<int> nodegraph,nodegraphindex;
	LOGPZ_DEBUG(logger,"before NodeToElGraph")
	NodeToElGraph(elgraph,elgraphindex,nodegraph,nodegraphindex);
	LOGPZ_DEBUG(logger,"after NodeToElGraph")
	int nelements = elgraphindex.NElements()-1;
	eltoelgraphindex.Resize(nelements+1);
	eltoelgraphindex[0] = 0;
	eltoelgraph.Resize(1000);
	eltoelweight.Resize(1000);
	int iel;
	for(iel=0; iel<nelements; iel++)
	{
		map<int,int> elset;
		int firstnodeindex = elgraphindex[iel];
		int lastnodeindex = elgraphindex[iel+1];
		int nodeindex;
		for(nodeindex = firstnodeindex; nodeindex< lastnodeindex; nodeindex++)
		{
			int node = elgraph[nodeindex];
			int firstelindex = nodegraphindex[node];
			int lastelindex = nodegraphindex[node+1];
			int elindex;
			for(elindex = firstelindex; elindex < lastelindex; elindex++)
			{
				int element = nodegraph[elindex];
				if(element != iel)
				{
					int nweight = 0;
					if(node < fNodeWeights.NElements()) nweight = fNodeWeights[node];
					elset[element] += nweight;
				}
			}
		}
		int eltoelsize = eltoelgraph.NElements();
		if(eltoelgraphindex[iel]+(int)elset.size() >= eltoelsize)
		{
			eltoelgraph.Resize(eltoelsize+elset.size()+1000);
			eltoelweight.Resize(eltoelgraph.NElements());
		}
		int count = eltoelgraphindex[iel];
		map<int,int>::iterator it;
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

void TPZRenumbering::SetElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex){
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
void TPZRenumbering::CornerEqs(int mincorners, int nelconsider, std::set<int> &cornernodes)
{
	
	TPZVec<int> nodtoelgraphindex;
	TPZVec<int> nodtoelgraph;
	int sub = 0;
	
	NodeToElGraph(fElementGraph,fElementGraphIndex,nodtoelgraph,nodtoelgraphindex);
	
	int nelem = fElementGraphIndex.NElements()-1;
    if (nelconsider > nelem) {
        DebugStop();
    }
	int element;
	for (element=0; element < nelconsider; element++) 
	{
		int firstindex = fElementGraphIndex[element];
		int lastindex = fElementGraphIndex[element+1];
		if (firstindex == lastindex) {
			continue;
		}
		TPZStack<int> corners;
		// a vector of sets of element connections for each node
		std::multimap<int,std::pair<int,std::set<int> > > connectivities;
		typedef std::multimap<int,std::pair<int,std::set<int> > > map_type;
		int ind;
		int maxelcon = 0;
		for (ind=firstindex; ind<lastindex; ind++) {
			int node = fElementGraph[ind];
			int firstelind = nodtoelgraphindex[node];
			int lastelind = nodtoelgraphindex[node+1];
			std::set<int> elcon;
			elcon.insert(&nodtoelgraph[firstelind],(&nodtoelgraph[lastelind-1])+1);
			maxelcon = maxelcon < (int)elcon.size() ? elcon.size() : maxelcon;
			connectivities.insert(map_type::value_type(elcon.size(), std::pair<int, std::set<int> >(node,elcon)));
		}
		
		map_type::reverse_iterator it = connectivities.rbegin();		
		// put all nodes with maximum connectivities on the stack
		int maxconnect = it->first;
		std::pair<map_type::const_iterator, map_type::const_iterator> p = connectivities.equal_range(it->first);
		map_type::const_iterator ithead;
		for (ithead = p.first; ithead != p.second; ithead++) 
		{
			corners.Push(ithead->second.first);
			cornernodes.insert(ithead->second.first);
		}
		
		// look the included sets
		map_type::iterator itf = connectivities.begin();		
		for (itf=connectivities.begin(); itf != connectivities.end(); itf++)
		{
			map_type::reverse_iterator it = connectivities.rbegin();
			while (it->first > itf->first && it!=connectivities.rend()) 
			{
				std::set<int>::iterator smallsetbeg, smallsetend, largesetbeg, largesetend;
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
				corners.Push(itf->second.first);
				cornernodes.insert(itf->second.first);
			}
		}
		if (corners.NElements() < mincorners) {
			it = connectivities.rbegin();
			int numelconnected = it->first-1;
			int ncorners = corners.NElements();
			while (numelconnected > 1 && ncorners < mincorners) {
				std::pair<map_type::const_iterator, map_type::const_iterator> p = connectivities.equal_range(numelconnected);
				map_type::const_iterator ithead;
				for (ithead = p.first; ithead != p.second; ithead++) 
				{
					corners.Push(ithead->second.first);
					cornernodes.insert(ithead->second.first);
					ncorners++;
				}
				numelconnected--;
			}
		}
		int ieq;
		for (ieq=0; ieq<corners.NElements(); ieq++) {
			if (cornernodes.find(corners[ieq]) == cornernodes.end()) {
				DebugStop();
			}
		}
		
		sub++;
	}
}

