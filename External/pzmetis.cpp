/**
 * @file
 * @brief Contains the implementation of the TPZMetis methods. 
 */

#include "pzmetis.h"

#ifdef USING_METIS
#include <math.h>
extern "C" {
#include "metis.h"
};
#endif

#include "pzlog.h"

#ifdef PZ_LOG
static TPZLogger logger("pz.metis");
#endif

#include <iostream>
using namespace std;

TPZMetis::TPZMetis() : TPZRenumbering()
{
    PZError<<"TPZMetis depends on the Metis library\n";
    PZError<<"Please reconfigure NeoPZ library using:\n";
    PZError<<"USING_METIS=ON"<<std::endl;
    DebugStop();
}

void TPZMetis::Print(std::ostream &out,char * title) {
	int nel = fElementGraphIndex.NElements()-1;
	int el;
	out << title;
	out << "\nTPZMetis::Print fNElements = " << fNElements << " fNNodes = " << fNNodes << endl;
	for (el=0;el<nel;el++) {
		int firstindex = fElementGraphIndex[el];
		int lastindex = fElementGraphIndex[el+1];
		int index;
		out << "Element number " << el << " : ";
		for (index=firstindex;index<lastindex;index++) {
			out << fElementGraph[index] << " ";
		}
		out << endl;
	}
	TPZManVector<int64_t> nodegraph(0),nodegraphindex(0);
	ConvertGraph(fElementGraph,fElementGraphIndex,nodegraph,nodegraphindex);
	int numelnodegraph = nodegraphindex[fNNodes];
	if (numelnodegraph == nodegraph.NElements() ) {
		nodegraph.Resize(numelnodegraph+1);
	}
	int nod;
	for (nod = numelnodegraph; nod>0; nod--) nodegraph[nod] = nodegraph[nod-1];
	for (el=0;el<fNNodes;el++) {
		int firstindex = nodegraphindex[el];
		int lastindex = nodegraphindex[el+1];
		int index;
		out << "Node number " << el << " : ";
		for (index=firstindex;index<lastindex;index++) {
			out << nodegraph[index+1] << " ";
		}
		out << endl;
	}
}

void TPZMetis::Print(std::ostream &out) {
	
	// 	int nel = fElementGraphIndex.NElements()-1;
	int el;
	out << fNNodes << std::endl;
	/*   for(el=0;el<nel;el++) {
	 int firstindex = fElementGraphIndex[el];
	 int lastindex = fElementGraphIndex[el+1];
	 int index;
	 for(index=firstindex;index<lastindex;index++) {
	 if(fElementGraph[index] == index) continue;//o atual n�o deve aparecer
	 out << (fElementGraph[index]+1) << " ";
	 }
	 out << endl;
	 } */
	TPZManVector<int64_t> nodegraph(0),nodegraphindex(0);
	ConvertGraph(fElementGraph,fElementGraphIndex,nodegraph,nodegraphindex);
	int numelnodegraph = nodegraphindex[fNNodes];
	if (numelnodegraph == nodegraph.NElements() ) {
		nodegraph.Resize(numelnodegraph+1);
	}
	int nod;
	for (nod = numelnodegraph; nod>0; nod--) nodegraph[nod] = nodegraph[nod-1];
	for (el=0;el<fNNodes;el++) {
		int firstindex = nodegraphindex[el];
		int lastindex = nodegraphindex[el+1];//come�o do pr�ximo
		int index;
		for (index=firstindex;index<lastindex;index++) {
			out << (nodegraph[index+1]+1) << " ";
		}
		out << endl;
	}
}

void TPZMetis::Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &inverseperm) {
	TPZManVector<int64_t> nodegraph(0),nodegraphindex(0);
	ConvertGraph(fElementGraph,fElementGraphIndex,nodegraph,nodegraphindex);
	int64_t numelnodegraph = nodegraphindex[fNNodes];
	if (numelnodegraph == nodegraph.NElements() )
	{
		nodegraph.Resize(numelnodegraph+1);
	}
	int64_t nod;
	for (nod = numelnodegraph; nod>0; nod--) nodegraph[nod] = nodegraph[nod-1];
	perm.Resize(fNNodes);
	inverseperm.Resize(fNNodes);
	for(nod=0;nod<fNNodes;nod++)
	{
		perm[nod] = inverseperm[nod] = nod;
	}

#ifdef USING_METIS
	TPZVec<int> nodegraphInt(0),nodegraphindexInt(0);
	int NNodes = (int) fNNodes;
	int64_t n, sz = nodegraph.NElements();
	nodegraph.Resize(sz);
	for(n=0;n<sz;n++)
		nodegraphInt[n] = (int)nodegraph[n];
	sz = nodegraphindex.NElements();
	for(n=0;n<sz;n++)
		nodegraphindexInt[n] = (int)nodegraphindex[n];

	// Using external library METIS 5
	int numflag = 0;
	int options = 0;
    int nperms = perm.NElements();
    int ninvers = inverseperm.NElements();
    int *permint = new int[nperms];
    int *inversepermint = new int[ninvers];
    if(!permint || !inverseperm.size()) {
        std::cout << "TPZMetis::Resequence memory is not enough.\n";
        return;
    }
    int64_t i;
    for(i=0L;i<nperms;i++)
        permint[i] = (int)perm[i];
    for(i=0L;i<nperms;i++)
        inversepermint[i] = (int)inverseperm[i];
    
//	METIS_NodeND(&fNNodes,&nodegraphindex[0],&nodegraph[1],&numflag,&options,&perm[0],&inverseperm[0]);
    METIS_NodeND(&NNodes,&nodegraphindexInt[0],&nodegraphInt[1],&numflag,&options,permint,inversepermint);
	fNNodes = (int64_t)NNodes;
#endif
}

void TPZMetis::Subdivide(int nParts, TPZVec < int > & Domains)
{
	TPZManVector<int64_t> Adjacency,AdjacencyIndex;
	TPZManVector<int> AdjacencyWeight;
	ConvertToElementoToElementGraph(fElementGraph,fElementGraphIndex,Adjacency,AdjacencyWeight,AdjacencyIndex);
	
#ifdef PZ_LOG
	{
		std::stringstream sout;
		TPZRenumbering::Print(Adjacency,AdjacencyIndex,"Element to element graph",sout);
		if (logger.isDebugEnabled())
		{
			LOGPZ_DEBUG(logger, sout.str())
		}
	}
#endif
	
#ifdef USING_METIS
	TPZManVector<int> AdjacencyInt,AdjacencyIndexInt;
	int64_t n, nVertices = AdjacencyIndex.NElements();
	AdjacencyIndexInt.Resize(nVertices,0);
	for(n=0;n<nVertices;n++)
		AdjacencyIndexInt[n] = (int)AdjacencyIndex[n];
	int64_t nEdges = Adjacency.NElements();
	AdjacencyInt.Resize(nEdges,0);
	for(n=0;n<nEdges;n++)
		AdjacencyInt[n] = (int)Adjacency[n];
	int nvertices = (int)nVertices-1;
//	int nedges = (int)nEdges;
	Domains.Resize(nvertices,0);
	// Upon successful completion, nEdgesCutted stores the edge-cut or the total communication volume of the partitioning solution.
	int nEdgesCutted = 0;

	TPZVec<int> Options(METIS_NOPTIONS);
	METIS_SetDefaultOptions(&Options[0]);
	
    int ncon = 2;
	if(METIS_PartGraphRecursive(&nvertices, &ncon, &AdjacencyIndexInt[0], &AdjacencyInt[0], NULL, NULL, NULL,   // &AdjacencyWeight[0],
					&nParts, NULL, NULL, &Options[0], &nEdgesCutted, &Domains[0]) != METIS_OK)
		DebugStop();
#else
    DebugStop();
#endif
	
}


