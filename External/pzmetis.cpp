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

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.metis"));
#endif

#include <iostream>
using namespace std;

TPZMetis::TPZMetis(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes)
{
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
	TPZVec<int> nodegraph(0),nodegraphindex(0);
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
	TPZVec<int> nodegraph(0),nodegraphindex(0);
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

void TPZMetis::Resequence(TPZVec<int> &perm, TPZVec<int> &inverseperm){
	TPZVec<int> nodegraph(0),nodegraphindex(0);
	ConvertGraph(fElementGraph,fElementGraphIndex,nodegraph,nodegraphindex);
	int numelnodegraph = nodegraphindex[fNNodes];
	if (numelnodegraph == nodegraph.NElements() )
	{
		nodegraph.Resize(numelnodegraph+1);
	}
	int nod;
	for (nod = numelnodegraph; nod>0; nod--) nodegraph[nod] = nodegraph[nod-1];
	perm.Resize(fNNodes);
	inverseperm.Resize(fNNodes);
	for(nod=0;nod<fNNodes;nod++) 
	{
		perm[nod] = inverseperm[nod] = nod;
	}
	
#ifdef USING_METIS
	int numflag = 0;
	int options = 0;
	METIS_NodeND(&fNNodes,&nodegraphindex[0],&nodegraph[1],&numflag,&options,&perm[0],&inverseperm[0]);
#endif
}

void TPZMetis::Subdivide(int nParts, TPZVec < int > & Domains)
{
	TPZManVector<int> Adjacency,AdjacencyIndex;
	TPZManVector<int> AdjacencyWeight;
	ConvertToElementoToElementGraph(fElementGraph,fElementGraphIndex,Adjacency,AdjacencyWeight,AdjacencyIndex);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		TPZRenumbering::Print(Adjacency,AdjacencyIndex,"Element to element graph",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
#ifdef USING_METIS
	int nVertices = AdjacencyIndex.NElements() -1;

	Domains.Resize(nVertices);
	// Upon successful completion, nEdgesCutted stores the edge-cut or the total communication volume of the partitioning solution.
	int nEdgesCutted = 0;

	TPZVec<int> Options(METIS_NOPTIONS);
	METIS_SetDefaultOptions(&Options[0]);
	
    int ncon = 2;
	if(METIS_PartGraphRecursive(&nVertices, &ncon, &AdjacencyIndex[0], &Adjacency[0], NULL, NULL, NULL,   // &AdjacencyWeight[0],
					&nParts, NULL, NULL, &Options[0], &nEdgesCutted, &Domains[0]) != METIS_OK)
		DebugStop();
#else
    DebugStop();
#endif
	
}


