/**
 * @file
 * @brief Contains the implementation of the TPZSloan methods. 
 */

#include "pzsloan.h"
#include "math.h"
#include <fstream>
using namespace std;

#include "pz_pthread.h"

TPZSloan::TPZSloan(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes)
{
	fNNodes = NNodes;
	fNElements = NElements;
        fMaxNodesElement = 27;
}

void TPZSloan::Resequence(TPZVec<long> &perm, TPZVec<long> &iperm)
{

	if(!fNNodes || !fNElements)
 	{
	   std::cout << __PRETTY_FUNCTION__ << " called with fNNodes = " << fNNodes
	       << " and fNElements = " << fNElements << std::endl;
	   return;
	}
	
	/**
	int nen = 2;                          // number of nodes per element
	int *npn = new int[nen * n_elements]; // npn[i..i+nen] contains adjacent nodes of element i
	int *xnpn = new int[n_elements+1];	  // xnpn[i] index of element i in npn

	Now npn and xnpn are fElementGraph and fElementGraphIndex respectively
	*/
	//feed npn and xnpn with data from jj and jk
	long i;
//	int k=0;
	int nnodes_per_element=fMaxNodesElement;
	
	TPZVec<long> elementgraph(fElementGraph.NElements()+1);
	TPZVec<long> elementgraphindex(fElementGraphIndex.NElements()+1);
	
	long n = fElementGraph.NElements();
	for (i=0; i<n; i++) {
		elementgraph[i+1] = fElementGraph[i]+1;
	}
	n = fElementGraphIndex.NElements();
	for (i=0; i<n; i++) {
		elementgraphindex[i+1] = fElementGraphIndex[i]+1;
	}

	long iadj = elementgraph.NElements()* nnodes_per_element* (nnodes_per_element-1);
	TPZVec<long> adj(iadj);
	TPZVec<long> xadj(fNNodes+2,-1);
	int nop = 0;
	long inpn = nnodes_per_element * fNElements;

	//adjacency list generation
	//ofstream sai("felementgraph.txt");
	//fElementGraph.Print(sai);
	//int *npn = new int [iadj];
	//npn=fElementGraph;
	//int *xnpn = new int[iadj];
	//xnpn=fElementGraphIndex;
	//int *adj = new int [iadj];
	//int *xadj = new int[iadj];
	// Print the element graph
#ifdef SLOANDEBUG
	int iel;
	for(iel=0; iel<fNElements; iel++)
	{
		int firstindex = elementgraphindex[iel+1];
		int lastindex = elementgraphindex[iel+2];
		int no;
		cout << "Element : " << iel+1 << " : ";
		for(no = firstindex; no<lastindex; no++)
		{
			cout << elementgraph[no] << " ";
		}
		cout << endl;
	}
#endif
  static pthread_mutex_t Lock_clindex = PTHREAD_MUTEX_INITIALIZER;
	PZ_PTHREAD_MUTEX_LOCK(&Lock_clindex,"TPZSloan::Resequence()");
//	int elementgraph_int = (int)elementgraph[1];
//	TPZVec<int> elementgraphindex_int = ((TPZVec<int>)elementgraphindex);
	gegra_(&fNNodes, &fNElements, &inpn, &elementgraph[1], &elementgraphindex[1], &iadj, &adj[0], &xadj[0], &nop);
//	elementgraph[1] = elementgraph_int;
//	elementgraphindex[1] = elementgraphindex_int;
	PZ_PTHREAD_MUTEX_UNLOCK(&Lock_clindex,"TPZSloan::Resequence()");
	//gegra_(&fNNodes, &fNElements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
#ifdef SLOANDEBUG
	cout << "node index vector ";
	long no;
	for(no=0; no<xadj.NElements(); no++)
	{
		cout << xadj[no] << " ";
		if(no && !(no%30)) cout << endl;

	}
	cout << endl;
	for(no=0; no<fNNodes; no++)
	{
		int firstindex = xadj[no]-1;
		int lastindex = xadj[no+1]-1;
		int el;
		cout << "Node : " << no+1 << " : ";
		for(el = firstindex; el<lastindex; el++)
		{
			cout << adj[el] << " ";
		}
		cout << endl;
	}
#endif
	TPZVec<long> iw(nnodes_per_element*(fNElements+1)*10);
	//int *iw = new int [nnodes_per_element*(fNElements+1)];
	long e2 = xadj[fNNodes]-1;

	//int *NowPw
	//Where to obtain n_nodes ?
	//rewriting with new order
	//array nnn contains the new enumeration (nnn = perm on Metis)
	long old_profile=0;
	long new_profile=0;
	perm.Resize(fNNodes+1);

	PZ_PTHREAD_MUTEX_LOCK(&Lock_clindex,"TPZSloan::Resequence()");
  label_(&fNNodes , &e2, &adj[0], &xadj[0], &perm[1], &iw[1], &old_profile, &new_profile);
	PZ_PTHREAD_MUTEX_UNLOCK(&Lock_clindex,"TPZSloan::Resequence()");

	//label_(&fNNodes , &e2, adj, xadj, NowPerm, iw, &old_profile, &new_profile);
	//cout << __PRETTY_FUNCTION__ << " oldprofile " << old_profile << " newprofile " << new_profile << endl;
	TPZVec <long> aux(perm.NElements());
	aux = perm;
	perm.Resize(aux.NElements()-1);
	for(i=0;i<perm.NElements();i++) perm[i]=aux[i+1]-1;
	iperm.Resize(perm.NElements());
	for(i=0;i<perm.NElements();i++) iperm[perm[i]]=i;
}

