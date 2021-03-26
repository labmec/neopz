/**
 * @file
 * @brief Contains the implementation of the TPZSloan methods. 
 */

#include <fstream>
#include <mutex>
#include "pzsloan.h"
#include "math.h"
using namespace std;

TPZSloan::TPZSloan(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes)
{
	fNNodes = NNodes;
	fNElements = NElements;
        fMaxNodesElement = 27;
}

void TPZSloan::Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &iperm)
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
	int64_t i;
//	int k=0;
	int nnodes_per_element=fMaxNodesElement;
	
	TPZVec<int64_t> elementgraph(fElementGraph.NElements()+1);
	TPZVec<int64_t> elementgraphindex(fElementGraphIndex.NElements()+1);
	
	int64_t n = fElementGraph.NElements();
	for (i=0; i<n; i++) {
		elementgraph[i+1] = fElementGraph[i]+1;
	}
	n = fElementGraphIndex.NElements();
	for (i=0; i<n; i++) {
		elementgraphindex[i+1] = fElementGraphIndex[i]+1;
	}

	int64_t iadj = elementgraph.NElements()* nnodes_per_element* (nnodes_per_element-1);
	TPZVec<int64_t> adj(iadj);
	TPZVec<int64_t> xadj(fNNodes+2,-1);
	int nop = 0;
	int64_t inpn = nnodes_per_element * fNElements;

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
  static std::mutex mutex_clindex;
  std::unique_lock<std::mutex> lck(mutex_clindex);
//	int elementgraph_int = (int)elementgraph[1];
//	TPZVec<int> elementgraphindex_int = ((TPZVec<int>)elementgraphindex);
	gegra_(&fNNodes, &fNElements, &inpn, &elementgraph[1], &elementgraphindex[1], &iadj, &adj[0], &xadj[0], &nop);
//	elementgraph[1] = elementgraph_int;
//	elementgraphindex[1] = elementgraphindex_int;
  lck.unlock();
	//gegra_(&fNNodes, &fNElements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
#ifdef SLOANDEBUG
	cout << "node index vector ";
	int64_t no;
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
	TPZVec<int64_t> iw(nnodes_per_element*(fNElements+1)*10);
	//int *iw = new int [nnodes_per_element*(fNElements+1)];
	int64_t e2 = xadj[fNNodes]-1;

	//int *NowPw
	//Where to obtain n_nodes ?
	//rewriting with new order
	//array nnn contains the new enumeration (nnn = perm on Metis)
	int64_t old_profile=0;
	int64_t new_profile=0;
	perm.Resize(fNNodes+1);
  
  lck.lock();
  label_(&fNNodes , &e2, &adj[0], &xadj[0], &perm[1], &iw[1], &old_profile, &new_profile);
  lck.unlock();

	//label_(&fNNodes , &e2, adj, xadj, NowPerm, iw, &old_profile, &new_profile);
	//cout << __PRETTY_FUNCTION__ << " oldprofile " << old_profile << " newprofile " << new_profile << endl;
	TPZVec <int64_t> aux(perm.NElements());
	aux = perm;
	perm.Resize(aux.NElements()-1);
	for(i=0;i<perm.NElements();i++) perm[i]=aux[i+1]-1;
	iperm.Resize(perm.NElements());
	for(i=0;i<perm.NElements();i++) iperm[perm[i]]=i;
}
