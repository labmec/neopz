/**
 * @file
 * @brief Contains the implementation of the TPZSloan methods. 
 */
#include "pzsloan.h"
#include "math.h"
#include <fstream>
using namespace std;

#include "pthread.h"

TPZSloan::TPZSloan(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes)
{
	fNNodes = NNodes;
	fNElements = NElements;
        fMaxNodesElement = 27;
}

// void TPZSloan::Resequence(int n_nodes, int n_elements, int *nnn,int *npn, int *xnpn, int old_profile, int new_profile)
// {
// 	//Sloan routine with npn and xnpn as parameters
// 	//Renumbering using the npn and xnpn vectors already built.
// 	int nen;
// //	int k=0;
// 
// 	/**
// 	*Must be updated to work in accordance to TPZMetis!
// 	*/
// 	nen = 27; //Assumed as the highest node amount per element
// 
// 	int iadj = n_elements * nen * (nen-1);
// 	int *adj = new int[iadj];
// 	int *xadj = new int[iadj];
// 	int nop = 0;
// 	int inpn = 2 * n_elements;
// 
// 	//adjacency list generation
// 	gegra_(&n_nodes, &n_elements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
// 
// 	int *iw = new int[2*(n_elements+1)];
// 	int e2 = xadj[n_nodes]-1;
// 
// 	//rewriting with new order
// 	//array nnn contains the new enumeration
// 	label_(&n_nodes , &e2, adj, xadj, nnn, iw, &old_profile, &new_profile);
// 	old_profile=0;
// 	new_profile=0;
// 
// }
// 
// void TPZSloan::Resequence(int * jj, int * jk, int n_nodes, int n_elements, int * nnn, int old_profile, int new_profile)
// {
// 	//Sloan renumbering routine passing jj and jk as parameter
// 	//use only with bar elements
// 
// 	int nen = 2;                          // number of nodes per element
// 	int *npn = new int[nen * n_elements]; // npn[i..i+nen] contains adjacent nodes of element i
// 	int *xnpn = new int[n_elements+1];	  // xnpn[i] index of element i in npn
// 
// 	//feed npn and xnpn with data from jj and jk
// 	int i;
// 	int k=0;
// 	for(i=1;i<=nen*n_elements;i=i+nen)
// 	{
// 		k=k+1;
// 		npn[i-1]=jj[k];
// 		npn[i]=jk[k];
// 	}
// 
// 	int j=1;
// 	for(i = 0;i <= n_elements; i++)
// 	{
// 		xnpn[i] = j;
// 		j=j+2;
// 	}
// 
// 
// 	int iadj = n_elements * nen * (nen-1);
// 	int *adj = new int[iadj];
// 	int *xadj = new int[iadj];
// 	int nop = 0;
// 	int inpn = 2 * n_elements;
// 
// 	//adjacency list generation
//   cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!\n" ;
// 	gegra_(&n_nodes, &n_elements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
// 
// 	int *iw = new int[2*(n_elements+1)];
// 	int e2 = xadj[n_nodes]-1;
// 
// 	//rewriting with new order
// 	//array nnn contains the new enumeration
// 	label_(&n_nodes , &e2, adj, xadj, nnn, iw, &old_profile, &new_profile);
// 	old_profile=0;
// 	new_profile=0;
// 
// 	for(i=1;i<=n_elements;i++)
// 	{
// 		if(abs(jj[i]-jk[i]) > old_profile)
// 		{
// 			old_profile = abs(jj[i]-jk[i]);
// 		}
// 		if(abs(nnn[jj[i]-1]-nnn[jk[i]-1]) > new_profile)
// 		{
// 			new_profile = abs(nnn[jj[i]-1]-nnn[jk[i]-1]);
// 		}
// 	}
// 
// 
// }



void TPZSloan::Resequence(TPZVec<int> &perm, TPZVec<int> &iperm)
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
	int i;
//	int k=0;
	int nnodes_per_element=fMaxNodesElement;
	
	TPZVec<int> elementgraph(fElementGraph.NElements()+1);
	TPZVec<int> elementgraphindex(fElementGraphIndex.NElements()+1);
	
	int n = fElementGraph.NElements();
	for (i=0; i<n; i++) {
		elementgraph[i+1] = fElementGraph[i]+1;
	}
	n = fElementGraphIndex.NElements();
	for (i=0; i<n; i++) {
		elementgraphindex[i+1] = fElementGraphIndex[i]+1;
	}

	int iadj = elementgraph.NElements()* nnodes_per_element* (nnodes_per_element-1);
	TPZVec<int> adj(iadj);
	TPZVec<int> xadj(fNNodes+2,-1);
	int nop = 0;
	int inpn = nnodes_per_element * fNElements;

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
	pthread_mutex_lock(&Lock_clindex);
	gegra_(&fNNodes, &fNElements, &inpn, &elementgraph[1], &elementgraphindex[1], &iadj, &adj[0], &xadj[0], &nop);
	pthread_mutex_unlock(&Lock_clindex);
	//gegra_(&fNNodes, &fNElements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
#ifdef SLOANDEBUG
	cout << "node index vector ";
	int no;
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
	TPZVec<int> iw(nnodes_per_element*(fNElements+1)*10);
	//int *iw = new int [nnodes_per_element*(fNElements+1)];
	int e2 = xadj[fNNodes]-1;

	//int *NowPw
	//Where to obtain n_nodes ?
	//rewriting with new order
	//array nnn contains the new enumeration (nnn = perm on Metis)
	int old_profile=0;
	int new_profile=0;
	perm.Resize(fNNodes+1);

	pthread_mutex_lock(&Lock_clindex);
  label_(&fNNodes , &e2, &adj[0], &xadj[0], &perm[1], &iw[1], &old_profile, &new_profile);
	pthread_mutex_unlock(&Lock_clindex);

	//label_(&fNNodes , &e2, adj, xadj, NowPerm, iw, &old_profile, &new_profile);
	//cout << __PRETTY_FUNCTION__ << " oldprofile " << old_profile << " newprofile " << new_profile << endl;
	TPZVec <int> aux(perm.NElements());
	aux = perm;
	perm.Resize(aux.NElements()-1);
	for(i=0;i<perm.NElements();i++) perm[i]=aux[i+1]-1;
	iperm.Resize(perm.NElements());
	for(i=0;i<perm.NElements();i++) iperm[perm[i]]=i;
}

