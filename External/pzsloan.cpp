#include "pzsloan.h"
#include "math.h"
#include <fstream>
using namespace std;

TPZSloan::TPZSloan(int NElements, int NNodes) : TPZRenumbering(NElements,NNodes),
	fNodeWeights(0), fElementGraph(0), fElementGraphIndex(0) 
{
	fNNodes = NNodes;
	fNElements = NElements;
}

void TPZSloan::Resequence(int n_nodes, int n_elements, int *nnn,int *npn, int *xnpn, int old_profile, int new_profile)
{
	//Sloan routine with npn and xnpn as parameters
	//Renumbering using the npn and xnpn vectors already built.
	int nen;
//	int k=0;

	/**
	*Must be updated to work in accordance to TPZMetis!
	*/
	nen = 4; //Assumed as the highest node amount per element 

	int iadj = n_elements * nen * (nen-1);
	int *adj = new int[iadj];
	int *xadj = new int[iadj];
	int nop = 0;
	int inpn = 2 * n_elements;
	
	//adjacency list generation
	gegra_(&n_nodes, &n_elements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
	
	int *iw = new int[2*(n_elements+1)];
	int e2 = xadj[n_nodes]-1;

	//rewriting with new order
	//array nnn contains the new enumeration
	label_(&n_nodes , &e2, adj, xadj, nnn, iw, &old_profile, &new_profile);
	old_profile=0;
	new_profile=0;

}

void TPZSloan::Resequence(int * jj, int * jk, int n_nodes, int n_elements, int * nnn, int old_profile, int new_profile)
{
	//Sloan renumbering routine passing jj and jk as parameter
	//use only with bar elements
	
	int nen = 2;                          // number of nodes per element
	int *npn = new int[nen * n_elements]; // npn[i..i+nen] contains adjacent nodes of element i
	int *xnpn = new int[n_elements+1];	  // xnpn[i] index of element i in npn

	//feed npn and xnpn with data from jj and jk
	int i;
	int k=0;
	for(i=1;i<=nen*n_elements;i=i+nen)
	{
		k=k+1;
		npn[i-1]=jj[k];
		npn[i]=jk[k];
	}

	int j=1;
	for(i = 0;i <= n_elements; i++)
	{
		xnpn[i] = j;
		j=j+2;
	}


	int iadj = n_elements * nen * (nen-1);
	int *adj = new int[iadj];
	int *xadj = new int[iadj];
	int nop = 0;
	int inpn = 2 * n_elements;
	
	//adjacency list generation
	gegra_(&n_nodes, &n_elements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
	
	int *iw = new int[2*(n_elements+1)];
	int e2 = xadj[n_nodes]-1;

	//rewriting with new order
	//array nnn contains the new enumeration
	label_(&n_nodes , &e2, adj, xadj, nnn, iw, &old_profile, &new_profile);
	old_profile=0;
	new_profile=0;

	for(i=1;i<=n_elements;i++)
	{
		if(abs(jj[i]-jk[i]) > old_profile)
		{
			old_profile = abs(jj[i]-jk[i]);
		}
		if(abs(nnn[jj[i]-1]-nnn[jk[i]-1]) > new_profile)
		{
			new_profile = abs(nnn[jj[i]-1]-nnn[jk[i]-1]);
		}
	}


}



void TPZSloan::Resequence(TPZVec<int> &perm)
{
	/**
	int nen = 2;                          // number of nodes per element
	int *npn = new int[nen * n_elements]; // npn[i..i+nen] contains adjacent nodes of element i
	int *xnpn = new int[n_elements+1];	  // xnpn[i] index of element i in npn

	Now npn and xnpn are fElementGraph and fElementGraphIndex respectively
	*/
	//feed npn and xnpn with data from jj and jk
	int i;
//	int k=0;
	int nnodes_per_element=27;
	
	int iadj = fElementGraph.NElements()* nnodes_per_element* (nnodes_per_element-1);
	TPZVec<int> adj(iadj);
	TPZVec<int> xadj(iadj);
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

	gegra_(&fNNodes, &fNElements, &inpn, &fElementGraph[1], &fElementGraphIndex[1], &iadj, &adj[0], &xadj[0], &nop);
	//gegra_(&fNNodes, &fNElements, &inpn, npn, xnpn, &iadj, adj, xadj, &nop);
	
	TPZVec<int> iw(nnodes_per_element*(fNElements+1));
	//int *iw = new int [nnodes_per_element*(fNElements+1)];
	int e2 = xadj[fNNodes]-1;

	//int *NowPw
	//Where to obtain n_nodes ?
	//rewriting with new order
	//array nnn contains the new enumeration (nnn = perm on Metis)
	int old_profile=0;
	int new_profile=0;
	perm.Resize(fNNodes+1);

	
	label_(&fNNodes , &e2, &adj[1], &xadj[1], &perm[1], &iw[1], &old_profile, &new_profile);
	//label_(&fNNodes , &e2, adj, xadj, NowPerm, iw, &old_profile, &new_profile);
	TPZVec <int> aux(perm.NElements());
	aux = perm;
	perm.Resize(aux.NElements()-1);
	for(i=0;i<perm.NElements();i++) perm[i]=aux[i+1]-1;
}
void TPZSloan::ClearDataStructures(){
	fNodeWeights.Resize(0);
	fElementGraph.Resize(0);
	fElementGraphIndex.Resize(0);
}

void TPZSloan::SetElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex){
	int i;
	fElementGraph.Resize(elgraph.NElements()+1);
	fElementGraphIndex.Resize(elgraphindex.NElements()+1);
	for(i=0;i<elgraph.NElements();i++) fElementGraph[i+1] = elgraph[i]+1;
	for(i=0;i<elgraphindex.NElements();i++) fElementGraphIndex[i+1] = elgraphindex[i]+1;
}
void TPZSloan::SetNodeWeights(TPZVec<int> &weights){
	fNodeWeights = weights;
}
