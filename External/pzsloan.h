#ifndef TPZSLOAN_H
#define TPZSLOAN_H

#include "pzrenumbering.h"
//#include "sloan\\sloan.h"


#include "sloan.h"


class TPZSloan : TPZRenumbering {
 public:
  void Resequence(int * jj, int * jk, int n_nodes, int n_elements, int * nnn, int old_profile, int new_profile);
  void Resequence(int n_nodes, int n_elements, int *nnn,int *npn, int *xnpn, int old_profile, int new_profile);
  void Resequence(TPZVec<int> &perm, TPZVec<int> &iperm);

  void SetElementGraph(TPZVec<int> &elgraph, TPZVec<int> &elgraphindex);
  /**Sets the number of equations associated with each node
     The derived class may or may not take this data into 
     consideration*/
  void SetNodeWeights(TPZVec<int> &weights);
  /**This will reset all datastructures the object may contain.
     Node resequencing algorithms may require a possibly large 
     amount of temporary data*/
  virtual void ClearDataStructures();



  
 private:

  //ConvertGraph
  /**Number of equations associated with each node*/
  TPZVec<int> fNodeWeights;
  

  /**Node number of each element*/
  TPZVec<int> fElementGraph;
  

  /**Indicates for each element the index of the first entry with
     fElementGraph for that element
     Size of the vector fNElements+1*/
  TPZVec<int> fElementGraphIndex;



  int fNNodes;

  int fNElements;
	
 public:
  TPZSloan(int NElements, int NNodes);

  virtual ~TPZSloan()
     {
     }

};

#endif

