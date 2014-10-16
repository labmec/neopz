//---------------------------------------------------------------------------

#ifndef TPZCutHillMcKeeH
#define TPZCutHillMcKeeH

#include "pzrenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <queue>
#include "pzstack.h"

class TPZCutHillMcKee : public TPZRenumbering {

  private:

  struct SGraph{

    TPZManVector<long> fnodegraph;

    TPZManVector<long> fnodegraphindex;

    int NNodes(){
      return fnodegraphindex.NElements()-1;
    }

    int Degree(int node){
      return fnodegraphindex[node+1]-fnodegraphindex[node];
    }

    int SmallestDegree(TPZVec<int> &ExploredNodes){
      std::cout << "TPZCutHillMcKee::SGraph::SmallestDegree\n";std::cout.flush();
      int mindegree = 1000000;
      int found = -1;
      const int n = this->NNodes();
      for(int i = 0; i < n; i++){
        if(ExploredNodes[i] == 1) continue;
        int degree = this->Degree(i);
        if(degree < mindegree){
          mindegree = degree;
          found = i;
        }
      }//for
      return found;
    }

    void AdjacentNodesOrdered(int parent, TPZVec<int> &adjNodes){
      std::multimap<int,int> order;

      const unsigned int n = this->Degree(parent);
      for(int i = 0; i < n; i++){
        const int adj = fnodegraph[fnodegraphindex[parent]+i];
        const int adjDegree = this->Degree(adj);
        order.insert(std::make_pair(adjDegree,adj));
      }//for
      adjNodes.Resize(n);
#ifdef DEBUG
      if(n!=order.size()){
        DebugStop();
      }
#endif
      std::multimap<int,int>::const_iterator mapIt;
      int i;
      for(i = 0, mapIt = order.begin(); mapIt != order.end(); i++, mapIt++){
        adjNodes[i] = mapIt->second;
      }
    }//void

  };//SGraph


  void ProcessParentNode(int Parent,
                         SGraph &graph,
                         TPZVec<int> &ExploredNodes,
                         TPZStack<int> &R,
                         std::queue<int> &Q,
                         TPZVec<int> &adjNodes);

  bool fReverse;

  public:

    virtual void Resequence(TPZVec<int> &perm, TPZVec<int> &iperm);

    TPZCutHillMcKee();

    TPZCutHillMcKee(int NElements, int NNodes, bool Reverse = true);

    virtual ~TPZCutHillMcKee();

};
//---------------------------------------------------------------------------
#endif
