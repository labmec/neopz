//---------------------------------------------------------------------------

#ifndef TPZCutHillMcKeeH
#define TPZCutHillMcKeeH

#include "pzrenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <queue>
#include "pzstack.h"

class TPZCutHillMcKee : public TPZRenumbering {

  public:

  struct SGraph{

    TPZManVector<long> fnodegraph;

    TPZManVector<long> fnodegraphindex;

    long NNodes(){
      return fnodegraphindex.NElements()-1;
    }

    long Degree(long node){
      return fnodegraphindex[node+1]-fnodegraphindex[node];
    }

  //accorging to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  //Previously, Sloan had suggested LastLevel.Resize( (LastLevel.NElements()+2)/2 );
    void ShrinkLastLevel(TPZVec<long> &LastLevel);

    void RootedLevelStructure(long rootNode, TPZStack< TPZStack<long> > &LevelStructure);

    long SmallestDegree(TPZVec<long> &ExploredNodes);

    void AdjacentNodes(long parent, TPZVec<long> &adjNodes);

    long * AdjacentNodesPtr(long parent, long &n){
        n = this->Degree(parent);
       return & ( fnodegraph[ fnodegraphindex[parent] ] );
    }//method

    void GetAdjacentNodes(const TPZVec<long> &parents,
                          const TPZVec< long > &exceptedNodes,
                          TPZStack<long> &adjNodes);

    void AdjacentNodesOrdered(long parent, const TPZVec<long> &exceptedNodes, TPZVec<long> &adjNodes);

    void SortNodes(TPZVec<long> &nodes);

    void PseudoPeripheralNodes(long &startNode, long &endNode);

    void Set2Vec(const std::set<long> &myset, TPZVec<long> &myVec) const;

  };//SGraph

  private:

  void ProcessParentNode(long Parent,
                         SGraph &graph,
                         TPZVec<long> &ExploredNodes,
                         TPZStack<long> &R,
                         std::queue<long> &Q,
                         TPZVec<long> &adjNodes);

  bool fReverse;

  virtual void Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter,
                          TPZVec<long> &permGatherReverse, TPZVec<long> &permScatterReverse);

  public:


    virtual void Resequence(TPZVec<long> &perm, TPZVec<long> &iperm);

    TPZCutHillMcKee();

    TPZCutHillMcKee(long NElements, long NNodes, bool Reverse = true);

    virtual ~TPZCutHillMcKee();

    bool fVerbose;

};
//---------------------------------------------------------------------------
#endif
