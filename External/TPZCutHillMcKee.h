//---------------------------------------------------------------------------

#ifndef TPZCutHillMcKeeH
#define TPZCutHillMcKeeH

#include "TPZRenumbering.h"
#include <map>
#include "pzmanvector.h"
#include <queue>
#include "pzstack.h"

class TPZCutHillMcKee : public TPZRenumbering {

  public:

  struct SGraph{

    TPZManVector<int64_t> fnodegraph;

    TPZManVector<int64_t> fnodegraphindex;

    int64_t NNodes(){
      return fnodegraphindex.NElements()-1;
    }

    int64_t Degree(int64_t node){
      return fnodegraphindex[node+1]-fnodegraphindex[node];
    }

  //according to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  //Previously, Sloan had suggested LastLevel.Resize( (LastLevel.NElements()+2)/2 );
    void ShrinkLastLevel(TPZVec<int64_t> &LastLevel);

    void RootedLevelStructure(int64_t rootNode, TPZStack< TPZStack<int64_t> > &LevelStructure);

    int64_t SmallestDegree(TPZVec<int64_t> &ExploredNodes);

    void AdjacentNodes(int64_t parent, TPZVec<int64_t> &adjNodes);

    int64_t * AdjacentNodesPtr(int64_t parent, int64_t &n){
        n = this->Degree(parent);
       return & ( fnodegraph[ fnodegraphindex[parent] ] );
    }//method

    void GetAdjacentNodes(const TPZVec<int64_t> &parents,
                          const TPZVec< int64_t > &exceptedNodes,
                          TPZStack<int64_t> &adjNodes);

    void AdjacentNodesOrdered(int64_t parent, const TPZVec<int64_t> &exceptedNodes, TPZVec<int64_t> &adjNodes);

    void SortNodes(TPZVec<int64_t> &nodes);

    void PseudoPeripheralNodes(int64_t &startNode, int64_t &endNode);

    void Set2Vec(const std::set<int64_t> &myset, TPZVec<int64_t> &myVec) const;

  };//SGraph

  private:

  void ProcessParentNode(int64_t Parent,
                         SGraph &graph,
                         TPZVec<int64_t> &ExploredNodes,
                         TPZStack<int64_t> &R,
                         std::queue<int64_t> &Q,
                         TPZVec<int64_t> &adjNodes);

  bool fReverse;

  virtual void Resequence(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter,
                          TPZVec<int64_t> &permGatherReverse, TPZVec<int64_t> &permScatterReverse);

  public:


    virtual void Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &iperm);

    TPZCutHillMcKee();

    TPZCutHillMcKee(int64_t NElements, int64_t NNodes, bool Reverse = true);

    virtual ~TPZCutHillMcKee();

    bool fVerbose;

};
//---------------------------------------------------------------------------
#endif
