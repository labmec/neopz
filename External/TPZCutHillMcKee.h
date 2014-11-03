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

    int NNodes(){
      return fnodegraphindex.NElements()-1;
    }

    int Degree(int node){
      return fnodegraphindex[node+1]-fnodegraphindex[node];
    }

  //accorging to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  //Previously, Sloan had suggested LastLevel.Resize( (LastLevel.NElements()+2)/2 );
    void ShrinkLastLevel(TPZVec<int> &LastLevel);

    void RootedLevelStructure(int rootNode, TPZStack< TPZVec<int> > &LevelStructure);

    int SmallestDegree(TPZVec<int> &ExploredNodes);

    void AdjacentNodes(int parent, TPZVec<int> &adjNodes);

    void GetAdjacentNodes(const TPZVec<int> &parents,
                          const TPZVec< int > &exceptedNodes,
                          std::set<int> &adjNodes);

    void AdjacentNodesOrdered(int parent, TPZVec<int> &adjNodes);

    void SortNodes(TPZVec<int> &nodes);

    void PseudoPeripheralNodes(int &startNode, int &endNode);

    void Set2Vec(const std::set<int> &myset, TPZVec<int> &myVec) const;

  };//SGraph

  private:

  void ProcessParentNode(int Parent,
                         SGraph &graph,
                         TPZVec<int> &ExploredNodes,
                         TPZStack<int> &R,
                         std::queue<int> &Q,
                         TPZVec<int> &adjNodes);

  bool fReverse;

  virtual void Resequence(TPZVec<int> &permGather, TPZVec<int> &permScatter,
                          TPZVec<int> &permGatherReverse, TPZVec<int> &permScatterReverse);

  public:


    virtual void Resequence(TPZVec<int> &perm, TPZVec<int> &iperm);

    TPZCutHillMcKee();

    TPZCutHillMcKee(int NElements, int NNodes, bool Reverse = true);

    virtual ~TPZCutHillMcKee();

    bool fVerbose;

};
//---------------------------------------------------------------------------
#endif
