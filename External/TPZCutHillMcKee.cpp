
#include "TPZCutHillMcKee.h"


TPZCutHillMcKee::TPZCutHillMcKee():TPZRenumbering(){
  fReverse = true;
}

TPZCutHillMcKee::TPZCutHillMcKee(int NElements, int NNodes, bool Reverse):TPZRenumbering(NElements,NNodes){
  fReverse = Reverse;
}

TPZCutHillMcKee::~TPZCutHillMcKee(){

}

void TPZCutHillMcKee::Resequence(TPZVec<int> &permGather, TPZVec<int> &permScatter){

  std::cout << "TPZCutHillMcKee ConvertGraph...";std::cout.flush();
  SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();
  std::cout << " done!\n";std::cout.flush();

        //CutHillMcKee
  std::queue<int> Q;
  TPZStack<int> R;
  TPZManVector<int> adjNodes;
  //reservando memoria em R
  const int nnodes = graph.NNodes();
  R.Resize(nnodes);
  R.Resize(0);
  TPZVec<int> ExploredNodes(nnodes,0);
  std::cout << "TPZCutHillMcKee process...\n";std::cout.flush();
  while(R.NElements() != graph.NNodes()){
    const int Parent = graph.SmallestDegree(ExploredNodes);
    if(Parent == -1){
      if(R.NElements() == graph.NNodes()) return;
      else DebugStop();
    }

    this->ProcessParentNode(Parent, graph, ExploredNodes, R, Q, adjNodes);
    while(Q.size()){
      const int Child = Q.front();
      Q.pop();
      this->ProcessParentNode(Child, graph, ExploredNodes, R, Q, adjNodes);
    }

  }//loop completo
  std::cout << " done!\n";std::cout.flush();

  const int n = R.NElements();
  if(n != nnodes) DebugStop();

#ifdef DEBUG
{//verificando se ha duplicados
  std::set<int> check;
  for(int i = 0; i < n; i++) check.insert(R[i]);
  if( ((int)(check.size())) != n) DebugStop();
}
#endif

  std::cout << "TPZCutHillMcKee Filling perm and iperm vectors...";std::cout.flush();
  if(fReverse){
    permScatter.Resize(n);
    for(int i = 0; i < n; i++) permScatter[n-1-i] = R[i];
  }
  else{
    permScatter = R;
  }
  permGather.Resize(n);
  for(int i = 0; i < n; i++) permGather[ permScatter[i] ] = i;
  std::cout << " done!\n";std::cout.flush();

}///void

void TPZCutHillMcKee::ProcessParentNode(int Parent,
                                        SGraph &graph,
                                        TPZVec<int> &ExploredNodes,
                                        TPZStack<int> &R,
                                        std::queue<int> &Q,
                                        TPZVec<int> &adjNodes){

    R.Push( Parent );
    ExploredNodes[Parent] = 1;
    graph.AdjacentNodesOrdered(Parent, adjNodes);
    const int nadj = adjNodes.NElements();
    for(int i = 0; i < nadj; i++){
      if(ExploredNodes[adjNodes[i]] == 0){//then a new element
        Q.push(adjNodes[i]);
        ExploredNodes[adjNodes[i]] = 1;
      }
    }//for


}//void

