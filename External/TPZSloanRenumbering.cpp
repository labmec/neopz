//---------------------------------------------------------------------------

#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

TPZSloanRenumbering::TPZSloanRenumbering(int NElements, int NNodes):TPZRenumbering(NElements,NNodes){

}

TPZSloanRenumbering::~TPZSloanRenumbering(){

}

void TPZSloanRenumbering::Resequence(TPZVec<int> &permGather, TPZVec<int> &permScatter){

  //computing graph
  TPZCutHillMcKee::SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();

  //finding pseudo peripheral nodes
  const int nnodes = graph.NNodes();
  int startNode = -1, endNode = -1;
  graph.PseudoPeripheralNodes(startNode, endNode);

  //computing distances to end node
  TPZStack< TPZVec<int> > LevelStructure;
  graph.RootedLevelStructure(endNode, LevelStructure);
  TPZVec<int> DistanceToEndNode( nnodes, -1 );
  for(int ilevel = 0; ilevel < LevelStructure.NElements(); ilevel++){
    for(int iel = 0; iel < LevelStructure[ ilevel ].NElements(); iel++){
      const int node = LevelStructure[ ilevel ][ iel ];
      DistanceToEndNode[ node ] = ilevel;
    }//for iel
  }//for ilevel
#ifdef DEBUG
  for(int i = 0; i < DistanceToEndNode.NElements(); i++){
    if(DistanceToEndNode[i] == -1) DebugStop();
  }
#endif

  //assigning initial status and priority
  TPZVec<int> priority( nnodes, -1);
  for(int i = 0; i < nnodes; i++){
    const int degree = graph.Degree(i);
    priority[i] = this->W1()*DistanceToEndNode[i] - this->W2()*(degree+1);
  }//for i

  TPZVec<int> Status(nnodes,EInactive);
  TPZStack<int> R;
  TPZManVector<int,10000> adjNodes, adjNodes2J;
  std::list<int> Q;
  Status[startNode] = EPreActive;
  Q.push_back( startNode );
  while(Q.size()){
    const int currnode = this->FindHighestPriority(Q,priority);
    Q.remove(currnode);
    if(Status[currnode] == EPreActive){
      graph.AdjacentNodes(currnode, adjNodes);
      for(int jadj = 0; jadj < adjNodes.NElements(); jadj++){
        const int adjNode = adjNodes[jadj];
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        if(Status[ adjNode] == EInactive){
          Status[ adjNode] = EPreActive;
          Q.push_back( adjNode );
        }
      }//jadj
    }//not a preactive node
    R.Push( currnode );
    Status[ currnode ] = EPostActive;

    //step 9: updating priorities and queue
    graph.AdjacentNodes(currnode, adjNodes);
    for(int jadj = 0; jadj < adjNodes.NElements(); jadj++){
      const int adjNode = adjNodes[jadj];
      if( Status[ adjNode ] == EPreActive ){
        Status[ adjNode ] = EActive;
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        graph.AdjacentNodes(adjNode, adjNodes2J);
        for(int k = 0; k < adjNodes2J.NElements(); k++){
          const int knode = adjNodes2J[k];
          if(Status[knode] != EPostActive){
            priority[knode] += this->W2();
          }//if ! post active
          if( Status[knode] == EInactive){
            Status[knode] = EPreActive;
            Q.push_back(knode);
          }// if inactive
        }//loop of k nodes
      }//if
    }//for jadj

  }//while

  if(R.NElements() != nnodes){
    std::cout << "TPZSloanRenumbering::Resequence error! Not all nodes were processed.\n"
              << "It may be an implementation bug or simply a graph with independent nodes\n"
              << "as if the mesh has subdomains that are not connected.\n"
              << "For the later, use CuttHillMcKee, which is prepared for this particular case\n";
    DebugStop();
  }

#ifdef DEBUG
{//verificando se ha duplicados
  std::set<int> check;
  for(int i = 0; i < nnodes; i++) check.insert(R[i]);
  if( ((int)(check.size())) != nnodes) DebugStop();
}
#endif

  permScatter = R;
  permGather.Resize(nnodes);
  for(int i = 0; i < nnodes; i++) permGather[ permScatter[i] ] = i;

}//void

int TPZSloanRenumbering::FindHighestPriority(const std::list<int> &Q, const TPZVec<int> &priority) const{
  if(Q.size() == 0) DebugStop();
  std::list<int>::const_iterator w = Q.begin(), e = Q.end();
  int result = *w;
  int nodepriority = priority[result];
  for( ; w != e; w++){
    int lcnode = *w;
    if(priority[lcnode] > nodepriority){
      result = lcnode;
      nodepriority = priority[lcnode];
    }
  }
  return result;

}//int



