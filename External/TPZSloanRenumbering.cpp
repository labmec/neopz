//---------------------------------------------------------------------------

#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

TPZSloanRenumbering::TPZSloanRenumbering(long NElements, long NNodes):TPZRenumbering(NElements,NNodes){

}

TPZSloanRenumbering::~TPZSloanRenumbering(){

}

void TPZSloanRenumbering::Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter){

  //computing graph
  TPZCutHillMcKee::SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();

  //finding pseudo peripheral nodes
  const long nnodes = graph.NNodes();
  long startNode = -1, endNode = -1;
  graph.PseudoPeripheralNodes(startNode, endNode);

  //computing distances to end node
  TPZStack< TPZVec<long> > LevelStructure;
  graph.RootedLevelStructure(endNode, LevelStructure);
  TPZVec<long> DistanceToEndNode( nnodes, -1 );
  for(long ilevel = 0; ilevel < LevelStructure.NElements(); ilevel++){
    for(long iel = 0; iel < LevelStructure[ ilevel ].NElements(); iel++){
      const long node = LevelStructure[ ilevel ][ iel ];
      DistanceToEndNode[ node ] = ilevel;
    }//for iel
  }//for ilevel
#ifdef DEBUG
  for(long i = 0; i < DistanceToEndNode.NElements(); i++){
    if(DistanceToEndNode[i] == -1) DebugStop();
  }
#endif

  //assigning initial status and priority
  TPZVec<long> priority( nnodes, -1);
  for(long i = 0; i < nnodes; i++){
    const long degree = graph.Degree(i);
    priority[i] = this->W1()*DistanceToEndNode[i] - this->W2()*(degree+1);
  }//for i

  TPZVec<long> Status(nnodes,EInactive);
  TPZStack<long> R;
  TPZManVector<long,10000> adjNodes, adjNodes2J;
  std::list<long> Q;
  TPZVec<long> HistoryOfQ(nnodes,0);
  Status[startNode] = EPreActive;
  Q.push_back( startNode );
  HistoryOfQ[ startNode ] = 1;
  while(Q.size()){
    const long currnode = this->FindHighestPriority(Q,priority);
    Q.remove(currnode);
    if(Status[currnode] == EPreActive){
      graph.AdjacentNodes(currnode, adjNodes);
      for(long jadj = 0; jadj < adjNodes.NElements(); jadj++){
        const long adjNode = adjNodes[jadj];
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        if(HistoryOfQ[adjNode] == 0 && Status[ adjNode] == EInactive){
          Status[ adjNode] = EPreActive;
          Q.push_back( adjNode );
          HistoryOfQ[ adjNode ] = 1;
        }
      }//jadj
    }//not a preactive node
    R.Push( currnode );
    Status[ currnode ] = EPostActive;

    //step 9: updating priorities and queue
    graph.AdjacentNodes(currnode, adjNodes);
    for(long jadj = 0; jadj < adjNodes.NElements(); jadj++){
      const long adjNode = adjNodes[jadj];
      if( Status[ adjNode ] == EPreActive ){
        Status[ adjNode ] = EActive;
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        graph.AdjacentNodes(adjNode, adjNodes2J);
        for(long k = 0; k < adjNodes2J.NElements(); k++){
          const long knode = adjNodes2J[k];
          if(Status[knode] != EPostActive){
            priority[knode] += this->W2();
          }//if ! post active
          if( HistoryOfQ[knode] == 0 && Status[knode] == EInactive){
            Status[knode] = EPreActive;
            Q.push_back( knode );
            HistoryOfQ[ knode ] = 1;
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
  std::set<long> check;
  for(long i = 0; i < nnodes; i++) check.insert(R[i]);
  if( ((long)(check.size())) != nnodes) DebugStop();
}
#endif

  permScatter = R;
  permGather.Resize(nnodes);
  for(long i = 0; i < nnodes; i++) permGather[ permScatter[i] ] = i;

}//void

long TPZSloanRenumbering::FindHighestPriority(const std::list<long> &Q, const TPZVec<long> &priority) const{
  if(Q.size() == 0) DebugStop();
  std::list<long>::const_iterator w = Q.begin(), e = Q.end();
  long result = *w;
  long nodepriority = priority[result];
  for( ; w != e; w++){
    long lcnode = *w;
    if(priority[lcnode] > nodepriority){
      result = lcnode;
      nodepriority = priority[lcnode];
    }
  }
  return result;

}//long



