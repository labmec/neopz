//---------------------------------------------------------------------------

#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

TPZSloanRenumbering::TPZSloanRenumbering(long NElements, long NNodes):TPZRenumbering(NElements,NNodes){

}

TPZSloanRenumbering::~TPZSloanRenumbering(){

}

//#define DEBUG_SLOAN_RENUMBERING
#ifdef DEBUG_SLOAN_RENUMBERING
#include <fstream>
std::ofstream myfile("c:\\Temp\\sloanrenumbering.txt");
#endif
void TPZSloanRenumbering::Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter){

#ifdef DEBUG_SLOAN_RENUMBERING
   	time_t StartTime = time(NULL);
#endif


  //computing graph
  TPZCutHillMcKee::SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();

#ifdef DEBUG_SLOAN_RENUMBERING
  {
   	time_t finalTime = time(NULL);
  	const double elapsedtime = difftime(finalTime, StartTime);
    myfile << "Tempo ConvertGraph = " << elapsedtime << "\n";
    myfile.flush();
    StartTime = time(NULL);
  }
#endif

  //finding pseudo peripheral nodes
  const long nnodes = graph.NNodes();
  long startNode = -1, endNode = -1;
  graph.PseudoPeripheralNodes(startNode, endNode);

  //computing distances to end node
  TPZStack< TPZStack<long> > LevelStructure;
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
    if(DistanceToEndNode[i] == -1){
      DebugStop();
    }
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
  R.Resize(nnodes);//pre allocating memory
  R.Resize(0);

//  TPZManVector<long,10000> adjNodes, adjNodes2J;
  SList Q(nnodes);
  Status[startNode] = EPreActive;
  Q.push_back( startNode );
  while(Q.fSize){
    long index;
    const long currnode = this->FindHighestPriority(Q,priority,index);
    Q.remove(index);
    if(Status[currnode] == EPreActive){
      long nadjNodes;
      long * adjNodePtr = graph.AdjacentNodesPtr(currnode, nadjNodes);
      for(long jadj = 0; jadj < nadjNodes; jadj++){
        const long adjNode = adjNodePtr[jadj];
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        if(Q.push_back(adjNode) == true && Status[adjNode] == EInactive){
          Status[ adjNode] = EPreActive;
          //ja foi inserido no if acima Q.push_back( adjNode );
        }
      }//jadj
    }//not a preactive node
    R.Push( currnode );
    Status[ currnode ] = EPostActive;

    //step 9: updating priorities and queue
    long nadjNodes;
    long * adjNodePtr = graph.AdjacentNodesPtr(currnode, nadjNodes);
    for(long jadj = 0; jadj < nadjNodes; jadj++){
      const long adjNode = adjNodePtr[jadj];
      if( Status[ adjNode ] == EPreActive ){
        Status[ adjNode ] = EActive;
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        long nadjNodes2J;
        long * adjNode2JPtr = graph.AdjacentNodesPtr(adjNode, nadjNodes2J);
        for(long k = 0; k < nadjNodes2J; k++){
          const long knode = adjNode2JPtr[k];
          if(Status[knode] != EPostActive){
            priority[knode] += this->W2();
          }//if ! post active
          if( Q.push_back(knode) == true && Status[knode] == EInactive){
            Status[knode] = EPreActive;
            //ja foi inserido no if acima Q.push_back( knode );
          }// if inactive
        }//loop of k nodes
      }//if
    }//for jadj

   /* if(R.NElements() % 50 == 0){
      myfile << R.NElements() << "\tQ.size = " << Q.size() << ", nnodes = " << nnodes << "\n";
      myfile.flush();
    }  */
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

#ifdef DEBUG_SLOAN_RENUMBERING
  {
   	time_t finalTime = time(NULL);
  	const double elapsedtime = difftime(finalTime, StartTime);
    myfile << "Tempo restante = " << elapsedtime << "\n\n\n";
    myfile.flush();
    StartTime = time(NULL);
  }
#endif


}//void

long TPZSloanRenumbering::FindHighestPriority(const SList &Q, const TPZVec<long> &priority, long &Qindex) const{
  if(Q.fSize == 0) DebugStop();
  long result = Q.fList[0];
  long nodepriority = priority[result];
  Qindex = 0;
  for(long i = 1; i < Q.fSize; i++){
    const long lcnode = Q.fList[i];
    if(priority[lcnode] > nodepriority){
      result = lcnode;
      nodepriority = priority[lcnode];
      Qindex = i;
    }
  }
  return result;

}//long



