//---------------------------------------------------------------------------

#include "TPZSloanRenumbering.h"
#include "TPZCutHillMcKee.h"

TPZSloanRenumbering::TPZSloanRenumbering(int64_t NElements, int64_t NNodes):TPZRenumbering(NElements,NNodes){

}

TPZSloanRenumbering::TPZSloanRenumbering():TPZRenumbering(){
}


TPZSloanRenumbering::~TPZSloanRenumbering(){

}

//#define DEBUG_SLOAN_RENUMBERING
#ifdef PZDEBUG_SLOAN_RENUMBERING
#include <fstream>
std::ofstream myfile("c:\\Temp\\sloanrenumbering.txt");
#endif
void TPZSloanRenumbering::Resequence(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter){

#ifdef PZDEBUG_SLOAN_RENUMBERING
   	time_t StartTime = time(NULL);
#endif


  //computing graph
  TPZCutHillMcKee::SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();

#ifdef PZDEBUG_SLOAN_RENUMBERING
  {
   	time_t finalTime = time(NULL);
  	const double elapsedtime = difftime(finalTime, StartTime);
    myfile << "Tempo ConvertGraph = " << elapsedtime << "\n";
    myfile.flush();
    StartTime = time(NULL);
  }
#endif

  //finding pseudo peripheral nodes
  const int64_t nnodes = graph.NNodes();
  int64_t startNode = -1, endNode = -1;
  graph.PseudoPeripheralNodes(startNode, endNode);

  //computing distances to end node
  TPZStack< TPZStack<int64_t> > LevelStructure;
  graph.RootedLevelStructure(endNode, LevelStructure);
  TPZManVector<int64_t> DistanceToEndNode( nnodes, -1 );
  int64_t nelslevel = LevelStructure.NElements();
  for(int64_t ilevel = 0; ilevel < nelslevel; ilevel++){
      int64_t nelsilevel= LevelStructure[ ilevel ].NElements();
    for(int64_t iel = 0; iel < nelsilevel; iel++){
      const int64_t node = LevelStructure[ ilevel ][ iel ];
      DistanceToEndNode[ node ] = ilevel;
    }//for iel
  }//for ilevel
#ifdef PZDEBUG
  for(int64_t i = 0; i < DistanceToEndNode.NElements(); i++){
    if(DistanceToEndNode[i] == -1){
      DebugStop();
    }
  }
#endif

//    std::cout << "Computed DistanceToEndNode\n";
  //assigning initial status and priority
    int64_t maxpriority = 0;
  TPZVec<int64_t> priority( nnodes, -1);
  for(int64_t i = 0; i < nnodes; i++){
    const int64_t degree = graph.Degree(i);
    priority[i] = this->W1()*DistanceToEndNode[i] - this->W2()*(degree+1);
      maxpriority = MAX(maxpriority,priority[i]);
  }//for i

  TPZVec<int64_t> Status(nnodes,EInactive);
  TPZStack<int64_t> R;
  R.Resize(nnodes);//pre allocating memory
  R.Resize(0);

//  TPZManVector<int64_t,10000> adjNodes, adjNodes2J;
  SList Q(nnodes);
  Status[startNode] = EPreActive;
  Q.push_back( startNode );
  while(Q.fSize){
    int64_t index;
      maxpriority = 0;
    const int64_t currnode = this->FindHighestPriority(Q,priority,index);
      maxpriority = MAX(maxpriority,priority[currnode]);
    Q.remove(index);
    if(Status[currnode] == EPreActive){
      int64_t nadjNodes;
      int64_t * adjNodePtr = graph.AdjacentNodesPtr(currnode, nadjNodes);
      for(int64_t jadj = 0; jadj < nadjNodes; jadj++){
        const int64_t adjNode = adjNodePtr[jadj];
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
        if(Q.push_back(adjNode) == true && Status[adjNode] == EInactive){
          Status[ adjNode] = EPreActive;
            maxpriority = MAX(maxpriority,priority[adjNode]);
          //ja foi inserido no if acima Q.push_back( adjNode );
        }
      }//jadj
    }//not a preactive node
    R.Push( currnode );
    Status[ currnode ] = EPostActive;

    //step 9: updating priorities and queue
    int64_t nadjNodes;
    int64_t * adjNodePtr = graph.AdjacentNodesPtr(currnode, nadjNodes);
    for(int64_t jadj = 0; jadj < nadjNodes; jadj++){
      const int64_t adjNode = adjNodePtr[jadj];
      if( Status[ adjNode ] == EPreActive ){
        Status[ adjNode ] = EActive;
        priority[ adjNode ] += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
          maxpriority = MAX(maxpriority,priority[adjNode]);
        int64_t nadjNodes2J;
        int64_t * adjNode2JPtr = graph.AdjacentNodesPtr(adjNode, nadjNodes2J);
        for(int64_t k = 0; k < nadjNodes2J; k++){
          const int64_t knode = adjNode2JPtr[k];
          if(Status[knode] != EPostActive){
              int pr = priority[knode];
            priority[knode] = pr+this->W2();
              maxpriority = MAX(maxpriority,priority[knode]);
          }//if ! post active
          if( Q.push_back(knode) == true && Status[knode] == EInactive){
            Status[knode] = EPreActive;
            //ja foi inserido no if acima Q.push_back( knode );
          }// if inactive
        }//loop of k nodes
      }//if
    }//for jadj

    if(R.NElements() % 100000 == 0){
        std::cout << R.NElements() << "\tQ.size = " << Q.fSize << ", %done = " << 100.*R.NElements()/nnodes << " maxpriority = " << maxpriority << "\n";
    }  
  }//while

  if(R.NElements() != nnodes){
    std::cout << "TPZSloanRenumbering::Resequence error! Not all nodes were processed.\n"
              << "It may be an implementation bug or simply a graph with independent nodes\n"
              << "as if the mesh has subdomains that are not connected.\n"
              << "For the later, use CuttHillMcKee, which is prepared for this particular case\n";
    DebugStop();
  }

#ifdef PZDEBUG
{//verificando se ha duplicados
  std::set<int64_t> check;
  for(int64_t i = 0; i < nnodes; i++) check.insert(R[i]);
  if( ((int64_t)(check.size())) != nnodes) DebugStop();
}
#endif

  permScatter = R;
  permGather.Resize(nnodes);
  for(int64_t i = 0; i < nnodes; i++) permGather[ permScatter[i] ] = i;

#ifdef PZDEBUG_SLOAN_RENUMBERING
  {
   	time_t finalTime = time(NULL);
  	const double elapsedtime = difftime(finalTime, StartTime);
    myfile << "Tempo restante = " << elapsedtime << "\n\n\n";
    myfile.flush();
    StartTime = time(NULL);
  }
#endif


}//void

int64_t TPZSloanRenumbering::FindHighestPriority(const SList &Q, const TPZVec<int64_t> &priority, int64_t &Qindex) const{
  if(Q.fSize == 0) DebugStop();
  int64_t result = Q.fList[0];
  int64_t nodepriority = priority[result];
  Qindex = 0;
  for(int64_t i = 1; i < Q.fSize; i++){
    const int64_t lcnode = Q.fList[i];
    if(priority[lcnode] > nodepriority){
      result = lcnode;
      nodepriority = priority[lcnode];
      Qindex = i;
    }
  }
  return result;

}//int64_t


void TPZSloanRenumbering::Resequence2(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter)
{
    //computing graph
    TPZCutHillMcKee::SGraph graph;
    ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
    graph.fnodegraph.Shrink();
    graph.fnodegraphindex.Shrink();
    
    //finding pseudo peripheral nodes
    const int64_t nnodes = graph.NNodes();
    
    fAllNodes.Resize(nnodes);
    for (int64_t i=0; i<nnodes; i++) {
        fAllNodes[i].fIndex = i;
    }
    
    
    int64_t startNode = -1, endNode = -1;
    graph.PseudoPeripheralNodes(startNode, endNode);
    
    //computing distances to end node
    TPZStack< TPZStack<int64_t> > LevelStructure;
    graph.RootedLevelStructure(endNode, LevelStructure);
    TPZVec<int> DistanceToEndNode( nnodes, -1 );
    for(int64_t ilevel = 0; ilevel < LevelStructure.NElements(); ilevel++){
        for(int64_t iel = 0; iel < LevelStructure[ ilevel ].NElements(); iel++){
            const int64_t node = LevelStructure[ ilevel ][ iel ];
            DistanceToEndNode[ node ] = ilevel;
        }//for iel
    }//for ilevel
#ifdef PZDEBUG
    for(int64_t i = 0; i < DistanceToEndNode.NElements(); i++){
        if(DistanceToEndNode[i] == -1){
            DebugStop();
        }
    }
#endif
    
//    std::cout << "Computed DistanceToEndNode\n";
    //assigning initial status and priority
    for(int64_t i = 0; i < nnodes; i++){
        const int degree = graph.Degree(i);
        fAllNodes[i].fPriority = this->W1()*DistanceToEndNode[i] - this->W2()*(degree+1);
    }//for i
    
    int64_t countnumbered = 0;
    int64_t numactive = 1;
    fAllNodes[startNode].fStatus = EPreActive;
    InsertNode(&fAllNodes[startNode]);

    while (fActive.size())
    {
        TNo *worknode = PopHighestPriorityNode();
        numactive--;
        
        if (worknode->fStatus == EPreActive)
        {
            int64_t nadjNodes;
            int64_t * adjNodePtr = graph.AdjacentNodesPtr(worknode->fIndex, nadjNodes);
            for(int64_t jadj = 0; jadj < nadjNodes; jadj++){
                const int64_t adjNode = adjNodePtr[jadj];
                if (fAllNodes[adjNode].fStatus == EInactive) {
                    fAllNodes[adjNode].fPriority += this->W2();
                    fAllNodes[adjNode].fStatus = EPreActive;
                    InsertNode(&fAllNodes[adjNode]);
                    numactive++;
                }
                else if(fAllNodes[adjNode].fStatus != EPostActive)
                {
                    int priority = fAllNodes[adjNode].fPriority;
                    priority += this->W2();
                    TransferPriority(&fAllNodes[adjNode],priority);
                }
            }
        }
        worknode->fStatus = EPostActive;
        worknode->fSequenceNumber = countnumbered;
        countnumbered++;

        //step 9: updating priorities and queue
        int64_t nadjNodes;
        int64_t * adjNodePtr = graph.AdjacentNodesPtr(worknode->fIndex, nadjNodes);
        for(int64_t jadj = 0; jadj < nadjNodes; jadj++){
            const int64_t adjNode = adjNodePtr[jadj];
            if( fAllNodes[ adjNode ].fStatus == EPreActive ){
                fAllNodes[ adjNode ].fStatus = EActive;
                int priority = fAllNodes[adjNode].fPriority;
                priority += this->W2();//no artigo de 86 o sloan soma W1, no de de 89 soma W2. To seguindo a notacao de 89
                TransferPriority(&fAllNodes[adjNode], priority);
                int64_t nadjNodes2J;
                int64_t * adjNode2JPtr = graph.AdjacentNodesPtr(adjNode, nadjNodes2J);
                for(int64_t k = 0; k < nadjNodes2J; k++){
                    const int64_t knode = adjNode2JPtr[k];
                    if(fAllNodes[knode].fStatus != EPostActive){
                        int priority = fAllNodes[knode].fPriority;
                        priority += this->W2();
                        if (fAllNodes[knode].fStatus == EInactive) {
                            fAllNodes[knode].fPriority = priority;
                            InsertNode(&fAllNodes[knode]);
                            fAllNodes[knode].fStatus = EPreActive;
                            numactive++;
                        }
                        else
                        {
                            TransferPriority(&fAllNodes[knode], priority);
                        }
                    }//if ! post active
                }//loop of k nodes
            }//if
        }//for jadj

        if(countnumbered % 1000000 == 0){
            std::cout << countnumbered << "\tActive size = " << numactive << ", %done = " << 100.*countnumbered/nnodes << "\n";
        }

    }
    std::cout << "number of elements sequenced " << countnumbered << std::endl;
}

