
#include "TPZCutHillMcKee.h"
#include "fstream"

void TPZCutHillMcKee::SGraph::Set2Vec(const std::set<int> &myset,
                                      TPZVec<int> &myVec) const{
  const int n = myset.size();
  myVec.Resize(n);
  std::set<int>::const_iterator w, e = myset.end();
  int i;
  for(w = myset.begin(), i = 0; w != e; w++, i++){
    myVec[i] = *w;
  }
}//void

int TPZCutHillMcKee::SGraph::SmallestDegree(TPZVec<int> &ExploredNodes){
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

void TPZCutHillMcKee::SGraph::GetAdjacentNodes(const TPZVec<int> &parents,
                                               const TPZVec< int > &exceptedNodes,
                                               std::set<int> &adjNodes){
  adjNodes.clear();
  TPZManVector<int,1000> localAdjNodes;
  for(int i = 0; i < parents.NElements(); i++){
    this->AdjacentNodes(parents[i],localAdjNodes);
    for(int j = 0; j < localAdjNodes.NElements(); j++){
      const int node = localAdjNodes[j];
      if(exceptedNodes[node] == 0){
        adjNodes.insert(node);
      }
    }//j
  }//i
}//void

void TPZCutHillMcKee::SGraph::AdjacentNodes(int parent, TPZVec<int> &adjNodes){
  const int n = this->Degree(parent);
  adjNodes.Resize(n);
  for(int i = 0; i < n; i++){
    const int adj = fnodegraph[fnodegraphindex[parent]+i];
    adjNodes[i] = adj;
  }//for
}//void

void TPZCutHillMcKee::SGraph::AdjacentNodesOrdered(int parent, TPZVec<int> &adjNodes){
  std::multimap<int,int> order;

  const int n = this->Degree(parent);
  for(int i = 0; i < n; i++){
    const int adj = fnodegraph[fnodegraphindex[parent]+i];
    const int adjDegree = this->Degree(adj);
    order.insert(std::make_pair(adjDegree,adj));
  }//for
  adjNodes.Resize(n);
#ifdef DEBUG
  if(n!= (int)(order.size()) ){
    DebugStop();
  }
#endif
  std::multimap<int,int>::const_iterator mapIt;
  int i;
  for(i = 0, mapIt = order.begin(); mapIt != order.end(); i++, mapIt++){
    adjNodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::RootedLevelStructure(int rootNode,
                                                   TPZStack< TPZVec<int> > &LevelStructure){
  LevelStructure.Resize(this->NNodes());//preallocating memory
  LevelStructure.Resize(0);

  TPZManVector<int,1000> adjNodes;
  TPZManVector<int,1000> thisLevel(1);
  std::set<int> SetOfAdjNodes;
  TPZVec<int> ProcessedNodes(this->NNodes(),0);
  thisLevel[0] = rootNode;
  for(int iLevel = 0; iLevel < this->NNodes(); iLevel++){

    const int nThisLevel = thisLevel.NElements();
    if(nThisLevel == 0) break;
    for(int i = 0; i < nThisLevel; i++){
      const int node = thisLevel[i];
      ProcessedNodes[node] = 1;
    }
    LevelStructure.Push( thisLevel);

    this->GetAdjacentNodes(thisLevel, ProcessedNodes, SetOfAdjNodes);
    this->Set2Vec(SetOfAdjNodes,thisLevel);

  }//for que de fato eh while

  LevelStructure.Shrink();//saving memory

}//void

void TPZCutHillMcKee::SGraph::SortNodes(TPZVec<int> &nodes){
  std::multimap<int,int> mymap;
  for(int i = 0; i < nodes.NElements(); i++){
    const int node = nodes[i];
    const int degree = this->Degree(node);
    mymap.insert(std::make_pair(degree,node));
  }//for
  std::multimap<int,int>::const_iterator mapIt;
  int i;
  for(i = 0, mapIt = mymap.begin(); mapIt != mymap.end(); i++, mapIt++){
    nodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::ShrinkLastLevel(TPZVec<int> &LastLevel){
  //accorging to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  //Previously, Sloan had suggested LastLevel.Resize( (LastLevel.NElements()+2)/2 );
  std::map< int, int > mymap;
  for(int i = 0; i < LastLevel.NElements(); i++){
    const int node = LastLevel[i];
    const int degree = this->Degree( node );
    mymap[ degree ] = node;
  }
  LastLevel.Resize(mymap.size());
  std::map< int, int >::const_iterator w, e = mymap.end();
  int i;
  for(i = 0, w = mymap.begin(); w != e; w++, i++){
    LastLevel[i] = w->second;
  }

}

void TPZCutHillMcKee::SGraph::PseudoPeripheralNodes(int &startNode, int &endNode){

  endNode = -1;

  TPZVec<int> emptyVec(this->NNodes(), 0);
  //step 1: first guess
  startNode = SmallestDegree(emptyVec);

  //step 2: rooted level structure
  TPZStack< TPZVec<int> > LevelStructure;
  this->RootedLevelStructure(startNode, LevelStructure);

#ifdef DEBUG_CM
  std::ofstream myfile("c:\\Temp\\PseudoPeripheralNodes.txt");
#endif

  int count = 0;
  while(1){
    count++;
    const int cpStart = startNode;
    const int cpEnd = endNode;
#ifdef DEBUG_CM
    myfile << "\nstart = " << startNode << "  end = " << endNode << "\n\n";
#endif
    const int nlevels = LevelStructure.NElements();
    TPZManVector<int,1000> LastLevel = LevelStructure[ nlevels-1 ];
    //step 3 - sort the last level
    this->SortNodes(LastLevel);
    //step 4 - shrink the last level
    ShrinkLastLevel(LastLevel);
    //step 5
    int we = 1e9;
    int hs = nlevels;
    for(int iQ = 0; iQ < LastLevel.NElements(); iQ++){
      TPZStack< TPZVec<int> > localLevelStructure;
      const int nodeAtQ = LastLevel[iQ];
#ifdef DEBUG_CM
      myfile << count << "\t" << iQ << "\t";myfile.flush();
#endif
      this->RootedLevelStructure(nodeAtQ, localLevelStructure);
#ifdef DEBUG_CM
      myfile << "rooted\t";myfile.flush();
#endif
      const int h = localLevelStructure.NElements();
      int w = 0;
      for(int j = 0; j < localLevelStructure.NElements(); j++){
        const int localW = localLevelStructure[j].NElements();
        if(localW > w) w = localW;
      }
      if(h > hs && w < we){
        startNode = nodeAtQ;
        LevelStructure = localLevelStructure;
#ifdef DEBUG_CM
        myfile << "break\n";myfile.flush();
#endif
        break;
      }
      else if(w < we){
        endNode = nodeAtQ;
        we = w;
      }
#ifdef DEBUG_CM
      myfile << "end\n";myfile.flush();
#endif
    }//for iStep
#ifdef DEBUG_CM
    myfile << "\n"; myfile.flush();
#endif

    if(cpStart == startNode && cpEnd == endNode){
      break;
    }
  }//while

}//void

////////////////////////////////////////////////////////////////////////////////

TPZCutHillMcKee::TPZCutHillMcKee():TPZRenumbering(){
  fReverse = true;
#ifdef DEBUG
  fVerbose = true;
#else
  fVerbose = false;
#endif
}

TPZCutHillMcKee::TPZCutHillMcKee(int NElements, int NNodes, bool Reverse):TPZRenumbering(NElements,NNodes){
  fReverse = Reverse;
#ifdef DEBUG
  fVerbose = true;
#else
  fVerbose = false;
#endif
}

TPZCutHillMcKee::~TPZCutHillMcKee(){

}


void TPZCutHillMcKee::Resequence(TPZVec<int> &perm, TPZVec<int> &iperm){
  TPZVec<int> permGather, permScatter;
  TPZVec<int> permGatherReverse, permScatterReverse;

  this->Resequence(permGather, permScatter, permGatherReverse, permScatterReverse);

  if(fReverse){
    perm = permGatherReverse;
    iperm = permScatterReverse;
  }
  else{
    perm = permGather;
    iperm = permScatter;
  }
}//void

void TPZCutHillMcKee::Resequence(TPZVec<int> &permGather, TPZVec<int> &permScatter,
                                 TPZVec<int> &permGatherReverse, TPZVec<int> &permScatterReverse){

  if(this->fVerbose) {
    std::cout << "TPZCutHillMcKee ConvertGraph...";
    std::cout.flush();
  }
  SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();


#ifdef DEBUG_CM
  if(0){
    int grafoindex[7] = {0,3,6,9,14,19,22};
    int grafo[22] = {3,4,5,  2,3,4,  1,3,4,  0,1,2,4,5,  0,1,2,3,5,  0,3,4};
    graph.fnodegraphindex.Resize(7);
    for(int i = 0; i < 7; i++) graph.fnodegraphindex[i] = grafoindex[i];
    graph.fnodegraph.Resize(22);
    for(int i = 0; i < 22; i++) graph.fnodegraph[i] = grafo[i];

    {
      TPZStack< TPZVec<int> > LevelStructure;
      graph.RootedLevelStructure(1, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructure.txt");
      for(int iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(int j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

    {
      int startNode = -1, endNode = -1;
      graph.PseudoPeripheralNodes(startNode, endNode);
      TPZStack< TPZVec<int> > LevelStructure;
      graph.RootedLevelStructure(startNode, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructureAfterPeripheral.txt");
      for(int iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(int j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

  }//if
#endif

        //CutHillMcKee
  std::queue<int> Q;
  TPZStack<int> R;
  TPZManVector<int> adjNodes;
  //reservando memoria em R
  const int nnodes = graph.NNodes();
  R.Resize(nnodes);
  R.Resize(0);
  TPZVec<int> ExploredNodes(nnodes,0);
  if(this->fVerbose) {
    std::cout << "TPZCutHillMcKee process...\n";std::cout.flush();
  }

  bool firstTime = true;
  while(R.NElements() != graph.NNodes()){
    if(this->fVerbose){
      std::cout << "Calling TPZCutHillMcKee::SGraph::SmallestDegree\n";
      std::cout.flush();
    }
    int Parent = -1;
    if(firstTime){
      int startNode = -1, endNode = -1;
      graph.PseudoPeripheralNodes(startNode, endNode);
      Parent = startNode;
      firstTime = false;
    }
    else Parent = graph.SmallestDegree(ExploredNodes);

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

  const int n = R.NElements();
  if(n != nnodes) DebugStop();

#ifdef DEBUG
{//verificando se ha duplicados
  std::set<int> check;
  for(int i = 0; i < n; i++) check.insert(R[i]);
  if( ((int)(check.size())) != n) DebugStop();
}
#endif

  if(this->fVerbose){
    std::cout << "TPZCutHillMcKee Filling perm and iperm vectors...";std::cout.flush();
  }


  permScatterReverse.Resize(n);
  for(int i = 0; i < n; i++) permScatterReverse[n-1-i] = R[i];
  permGatherReverse.Resize(n);
  for(int i = 0; i < n; i++) permGatherReverse[ permScatterReverse[i] ] = i;


  permScatter = R;
  permGather.Resize(n);
  for(int i = 0; i < n; i++) permGather[ permScatter[i] ] = i;


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

