
#include "TPZCutHillMcKee.h"
#include "fstream"

void TPZCutHillMcKee::SGraph::Set2Vec(const std::set<int64_t> &myset,
                                      TPZVec<int64_t> &myVec) const{
  const int64_t n = myset.size();
  myVec.Resize(n);
  std::set<int64_t>::const_iterator w, e = myset.end();
  int64_t i;
  for(w = myset.begin(), i = 0; w != e; w++, i++){
    myVec[i] = *w;
  }
}//void

int64_t TPZCutHillMcKee::SGraph::SmallestDegree(TPZVec<int64_t> &ExploredNodes){
//    std::cout << "TPZCutHillMcKee::SGraph::SmallestDegree size = " << this->NNodes() << std::endl;
//    std::cout.flush();
  int64_t mindegree = 1000000;
  int64_t found = -1;
  const int64_t n = this->NNodes();
  for(int64_t i = 0; i < n; i++){
    if(ExploredNodes[i] == 1) continue;
    int64_t degree = this->Degree(i);
    if(degree < mindegree){
      mindegree = degree;
      found = i;
    }
  }//for
  return found;
}


void TPZCutHillMcKee::SGraph::GetAdjacentNodes(const TPZVec<int64_t> &parents,
                                               const TPZVec< int64_t > &exceptedNodes,
                                               TPZStack<int64_t> &adjNodes){
  adjNodes.Resize(this->NNodes());//pre allocating memory
  int64_t count = 0;
  TPZVec<int64_t> History(this->NNodes(),0); //avoiding duplicates in adjNodes. it could be a std::set. set is better when few elements
  for(int64_t i = 0; i < parents.NElements(); i++){
    int64_t nlocal;
    int64_t * localAdjNodes = AdjacentNodesPtr(parents[i], nlocal);
    for(int64_t j = 0; j < nlocal; j++){
      const int64_t node = localAdjNodes[j];
      if(exceptedNodes[node] == 0 && History[node] == 0){
        adjNodes[count] = node;
        count++;
        History[node] = 1;
      }
    }//j
  }//i
  adjNodes.Resize(count);
  adjNodes.Shrink(); //saving memory
}//void

void TPZCutHillMcKee::SGraph::AdjacentNodes(int64_t parent, TPZVec<int64_t> &adjNodes){
  const int64_t n = this->Degree(parent);
  adjNodes.Resize(n);
  for(int64_t i = 0; i < n; i++){
    const int64_t adj = fnodegraph[fnodegraphindex[parent]+i];
    adjNodes[i] = adj;
  }//for
}//void

void TPZCutHillMcKee::SGraph::AdjacentNodesOrdered(int64_t parent,
                                                   const TPZVec<int64_t> &exceptedNodes,
                                                   TPZVec<int64_t> &adjNodes){
  std::multimap<int64_t,int64_t> order;

  const int64_t ndegree = this->Degree(parent);
  int64_t count = 0;
  for(int64_t i = 0; i < ndegree; i++){
    const int64_t adj = fnodegraph[fnodegraphindex[parent]+i];
    if(exceptedNodes[adj] == 1) continue;
    const int64_t adjDegree = this->Degree(adj);
    order.insert(std::make_pair(adjDegree,adj));
    count++;
  }//for

  adjNodes.Resize(count);
#ifdef PZDEBUG
  if(count != (int64_t)(order.size()) ){
    DebugStop();
  }
#endif
  std::multimap<int64_t,int64_t>::const_iterator mapIt;
  int64_t i;
  for(i = 0, mapIt = order.begin(); mapIt != order.end(); i++, mapIt++){
    adjNodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::RootedLevelStructure(int64_t rootNode,
                                                   TPZStack< TPZStack<int64_t> > &LevelStructure){
  LevelStructure.Resize(0);

  TPZManVector<int64_t,1000> adjNodes;

  TPZStack<int64_t> thisLevel;
  thisLevel.Resize(1);
  thisLevel[0] = rootNode;

  TPZVec<int64_t> ProcessedNodes(this->NNodes(),0);

  for(int64_t iLevel = 0; iLevel < this->NNodes(); iLevel++){

    const int64_t nThisLevel = thisLevel.NElements();
    if(nThisLevel == 0) break;
    for(int64_t i = 0; i < nThisLevel; i++){
      const int64_t node = thisLevel[i];
      ProcessedNodes[node] = 1;
    }
    LevelStructure.Push( thisLevel);

    this->GetAdjacentNodes(LevelStructure[LevelStructure.NElements() - 1], ProcessedNodes, thisLevel);

  }//for que de fato eh while

  LevelStructure.Shrink();//saving memory

}//void

void TPZCutHillMcKee::SGraph::SortNodes(TPZVec<int64_t> &nodes){
  std::multimap<int64_t,int64_t> mymap;
  for(int64_t i = 0; i < nodes.NElements(); i++){
    const int64_t node = nodes[i];
    const int64_t degree = this->Degree(node);
    mymap.insert(std::make_pair(degree,node));
  }//for
  std::multimap<int64_t,int64_t>::const_iterator mapIt;
  int64_t i;
  for(i = 0, mapIt = mymap.begin(); mapIt != mymap.end(); i++, mapIt++){
    nodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::ShrinkLastLevel(TPZVec<int64_t> &LastLevel){

  const int64_t nelsOrig = LastLevel.NElements();

  //accorging to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  std::map< int64_t, int64_t > mymap;
  for(int64_t i = 0; i < LastLevel.NElements(); i++){
    const int64_t node = LastLevel[i];
    const int64_t degree = this->Degree( node );
    mymap[ degree ] = node;
  }
  LastLevel.Resize(mymap.size());
  std::map< int64_t, int64_t >::const_iterator w, e = mymap.end();
  int64_t i;
  for(i = 0, w = mymap.begin(); w != e; w++, i++){
    LastLevel[i] = w->second;
  }

  //Previously, Sloan had suggested, in the article of 1986, LastLevel.Resize( (LastLevel.NElements()+2)/2 );
  const int64_t newsize = (nelsOrig+2)/2;
  if(newsize < LastLevel.NElements()){
    LastLevel.Resize( newsize );
  }

}

void TPZCutHillMcKee::SGraph::PseudoPeripheralNodes(int64_t &startNode, int64_t &endNode){

  endNode = -1;

  TPZManVector<int64_t> emptyVec(this->NNodes(), 0);
  //step 1: first guess
  startNode = SmallestDegree(emptyVec);

  //step 2: rooted level structure
  TPZStack< TPZStack<int64_t> > LevelStructure;
  this->RootedLevelStructure(startNode, LevelStructure);

#ifdef PZDEBUG_CM
  std::ofstream myfile("c:\\Temp\\PseudoPeripheralNodes.txt");
#endif

  int64_t count = 0;
  const int64_t maxCount = 10;
  TPZStack< TPZStack<int64_t> > localLevelStructure;
  TPZStack<int64_t> LastLevel;
  LastLevel.Resize( this->NNodes() );//pre allocating memory
  LastLevel.Resize(0);
  while(count < maxCount || endNode == -1){ //maxCount pra não ficar aqui eternamente

    count++;
    const int64_t cpStart = startNode;
    const int64_t cpEnd = endNode;
#ifdef PZDEBUG_CM
    myfile << "\nstart = " << startNode << "  end = " << endNode << "\n\n";
#endif
    const int64_t nlevels = LevelStructure.NElements();
      if(!nlevels) break;
    LastLevel = LevelStructure[ nlevels-1 ];
    //step 3 - sort the last level
    this->SortNodes(LastLevel);
    //step 4 - shrink the last level
    ShrinkLastLevel(LastLevel);

    //step 5
    int64_t we = 1000000000L;
    int64_t hs = nlevels;
    const int64_t nelLastLevel = LastLevel.NElements();
    for(int64_t iQ = 0; iQ < nelLastLevel; iQ++){

      const int64_t nodeAtQ = LastLevel[iQ];
#ifdef PZDEBUG_CM
      myfile << count << "\t" << iQ << "\t";myfile.flush();
#endif
      this->RootedLevelStructure(nodeAtQ, localLevelStructure);
#ifdef PZDEBUG_CM
      myfile << "rooted\t";myfile.flush();
#endif

      const int64_t h = localLevelStructure.NElements();
      int64_t w = 0;
      for(int64_t j = 0; j < h; j++){
        const int64_t localW = localLevelStructure[j].NElements();
        if(localW > w) w = localW;
      }
      if(h > hs && w < we){
        startNode = nodeAtQ;
        LevelStructure = localLevelStructure;
#ifdef PZDEBUG_CM
        myfile << "break\n";myfile.flush();
#endif
        break;
      }
      else if(w < we){
        endNode = nodeAtQ;
        we = w;
      }
#ifdef PZDEBUG_CM
      myfile << "end\n";myfile.flush();
#endif

    }//for iStep
#ifdef PZDEBUG_CM
    myfile << "\n"; myfile.flush();
#endif

    if(cpStart == startNode && cpEnd == endNode){
      break;
    }
  }//while

  if(startNode == -1 || endNode == -1) DebugStop();

}//void

////////////////////////////////////////////////////////////////////////////////

TPZCutHillMcKee::TPZCutHillMcKee():TPZRenumbering(){
  fReverse = true;
#ifdef PZDEBUG
  fVerbose = true;
#else
  fVerbose = false;
#endif
}

TPZCutHillMcKee::TPZCutHillMcKee(int64_t NElements, int64_t NNodes, bool Reverse):TPZRenumbering(NElements,NNodes){
  fReverse = Reverse;
#ifdef PZDEBUG
  fVerbose = true;
#else
  fVerbose = false;
#endif
}

TPZCutHillMcKee::~TPZCutHillMcKee(){

}


void TPZCutHillMcKee::Resequence(TPZVec<int64_t> &perm, TPZVec<int64_t> &iperm){

  TPZVec<int64_t> permGather, permScatter;
  TPZVec<int64_t> permGatherReverse, permScatterReverse;

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

void TPZCutHillMcKee::Resequence(TPZVec<int64_t> &permGather, TPZVec<int64_t> &permScatter,
                                 TPZVec<int64_t> &permGatherReverse, TPZVec<int64_t> &permScatterReverse){

  if(this->fVerbose) {
    std::cout << "TPZCutHillMcKee ConvertGraph...";
    std::cout.flush();
  }
  SGraph graph;
  ConvertGraph(fElementGraph, fElementGraphIndex, graph.fnodegraph, graph.fnodegraphindex);
  graph.fnodegraph.Shrink();
  graph.fnodegraphindex.Shrink();


#ifdef PZDEBUG_CM
  if(0){
    int64_t grafoindex[7] = {0,3,6,9,14,19,22};
    int64_t grafo[22] = {3,4,5,  2,3,4,  1,3,4,  0,1,2,4,5,  0,1,2,3,5,  0,3,4};
    graph.fnodegraphindex.Resize(7);
    for(int64_t i = 0; i < 7; i++) graph.fnodegraphindex[i] = grafoindex[i];
    graph.fnodegraph.Resize(22);
    for(int64_t i = 0; i < 22; i++) graph.fnodegraph[i] = grafo[i];

    {
      TPZStack< TPZStack<int64_t> > LevelStructure;
      graph.RootedLevelStructure(1, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructure.txt");
      for(int64_t iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(int64_t j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

    {
      int64_t startNode = -1, endNode = -1;
      graph.PseudoPeripheralNodes(startNode, endNode);
      TPZStack< TPZStack<int64_t> > LevelStructure;
      graph.RootedLevelStructure(startNode, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructureAfterPeripheral.txt");
      for(int64_t iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(int64_t j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

  }//if
#endif

        //CutHillMcKee
  std::queue<int64_t> Q;
  TPZStack<int64_t> R;
  TPZManVector<int64_t> adjNodes;
  //reservando memoria em R
  const int64_t nnodes = graph.NNodes();
  R.Resize(nnodes);
  R.Resize(0);
  TPZVec<int64_t> ExploredNodes(nnodes,0);
  if(this->fVerbose) {
    std::cout << "TPZCutHillMcKee process...\n";std::cout.flush();
  }

  bool firstTime = true;
  while(R.NElements() != graph.NNodes()){
    if(this->fVerbose){
      std::cout << "Calling TPZCutHillMcKee::SGraph::SmallestDegree\n";
      std::cout.flush();
    }
    int64_t Parent = -1;
    if(firstTime){
      int64_t startNode = -1, endNode = -1;
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
      const int64_t Child = Q.front();
      Q.pop();
      this->ProcessParentNode(Child, graph, ExploredNodes, R, Q, adjNodes);
    }

  }//loop completo

  const int64_t n = R.NElements();
  if(n != nnodes) DebugStop();

#ifdef PZDEBUG
{//verificando se ha duplicados
  std::set<int64_t> check;
  for(int64_t i = 0; i < n; i++) check.insert(R[i]);
  if( ((int64_t)(check.size())) != n) DebugStop();
}
#endif

  if(this->fVerbose){
    std::cout << "TPZCutHillMcKee Filling perm and iperm vectors...";std::cout.flush();
  }


  permScatterReverse.Resize(n);
  for(int64_t i = 0; i < n; i++) permScatterReverse[n-1-i] = R[i];
  permGatherReverse.Resize(n);
  for(int64_t i = 0; i < n; i++) permGatherReverse[ permScatterReverse[i] ] = i;


  permScatter = R;
  permGather.Resize(n);
  for(int64_t i = 0; i < n; i++) permGather[ permScatter[i] ] = i;


}///void

void TPZCutHillMcKee::ProcessParentNode(int64_t Parent,
                                        SGraph &graph,
                                        TPZVec<int64_t> &ExploredNodes,
                                        TPZStack<int64_t> &R,
                                        std::queue<int64_t> &Q,
                                        TPZVec<int64_t> &adjNodes){

    R.Push( Parent );
    ExploredNodes[Parent] = 1;
    graph.AdjacentNodesOrdered(Parent, ExploredNodes, adjNodes);
    const int64_t nadj = adjNodes.NElements();
    for(int64_t i = 0; i < nadj; i++){
      if(ExploredNodes[adjNodes[i]] == 0){//then a new element
        Q.push(adjNodes[i]);
        ExploredNodes[adjNodes[i]] = 1;
      }
    }//for


}//void

