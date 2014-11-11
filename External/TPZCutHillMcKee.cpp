
#include "TPZCutHillMcKee.h"
#include "fstream"

void TPZCutHillMcKee::SGraph::Set2Vec(const std::set<long> &myset,
                                      TPZVec<long> &myVec) const{
  const long n = myset.size();
  myVec.Resize(n);
  std::set<long>::const_iterator w, e = myset.end();
  long i;
  for(w = myset.begin(), i = 0; w != e; w++, i++){
    myVec[i] = *w;
  }
}//void

long TPZCutHillMcKee::SGraph::SmallestDegree(TPZVec<long> &ExploredNodes){
    std::cout << "TPZCutHillMcKee::SGraph::SmallestDegree size = " << this->NNodes() << std::endl;
    std::cout.flush();
  long mindegree = 1000000;
  long found = -1;
  const long n = this->NNodes();
  for(long i = 0; i < n; i++){
    if(ExploredNodes[i] == 1) continue;
    long degree = this->Degree(i);
    if(degree < mindegree){
      mindegree = degree;
      found = i;
    }
  }//for
  return found;
}


void TPZCutHillMcKee::SGraph::GetAdjacentNodes(const TPZVec<long> &parents,
                                               const TPZVec< long > &exceptedNodes,
                                               TPZStack<long> &adjNodes){
  adjNodes.Resize(this->NNodes());//pre allocating memory
  long count = 0;
  TPZVec<long> History(this->NNodes(),0); //avoiding duplicates in adjNodes. it could be a std::set. set is better when few elements
  for(long i = 0; i < parents.NElements(); i++){
    long nlocal;
    long * localAdjNodes = AdjacentNodesPtr(parents[i], nlocal);
    for(long j = 0; j < nlocal; j++){
      const long node = localAdjNodes[j];
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

void TPZCutHillMcKee::SGraph::AdjacentNodes(long parent, TPZVec<long> &adjNodes){
  const long n = this->Degree(parent);
  adjNodes.Resize(n);
  for(long i = 0; i < n; i++){
    const long adj = fnodegraph[fnodegraphindex[parent]+i];
    adjNodes[i] = adj;
  }//for
}//void

void TPZCutHillMcKee::SGraph::AdjacentNodesOrdered(long parent,
                                                   const TPZVec<long> &exceptedNodes,
                                                   TPZVec<long> &adjNodes){
  std::multimap<long,long> order;

  const long ndegree = this->Degree(parent);
  long count = 0;
  for(long i = 0; i < ndegree; i++){
    const long adj = fnodegraph[fnodegraphindex[parent]+i];
    if(exceptedNodes[adj] == 1) continue;
    const long adjDegree = this->Degree(adj);
    order.insert(std::make_pair(adjDegree,adj));
    count++;
  }//for

  adjNodes.Resize(count);
#ifdef DEBUG
  if(count != (long)(order.size()) ){
    DebugStop();
  }
#endif
  std::multimap<long,long>::const_iterator mapIt;
  long i;
  for(i = 0, mapIt = order.begin(); mapIt != order.end(); i++, mapIt++){
    adjNodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::RootedLevelStructure(long rootNode,
                                                   TPZStack< TPZStack<long> > &LevelStructure){
  LevelStructure.Resize(0);

  TPZManVector<long,1000> adjNodes;

  TPZStack<long> thisLevel;
  thisLevel.Resize(1);
  thisLevel[0] = rootNode;

  TPZVec<long> ProcessedNodes(this->NNodes(),0);

  for(long iLevel = 0; iLevel < this->NNodes(); iLevel++){

    const long nThisLevel = thisLevel.NElements();
    if(nThisLevel == 0) break;
    for(long i = 0; i < nThisLevel; i++){
      const long node = thisLevel[i];
      ProcessedNodes[node] = 1;
    }
    LevelStructure.Push( thisLevel);

    this->GetAdjacentNodes(LevelStructure[LevelStructure.NElements() - 1], ProcessedNodes, thisLevel);

  }//for que de fato eh while

  LevelStructure.Shrink();//saving memory

}//void

void TPZCutHillMcKee::SGraph::SortNodes(TPZVec<long> &nodes){
  std::multimap<long,long> mymap;
  for(long i = 0; i < nodes.NElements(); i++){
    const long node = nodes[i];
    const long degree = this->Degree(node);
    mymap.insert(std::make_pair(degree,node));
  }//for
  std::multimap<long,long>::const_iterator mapIt;
  long i;
  for(i = 0, mapIt = mymap.begin(); mapIt != mymap.end(); i++, mapIt++){
    nodes[i] = mapIt->second;
  }
}//void

void TPZCutHillMcKee::SGraph::ShrinkLastLevel(TPZVec<long> &LastLevel){

  const long nelsOrig = LastLevel.NElements();

  //accorging to suggestion of INTERNATIONAL JOURNAL FOR NUMERICAL METHODS IN ENGINEERING, VOL. 28,2651-2679 (1989)
  //A FORTRAN PROGRAM FOR PROFILE AND WAVEFRONT REDUCTION by S. W. SLOAN
  //removing nodes with same degree
  std::map< long, long > mymap;
  for(long i = 0; i < LastLevel.NElements(); i++){
    const long node = LastLevel[i];
    const long degree = this->Degree( node );
    mymap[ degree ] = node;
  }
  LastLevel.Resize(mymap.size());
  std::map< long, long >::const_iterator w, e = mymap.end();
  long i;
  for(i = 0, w = mymap.begin(); w != e; w++, i++){
    LastLevel[i] = w->second;
  }

  //Previously, Sloan had suggested, in the article of 1986, LastLevel.Resize( (LastLevel.NElements()+2)/2 );
  const long newsize = (nelsOrig+2)/2;
  if(newsize < LastLevel.NElements()){
    LastLevel.Resize( newsize );
  }

}

void TPZCutHillMcKee::SGraph::PseudoPeripheralNodes(long &startNode, long &endNode){

  endNode = -1;

  TPZVec<long> emptyVec(this->NNodes(), 0);
  //step 1: first guess
  startNode = SmallestDegree(emptyVec);

  //step 2: rooted level structure
  TPZStack< TPZStack<long> > LevelStructure;
  this->RootedLevelStructure(startNode, LevelStructure);

#ifdef DEBUG_CM
  std::ofstream myfile("c:\\Temp\\PseudoPeripheralNodes.txt");
#endif

  long count = 0;
  const long maxCount = 10;
  TPZStack< TPZStack<long> > localLevelStructure;
  TPZStack<long> LastLevel;
  LastLevel.Resize( this->NNodes() );//pre allocating memory
  LastLevel.Resize(0);
  while(count < maxCount || endNode == -1){ //maxCount pra não ficar aqui eternamente

    count++;
    const long cpStart = startNode;
    const long cpEnd = endNode;
#ifdef DEBUG_CM
    myfile << "\nstart = " << startNode << "  end = " << endNode << "\n\n";
#endif
    const long nlevels = LevelStructure.NElements();
    LastLevel = LevelStructure[ nlevels-1 ];
    //step 3 - sort the last level
    this->SortNodes(LastLevel);
    //step 4 - shrink the last level
    ShrinkLastLevel(LastLevel);

    //step 5
    long we = 1e9;
    long hs = nlevels;
    const long nelLastLevel = LastLevel.NElements();
    for(long iQ = 0; iQ < nelLastLevel; iQ++){

      const long nodeAtQ = LastLevel[iQ];
#ifdef DEBUG_CM
      myfile << count << "\t" << iQ << "\t";myfile.flush();
#endif
      this->RootedLevelStructure(nodeAtQ, localLevelStructure);
#ifdef DEBUG_CM
      myfile << "rooted\t";myfile.flush();
#endif

      const long h = localLevelStructure.NElements();
      long w = 0;
      for(long j = 0; j < h; j++){
        const long localW = localLevelStructure[j].NElements();
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

  if(startNode == -1 || endNode == -1) DebugStop();

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

TPZCutHillMcKee::TPZCutHillMcKee(long NElements, long NNodes, bool Reverse):TPZRenumbering(NElements,NNodes){
  fReverse = Reverse;
#ifdef DEBUG
  fVerbose = true;
#else
  fVerbose = false;
#endif
}

TPZCutHillMcKee::~TPZCutHillMcKee(){

}


void TPZCutHillMcKee::Resequence(TPZVec<long> &perm, TPZVec<long> &iperm){

  TPZVec<long> permGather, permScatter;
  TPZVec<long> permGatherReverse, permScatterReverse;

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

void TPZCutHillMcKee::Resequence(TPZVec<long> &permGather, TPZVec<long> &permScatter,
                                 TPZVec<long> &permGatherReverse, TPZVec<long> &permScatterReverse){

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
    long grafoindex[7] = {0,3,6,9,14,19,22};
    long grafo[22] = {3,4,5,  2,3,4,  1,3,4,  0,1,2,4,5,  0,1,2,3,5,  0,3,4};
    graph.fnodegraphindex.Resize(7);
    for(long i = 0; i < 7; i++) graph.fnodegraphindex[i] = grafoindex[i];
    graph.fnodegraph.Resize(22);
    for(long i = 0; i < 22; i++) graph.fnodegraph[i] = grafo[i];

    {
      TPZStack< TPZStack<long> > LevelStructure;
      graph.RootedLevelStructure(1, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructure.txt");
      for(long iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(long j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

    {
      long startNode = -1, endNode = -1;
      graph.PseudoPeripheralNodes(startNode, endNode);
      TPZStack< TPZStack<long> > LevelStructure;
      graph.RootedLevelStructure(startNode, LevelStructure);
      std::ofstream myfile("c:\\Temp\\LevelStructureAfterPeripheral.txt");
      for(long iLevel = 0; iLevel < LevelStructure.NElements(); iLevel++){
        myfile << "Level " << iLevel << ": ";
        for(long j = 0; j < LevelStructure[iLevel].NElements(); j++){
          myfile << LevelStructure[iLevel][j] << "  ";
        }
        myfile << "\n";
      }
    }

  }//if
#endif

        //CutHillMcKee
  std::queue<long> Q;
  TPZStack<long> R;
  TPZManVector<long> adjNodes;
  //reservando memoria em R
  const long nnodes = graph.NNodes();
  R.Resize(nnodes);
  R.Resize(0);
  TPZVec<long> ExploredNodes(nnodes,0);
  if(this->fVerbose) {
    std::cout << "TPZCutHillMcKee process...\n";std::cout.flush();
  }

  bool firstTime = true;
  while(R.NElements() != graph.NNodes()){
    if(this->fVerbose){
      std::cout << "Calling TPZCutHillMcKee::SGraph::SmallestDegree\n";
      std::cout.flush();
    }
    long Parent = -1;
    if(firstTime){
      long startNode = -1, endNode = -1;
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
      const long Child = Q.front();
      Q.pop();
      this->ProcessParentNode(Child, graph, ExploredNodes, R, Q, adjNodes);
    }

  }//loop completo

  const long n = R.NElements();
  if(n != nnodes) DebugStop();

#ifdef DEBUG
{//verificando se ha duplicados
  std::set<long> check;
  for(long i = 0; i < n; i++) check.insert(R[i]);
  if( ((long)(check.size())) != n) DebugStop();
}
#endif

  if(this->fVerbose){
    std::cout << "TPZCutHillMcKee Filling perm and iperm vectors...";std::cout.flush();
  }


  permScatterReverse.Resize(n);
  for(long i = 0; i < n; i++) permScatterReverse[n-1-i] = R[i];
  permGatherReverse.Resize(n);
  for(long i = 0; i < n; i++) permGatherReverse[ permScatterReverse[i] ] = i;


  permScatter = R;
  permGather.Resize(n);
  for(long i = 0; i < n; i++) permGather[ permScatter[i] ] = i;


}///void

void TPZCutHillMcKee::ProcessParentNode(long Parent,
                                        SGraph &graph,
                                        TPZVec<long> &ExploredNodes,
                                        TPZStack<long> &R,
                                        std::queue<long> &Q,
                                        TPZVec<long> &adjNodes){

    R.Push( Parent );
    ExploredNodes[Parent] = 1;
    graph.AdjacentNodesOrdered(Parent, ExploredNodes, adjNodes);
    const long nadj = adjNodes.NElements();
    for(long i = 0; i < nadj; i++){
      if(ExploredNodes[adjNodes[i]] == 0){//then a new element
        Q.push(adjNodes[i]);
        ExploredNodes[adjNodes[i]] = 1;
      }
    }//for


}//void

