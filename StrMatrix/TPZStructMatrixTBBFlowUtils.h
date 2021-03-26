#ifndef TPZSRTMATRIXTBBFLOWUTILS_H
#define TPZSRTMATRIXTBBFLOWUTILS_H
/*This header file contains auxiliary classes for TPZStructMatrixTBBFlow.
It is only included if USING_TBB=ON, and it is not installed.*/

#ifdef USING_TBB
#include "tbb/tbb.h"
#include "tbb/flow_graph.h"

#include "pzmanvector.h"
#include "tpzautopointer.h"
template<class T>
class TPZFMatrix;
template <class T>
class TPZMatrix;
class TPZCompMesh;
class TPZGuiInterface;
class TPZStructMatrixTBBFlow;


class TPZFlowGraph {
public:
  TPZFlowGraph(TPZStructMatrixTBBFlow *strmat);
  ~TPZFlowGraph();
  TPZFlowGraph(TPZFlowGraph const &copy);

  // vectors for mesh coloring
  TPZVec<int> fnextBlocked, felSequenceColor, felSequenceColorInv;
  TPZManVector<int> fElementOrder;

  TPZCompMesh *cmesh;

  tbb::flow::graph fGraph;
  tbb::flow::broadcast_node<tbb::flow::continue_msg> fStartNode;
  std::vector<tbb::flow::continue_node<tbb::flow::continue_msg> *> fNodes;

  void OrderElements();
  void ElementColoring();
  void CreateGraph();
  void ExecuteGraph(TPZFMatrix<STATE> *rhs, TPZMatrix<STATE> *matrix = 0);

  /// current structmatrix object
  TPZStructMatrixTBBFlow *fStruct;
  /// gui interface object
  TPZAutoPointer<TPZGuiInterface> fGuiInterface;
  /// global matrix
  TPZMatrix<STATE> *fGlobMatrix;
  /// global rhs vector
  TPZFMatrix<STATE> *fGlobRhs;
};

class TPZFlowNode {
public:
  TPZFlowNode(TPZFlowGraph *graph, int el) : myGraph(graph), iel(el){};

  ~TPZFlowNode(){};

  void operator()(tbb::flow::continue_msg) const;

  TPZFlowGraph *myGraph;

  /// element to be processed by this node
  int iel;
};

struct TPZGraphThreadData {
  // copy constructor
  TPZGraphThreadData(TPZCompMesh *cmesh, TPZVec<int> &fnextBlocked,
                     TPZVec<int> &felSequenceColor,
                     TPZVec<int> &felSequenceColorInv);
  // destructor
  ~TPZGraphThreadData();
  // tbb tasks graph
  tbb::flow::graph fAssembleGraph;
  // initial node
  tbb::flow::broadcast_node<tbb::flow::continue_msg> fStart;
  // store all the nodes
  std::vector<tbb::flow::continue_node<tbb::flow::continue_msg> *> fGraphNodes;
  // vector for coloring mesh
  TPZVec<int> felSequenceColor;

  /// current structmatrix object
  TPZStructMatrixTBBFlow *fStruct;
  /// gui interface object
  TPZAutoPointer<TPZGuiInterface> fGuiInterface;
  /// global matrix
  TPZMatrix<STATE> *fGlobMatrix;
  /// global rhs vector
  TPZFMatrix<STATE> *fGlobRhs;
};

struct TPZGraphThreadNode {
  TPZGraphThreadData *data;
  int iel;
  TPZGraphThreadNode(TPZGraphThreadData *data, int el) : data(data), iel(el) {}
  void operator()(tbb::flow::continue_msg) const;
};
#else
class TPZFlowGraph{};
#endif
#endif