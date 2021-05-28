#ifndef TPZSRTMATRIXTBBFLOWUTILS_H
#define TPZSRTMATRIXTBBFLOWUTILS_H
/*This header file contains auxiliary classes for TPZStructMatrixTBBFlow.
It is only included if USING_TBB=ON, and it is not installed.*/

#ifdef USING_TBB
#include "tbb/tbb.h"
#include "tbb/flow_graph.h"

#include "pzmanvector.h"
#include "tpzautopointer.h"
class TPZBaseMatrix;
template<class T>
class TPZFMatrix;
template <class T>
class TPZMatrix;
class TPZCompMesh;
class TPZGuiInterface;
template<class T>
class TPZStructMatrixTBBFlow;

template<class TVar = STATE>
class TPZFlowGraph {
public:
  TPZFlowGraph(TPZStructMatrixTBBFlow<TVar> *strmat, bool computeRhs);
  ~TPZFlowGraph();
  TPZFlowGraph(TPZFlowGraph const &copy);

  // vectors for mesh coloring
  TPZVec<int> fnextBlocked, felSequenceColor, felSequenceColorInv;
  TPZManVector<int> fElementOrder;

  TPZCompMesh *fMesh;

  tbb::flow::graph fGraph;
  tbb::flow::broadcast_node<tbb::flow::continue_msg> fStartNode;
  std::vector<tbb::flow::continue_node<tbb::flow::continue_msg> *> fNodes;

  void OrderElements();
  void ElementColoring();
  void CreateGraph();
  void ExecuteGraph(TPZBaseMatrix *rhs, TPZBaseMatrix *matrix = 0);

  /// current structmatrix object
  TPZStructMatrixTBBFlow<TVar> *fStruct;
  /// gui interface object
  TPZAutoPointer<TPZGuiInterface> fGuiInterface;
  /// global matrix
  TPZMatrix<TVar> *fGlobMatrix;
  /// global rhs vector
  TPZFMatrix<TVar> *fGlobRhs;
  /// whether to compute rhs
  bool fComputeRhs{true};
};

template<class TVar = STATE>
class TPZFlowNode {
public:
  TPZFlowNode(TPZFlowGraph<TVar> *graph, int el, bool computeRhs)
    : myGraph(graph), iel(el), fComputeRhs(computeRhs){};

  ~TPZFlowNode(){};

  void operator()(tbb::flow::continue_msg) const;

  TPZFlowGraph<TVar> *myGraph;

  /// element to be processed by this node
  int iel;
  /// whether to compute rhs
  bool fComputeRhs{true};
};

template<class TVar>
struct TPZGraphThreadData {
  // copy constructor
  TPZGraphThreadData(TPZCompMesh *cmesh, TPZVec<int> &fnextBlocked,
                     TPZVec<int> &felSequenceColor,
                     TPZVec<int> &felSequenceColorInv,
                     bool computeRhs);
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
  TPZStructMatrixTBBFlow<TVar> *fStruct;
  /// gui interface object
  TPZAutoPointer<TPZGuiInterface> fGuiInterface;
  /// global matrix
  TPZMatrix<STATE> *fGlobMatrix;
  /// global rhs vector
  TPZFMatrix<STATE> *fGlobRhs;
  /// whether to compute rhs
  bool fComputeRhs{true};
};

template<class TVar>
struct TPZGraphThreadNode {
  TPZGraphThreadData<TVar> *data;
  int iel;
  TPZGraphThreadNode(TPZGraphThreadData<TVar> *data, int el) : data(data), iel(el) {}
  void operator()(tbb::flow::continue_msg) const;
};
#else
template<class TVar>
class TPZFlowGraph{};
#endif
#endif
