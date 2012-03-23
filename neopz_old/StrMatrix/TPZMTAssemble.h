//$Id: TPZMTAssemble.h,v 1.1 2008-03-18 12:15:35 cesar Exp $

#ifndef TPZMTASSEMBLE_H
#define TPZMTASSEMBLE_H

#include <set>
#include <vector>
#include "pzvec.h"
#include "tpzautopointer.h"
class TPZCompMesh;
class TPZMatrix;
class TPZFMatrix;
class TPZCompEl;
class TPZElementMatrix;

/** Auxiliar structure
 */
struct SMTAssembleResidual{
  TPZCompEl * compel;
  TPZFMatrix * rhs;
  std::set<int> *MaterialIds;
  int mineq, maxeq, index;
  SMTAssembleResidual( TPZCompEl * _compel, TPZFMatrix * _rhs, int _mineq, int _maxeq, std::set<int> * _MaterialIds, int _index)
  {
    compel = _compel;
    rhs = _rhs;
    mineq = _mineq;
    maxeq = _maxeq;
    MaterialIds = _MaterialIds;
    index = _index;
  }
  ~SMTAssembleResidual(){
    ///nothing to be done here
  }

};

/**
 * Class of static methods to replace TPZStructMatrix::Assemble in a multi-threading execution
 */
class TPZMTAssemble{
public:

  /** Default constructor
   */
  TPZMTAssemble();

  /** Destructor
   */
  ~TPZMTAssemble();
  
  /** Multi-threading assemblage process
   */
  static void AssembleMT(TPZFMatrix & rhs, TPZCompMesh &mesh, int mineq, int maxeq, std::set<int> *MaterialIds);

protected:

  /** Stack of computed element residuals waiting to be assembled.
   */
  static std::vector< std::pair< TPZElementMatrix *, SMTAssembleResidual * > > gComputedEF;

  /** Assemble computed element residuals.
   */
  static void ContributeEFs();
  
  /** The thread execution
   */
  static void * ExecuteAssembleResidualMT(void * data);

};

#endif
