/**
 * @file
 * @brief Contains the TPZMTAssemble class which replaces TPZStructMatrix::Assemble in a multi-threading execution.\n
 * Also contains SMTAssembleResidual structure.
 */

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
struct TPZElementMatrix;

/** 
 * @ingroup structural
 * @brief Auxiliar structure to assembling residual vector. \ref structural "Structural Matrix"
 */
struct SMTAssembleResidual{
	/// Pointer to computational element
	TPZCompEl * compel;
	/// Residual matrix
	TPZFMatrix * rhs;
	/// Material identifiers
	std::set<int> *MaterialIds;
	/// Index of the first equation
	int mineq;
	/// Index of the last equation
	int maxeq;
	/// Index to assemble
	int index;

	/// Constructor
	SMTAssembleResidual( TPZCompEl * _compel, TPZFMatrix * _rhs, int _mineq, int _maxeq, std::set<int> * _MaterialIds, int _index)
	{
		compel = _compel;
		rhs = _rhs;
		mineq = _mineq;
		maxeq = _maxeq;
		MaterialIds = _MaterialIds;
		index = _index;
	}
	/// Destructor
	~SMTAssembleResidual(){
		//nothing to be done here
	}
};

/**
 * @ingroup structural
 * @brief Class of static methods to replace TPZStructMatrix::Assemble in a multi-threading execution. \ref structural "Structural Matrix"
 */
class TPZMTAssemble{
public:
	
	/** @brief Default constructor */
	TPZMTAssemble();
	
	/** @brief Destructor */
	~TPZMTAssemble();
	
	/** @brief Multi-threading assemblage process */
	static void AssembleMT(TPZFMatrix & rhs, TPZCompMesh &mesh, int mineq, int maxeq, std::set<int> *MaterialIds);
	
protected:
	
	/** @brief Stack of computed element residuals waiting to be assembled */
	static std::vector< std::pair< TPZElementMatrix *, SMTAssembleResidual * > > gComputedEF;
	
	/** @brief Assembles computed element residuals */
	static void ContributeEFs();
	
	/** @brief The thread execution */
	static void * ExecuteAssembleResidualMT(void * data);
};

#endif
