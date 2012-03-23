// TPZSubMeshAnalysis.h: interface for the TPZSubMeshAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZSUBMESHANALYSIS_H
#define TPZSUBMESHANALYSIS_H


#include "pzanalysis.h"
#include "pzmatred.h"
class TPZSubCompMesh;

#include "pzfmatrix.h"	// Added by ClassView
class TPZSubMeshAnalysis : public TPZAnalysis  
{
private:
	TPZFMatrix fReferenceSolution;
	
	TPZAutoPointer<TPZMatrix > fReducableStiff;
	TPZSubCompMesh *fMesh;

public:
	virtual void LoadSolution(TPZFMatrix &sol);
  /**
    *Constructor: create an object analysis from one mesh
	**/
	TPZSubMeshAnalysis(TPZSubCompMesh *mesh);

  /**
    *Destructor
	**/
	virtual ~TPZSubMeshAnalysis();


  /**
    *Run: assemble the stiffness matrix
	**/
	void Run(std::ostream &out);

  /**
    *CondensedSolution: returns the condensed stiffness
    *matrix - ek - and the condensed solution vector - ef
    */
  void CondensedSolution(TPZFMatrix &ek, TPZFMatrix &ef);

  /**
   * Assemble the global stiffness matrix and put it into the 
   * reducable stiffness matrix
   */
  virtual void Assemble();
};

#endif 
