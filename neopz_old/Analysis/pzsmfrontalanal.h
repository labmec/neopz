// TPZSubMeshFrontalAnalysis.h: interface for the TPZSubMeshFrontalAnalysis class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZSUBMESHFRONTANALYSIS_H
#define TPZSUBMESHFRONTANALYSIS_H


#include "pzanalysis.h"
#include "pzmatred.h"
#include "TPZFrontMatrix.h"
class TPZSubCompMesh;
class TPZFront;

#include "pzfmatrix.h"	// Added by ClassView
/**
 * Analysis for substructuring
 */
class TPZSubMeshFrontalAnalysis : public TPZAnalysis  
{
private:
  TPZFMatrix fReferenceSolution;
	
  //	TPZMatRed fReducableStiff;
//  TPZFrontMatrix<TPZFileEqnStorage, TPZFrontSym> fFrontalStiff;
  
  TPZSubCompMesh *fMesh;

  TPZFront *fFront;

public:
  virtual void LoadSolution(TPZFMatrix &sol);
  /**
   *Constructor: create an object analysis from one mesh
   **/
  TPZSubMeshFrontalAnalysis(TPZSubCompMesh *mesh);

  /**
   *Destructor
   **/
  virtual ~TPZSubMeshFrontalAnalysis();


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
   *Set the front matrix
   */
  void SetFront(TPZFront &front) { fFront = &front;}

};

#endif 
