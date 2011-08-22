/**
 * @file
 * @brief Contains TPZAnalysisError class which implements analysis procedures with hp adaptivity.
 */

#ifndef TPZANALYSISERRORH
#define TPZANALYSISERRORH

#include <fstream>
#include <iostream>

#include "pzanalysis.h"
#include "pzstack.h"
#include "pzcompel.h"
class TPZCompMesh;
class TPZCompElSide;

/**
 * @brief Implements analysis procedures with hp adaptivity. \ref analysis "Analysis"
 * @ingroup analysis
 * @note This class implements : "A Fast hp Adaptive Finite Element Mesh Design"
 * for : A. A. Novotny et al.
 */
class TPZAnalysisError : public TPZAnalysis {
	
	TPZManVector<int> fElIndexes;
	TPZManVector<REAL> fElErrors;
	TPZStack<TPZCompElSide> fSingular;
	REAL fTotalError;
	REAL fAdmissibleError;
	
	REAL fEtaAdmissible;
	int fNIterations;
	
	public :
	/** @brief Object constructors*/
	//TPZAnalysisError(TPZAnalysis &an);
	TPZAnalysisError(TPZCompMesh *mesh,std::ostream &out);
	/** @brief Delete objects*/
	~TPZAnalysisError() {};
	
	/** @brief Set the parameters which will govern the adaptive process*/
	void SetAdaptivityParameters(REAL EtaAdmissible, int NIterations);
	
	/** @brief Run the algorithm of the fast hp adaptive finite element mesh design*/
	void hp_Adaptive_Mesh_Design(std::ostream &out,REAL &EtaAdmissible);
	
	/** @brief Search the element whith contain this point*/
	void GetSingularElements(TPZStack<TPZCompElSide> &listel);
	
	/**
	 * @brief Calculate pn and hn parameters for the elements neighbours to the element with contain the singular point
	 * @param elside Computational element and side to zoom
	 * @param csi determines the number of layers of refinement
	 * @param singularity_order determines the strength of the singularity 
	 */
	void ZoomInSingularity(REAL csi, TPZCompElSide elside, REAL singularity_order = 0.9);
	
	/** @brief Run one iteration of HP adaptivity*/
	void HPAdapt(REAL CurrentEtaAdmissible, std::ostream &out);
	
	/** @brief Returns the maximal local error of the elements of the mesh*/
	REAL MaximLocalError();
	
	/** @brief Calculate the h parameter of the element*/
	REAL h_Parameter(TPZCompEl *cel);
	
	/** @brief Plot to the aproximated solution of the FEM with Mathematica package*/
	void MathematicaPlot();
	
	/** @brief Compute the list of errors of all elements and also the admissible error
	 for any element in the grid
	 
	 Is called from HPAdapt()
	 */
	void EvaluateError(REAL CurrentEtaAdmissible, std::ostream &out);
	
	/** @brief Postprocess the intermediate solutions*/
private:
	void PlotLocal(int iter, REAL CurrentEtaAdmissible, std::ostream &out);
	
	void ExpandConnected(TPZStack<TPZCompElSide> &singel);
	
};

#endif
