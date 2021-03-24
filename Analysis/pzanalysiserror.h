/**
 * @file
 * @brief Contains TPZAnalysisError class which implements analysis procedures with hp adaptivity.
 */

#ifndef TPZANALYSISERRORH
#define TPZANALYSISERRORH

#include <fstream>        // for ostream, operator<<
#include "pzanalysis.h"   // for TPZAnalysis
#include "pzcompel.h"     // for TPZCompElSide, TPZCompEl (ptr only)
#include "pzmanvector.h"  // for TPZManVector
#include "pzreal.h"       // for REAL
#include "pzstack.h"      // for TPZStack
class TPZCompMesh;

/**
 * @brief Implements analysis procedures with hp adaptivity. \ref analysis "Analysis"
 * @ingroup analysis
 * @note This class implements : "A Fast hp Adaptive Finite Element Mesh Design" for : A. A. Novotny et al.
 */
class TPZAnalysisError : public TPZAnalysis {
	
	/** @brief Indexes of the elements vector */
	TPZManVector<int64_t> fElIndexes;
	/** @brief Vector with error values by elements */
	TPZManVector<REAL> fElErrors;
	TPZStack<TPZCompElSide> fSingular;
	/** Total error computed */
	REAL fTotalError;
	REAL fAdmissibleError;
	
	REAL fEtaAdmissible;
	/** @brief Number of iterations */
	int64_t fNIterations;
	
	public :
	/** @brief Object constructors*/
	TPZAnalysisError(TPZCompMesh *mesh,std::ostream &out);
	/** @brief Delete objects*/
	~TPZAnalysisError() {};
	
	/** @brief Set the parameters which will govern the adaptive process*/
	void SetAdaptivityParameters(REAL EtaAdmissible, int64_t NIterations);
	
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
	
	/** @brief Compute the list of errors of all elements and also the admissible error for any element in the grid */
	/** Is called from HPAdapt() */
	void EvaluateError(REAL CurrentEtaAdmissible, bool store_error, std::ostream &out);
	
private:
	/** @brief Postprocess the intermediate solutions*/
	void PlotLocal(int64_t iter, REAL CurrentEtaAdmissible, std::ostream &out);
	
	void ExpandConnected(TPZStack<TPZCompElSide> &singel);	
};

#endif
