/**This class implementent : "A Fast hp Adaptive Finite Element Mesh Design"
                       for : A. A. Novotny et al.*/
// -*- c++ -*-
                       
#ifndef TPZANALYSISERRORH
#define TPZANALYSISERRORH

#include <fstream>
#include <iostream>
using namespace std;

#include "pzanalysis.h"
#include "pzstack.h"
#include "pzcompel.h"
class TPZCompMesh;
class TPZCompElSide;

class TPZAnalysisError : public TPZAnalysis {

   TPZManVector<int> fElIndexes;
   TPZManVector<REAL> fElErrors;
   TPZStack<TPZCompElSide> fSingular;
   REAL fTotalError;
   REAL fAdmissibleError;

   REAL fEtaAdmissible;
   int fNIterations;

public :
   /**Object constructors*/
	//TPZAnalysisError(TPZAnalysis &an);
   TPZAnalysisError(TPZCompMesh *mesh,ostream &out);
   /**Delete objects*/
   ~TPZAnalysisError() {};

   /**Set the parameters which will govern the adaptive process*/
   void SetAdaptivityParameters(REAL EtaAdmissible, int NIterations);

   /**Run the algorithm of the fast hp adaptive finite element mesh design*/
   void hp_Adaptive_Mesh_Design(ostream &out,REAL &EtaAdmissible);

   /**Search the element whith contain this point*/
   void GetSingularElements(TPZStack<TPZCompElSide> &listel);

   /**Calculate pn and hn parameters for the elements neighbours to the element
      with contain the singular point
	the parameter csi determines the number of layers of refinement
	singularity_order determines the strength of the singularity */
   void ZoomInSingularity(REAL csi, TPZCompElSide elside, REAL singularity_order = 0.9);

   /**Run one iteration of HP adaptivity*/
   void HPAdapt(REAL CurrentEtaAdmissible, ostream &out);

   /**Return the maximal local error of the elements of the mesh*/
   REAL MaximLocalError();

   /**Calculate the h parameter of the element*/
   REAL h_Parameter(TPZCompEl *cel);

   /**Plot to the aproximated solution of the FEM with Mathematica package*/
	void MathematicaPlot();

   /** Compute the list of errors of all elements and also the admissible error
   for any element in the grid
   Is called from HPAdapt()
   */
   void EvaluateError(REAL CurrentEtaAdmissible, ostream &out);

   /**Postprocess the intermediate solutions*/
private:
	void PlotLocal(int iter, REAL CurrentEtaAdmissible, ostream &out);

	void ExpandConnected(TPZStack<TPZCompElSide> &singel);

};

#endif
