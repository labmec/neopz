
#include <stdio.h>
#include <cstdlib>
#include <fstream>
#include <cmath>
//#include "pztempmat.h"

template <class TConv>
void CheckConvergence(TConv &obj, TPZFMatrix<REAL> &state, TPZFMatrix<REAL> &range, TPZVec<REAL> &coefs){

   TPZFMatrix<REAL> incval(range);
   int i,j,numrows;
   numrows = state.Rows();
   for(i=0; i<numrows; i++) {
     double randnum = (rand()&1000)/999.;
      incval(i,0) = range(i,0)*randnum;
   }
   std::ofstream logfile("conv.log");
   int icase;
   int numcases = obj.NumCases();
   int ncoefs = coefs.NElements();
   for(icase = 0; icase < numcases; icase++) {
      obj.LoadSolution(state);
      TPZFMatrix<REAL> ReferenceResidual, Tangent, EstimateRes;
      obj.ComputeTangent(Tangent,coefs,icase);
      obj.Residual(ReferenceResidual,icase);
      EstimateRes.Redim(ReferenceResidual.Rows(),ReferenceResidual.Cols());
      Tangent.Multiply(incval,EstimateRes);
      int interval;
      double difnorm[10] = {0.};
      for(interval = 1; interval < 10; interval++) {
         TPZFMatrix<REAL> actualstate(state);
         TPZFMatrix<REAL> residual;
         for(i=0; i<numrows; i++) {
            for(j=0; j<ncoefs; j++) {
               actualstate(i,j) += (interval/10.)*incval(i,0)*coefs[j];
            }
         }
         obj.LoadSolution(actualstate);
         obj.Residual(residual,icase);
         residual -= ReferenceResidual;
         residual = residual - EstimateRes*(interval/10.);
         difnorm[interval] = Norm(residual);
      }
      std::cout << "icase = " << icase << std::endl;
      logfile << "icase = " << icase << std::endl;
      for(interval = 2; interval<10; interval++) {
         if(fabs(difnorm[interval]) < 1.e-12 || fabs(difnorm[interval-1]) <1.e-12) {
	   std::cout << "residual too small\n";
           logfile << "residual too small\n";
	   break;
	 }
         REAL val = ( std::log10(difnorm[interval]) - std::log10(difnorm[interval-1]));
         REAL intval = ( std::log10((REAL)interval) - std::log10((REAL)(interval-1)));
         val = val / intval;
         std::cout << val << std::endl;
         logfile << val << std::endl;
      }
   }
   logfile.flush();
}
