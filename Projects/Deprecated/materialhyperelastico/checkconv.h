

#include <stdio.h>

#include <stdlib.h>

#include <fstream.h>

//#include "pztempmat.h"



template <class TConv>

void CheckConvergence(TConv &obj, TPZFMatrix &state, TPZFMatrix &range, TPZVec<REAL> &coefs){



   TPZFMatrix incval(range);

   int i,j,numrows;

   numrows = state.Rows();

   for(i=0; i<numrows; i++) {

     double randnum = (rand()&1000)/999.;

      incval(i,0) = range(i,0)*randnum;

   }

   ofstream log("conv.log");

   int icase;

   int numcases = obj.NumCases();

   int ncoefs = coefs.NElements();

   for(icase = 0; icase < numcases; icase++) {

      obj.LoadSolution(state);

      TPZFMatrix ReferenceResidual, Tangent, EstimateRes;

      obj.ComputeTangent(Tangent,coefs,icase);

      obj.Residual(ReferenceResidual,icase);

      EstimateRes.Redim(ReferenceResidual.Rows(),ReferenceResidual.Cols());

      Tangent.Multiply(incval,EstimateRes);

      int interval;

      double difnorm[10] = {0.};

//	  double resnorm[10] ={0.};

      for(interval = 1; interval < 10; interval++) {

         TPZFMatrix actualstate(state);

         TPZFMatrix residual;

         for(i=0; i<numrows; i++) {

            for(j=0; j<ncoefs; j++) {

               actualstate(i,j) += (interval/10.)*incval(i,0)*coefs[j];

            }

         }

         obj.LoadSolution(actualstate);

         obj.Residual(residual,icase);

         residual -= ReferenceResidual;

//		 resnorm[interval] = Norm(residual);

         residual = residual - EstimateRes*(interval/10.);

         difnorm[interval] = Norm(residual);

      }

      cout << "icase = " << icase << endl;

      log << "icase = " << icase << endl;

      for(interval = 2; interval<10; interval++) {

         if(fabs(difnorm[interval]) < 1.e-12 || fabs(difnorm[interval-1]) <1.e-12) {

	   cout << "residual too small\n";

           log << "residual too small\n";

	   break;

	 }

         cout << (log10(difnorm[interval])-log10(difnorm[interval-1]))/

                 (log10(interval)-log10(interval-1)) << endl;

         log << (log10(difnorm[interval])-log10(difnorm[interval-1]))/

                (log10(interval)-log10(interval-1)) << endl;

      }

   }

   log.flush();

}

