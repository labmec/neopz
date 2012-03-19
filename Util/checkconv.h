/**
 * @file
 * @brief Contains the implementation of the CheckConvergence function.
 */

#include <stdio.h>

#include <stdlib.h>

#include <fstream>


template <class TConv>

/**
 * @brief Implements a general procedure to check whether the class TConv implements a consistente tangent matrix for the 
 * object obj.\n
 */
/**
 * The following methods need to be implemented for TConv \n
 * -# LoadState : to load the state variables in local storage \n
 * -# ComputeTangent : which will compute a tangent matrix \n
 * -# Residual : which will compute a residual vector
 */
void CheckConvergence(TConv &obj, TPZFMatrix<REAL> &state, TPZFMatrix<REAL> &range, TPZVec<REAL> &coefs) {

	TPZFMatrix<REAL> incval(range);
	int i,numrows;
	
	numrows = state.Rows();
	for(i=0; i<numrows; i++) {
		REAL randnum = (rand()&1000)/999.;
		incval(i,0) = range(i,0)*randnum;
	}
	
	std::ofstream log("conv.log");
	
#ifdef _AUTODIFF
	int icase, j;
	
	int numcases = obj.NumCases();
	
	int ncoefs = coefs.NElements();
	
	
	for(icase = 0; icase < numcases; icase++) {
		
		obj.LoadState(state);
		TPZFMatrix<REAL> ReferenceResidual, Tangent, EstimateRes;
		
		obj.ComputeTangent(Tangent,coefs,icase);
		
		// Tangent.Print("tangent matrix");
		obj.Residual(ReferenceResidual,icase);
		
		// ReferenceResidual.Print("reference residual");
		EstimateRes.Redim(ReferenceResidual.Rows(),ReferenceResidual.Cols());
		Tangent.Multiply(incval,EstimateRes);
		
		int interval;
		
		REAL difnorm[10] = {0.};
		
		// REAL resnorm[10] ={0.};
		
		for(interval = 1; interval < 10; interval++) {
			
			TPZFMatrix<REAL> actualstate(state);
			TPZFMatrix<REAL> residual;
			
			for(i=0; i<numrows; i++) {
				for(j=0; j<ncoefs; j++) {
					actualstate(i,j) += (interval/10.)*incval(i,0)*coefs[j];
				}
			}
			
			obj.LoadState(actualstate);
			obj.Residual(residual,icase);
			// residual.Print("residual");
			residual -= ReferenceResidual;
			// resnorm[interval] = Norm(residual);
			residual = residual - EstimateRes*(REAL(interval/10.));
			// residual.Print("residual adaptado");
			difnorm[interval] = Norm(residual);
			// cout << "difnorm = " << difnorm[interval];
			
		}

		std::cout << "icase = " << icase << std::endl;
		log << "icase = " << icase << std::endl;
		
		for(interval = 2; interval<10; interval++) {
			
			if(fabs(difnorm[interval]) < REAL(1.e-12) || fabs(difnorm[interval-1]) <REAL(1.e-12)) {
				std::cout << "residual too small\n";
				log << "residual too small\n";
				break;
			}
			std::cout << (log10(difnorm[interval])-log10(difnorm[interval-1]))/
			
			(log10((float)interval)-log10(interval-1.0)) << std::endl;
			log << (log10(difnorm[interval])-log10(difnorm[interval-1]))/
			(log10((float)interval)-log10(interval-1.0)) << std::endl;
		}
		
	}

#endif
	
	log.flush();

}

