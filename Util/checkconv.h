/**
 * @file
 * @brief Contains the implementation of the CheckConvergence function.
 */

#include <stdio.h>
#include <stdlib.h>

#include <fstream>
#include <pzfmatrix.h>

template <class TConv>

/**
 * @brief Implements a general procedure to check whether the class TConv implements a consistente tangent matrix for the object obj.
 */
/**
 * The following methods need to be implemented for TConv 
 * \li LoadState : to load the state variables in local storage \n
 * \li ComputeTangent : which will compute a tangent matrix \n
 * \li Residual : which will compute a residual vector
 */
void CheckConvergence(TConv &obj, TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &range, TPZVec<REAL> &coefs) {
    
	TPZFMatrix<STATE> incval(range);
	int i,numrows;
	
	numrows = state.Rows();
	for(i=0; i<numrows; i++) {
		REAL randnum = (rand()&1000)/999.;
		incval(i,0) = range(i,0)*(STATE)randnum;
	}
	
	std::ofstream log("conv.log");
	
	int icase, j;
	
	int numcases = obj.NumCases();
	
	int ncoefs = coefs.NElements();
	
	
	for(icase = 0; icase < numcases; icase++) {
		
		obj.LoadState(state);
		TPZFMatrix<STATE> ReferenceResidual, Tangent, EstimateRes;
		
		obj.ComputeTangent(Tangent,coefs,icase);
		
		// Tangent.Print("tangent matrix");
		obj.Residual(ReferenceResidual,icase);
		
      //  ReferenceResidual.Print("reference residual");
		EstimateRes.Redim(ReferenceResidual.Rows(),ReferenceResidual.Cols());
		Tangent.Multiply(incval,EstimateRes);
		
		int interval;
		
		REAL difnorm[10] = {0.};
		
		for(interval = 1; interval < 10; interval++) {
			
			TPZFMatrix<STATE> actualstate(state);
			TPZFMatrix<STATE> residual,temp;
			
			for(i=0; i<numrows; i++) {
				for(j=0; j<ncoefs; j++) {
					actualstate(i,j) += (STATE)(interval/10.)*(STATE)incval(i,0)*(STATE)coefs[j];
				}
			}
			
			obj.LoadState(actualstate);
			obj.Residual(residual,icase);
			
			residual -= ReferenceResidual;
         //   residual.Print("residual-referenceresidual");
			// resnorm[interval] = Norm(residual);
			residual = residual - EstimateRes*(STATE(interval/10.));
            temp = EstimateRes*(STATE(interval/10.));
        //    std::cout << "\n EstimateRes*(STATE(interval/10.))  = " <<temp << std::endl;
			// residual.Print("residual adaptado");
			difnorm[interval] = Norm(residual);
       //     std::cout << "\n difnorm = " << difnorm[interval]<< std::endl;
			
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
	
	log.flush();
}

/*
 
 templates for the interface a TPZCheckConvergence class must specify
 /// Number of test cases implemented by this class
 int NumCasesCheckConv() const;
 /// Compute the tangent matrix for a particular case
 void TangentCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &tangent, int icase);
 /// Compute the residual for a particular state
 void ResidualCheckConv(TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &residual, int icase);
 */

/**
 * The following methods need to be implemented for TConv 
 * \li LoadStateCheckConv : to load the state variables in local storage \n
 * \li TangentCheckConv : which will compute a tangent matrix \n
 * \li ResidualCheckConv : which will compute a residual vector \n
 * \li NumCasesCheckConv() : number of test cases implemented
 */
template<class TConv>
void TPZCheckConvergence(TConv &obj, TPZFMatrix<STATE> &state, TPZFMatrix<STATE> &range, bool direction_range = false) {
    
	TPZFMatrix<STATE> incval(range);
	int numrows;
	
	numrows = state.Rows();
    REAL randnum = (rand()&1000)/999.;
	for(int i=0; i<numrows; i++) {
		incval(i,0) = range(i,0)*(STATE)randnum;
		if(!direction_range) randnum = (rand()&1000)/999.;
	}
	
	std::ofstream log("conv.log");
	
	int icase;
	
	int numcases = obj.NumCasesCheckConv();
	
	
	
	for(icase = 0; icase < numcases; icase++) {
		
		TPZFMatrix<STATE> ReferenceResidual, Tangent, EstimateRes;
		
		obj.TangentCheckConv(state, Tangent,icase);
		
		// Tangent.Print("tangent matrix");
		obj.ResidualCheckConv(state, ReferenceResidual,icase);
		
		// ReferenceResidual.Print("reference residual");
		EstimateRes.Redim(ReferenceResidual.Rows(),ReferenceResidual.Cols());
		Tangent.Multiply(incval,EstimateRes);
				
		REAL difnorm[10] = {0.};
		
		for(int interval = 1; interval < 10; interval++) 
        {
			
			TPZFMatrix<STATE> actualstate(state);
			TPZFMatrix<STATE> residual;
			
			for(int i=0; i<numrows; i++) 
            {
                actualstate(i) += (STATE)(interval/10.)*(STATE)incval(i,0);
			}
			
			obj.ResidualCheckConv(actualstate, residual,icase);
			// residual.Print("residual");
			residual -= ReferenceResidual;
			// resnorm[interval] = Norm(residual);
			residual = residual - EstimateRes*(STATE(interval/10.));
			// residual.Print("residual adaptado");
			difnorm[interval] = Norm(residual);
            std::cout << "difnorm = " << difnorm[interval];
			
		}
        
		std::cout << "icase = " << icase << std::endl;
		log << "icase = " << icase << std::endl;
		
		for(int interval = 2; interval<10; interval++) {
			
			if(fabs(difnorm[interval]) < REAL(1.e-12) || fabs(difnorm[interval-1]) <REAL(1.e-12)) {
				std::cout << "residual too small\n";
				log << "residual too small\n";
				break;
			}
			std::cout << (log10(difnorm[interval])-log10(difnorm[interval-1]))/
			
			(log10((REAL)interval)-log10(interval-1.0)) << std::endl;
			log << (log10(difnorm[interval])-log10(difnorm[interval-1]))/
			(log10((REAL)interval)-log10(interval-1.0)) << std::endl;
		}
		
	}	
	log.flush();
}
