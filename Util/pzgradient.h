/**
 * @file
 * @brief Contains the TPZGradient class which implements the methods to reconstruction gradient 
 */

#ifndef PZGRADIENTH 
#define PZGRADIENTH

#include "pzfunction.h"

#include <iostream>


class TPZGradient : public TPZFunction<STATE>
{
    TPZManVector<REAL,3> fCenter;
    TPZFNMatrix<3,STATE> fGradient;
    STATE fUc;
    
public:
    TPZGradient();
    
    TPZGradient(const TPZGradient &cp);
    
    ~TPZGradient() {
    }
    
    TPZGradient &operator=(const TPZGradient &copy) {
        fCenter = copy.fCenter;
        fGradient = copy.fGradient;
        fUc = copy.fUc;
        return *this;
	}
    
    /*
     *@brief Enter the data of the reconstruction gradient
     *@param center value of the coordinate at the center of the element
     *@param grad gradient reconstructed
     *@param u0 value of the approximate solution (by FEM) at the center point
     */
    void SetData(TPZManVector<REAL,3> &center, TPZFMatrix<STATE> &grad, STATE u0){
     
        fCenter = center;
        fGradient = grad;
        fUc = u0;
    }
    
    /*
     *@brief fill the function used in the l2 projection
     *@param pt points where the function is calculated
     *@param f function projected into finite element space
     */
    virtual void Execute(const TPZVec<REAL> &pt, TPZVec<STATE> &f){
        
        int nr = fGradient.Rows();
        int nc = fGradient.Cols();
        for(int j=0; j<nc; j++){
            f[j] = fUc;
            for(int i = 0; i<nr; i++){
                f[j] += (pt[i] - fCenter[i])*fGradient(i,j);
            }
        }
    }
    
    virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
        DebugStop();
    }
    
    /** @brief Execute method receiving axes. It is used in shape functions */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<STATE> &f, TPZFMatrix<STATE> &df){
        DebugStop();
    }
    
    /**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) {
        DebugStop();
    }
    
    
	/** @brief Returns number of functions. */
	virtual int NFunctions(){return 1;}
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder() {return -1;}
    
    };
    
#endif
