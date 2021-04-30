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
    //center of the element K
    TPZManVector<REAL,3> fCenter;
    
    //gradient estimate on the element K
    TPZManVector<STATE,3> fGradient;
    
    //cell averaged value
    STATE fUc;
    
    //Slope limiter
    STATE falphaK;
        
public:
    TPZGradient();
    
    TPZGradient(const TPZGradient &cp);
    
    ~TPZGradient() {
    }
    
    TPZGradient &operator=(const TPZGradient &copy){
        fCenter = copy.fCenter;
        fGradient = copy.fGradient;
        fUc = copy.fUc;
        falphaK = copy.falphaK;
        return *this;
	}
    
    /*
     *@brief Enter the data of the reconstruction gradient
     *@param center value of the coordinate at the center of the element
     *@param grad gradient reconstructed
     *@param u0 value of the approximate solution (by FEM) at the center point (cell averaged)
     *@param alphak: value of the slope limiter
     */
    void SetData(TPZManVector<REAL,3> &center, TPZManVector<STATE,3> &grad, STATE u0, STATE alphak){
     
        fCenter = center;
        fGradient = grad;
        fUc = u0;
        falphaK = alphak;
    }
    
    /*
     *@brief fill the function used in the l2 projection
     *@param pt points where the function is calculated
     *@param f function projected into finite element space
     */
    virtual void Execute(const TPZVec<REAL> &pt, TPZVec<STATE> &f) override {
        
        int nr = fGradient.size();
        
        f[0] = fUc;
        for(int i = 0; i<nr; i++)
        {
            f[0] += (STATE)(pt[i] - fCenter[i])*falphaK*fGradient[i];
        }
    }
    
    virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf) override{
        DebugStop();
    }
    
    /** @brief Execute method receiving axes. It is used in shape functions */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) override{
        DebugStop();
    }
    
    /**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) override {
        DebugStop();
    }
    
    
	/** @brief Returns number of functions. */
	virtual int NFunctions()const override {return 1;}
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder() const override {return -1;}
    
    };
    
#endif
