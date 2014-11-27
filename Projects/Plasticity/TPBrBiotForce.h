//
//  TPBrBiotForce.h
//  PZ
//
//  Created by Philippe Devloo on 5/13/14.
//
//

#ifndef __PZ__TPBrBiotForce__
#define __PZ__TPBrBiotForce__

#include <iostream>
#include "pzfunction.h"


class TPBrBiotForce : public TPZFunction<STATE>
{
    /// raio do poco
    STATE fRwell;
    /// raio do reservatorio
    STATE fRreservoir;
    /// pressao do poco
    STATE fPwell;
    /// pressao do reservatorio
    STATE fPreservoir;
    /// constante de Biot
    STATE fBiot;
    /// constant combination of above parameters
    STATE fConstant;
    
    /// auxiliary method
    void ComputeConstant()
    {
        fConstant = fBiot*(fPreservoir-fPwell)/log(fRreservoir/fRwell);
    }
    
public:
    /** @brief Class constructor */
    TPBrBiotForce() : fRwell(1.), fRreservoir(1.), fPwell(1.), fPreservoir(10.), fBiot(1.)
    {
        ComputeConstant();
    }

  /** @brief Class copy constructor */
  TPBrBiotForce(const TPBrBiotForce &cp) : TPZFunction<STATE>(cp), fRwell(cp.fRwell), fRreservoir(cp.fRreservoir), fPwell(cp.fPwell), fPreservoir(cp.fPreservoir), fBiot(cp.fBiot), fConstant(cp.fConstant)
  {

  }
  
    /** @brief Class copy constructor */
    TPBrBiotForce & operator=(const TPBrBiotForce &cp)
    {
        TPZFunction<STATE>::operator=(cp);
        fRwell = cp.fRwell;
        fRreservoir = cp.fRreservoir;
        fPwell = cp.fPwell;
        fPreservoir = cp.fPreservoir;
        fBiot = cp.fBiot;
        fConstant = cp.fConstant;
        return *this;
    }
  
	/** @brief Class destructor */
	virtual ~TPBrBiotForce()
    {
        
    }
    
    void SetConstants(STATE Rwell, STATE Rreservoir, STATE Pwell, STATE Preservoir, STATE biot)
    {
        fRwell = Rwell;
        fRreservoir = Rreservoir;
        fPwell = Pwell;
        fPreservoir = Preservoir;
        fBiot = biot;
        ComputeConstant();
    }
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values returns the gradient of the pressure multiplied by biot
	 * @param df pressure function multiplied by biot
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df)
    {
        REAL r2 = x[0]*x[0]+x[1]*x[1];
        f[0] = -fConstant*x[0]/r2;
        f[1] = -fConstant*x[1]/r2;
        
        df.Redim(1, 1);
        df(0,0) = fBiot*(fPwell+fConstant*log(r2/fRwell/fRwell)/2.);
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate
	 * @param f function values
	 * @param gradf function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<STATE> &f, TPZFMatrix<STATE> &gradf){
        DebugStop();
    }
    
    /** @brief Execute method receiving axes. It is used in shape functions */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<STATE> &f, TPZFMatrix<STATE> &df){
        DebugStop();
    }
    
    /** @brief Simpler version of Execute method which does not compute function derivatives */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f){
        REAL r2 = x[0]*x[0]+x[1]*x[1];
        f[0] = -fConstant*x[0]/r2;
        f[1] = -fConstant*x[1]/r2;
    }
    
	/** @brief Returns number of functions. */
	virtual int NFunctions(){
        return 2;
    }
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder()
    {
        return 4;
    }
    
    /** @brief Print a brief statement */
    virtual void Print(std::ostream &out)
    {
        out << __PRETTY_FUNCTION__ << (void *) this << std::endl;
        /// pressao do poco
        out << "Pressao do poco " <<  fPwell << std::endl;
        /// pressao do reservatorio
        out << "Pressao dos poros " << fPreservoir << std::endl;
        /// raio do poco
        out << "Raio do poco " <<  fRwell << std::endl;
        /// raio do reservatorio
        out << "Raio do reservatorio " <<  fRreservoir << std::endl;
        /// constante de Biot
        out << "Constante de biot " <<  fBiot << std::endl;
        /// constant combination of above parameters
        out << "Constante " << fConstant << std::endl;
        out << "NFunctions = " << NFunctions() << std::endl;
        out << "Polynomial Order = " << PolynomialOrder() << std::endl;
    }
    
    /** @brief Define the class id associated with the class */
	/**
	 * This id has to be unique for all classes
	 * A non unique id is flagged at the startup of the program
	 */
	virtual int ClassId() const ;
	
	/** @brief Writes this object to the TPZStream buffer. Include the classid if withclassid = true */
	virtual void Write(TPZStream &buf, int withclassid) ;
	
	/** @brief read objects from the stream */
	virtual void Read(TPZStream &buf, void *context);
	


};

#endif /* defined(__PZ__TPBrBiotForce__) */
