//$Id: pzfunction.h,v 1.3 2008-02-05 22:23:27 tiago Exp $

#ifndef PZFUNCTION_H
#define PZFUNCTION_H

#include "pzsave.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/** @ingroup util */
/** @brief TPZFunction class Id */
const int TPZFUNCTIONID = 9000;

/**
 * @ingroup util
 * @brief Implements a function. \ref util "Utility"
 * @note Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
 * @since August 01, 2007
 */
class TPZFunction : public TPZSaveable{
public:
	
	/** @brief Class constructor */
	TPZFunction();
	
	/** @brief Class destructor */
	~TPZFunction();
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df) = 0;

    /** Versao do Execute recebendo axes. Utilizado em shape functions.
     */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<REAL> &f, TPZFMatrix<REAL> &df){
        DebugStop();
    }
    
    /** Simpler version of Execute method which does not compute function derivatives
     */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f){
        DebugStop();
    }
    
	/** @brief Returns number of functions. */ 
	virtual int NFunctions() = 0;
	
	/** @brief Polynomial order of this function. \n 
	 * In case of non-polynomial function it can be a reasonable approximation order.
	 */
	virtual int PolynomialOrder() = 0;
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const;
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
};

class TPZDummyFunction : public TPZFunction
{

    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<REAL> &f);
public:
	
	/** @brief Class constructor */
	TPZDummyFunction()
    {
        fFunc = 0;
    }
	
	/** @brief Class destructor */
	virtual ~TPZDummyFunction()
    {
        
    }
    
    TPZDummyFunction(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<REAL> &val))
    {
        fFunc = FuncPtr;
    }
    
    TPZDummyFunction(const TPZDummyFunction &cp) : fFunc(cp.fFunc)
    {
        
    }
    
    TPZDummyFunction &operator=(const TPZDummyFunction &cp)
    {
        fFunc = cp.fFunc;
        return *this;
    }
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix<REAL> &df)
    {
        DebugStop();
    }
    
    /** Versao do Execute recebendo axes. Utilizado em shape functions.
     */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<REAL> &f, TPZFMatrix<REAL> &df){
        DebugStop();
    }
    
    /** Simpler version of Execute method which does not compute function derivatives
     */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f){
        fFunc(x,f);
    }
    
	/** @brief Returns number of functions. */ 
	virtual int NFunctions()
    {
        return 1;
    }
	
	/** @brief Polynomial order of this function. \n 
	 * In case of non-polynomial function it can be a reasonable approximation order.
	 */
	virtual int PolynomialOrder() 
    {
        return -1;
    }
	
	/** @brief Unique identifier for serialization purposes */
	virtual int ClassId() const
    {
        return -1;
    }
	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid)
    {
        DebugStop();
    }
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }

    
};

#endif
