/** @file */

#ifndef PZFUNCTION_H
#define PZFUNCTION_H

#include "TPZSavable.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/** @ingroup util */

/**
 * @ingroup util
 * @brief Implements a function. \ref util "Utility"
 * @note Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
 * @since August 01, 2007
 */
template<class TVar>
class TPZFunction : public virtual TPZSavable {
public:
	
	/** @brief Class constructor */
	TPZFunction(){
            
        }
	
	/** @brief Class destructor */
	~TPZFunction(){
            
        }
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
        DebugStop();
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate	 
	 * @param f function values
	 * @param gradf function derivatives
	 */	
	virtual void Execute(const TPZVec<REAL> &x, REAL time, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf){
        DebugStop();
    }

    /** @brief Execute method receiving axes. It is used in shape functions */
    virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
        DebugStop();
    }
    
    /** @brief Simpler version of Execute method which does not compute function derivatives */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f){
        DebugStop();
    }

    /** @brief Polynomial order of this function.
      *  In case of non-polynomial function it can be a reasonable approximation order. */
    virtual int PolynomialOrder() const {
        return 1;
    }
    
    /// number of values returned by this function
    virtual int NFunctions() const
    {
        return 1;
    }
    
    /** @brief Print a brief statement */
    virtual void Print(std::ostream &out)
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        out << "Polynomial Order = " << PolynomialOrder() << std::endl;
    }
    
public:
    virtual int ClassId() const;

    virtual void Write(TPZStream &buf, int withclassid) const { //ok
    }

    virtual void Read(TPZStream &buf, void *context) { //ok
    }
	
};

template<class TVar>
int TPZFunction<TVar>::ClassId() const{
    return Hash("TPZFunction") ^ ClassIdOrHash<TVar>()<<1;
}

template<class TVar>
class TPZDummyFunction : public TPZFunction<TVar>
{
	
    void (*fFunc)(const TPZVec<REAL> &x, TPZVec<TVar> &f);
    void (*fFunc2)(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    void (*fFunc3)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf);
    
    int fPorder;
public:
	
	/** @brief Class constructor */
	TPZDummyFunction() : TPZRegisterClassId(&TPZDummyFunction::ClassId), TPZFunction<TVar>(), fPorder(-1)
    {
        fFunc = 0;
		fFunc2 = 0;
		fFunc3 = 0;		
    }
	
	/** @brief Class destructor */
	virtual ~TPZDummyFunction()
    {
        
    }
    
    TPZDummyFunction(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val), int polynomialorder )
    : TPZRegisterClassId(&TPZDummyFunction::ClassId)
    {
        fFunc = FuncPtr;
		fFunc2 = 0;
		fFunc3 = 0;
		fPorder = polynomialorder;
    }
    
    TPZDummyFunction(void (*FuncPtr)(const TPZVec<REAL> &x, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf), int polynomialorder )
    : TPZRegisterClassId(&TPZDummyFunction::ClassId)
    {
		fFunc = 0;
        fFunc2 = FuncPtr;		
		fFunc3 = 0;
        fPorder = polynomialorder;
    }
	
    TPZDummyFunction(void (*FuncPtr)(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &val, TPZFMatrix<TVar> &gradf), int polynomialorder ) : TPZRegisterClassId(&TPZDummyFunction::ClassId)
    {
		fFunc = 0;
        fFunc2 = 0;		
		fFunc3 = FuncPtr;
        fPorder = polynomialorder;
    }	
    
    TPZDummyFunction(const TPZDummyFunction &cp) : TPZRegisterClassId(&TPZDummyFunction::ClassId), 
    fFunc(cp.fFunc), fFunc2(cp.fFunc2), fFunc3(cp.fFunc3), fPorder(cp.fPorder)
    {
        
    }
    

    TPZDummyFunction &operator=(const TPZDummyFunction &cp)
    {
        fFunc = cp.fFunc;
		fFunc2 = cp.fFunc2;
		fFunc3 = cp.fFunc3;
        fPorder = cp.fPorder;
        return *this;
    }
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
    {
        if (!fFunc2) {
			if (!fFunc)
				DebugStop();
			else {
				df.Zero();
				fFunc(x, f);
			}
		}
        else
            fFunc2(x, f, df);
    }
	
	/**
	 * @brief Performs time dependent function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param ftime  time to evaluate	 
	 * @param f function values
	 * @param gradf function derivatives
	 */	
	virtual void Execute(const TPZVec<REAL> &x, REAL ftime, TPZVec<TVar> &f, TPZFMatrix<TVar> &gradf)
    {
        if (!fFunc3) {
			DebugStop();
		}
        fFunc3(x, ftime, f, gradf);
    }	
    
	/**
	 * @brief Execute method receiving axes. It is used in shape functions
	 * @note NOT IMPLEMENTED
	 */
	virtual void Execute(const TPZVec<REAL> &x, const TPZFMatrix<REAL> &axes, TPZVec<TVar> &f, TPZFMatrix<TVar> &df){
        DebugStop();
    }
    
    /**
	 * @brief Simpler version of Execute method which does not compute function derivatives 
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 */
    virtual void Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f){
		if (!fFunc) {
			DebugStop();
		}
        fFunc(x,f);
    }
    
	/** @brief Returns number of functions. */ 
	virtual int NFunctions() const
    {
        return 1;
    }
    
    void SetPolynomialOrder(int porder)
    {
        fPorder = porder;
    }
	
	/** @brief Polynomial order of this function. */
	/** In case of non-polynomial function it can be a reasonable approximation order. */
	virtual int PolynomialOrder() const
    {
#ifdef PZDEBUG
        if (fPorder == -1) {
            DebugStop();
        }
#endif
        return fPorder;
    }
	
	/** @brief Unique identifier for serialization purposes */
	public:
virtual int ClassId() const;

	
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const
    {
        DebugStop();
    }
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context)
    {
        DebugStop();
    }
};

template<class TVar>
int TPZDummyFunction<TVar>::ClassId() const{
    return TPZFunction<TVar>::ClassId() ^ Hash("TPZDummyFunction");
}

#endif
