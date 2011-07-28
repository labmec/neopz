//$Id: pzfunction.h,v 1.3 2008-02-05 22:23:27 tiago Exp $

#ifndef PZFUNCTION_H
#define PZFUNCTION_H

#include "pzsave.h"
#include "pzvec.h"
#include "pzfmatrix.h"

/** @ingroup util */
const int TPZFUNCTIONID = 9000;

/**
 * @ingroup util
 * @brief Implements a function. 
 *
 * Its purpose is to replace void* fp() instances which are difficult to be transmitted in parallel executions.
 * @since August 01, 2007
 */
class TPZFunction : public TPZSaveable{
public:
	
	/**
	 * @brief Class constructor
	 */
	TPZFunction();
	
	/**
	 * @brief Class destructor
	 */
	~TPZFunction();
    
	/**
	 * @brief Performs function computation
	 * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
	 * @param f function values
	 * @param df function derivatives
	 */
	virtual void Execute(const TPZVec<REAL> &x, TPZVec<REAL> &f, TPZFMatrix &df) = 0;
	
	/** @brief Returns number of functions.
	 */ 
	virtual int NFunctions() = 0;
	
	/** @brief Polynomial order of this function. In case of non-polynomial
	 * function it can be a reasonable approximation order.
	 */
	virtual int PolynomialOrder() = 0;
	
	/**
	 * @brief Unique identifier for serialization purposes
	 */
	virtual int ClassId() const;
	
	/**
	 * @brief Save the element data to a stream
	 */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/**
	 * @brief Read the element data from a stream
	 */
	virtual void Read(TPZStream &buf, void *context);
	
};

#endif
