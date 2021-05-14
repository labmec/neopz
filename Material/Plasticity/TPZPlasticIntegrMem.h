// $Id: TPZPlasticIntegrMem.h,v 1.4 2009-06-29 22:54:00 erick Exp $

#ifndef TPZPLASTICINTEGRMEM_H
#define TPZPLASTICINTEGRMEM_H

#include "fadType.h" 
#include "TPZTensor.h"
#include "pzmanvector.h"
#include <iostream>


/**
  This class stores the whole content of a Plastic Integration Step in order to allow its re-evaluation
 */
template <class T, int N>
class TPZPlasticIntegrMem 
{

public:
	
	/**
	 * Default constructor - all values set to zero
	 */
    TPZPlasticIntegrMem(): m_elastoplastic_state(), fK(0.), fLambda(0. ), fDelGamma(N,T(0.)), fValidEqs(N,0), fForceYield(0)
	{ }

	/**
	 * Full constructor - all values set to zero
	 */
    TPZPlasticIntegrMem(const TPZPlasticState<T> & state,
						const REAL & k,
						const REAL & lambda,
						const TPZManVector<T,N> & delGamma,
						const TPZManVector<int, N> & validEqs,
						const int forceYield):
		    m_elastoplastic_state(state), fK(k), fLambda(lambda), fDelGamma(delGamma), fValidEqs(validEqs), fForceYield(forceYield)
	{ }
	
	/**
	 * Default destructor
	 */
	~TPZPlasticIntegrMem(){ }
	
	/**
	 * Copy constructor
	 */
	TPZPlasticIntegrMem(const TPZPlasticIntegrMem<T, N> & source):
		m_elastoplastic_state(source.m_elastoplastic_state),
	    fK(source.fK), fLambda(source.fLambda), fDelGamma(source.fDelGamma),
	    fValidEqs(source.fValidEqs), fForceYield(source.fForceYield)
	{ }
	
	/**
	 * Operator =
	 */
	TPZPlasticIntegrMem<T, N> & operator=(const TPZPlasticIntegrMem<T, N> & source)
	{
		m_elastoplastic_state = source.m_elastoplastic_state;
		fK = source.fK;
		fLambda = source.fLambda;
		fDelGamma = source.fDelGamma;
		fValidEqs = source.fValidEqs;
		fForceYield = source.fForceYield;
		
		return *this;
	}

	
public:	
	// Tensors to hold the total and plastic strain tensors
	TPZPlasticState<T> m_elastoplastic_state;
	
	// DeltaEpsT multiplier in the integration step, Plastic damage variable
	REAL fK, fLambda;
	
	// Plastic flow multiplier
	TPZManVector<T,N> fDelGamma;
	
	// set of valid Plastic flow equations
	TPZManVector<int,N> fValidEqs;
	
	// whether to force post peak yield behavior
	int fForceYield;
	
};
#endif
