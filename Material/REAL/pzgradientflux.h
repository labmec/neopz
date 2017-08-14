/**
 * \file
 * @brief Contains the TPZGradientFlux class.
 */

#ifndef TPZGRADIENTFLUX_H
#define TPZGRADIENTFLUX_H

#include "pzreal.h"  // for STATE, REAL
#include "pzvec.h"   // for TPZVec
class TPZMaterialData;

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZGradientFlux{
	
public:
	/** @brief Default constructor */
	TPZGradientFlux();
	/** @brief Copy constructor */
	TPZGradientFlux(const TPZGradientFlux &cp);
	/** @brief Destructor */
	~TPZGradientFlux();
	
	/** @brief Computes numerical flux */
	void ComputeFlux(TPZVec<STATE> &solL, TPZVec<STATE> &solR, const TPZVec<REAL> &normal, TPZVec<STATE> & F);
	
	/** @brief Apply limiter */
	void ApplyLimiter(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright);
	
private:
	
	void ApplyVanAlbadaLimiter(STATE &soll, STATE &solr,
							   const TPZVec<STATE>& gradL, const TPZVec<STATE> &gradR, 
							   const TPZVec<STATE> &normal, 
							   const TPZVec<STATE> &dL, const TPZVec<STATE> & dR);
	
	/** @brief It corrects \f$ soll \f$ and \f$ solr \f$ values */
	void ApplyMinModLimiter(STATE &soll, STATE &solr,
							const TPZVec<STATE>& gradL, const TPZVec<STATE> &gradR, 
							const TPZVec<STATE> &normal, 
							const TPZVec<STATE> &dL, const TPZVec<STATE> & dR);

	/** @brief Computes the dot product (scalar) */
	STATE Dot(const TPZVec<STATE> &A, const TPZVec<STATE> &B) {
		double result = 0.;
		int n = A.NElements();
		for(int i = 0; i < n; i++) result += A[i]*B[i];
		return result;
	}
	
};

#endif
