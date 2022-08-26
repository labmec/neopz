/**
 * @file TPZL2ProjectionHCurl.h
 * @brief Contains the TPZL2Projection class which implements an L2 projection of a given scalar solution
 */

#ifndef TPZL2PROJECTIONHCURL_H
#define TPZL2PROJECTIONHCURL_H

#include "Projection/TPZL2Projection.h"


template<class TVar=STATE>
class TPZL2ProjectionHCurl : public TPZL2Projection<TVar> {

public:

    TPZL2ProjectionHCurl() = default;

    TPZL2ProjectionHCurl(int id, int dim, int nstate=1);
	
	TPZL2ProjectionHCurl(int id, int dim, int nstate, const TPZVec<TVar> &sol);

    ~TPZL2ProjectionHCurl(){};


    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    /** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
	int NSolutionVariables(int var) const override;

};

#endif
