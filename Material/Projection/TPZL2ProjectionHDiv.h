/**
 * @file TPZL2ProjectionHDiv.h
 * @brief Contains the TPZL2Projection class which implements an L2 projection of a given scalar solution
 */

#ifndef TPZL2PROJECTIONHDIV_H
#define TPZL2PROJECTIONHDIV_H

#include "Projection/TPZL2Projection.h"


template<class TVar=STATE>
class TPZL2ProjectionHDiv : public TPZL2Projection<TVar> {

public:

    TPZL2ProjectionHDiv() = default;

    TPZL2ProjectionHDiv(int id, int dim, int nstate=1);
	
	TPZL2ProjectionHDiv(int id, int dim, int nstate, const TPZVec<TVar> &sol);

    ~TPZL2ProjectionHDiv(){};


    void Contribute(const TPZMaterialDataT<TVar> &data,
                    REAL weight,
                    TPZFMatrix<TVar> &ek, TPZFMatrix<TVar> &ef) override;

    /** @brief It returns the variable index associated with the name */
	int VariableIndex(const std::string &name) const override;
	
	int NSolutionVariables(int var) const override;

};

#endif
