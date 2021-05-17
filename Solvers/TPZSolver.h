/**
 * @file
 * @brief Contains TPZSolver class which defines a abstract class of solvers which will be used by matrix classes.
 */

#ifndef TPZSOLVER_H
#define TPZSOLVER_H

#include "TPZSavable.h"
/**
 * @ingroup solver
 * @brief Defines a abstract class of solvers which will be used by matrix classes.
 */
class TPZSolver: public TPZSavable
{
public:
    //! Defines the ClassId of TPZSolver class
    int ClassId() const override;
	
	//! Clones the current object returning a pointer of type TPZSolver
	virtual TPZSolver *Clone() const = 0;

	//! Destructor
	virtual ~TPZSolver() = default;
	
	/*! Resets the matrix associated with the solver. This is useful when the matrix needs to be recomputed in a non linear problem*/
	virtual void ResetMatrix() = 0;
};

#endif