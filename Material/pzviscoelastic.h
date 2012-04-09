/**
 * \file
 * @brief Contains the TPZViscoelastic class which implements an isotropic viscoelasticity material.
 */
#ifndef TPZVISCOELASTIC_H
#define TPZVISCOELASTIC_H

/** @brief First index to qsi vector for contribute method of viscoelasticity material */
const int _XX_ = 0;
/** @brief Second index to qsi vector for contribute method of viscoelasticity material */
const int _XY_ = 1;
/** @brief Third index to qsi vector for contribute method of viscoelasticity material */
const int _XZ_ = 2;
/** @brief Fourth index to qsi vector for contribute method of viscoelasticity material */
const int _YY_ = 3;
/** @brief Fifth index to qsi vector for contribute method of viscoelasticity material */
const int _YZ_ = 4;
/** @brief Sixth index to qsi vector for contribute method of viscoelasticity material */
const int _ZZ_ = 5;


#include <iostream>
#include "pzfmatrix.h"
#include "pzvec.h"
#include <vector>

#include "pzelast3d.h"
#include "pzmatwithmem.h"

/**
 * @ingroup material
 * @brief This class implements an isotropic viscoelasticity material.
 * @author Nathan Shauer
 * @since 7/16/2010.
 */
class TPZViscoelastic : public TPZMatWithMem<TPZFMatrix<REAL>, TPZElasticity3D>
{
		
public:
	enum SOLUTIONVARS{ENone = -1, EViscoStressX = 30, EViscoStressY = 31, EViscoStressZ = 32};
	
	TPZViscoelastic(int id,REAL ElaE, REAL poissonE, REAL lambdaV, REAL muV, REAL alphaT, TPZVec <REAL> &force);
	
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ek,
							TPZFMatrix<REAL> &ef);
	
	void UpdateQsi(TPZMaterialData &data);
	
	/** @brief Returns index of post-processing variable */
	virtual int VariableIndex(const std::string &name);
	
	/** @brief Number of data of variable var */
	virtual int NSolutionVariables(int var);
	
	/*
	 Computes the stress. 
	 Remember you cant update qsi if you want to calculate stress 
	 */
	virtual void ComputeStressTensor(TPZFMatrix<REAL> &Stress, TPZMaterialData &data) const;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);


	
protected:
	
	REAL flambdaE,fmuE,flambdaV,fmuV,falphaT,fElaVE,fPoissonVE;  
	//int fMemory;
	
};
#endif