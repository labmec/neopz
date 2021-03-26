/**
 * @file
 * @brief Contains the TPZViscoelastic class which implements an isotropic viscoelasticity material.
 * @author Nathan Shauer
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
#include "TPZMatWithMem.h"

/**
 * @ingroup material
 * @brief This class implements an isotropic viscoelasticity material.
 * @author Nathan Shauer
 * @since 7/16/2010.
 */
class TPZViscoelastic : public TPZMatWithMem<TPZFMatrix<STATE>, TPZElasticity3D>
{
		
public:
	enum SOLUTIONVARS{ENone = -1, EViscoStressX = 30, EViscoStressY = 31, EViscoStressZ = 32};
  
  /**
	 * @brief Empty constructor
	 */
	TPZViscoelastic();
	
	/**
	 * @brief id constructor
	 */
	TPZViscoelastic(int id);

  /**
	 * @brief Not a good variables initialization constructor. Uses Hooke for elastic part and lame for viscous part
	 */
	
	TPZViscoelastic(int id,STATE ElaE, STATE poissonE, STATE lambdaV, STATE muV, STATE alpha, STATE deltaT, TPZVec <STATE> &force);

  /**
	 * @brief Set material Data with hooke constants
	 */
	
	void SetMaterialDataHooke(STATE ElaE, STATE poissonE, STATE ElaV, STATE poissonV, STATE alpha, STATE deltaT, TPZVec <STATE> &force);
	
	
	virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> &ek,
							TPZFMatrix<STATE> &ef) override;
	
	void UpdateQsi(TPZMaterialData &data);
	
	/** @brief Returns index of post-processing variable */
	virtual int VariableIndex(const std::string &name) override;
	
	/** @brief Number of data of variable var */
	virtual int NSolutionVariables(int var) override;
	
	/*
	 Computes the stress. 
	 Remember you cant update qsi if you want to calculate stress 
	 */
	virtual void ComputeStressTensor(TPZFMatrix<STATE> &Stress, TPZMaterialData &data) const override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;

	/** @brief Fill material data parameter with necessary requirements for the Contribute method. */
	virtual void FillDataRequirements(TPZMaterialData &data) override;
	
	/** @brief Saves the element data to a stream */
	void Write(TPZStream &buf, int withclassid) const override;
	
	/** @brief Reads the element data from a stream */
	void Read(TPZStream &buf, void *context) override;

	public:
int ClassId() const override;

	
protected:
	
	/**
	 * @brief the stability coeficient alpha
	 */
	STATE fAlpha;
	
	/**
	 * @brief the stability coeficient alpha
	 */
	STATE fDeltaT;

	/**
	 * @brief viscoelasticity coeficients
	 */
	STATE fLambdaV;
	STATE fmuV;
};

#endif