/**
 * \file
 * @brief Contains the TPZNLMat1dRotatedEngStrain class which implements a non linear 1d material.
 */
//$Id: pznlmat1drotatedengstrain.h,v 1.7 2009-11-16 18:41:59 diogo Exp $
// -*- c++ -*-

#ifndef TPZNLMAT1DROTATEDENGSTRAIN_H
#define TPZNLMAT1DROTATEDENGSTRAIN_H

#include "pznlmat1d.h"
#include <fstream>
#include <iostream>

/**
 * @ingroup material
 * @brief Implements a non linear 1d material based on a rotated engineering strain measurement.
 * @author Edimar Cesar Rylo
 * @since May, 2006
 */
/**
 * The implementation is based on the section 3.1.1 of the book Non-linear \n
 * Finite Element Analysis of Solids and Structures - Volume 1: Essentials of M. A. Crisfield
 */
class TPZNLMat1dRotatedEngStrain : public TPZNLMat1d
{
public:
    /** @brief Simple Constructor */
    TPZNLMat1dRotatedEngStrain(int id);
	
    /** @brief Default destructor */
    virtual ~TPZNLMat1dRotatedEngStrain();
	
    virtual std::string Name() { return "nonlinear_RotatedEnginneringStrain"; }
	
    /** @brief Returns the integrable dimension of the material: Material is 1d */
    virtual int Dimension() {return  1;}
	
    /** @brief Returns the number of state variables associated with the material: Only w? */
    virtual int NStateVariables() {return  1;}
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout);
	
    /** @brief Returns the variable index associated with the name */
	virtual int VariableIndex(const std::string &name);
	
    virtual int NSolutionVariables(int var);
	
protected:
	virtual void Solution(TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol,
						  TPZFMatrix<REAL> &axes, int var, TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
    }
	
	
    virtual void Contribute(TPZMaterialData &data,
                            REAL weight,
                            TPZFMatrix<REAL> &ek,
							TPZFMatrix<REAL> &ef);
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<REAL> &ef);
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ek,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc);
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<REAL> &ef,
							  TPZBndCond &bc)
	{
		TPZNLMat1d::ContributeBC(data,weight,ef,bc);
	}
	
	/** @brief To create another material of the same type */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief Reads data of the material from a istream (file data) */
	virtual void SetData(std::istream &data);

    virtual void Write(TPZStream &buf, int withclassid);
	
    virtual void Read(TPZStream &buf, void *context);
	
    virtual REAL Eps(TPZVec<REAL> &sol,TPZFMatrix<REAL> &axes,TPZFMatrix<REAL> &dphi) = 0;
};

#endif
