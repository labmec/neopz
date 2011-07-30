/**
 * \file
 * @brief Contains the TPZMultCamada class.
 */
/* Generated by Together */

class TPZMatPlaca2;
class TPZBndCond;

#ifndef TPZMULTCAMADA_H
#define TPZMULTCAMADA_H
#include "pzmaterial.h"
#include "pzstack.h"

/**
 * @ingroup material
 * @brief 
 */
class TPZMultCamada : public TPZMaterial {
public:
	
	TPZMultCamada(int matindex) : TPZMaterial(matindex), fCamadas() {}
	
    void AddLayer(TPZMatPlaca2 * l) { fCamadas.Push(l); }
	
    /** Compute contribution to the stiffness matrix and right hand side at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix & ek,
							TPZFMatrix & ef);
	
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix & ek,
							  TPZFMatrix & ef,
							  TPZBndCond & bc);
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix & ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
    }
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix & ef,
							  TPZBndCond & bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }
	
	/** returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
	
    /** returns the number of variables associated with the variable indexed by var. var is obtained by calling VariableIndex */
    virtual int NSolutionVariables(int var);
protected:
	/** returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec < REAL > & Sol, TPZFMatrix & DSol, TPZFMatrix & axes, int var, TPZVec < REAL > & Solout);
public:
	/** returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec < REAL > & Solout)
	{
        Solution(data.sol,data.dsol,data.axes,var,Solout);
    }
	
    /** returns the number of state variables associated with the material*/
    virtual int NStateVariables();
	
	/**returns the integrable dimension of the material*/
    virtual int Dimension() {return 2;}
	
	
private:
    TPZStack < TPZMatPlaca2 * > fCamadas;
};


#endif  //TPZMULTCAMADA_H
