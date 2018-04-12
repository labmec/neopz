/**
 * \file
 * @brief Contains the TPZMultCamada class.
 */

class TPZMatPlaca2;
class TPZBndCond;

#ifndef TPZMULTCAMADA_H
#define TPZMULTCAMADA_H
#include "TPZMaterial.h"
#include "pzstack.h"

/**
 * @ingroup material
 * @brief DESCRIBE PLEASE
 */
class TPZMultCamada : public TPZMaterial {
public:
	/** @brief Default constructor */	
	TPZMultCamada(int matindex) : 
    TPZRegisterClassId(&TPZMultCamada::ClassId),
    TPZMaterial(matindex), fCamadas() {}

	/** @brief Add layer */
    void AddLayer(TPZMatPlaca2 * l) { fCamadas.Push(l); }

	/**
	 * @name Contribute methods (weak formulation)
	 * @{
	 */
	 
    /** @brief Computes contribution to the stiffness matrix and right hand side at an integration point */
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ek,
							TPZFMatrix<STATE> & ef);
	
    virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ek,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc);
	
    virtual void Contribute(TPZMaterialData &data,
							REAL weight,
							TPZFMatrix<STATE> & ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
    }
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> & ef,
							  TPZBndCond & bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
    }

	/** @} */
	
    virtual int VariableIndex(const std::string &name);
	
    virtual int NSolutionVariables(int var);

protected:
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
	virtual void Solution(TPZVec < STATE > & Sol, TPZFMatrix<STATE> & DSol, TPZFMatrix<REAL> & axes, int var, TPZVec < STATE > & Solout);
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec < STATE > & Solout)
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

        Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
    }
	
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables();
	
	/** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const {return 2;}
    public:
virtual int ClassId() const;

private:
	/** @brief Vector of layers */
    TPZStack < TPZMatPlaca2 * > fCamadas;

};

#endif  //TPZMULTCAMADA_H
