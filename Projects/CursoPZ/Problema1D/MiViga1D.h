/**
 * @file
 */
#ifndef TPZMIVIGAHH
#define TPZMIVIGAHH

#include "TPZMaterial.h"

/// @brief Clase para un material viga uni-dimensional
class TPZMiViga1D  : public TPZMaterial
{
	/** @brief  Elasticity coefficient */
	REAL fE;
	/** @brief Modulo de Poisson */
	REAL fNi;
	REAL fKsi;
	REAL fArea;
	/** @brief Momento de inercia */
	REAL fI;

	REAL fG;

	/** @brief Carga aplicada sobre la viga */
	REAL fCarga;
    
public:
	/** @brief Simple constructor with material id and dimension of the spatial domain */
	TPZMiViga1D(int nummat,REAL E,REAL Poisson,REAL Inercia,REAL Ksi = 5./6.,REAL Area = 1);
	
	/** @brief Default destructor */
	~TPZMiViga1D() { }

    /** @brief Returns the name of the material */
    virtual std::string Name() { return "TPZMiViga1D"; }

	void SetLoad(double carga) { fCarga = carga; }
	
	// MUST TO BE IMPLEMENTED
	
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const { return 1; }

    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const { return 2; }

	virtual int NSolutionVariables(int var);


    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param[in] data stores all input data
     * @param[in] weight is the weight of the integration rule
     * @param[out] ek is the stiffness matrix
     * @param[out] ef is the load vector
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data stores all input data
     * @param weight is the weight of the integration rule
     * @param ek is the stiffness matrix
     * @param ef is the load vector
	 * @param bc boundary condition
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);

	// @brief Para escribir la solucion en la forma que interesa
	void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
};

#endif
