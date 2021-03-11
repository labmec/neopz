//
//  PZMatPoissonD3.h
//  PZ
//
//  Created by Douglas Castro on 5/23/14.
//
//

#ifndef __PZ__PZMatPoissonD3__
#define __PZ__PZMatPoissonD3__

#include <iostream>

#include "TPZMaterial.h"


/** Material para problema de Poisson 3D */
/**  Div.(-fK(x,y,z) Grad u) = fF    em \Omega
 *   u  = uD                       em \partial\Omega
 *   du/dn = g                     em \partial\Omega
 */

/**
 *   q = -fK(x,y,z) Grad u          em \Omega
 *   Div.(q) = fF                    em \Omega
 *   u  = uD                       em \partial\Omega
 *   du/dn = g                     em \partial\Omega
 */


class TPZMatPoissonD3 : public TPZMaterial {
    
protected:
	/** Material Id */
    int fMatId;
    
    /** Valor da funcao de carga */
    REAL fF; //fF
    
    /** Dimensao do dominio */
    int fDim;
    
    /** Coeficiente que multiplica o gradiente */
    REAL fK;
    
    /** @brief fluid viscosity*/
	REAL fvisc;
    
    /** @brief permeability tensor. Coeficient which multiplies the gradient operator*/
	TPZFMatrix<REAL> fTensorK;
    
    /** @brief inverse of the permeability tensor.*/
	TPZFMatrix<REAL> fInvK;
    
    /** @brief Pointer to forcing function, it is the Permeability and its inverse */
    TPZAutoPointer<TPZFunction<STATE> > fPermeabilityFunction;
	
public:
    
	TPZMatPoissonD3();
	
	TPZMatPoissonD3(int matid, int dim);
    
	virtual ~TPZMatPoissonD3();
    
    /** @brief copy constructor */
    TPZMatPoissonD3(const TPZMatPoissonD3 &copy);
    
    TPZMatPoissonD3 &operator=(const TPZMatPoissonD3 &copy);
	
	virtual std::string Name() { return "TPZMatPoissonD3"; }
	
	int Dimension() const {return fDim;}
    
    int MatId()
    {
        return fMatId;
    }
	    
    void SetPermeability(REAL perm) {
		fK = perm;
	}
	virtual int NStateVariables(void) const { return 1; }

    //Set the permeability tensor and inverser tensor
    void SetPermeabilityTensor(TPZFMatrix<REAL> K, TPZFMatrix<REAL> invK){
        
//        if(K.Rows() != fDim || K.Cols() != fDim) DebugStop();
//        if(K.Rows()!=invK.Rows() || K.Cols()!=invK.Cols()) DebugStop();
        
        fTensorK = K;
        fInvK = invK;
    }
	
    void SetViscosity(REAL visc) {
		fvisc = visc;
	}
    
	void GetPermeability(REAL &perm) {
		perm = fK;
	}
	
	void SetInternalFlux(REAL flux) {
		fF = flux;
	}

    
    void Print(std::ostream &out);
    
	/** @name Contribute methods
	 * @{
	 */
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since June 2, 2014
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
		DebugStop();
	}
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since June 2, 2014
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
		DebugStop();
	}
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since June 2, 2014
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition object
	 * @since June 2, 2014
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    void         ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @param datavec [in]
	 * @param dataleft [in]
     * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 2, 2014
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @param data [in]
	 * @param dataleft [in]
	 * @param dataright [in]
	 * @param weight [in]
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since June 2, 2014
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {DebugStop();}
    
	
	//virtual int VariableIndex(const std::string &name);
	
	//virtual int NSolutionVariables(int var);
	
	/**
     * @brief It return a solution to multiphysics simulation.
	 * @param datavec [in] Data material vector
     * @param var [in] number of solution variables. See  NSolutionVariables() method
     * @param Solout [out] is the solution vector
     */
	//virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    virtual void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec);
    
    int VariableIndex(const std::string &name);
    
    int NSolutionVariables(int var);
    
    // metodo para gerar vtk
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout);
    // metodo para computar erros Pressao
    void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout);
    // metodo para computar erros Hdiv
    void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);
    
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                                 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
                TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
    
    void ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);
    
};



#endif /* defined(__PZ__PZMatPoissonD3__) */
