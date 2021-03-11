/*
 *  TPZDarcyPMaterial.h
 *  PZ
 *
 *  Created by Pablo G S Carvalho on 08/09/16.
 *  Copyright 2016 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMatWithMem.h"

#include "pzfmatrix.h"
#include "pzbndcond.h"
#include "pzlog.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"


#ifndef TPZDARCYPMATERIAL
#define TPZDARCYPMATERIAL



class TPZDarcyPMaterial : public TPZMatWithMem<TPZFMatrix<STATE>, TPZMaterial >  {
    
private:
    
    /// dimension of the material
    int fDimension;
    
    /// Aproximation Space for velocity
    int fSpace;
    
    /// viscosidade
    STATE fViscosity;
    
    /** @brief Medium permeability. Coeficient which multiplies the gradient operator*/
    REAL fk;
    
    /// termo contrario a beta na sua formulacao (para ser conforme a literatura)
    STATE fTheta;
    
    
public:
    
    
    /**
     * Empty Constructor
     */
    TPZDarcyPMaterial();
    
    /** Creates a material object and inserts it in the vector of
     *  material pointers of the mesh.
     */
    TPZDarcyPMaterial(int matid, int dimension, int space, STATE viscosity, STATE permeability, STATE theta);
    
    
    /** Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
     */
    TPZDarcyPMaterial(const TPZDarcyPMaterial &mat);
    
    /**
     * Destructor
     */
    ~TPZDarcyPMaterial();
    
    /** Fill material data parameter with necessary requirements for the
     * Contribute method. Here, in base class, all requirements are considered
     * as necessary. Each derived class may optimize performance by selecting
     * only the necessary data.
     * @since April 10, 2007
     */
    
    void FillDataRequirements(TPZMaterialData &data) override;
    
    void FillDataRequirements(TPZVec<TPZMaterialData> &datavec) override;
    
    virtual void FillDataRequirementsInterface(TPZMaterialData &data) override;
    
    void FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData> &datavec) override;
    
    void SetPermeability(REAL perm) {
        fk = perm;
    }
    
    /** returns the name of the material */
    std::string Name()  override {
        return "TPZDarcyPMaterial";
    }
    
    /** returns the integrable dimension of the material */
    int Dimension() const  override {return 2;}
    
    /** returns the number of state variables associated with the material */
    virtual int NStateVariables() const  override {return 4;} // for hdiv are 3, plus pressure, so 3 + 1 = 4 itapopo
    
    /** print out the data associated with the material */
    void Print(std::ostream &out = std::cout) override;
    
    /** returns the variable index associated with the name */
    int VariableIndex(const std::string &name) override;
    
    /** returns the number of variables associated with the variable
     indexed by var.  var is obtained by calling VariableIndex */
    int NSolutionVariables(int var) override;
    
    /** Computes the divergence over the parametric space */
    void ComputeDivergenceOnDeformed(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi);
    
    void ComputeDivergenceOnMaster(TPZVec<TPZMaterialData> &datavec, TPZFMatrix<STATE> &DivergenceofPhi, STATE &DivergenceofU);
    
    /** returns the solution associated with the var index based on
     * the finite element approximation */
    void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
    
    /** index of velocity */
    int VIndex(){ return 0; }
    
    /** index of pressure */
    int PIndex(){ return 1; }
    
    /** inner product of two tensors. See Gurtin (2003), p. 5. */
    template <typename TVar>
    TVar Inner(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
    /** inner product of two vectors. See Gurtin (2003), p. 5. */
    template <typename TVar>
    TVar InnerVec(TPZFMatrix<TVar> &S, TPZFMatrix<TVar> &T);
    
    /** trace of the tensor GradU = Div(U)*/
    STATE Tr(TPZFMatrix<REAL> &GradU );
    
    /** transpose of the tensor GradU = Div(U)*/
    STATE Transpose(TPZFMatrix<REAL> &GradU );
    
    /** Fill the vector of gradient for each phi */
    void FillGradPhi(TPZMaterialData &dataV, TPZVec< TPZFMatrix<REAL> > &GradPhi);
    
    /// transform a H1 data structure to a vector data structure
    void FillVecShapeIndex(TPZMaterialData &data);
    
    /** @brief Not used contribute methods */
    virtual void Contribute(TPZMaterialData &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {DebugStop();}
    
    // Contribute Methods being used - Multiphysics
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override
    {
        DebugStop();
    }
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        DebugStop();
    }
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override {
        DebugStop();
    }
    
    
    // Contribute Methods being used
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @since April 16, 2007
     */
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one BC interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ek[out] is the stiffness matrix
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
        DebugStop();
    }
    
    /**
     * It computes a contribution to the stiffness matrix and load vector at one internal interface integration point.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     * @since April 16, 2007
     */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &datavecleft, TPZVec<TPZMaterialData> &datavecright, REAL weight,TPZFMatrix<STATE> &ef) override {
        DebugStop();
    }
    
    /**
     * Save the element data to a stream
     */
    void Write(TPZStream &buf, int withclassid) const override;
    
    /**
     * Read the element data from a stream
     */
    void Read(TPZStream &buf, void *context) override;
    
    
    virtual int NEvalErrors()  override {return 6;}
    
    /**
     * It computes errors contribution in differents spaces.
     * @param data[in] stores all input data
     * @param weight[in] is the weight of the integration rule
     * @param ef[out] is the load vector
     * @param bc[in] is the boundary condition material
     */
    virtual void Errors(TPZVec<TPZMaterialData> &data, TPZVec<STATE> &u_exact, TPZFMatrix<STATE> &du_exact, TPZVec<REAL> &errors) override;
    
};

#endif
