//---------------------------------------------------------------------------

/**
 * @file
 * @brief Contains the TPZMatLaplacian class.
 */
#ifndef TPZMatDualHybridPoissonH
#define TPZMatDualHybridPoissonH


#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief \f$ -Laplac(u) = fXf  \f$
 */
/**
 * \f$ -Laplac(u) = fXf  \f$
 */

class TPZMatDualHybridPoisson : public TPZMaterial {
    
    protected :
    
    /** @brief Forcing function value */
    REAL fXf;
    
    REAL fBetaZero;

    int mydim;
    
public:
    

    
    TPZMatDualHybridPoisson(int nummat, REAL f, REAL betaZero);
    TPZMatDualHybridPoisson(int matid);
    
    TPZMatDualHybridPoisson();
    
    TPZMatDualHybridPoisson(const TPZMatDualHybridPoisson &copy);
    
    void SetDimension(int dim){
        mydim = dim;
    }
    
    virtual ~TPZMatDualHybridPoisson() override;
    
    REAL Beta(int p, REAL size) const{
        return p*p*this->fBetaZero/size;
    }
    
    virtual TPZMaterial * NewMaterial() override {
        return new TPZMatDualHybridPoisson(*this);
    }
    
    virtual int Dimension() const override { return mydim;}
    
    virtual int NStateVariables() const override{
        return 1;
    }
    
    virtual void Print(std::ostream & out) override;
    
    virtual std::string Name()  override { return "TPZMatDualHybridPoisson"; }
    
    /**
     * @name Contribute methods (weak formulation)
     * @{
     */
    
    virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZMaterialData &data,REAL weight,
                              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight,
                                     TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;
    
    virtual int VariableIndex(const std::string &name) override;
    
    virtual int NSolutionVariables(int var) override;
    
public:
    
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
protected:
    void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
                TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
                TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
public:
    
    virtual int NEvalErrors()  override {return 3;}
    
    public:
virtual int ClassId() const override;
    
    virtual void Write(TPZStream &buf, int withclassid) const override {
        DebugStop();
    }
    
    virtual void Read(TPZStream &buf, void *context) override {
        DebugStop();
    }
    
};

#endif

