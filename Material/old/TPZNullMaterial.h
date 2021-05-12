 /**
 * @file pzmaterial.h
 * @brief Header file for abstract class TPZMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef TPZNULLMATERIALH
#define TPZNULLMATERIALH

#include "TPZMaterial.h"

class  TPZNullMaterial : public TPZMaterial
{
    
public:
    TPZNullMaterial(int num) : TPZRegisterClassId(&TPZNullMaterial::ClassId), TPZMaterial(num) {
        fDim = 1;
        fNState = 1;
    }
    
    TPZNullMaterial(int num, int dimension, int nstate) : TPZRegisterClassId(&TPZNullMaterial::ClassId), TPZMaterial(num) {
        fDim = dimension;
        fNState = nstate;
    }
    
    TPZNullMaterial(const TPZNullMaterial &copy) : TPZRegisterClassId(&TPZNullMaterial::ClassId),
    TPZMaterial(copy)
    {
        fDim = copy.fDim;
        fNState = copy.fNState;
    }
    
    TPZNullMaterial &operator=(const TPZNullMaterial &copy)
    {
        fDim = copy.fDim;
        fNState = copy.fNState;
        TPZMaterial::operator=(copy);
        return *this;
    }
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** Upon return vectorindex contains the index of the material object within the vector */

    
    /** @brief Default constructor */
    TPZNullMaterial();
    
    /** @brief Default destructor */
    virtual ~TPZNullMaterial();
    
    /** @brief To create another material of the same type*/
    virtual TPZMaterial * NewMaterial() override
    {
        return new TPZNullMaterial(*this);
    }


    /** @brief Returns the name of the material */
    virtual std::string Name() override { return "TPZNullMaterial"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() const override
    {
        return fDim;
    }
    
    void SetDimension(int dim){
        fDim = dim;
    }
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() const override
    {
        return fNState;
    }
    
    void SetNStateVariables(int nstate)
    {
        fNState = nstate;
    }
    
	
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout)override;
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name)override;
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var) override;
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout) override;
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation around one interface element */
    virtual void Solution(TPZMaterialData &data, std::map<int, TPZMaterialData> &dataleftvec, std::map<int, TPZMaterialData> &datarightvec, int var, TPZVec<STATE> &Solout, TPZCompEl * left, TPZCompEl * ritgh) override;
    
protected:
    
    /** @deprecated Deprecated interface for Solution method which must use material data. */
    virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
    
    /** @brief Problem dimension */
    int fDim;
    
    /// Number of state variables
    int fNState = 1;
    
public:
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
    
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) override;
	
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)override
    {
        TPZFMatrix<STATE> ek(ef.Rows(),ef.Rows(),0.);
        Contribute(datavec, weight, ek, ef);
    }
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)override;
  

    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)override
    {
        
    }
    
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef) override;
    
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc) override;
	
    
    
public:

    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const override;

    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid) const override;
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context) override;
	
};


#endif
