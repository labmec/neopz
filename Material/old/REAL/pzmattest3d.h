/**
 * @file
 * @brief Contains the TPZMaterialTest3D class. Three-dimensional test.
 */

#ifndef MATTEST3DHPP
#define MATTEST3DHPP

#include "TPZMaterial.h"
#include "pzfmatrix.h"

/**
 * @ingroup material
 * @brief Implements a three dimensional linear material for test.
 */
class TPZMaterialTest3D : public TPZMaterial
{
private:
	
	/** @brief Source */
	TPZFMatrix<STATE> fXf;
	
	public :	
	/** @brief Default empty constructor */
	TPZMaterialTest3D();
	
	/** @brief Full data constructor */
	TPZMaterialTest3D(int nummat);
	
	/** @brief Destructor */
	virtual ~TPZMaterialTest3D();
	
	/** @brief Cedric : para testes no programa main 3dmaterial.c */
	static int geq3;
	
	/** @brief Set the flow */
	void SetMaterial(TPZFMatrix<STATE> &xfin);
	
	virtual int Dimension() const override;
	
	virtual int NStateVariables() const override;
	
	/** @brief Prints the object data structure */
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name()  override { return "TPZMaterialTest3D"; }
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef ) override;
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc ) override;
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ef ) override
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc ) override
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution( TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
						  int var,TPZVec<STATE> &Solout ) override;
public:
	/**
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation
	 */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	virtual TPZMaterial * NewMaterial() override;
	
protected:
	void Errors( TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx,
						TPZFMatrix<REAL> &axes, TPZVec<STATE> &u_exact,
						 TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values ) override;
	
	public:
virtual int ClassId() const override;

	
	virtual void Read(TPZStream &buf, void *context) override;
	
	virtual void Write(TPZStream &buf, int withclassid) const override;
};

#endif
