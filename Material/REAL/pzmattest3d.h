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
	
	virtual int Dimension() const;
	
	virtual int NStateVariables();
	
	/** @brief Prints the object data structure */
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMaterialTest3D"; }
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef );
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc );
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ef )
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<STATE> &ef,TPZBndCond &bc )
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
	
protected:
	virtual void Solution( TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,
						  int var,TPZVec<STATE> &Solout );
public:
	/**
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation
	 */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	virtual TPZMaterial * NewMaterial();
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux( TPZVec<REAL> &x, TPZVec<STATE> &Sol,
					  TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux );
	
	virtual void Errors( TPZVec<REAL> &x,TPZVec<STATE> &u,TPZFMatrix<STATE> &dudx,
						TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,TPZVec<STATE> &u_exact,
						TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values );
	
	public:
virtual int ClassId() const;

	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid) const;
};

#endif
