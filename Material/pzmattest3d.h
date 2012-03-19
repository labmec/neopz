/**
 * \file
 * @brief Contains the TPZMaterialTest3D class. Three-dimensional test.
 */
#ifndef MATTEST3DHPP
#define MATTEST3DHPP


#include "pzmaterial.h"
#include "pzfmatrix.h"


/**
 * @ingroup material
 * @brief Implements a three dimensional linear material for test.
 */
class TPZMaterialTest3D : public TPZMaterial
{
private:
	
	/** @brief Source */
	TPZFMatrix<REAL> fXf;
	
	public :
	
	/** @brief Default empty constructor */
	TPZMaterialTest3D();
	
	/** @brief Full data constructor */
	TPZMaterialTest3D(int nummat);
	
	/** @brief Destructor */
	virtual ~TPZMaterialTest3D();
	
public:
	/** @brief Cedric : para testes no programa main 3dmaterial.c */
	static int geq3;
	
	/** @brief Set the flow */
	void SetMaterial(TPZFMatrix<REAL> &xfin);
	
	virtual int Dimension();
	
	virtual int NStateVariables();
	
	/** @brief Prints the object data structure */
	virtual void Print(std::ostream & out);
	
	virtual std::string Name() { return "TPZMaterialTest3D"; }
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef );
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef,TPZBndCond &bc );
	
	virtual void Contribute( TPZMaterialData &data,REAL weight,
							TPZFMatrix<REAL> &ef )
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							  TPZFMatrix<REAL> &ef,TPZBndCond &bc )
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int var);
	
	virtual int NFluxes(){ return 3;}
	
protected:
	virtual void Solution( TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,
						  int var,TPZVec<REAL> &Solout );
public:
	/**
	 * @brief Returns the solution associated with the var index based on
	 * the finite element approximation
	 */
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	virtual void Flux( TPZVec<REAL> &x, TPZVec<REAL> &Sol,
					  TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux );
	
	virtual void Errors( TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix<REAL> &dudx,
						TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,TPZVec<REAL> &u_exact,
						TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values );
	
	virtual int ClassId() const;
	
	virtual void Read(TPZStream &buf, void *context);
	
	virtual void Write(TPZStream &buf, int withclassid);
};

#endif
