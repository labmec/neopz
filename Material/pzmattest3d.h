#ifndef MATTEST3DHPP
#define MATTEST3DHPP


#include "pzmaterial.h"
#include "pzfmatrix.h"

//#include "pzmanvector.h"


class TPZMaterialTest3D : public TPZMaterial
{
private:

  /**
   * Fonte
   */
  TPZFMatrix fXf;

public :

  /**
   * Default empty constructor
   */
  TPZMaterialTest3D();

  /**
   * Full data constructor
   */
  TPZMaterialTest3D(int nummat);

  /**
   * Destructor
   */
  virtual ~TPZMaterialTest3D();

public:
  /**
   *Cedric : para testes no programa main 3dmaterial.c
   */
  static int geq3;

  /**
   * Set the flow
   */
  void SetMaterial(TPZFMatrix &xfin);

  /**
   * @see TPZMaterial
   */
  virtual int Dimension();

  /**
   * @see TPZMaterial
   */
  virtual int NStateVariables();

  /**
   * Print the object data structure
   */
  virtual void Print(std::ostream & out);

  /**
   * @see TPZMaterial
   */
  virtual std::string Name() { return "TPZMaterialTest3D"; }

  /**
   * @see TPZMaterial
   */
  virtual void Contribute( TPZMaterialData &data,REAL weight,
                           TPZFMatrix &ek,TPZFMatrix &ef );

  /**
   * @see TPZMaterial
   */
  virtual void ContributeBC( TPZMaterialData &data,REAL weight,
                             TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc );

  /**
   * @see TPZMaterial
   */
  virtual void Contribute( TPZMaterialData &data,REAL weight,
						   TPZFMatrix &ef )
  {
		TPZMaterial::Contribute(data,weight,ef);
  }

  /**
   * @see TPZMaterial
   */
  virtual void ContributeBC( TPZMaterialData &data,REAL weight,
							 TPZFMatrix &ef,TPZBndCond &bc )
  {
		TPZMaterial::ContributeBC(data,weight,ef,bc);
  }

  /**
   * @see TPZMaterial
   */
  virtual int VariableIndex(const std::string &name);

  /**
   * @see TPZMaterial
   */
  virtual int NSolutionVariables(int var);

  /**
   * @see TPZMaterial
   */
  virtual int NFluxes(){ return 3;}

protected:
  /**
   * @see TPZMaterial
   */
  virtual void Solution( TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,
						 int var,TPZVec<REAL> &Solout );
public:
      /**returns the solution associated with the var index based on
       * the finite element approximation*/
virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
    TPZMaterial::Solution(data,var,Solout);
}


  /**
   * @see TPZMaterial
   */
  virtual TPZAutoPointer<TPZMaterial> NewMaterial();


  /**
   * Compute the value of the flux function to be used by ZZ error estimator
   */
  virtual void Flux( TPZVec<REAL> &x, TPZVec<REAL> &Sol,
                     TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux );

  /**
   * @see TPZMaterial
   */
  virtual void Errors( TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx,
                       TPZFMatrix &axes, TPZVec<REAL> &flux,TPZVec<REAL> &u_exact,
                       TPZFMatrix &du_exact,TPZVec<REAL> &values );

  /**
   * @see TPZSaveable
   */
  virtual int ClassId() const;

  /**
   * @see TPZSaveable
   */
  virtual void Read(TPZStream &buf, void *context);

  /**
   * @see TPZSaveable
   */
  virtual void Write(TPZStream &buf, int withclassid);
};

#endif
