//---------------------------------------------------------------------------

#ifndef TPZLinearConvecDiffH
#define TPZLinearConvecDiffH
#include <iostream>
#include "pzfmatrix.h"
#include "TPZMaterial.h"


/**
 * @ingroup material
 * @brief Convecção-difusão linear 2D
 */
/**
 * \f$ -fK Laplac(u) + div(fConv*u) = fXf  \f$
 */
class TPZLinearConvecDiff : public TPZMaterial {

	protected :

	/** @brief Forcing function value */
	STATE fXf;

	/** @brief Coeficient which multiplies the Laplacian operator. */
	STATE fK;

	/** @brief Convection vector */
	REAL fConvDir[2];

  /** @brief Multiplication value for the streamline diffusion term */
	STATE fSD;

public:

	TPZLinearConvecDiff(int nummat, REAL k, const TPZVec<REAL> &conv, REAL f, REAL SD);

  TPZLinearConvecDiff(int matid);

	TPZLinearConvecDiff();
	
	TPZLinearConvecDiff(const TPZLinearConvecDiff &c);

	virtual ~TPZLinearConvecDiff();

	virtual TPZMaterial * NewMaterial(){
		return new TPZLinearConvecDiff(*this);
	}

	virtual int Dimension() const { return 2;}

	int NStateVariables(){ return 1; }

	virtual void Print(std::ostream & out);

	virtual std::string Name() { return "TPZLinearConvecDiff"; }

	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef);

	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);

	virtual int VariableIndex(const std::string &name);

	virtual int NSolutionVariables(int var);

	virtual int NFluxes(){ return 2;}

	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout);

	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values);

	virtual int NEvalErrors() {return 3;}

        public:
virtual int ClassId() const;

	virtual void Write(TPZStream &buf, int withclassid) const{
    DebugStop();///implementar
  }

	virtual void Read(TPZStream &buf, void *context){
    DebugStop();///implementar
  }

};


#endif
