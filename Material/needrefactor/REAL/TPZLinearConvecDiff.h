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

	virtual TPZMaterial * NewMaterial() override {
		return new TPZLinearConvecDiff(*this);
	}

	virtual int Dimension() const  override { return 2;}

	virtual int NStateVariables() const override { return 1; }

	virtual void Print(std::ostream & out) override;

	virtual std::string Name()  override { return "TPZLinearConvecDiff"; }

	virtual void Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;

	virtual void ContributeBC(TPZMaterialData &data,REAL weight,
							              TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) override;

	virtual int VariableIndex(const std::string &name) override;

	virtual int NSolutionVariables(int var) override;

	virtual void Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout) override;

protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
public:

	virtual int NEvalErrors() override {return 3;}

        public:
virtual int ClassId() const override;

	virtual void Write(TPZStream &buf, int withclassid) const override {
    DebugStop();///implementar
  }

	virtual void Read(TPZStream &buf, void *context) override {
    DebugStop();///implementar
  }

};


#endif
