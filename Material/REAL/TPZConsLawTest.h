/**
 * @file
 * @brief Contains the TPZConsLawTest class for test. Material as conservation law
 */
#ifndef CONSLAWTESTHPP
#define CONSLAWTESTHPP

#include <iostream>
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzconslaw.h"


/**
 * @brief Only to test a material as conservation law. It was used for testing purposes
 * @ingroup material
 */
class TPZConsLawTest  : public TPZConservationLaw {
	
	TPZFMatrix<STATE> fXf;//fonte
	TPZVec<STATE> fB;
	int fArtificialDiffusion;

	/// Integer for integration degree of the initial solution
	int fTest;
	
	STATE fDelta;

	public :
	
	TPZConsLawTest(int nummat, TPZVec<STATE> B,int artdiff,STATE delta_t,int dim,STATE delta,int test=0);
	
	virtual ~TPZConsLawTest();
	
	void SetMaterial(TPZFMatrix<STATE> &xfin) {
		fXf = xfin;
	}
	
	STATE DeltaOtimo();
	
	STATE CFL(int degree);
	
	STATE B(int i,TPZVec<REAL> &x);
		
	REAL Delta();
	
	STATE T(int jn,TPZVec<REAL> &x);
	
	virtual int NStateVariables() const override;
	
	virtual void Print(std::ostream & out) override;
	
	virtual std::string Name()  override { return "TPZConsLawTest"; }
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) override;
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix<STATE> &ef) override {
		TPZConservationLaw::Contribute(data,weight,ef);
	}
	
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                     REAL weight,
                                     TPZFMatrix<STATE> &ek,
                                     TPZFMatrix<STATE> &ef) override;
	
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ek,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override;
	virtual void ContributeBC(TPZMaterialData &data,
							  REAL weight,
							  TPZFMatrix<STATE> &ef,
							  TPZBndCond &bc) override
	{
		TPZConservationLaw::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int VariableIndex(const std::string &name) override;
	
	virtual int NSolutionVariables(int var) override;
	
protected:
	virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout) override;
public:
	virtual void Solution(TPZMaterialData &data,int var,TPZVec<STATE> &Solout) override
	{
        int numbersol = data.sol.size();
        if (numbersol != 1) {
            DebugStop();
        }

		Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
	}
protected:
	void Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
				TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
				TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) override;
	
	void ComputeSolRight(TPZVec<STATE> &solr,TPZVec<STATE> &soll,TPZVec<REAL> &normal,TPZBndCond *bcright);
	void ComputeSolLeft(TPZVec<STATE> &solr,TPZVec<STATE> &soll,TPZVec<REAL> &normal,TPZBndCond *bcleft);

};

#endif
