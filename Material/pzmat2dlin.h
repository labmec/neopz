/**
 * \file
 * @brief Contains the TPZMat2dLin class which implements a bi-dimensional linear problem.
 */
#ifndef MAT2DLINHPP
#define MAT2DLINHPP

#include <iostream>
#include "pzmaterial.h"
#include "pzfmatrix.h"
//#include "pzreal.h"
struct TPZElementMatrix;
class TPZBndCond;
template<class T>
class TPZVec;

/**
 * @ingroup material
 * @brief Implements a bi-dimensional linear problem.
 */
class TPZMat2dLin : public TPZMaterial{
	
	TPZFMatrix    fKxx, fKxy, fKyx, fKyy, fKx0, fK0x, fKy0, fK0y, fK00, fXf;
	public :
	
    TPZMat2dLin(int num = 0) : TPZMaterial(num), fKxx(), fKxy(),
	fKyx() , fKyy(), fKx0(), fK0x(), fKy0(), fK0y(), fK00(), fXf() {
    }
	
	TPZMat2dLin(TPZMat2dLin &copy) : TPZMaterial(copy),
	fKxx(copy.fKxx), fKxy(copy.fKxy), fKyx(copy.fKyx), fKyy(copy.fKyy),
	fKx0(copy.fKx0), fK0x(copy.fK0x), fKy0(copy.fKy0),
	fK0y(copy.fK0y), fK00(copy.fK00), fXf(copy.fXf){}
	
	virtual int NStateVariables() { return fKxx.Rows(); }
	
	//  int NFluxes() { return NStateVariables(); }
	
	int Dimension() { return 2; }
	
	void Print(std::ostream & out = std::cout);
	
	void SetMaterial(TPZFMatrix &xkin,TPZFMatrix &xcin,TPZFMatrix &xfin){
		int r = xkin.Rows();
		fKxx = xkin;
		fKyy = xkin;
		fK00 = xcin;
		fXf = xfin;
		fKxy.Redim(r,r);
		fKyx.Redim(r,r);
		fKx0.Redim(r,r);
		fK0x.Redim(r,r);
		fKy0.Redim(r,r);
		fK0y.Redim(r,r);
	}
	
	void ConvectionDiffusion(REAL angle,REAL diff);
	
	TPZFMatrix &Xk() {return fKxx;}
	TPZFMatrix &Ck() {return fK00;}
	TPZFMatrix &Xf() {return fXf;}
	
	virtual std::string Name() { return "TPZMat2dLin"; }
	
	//  virtual TPZBndCond *CreateBC(int num,int typ,TPZFMatrix &val1,TPZFMatrix &val2);
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix &ek,TPZFMatrix &ef);
	
	virtual void Contribute(TPZMaterialData &data,REAL weight,
							TPZFMatrix &ef)
	{
		TPZMaterial::Contribute(data,weight,ef);
	}
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc);
	
	virtual void ContributeBC(TPZMaterialData &data, REAL weight,
							  TPZFMatrix &ef,TPZBndCond &bc)
	{
		TPZMaterial::ContributeBC(data,weight,ef,bc);
	}
	
	virtual int NFluxes();
	
	virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &fl);
	
	void Errors(TPZVec<REAL> &x,TPZVec<REAL> &u,TPZFMatrix &dudx,TPZFMatrix &axes,TPZVec<REAL> &flux,
				TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values);
	
	virtual int VariableIndex(const std::string &name);
	
	virtual int NSolutionVariables(int index);
	
protected:
	void Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes, int var,TPZVec<REAL> &Solout);
public:
	virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
	{
		TPZMaterial::Solution(data,var,Solout);
	}
	/** @brief Creates a copy of the material object */
	virtual TPZAutoPointer<TPZMaterial> NewMaterial();
	
	TPZBndCond *OutflowFlux(TPZAutoPointer<TPZMaterial> &reference, int bc);
	
	/** @brief returns the unique identifier for reading/writing objects to streams */
	virtual int ClassId() const;
	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid);
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);

};

#endif


