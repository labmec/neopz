/**
 * \file
 * @brief Contains the TPZBndCond class which implements a boundary condition for TPZMaterial objects.
 */

#ifndef BNDCONDHPP
#define BNDCONDHPP

#include <iostream>
#include <map>

#include "pzreal.h"
#include "pzdiscgal.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzfmatrix.h"
#include "tpzautopointer.h"

template <class T, int N>
class TPZManVector;

/**
 * @ingroup material
 * @brief This class defines the boundary condition for TPZMaterial objects
 */
/**
 * This class redirects the call to Contribute to calls ContributeBC of the referring material object
 */
class TPZBndCond : public TPZDiscontinuousGalerkin {
	
	friend class TPZMaterial;
    
protected:
    
    class TPZ_BCDefine : public TPZSavable {
       
        public :
        /** @brief second value of boundary condition */
        TPZFNMatrix<6,STATE> fBCVal2;

        /** @brief Pointer to forcing function, it is the right member at differential equation */
        TPZAutoPointer<TPZFunction<STATE> > fForcingFunction;
        
        /** @brief Pointer to exact solution function, needed to calculate exact error */
        TPZAutoPointer<TPZFunction<STATE> > fForcingFunctionExact;
        
        /** @brief Pointer to time dependent forcing function, it is the right member at differential equation */
        TPZAutoPointer<TPZFunction<STATE> > fTimeDependentForcingFunction;
        
        /** @brief Pointer to time dependent exact solution function, needed to calculate exact error */
        TPZAutoPointer<TPZFunction<STATE> > fTimedependentFunctionExact;
        
        /** @brief Pointer to bc forcing function, it is a variable boundary condition at differential equation */
        TPZAutoPointer<TPZFunction<STATE> > fBCForcingFunction;
        
        /** @brief Pointer to time dependent bc forcing function, it is a variable boundary condition at differential equation */
        TPZAutoPointer<TPZFunction<STATE> > fTimedependentBCForcingFunction;            

        TPZ_BCDefine() : fBCVal2(), fForcingFunction(NULL), fForcingFunctionExact(NULL),fTimeDependentForcingFunction(NULL),
                        fTimedependentFunctionExact(NULL), fBCForcingFunction(NULL),fTimedependentBCForcingFunction(NULL)
        {
            
        }
        TPZ_BCDefine(TPZFMatrix<STATE> Val2) : fBCVal2(), fForcingFunction(NULL), fForcingFunctionExact(NULL),fTimeDependentForcingFunction(NULL),
                fTimedependentFunctionExact(NULL), fBCForcingFunction(NULL),fTimedependentBCForcingFunction(NULL)
        {
            
        }
        TPZ_BCDefine(const TPZ_BCDefine &cp) : fBCVal2(cp.fBCVal2), fForcingFunction(cp.fForcingFunction),fForcingFunctionExact(cp.fForcingFunctionExact),
                fTimeDependentForcingFunction(cp.fTimeDependentForcingFunction), fTimedependentFunctionExact(cp.fTimedependentFunctionExact),
                fBCForcingFunction(cp.fBCForcingFunction),fTimedependentBCForcingFunction(cp.fTimedependentBCForcingFunction)
        {
            
        }
        TPZ_BCDefine &operator=(const TPZ_BCDefine &cp)
        {
            fBCVal2 = cp.fBCVal2;
            fForcingFunction = cp.fForcingFunction;
            fForcingFunctionExact = cp.fForcingFunctionExact;
            fTimeDependentForcingFunction = cp.fTimeDependentForcingFunction;
            fTimedependentFunctionExact = cp.fTimedependentFunctionExact;
            fBCForcingFunction = cp.fBCForcingFunction;
            fTimedependentBCForcingFunction = cp.fTimedependentBCForcingFunction;
            return *this;
        }
        ~TPZ_BCDefine()
        {
            
        }
        int ClassId() const;
        void Read(TPZStream& buf, void* context);
        void Write(TPZStream& buf, int withclassid) const;

    };
    
    TPZVec<TPZ_BCDefine> fBCs;
    
	/** @brief boundary condition type */
	int 		fType;
	/** @brief first value of boundary condition */
 	TPZFNMatrix<9,STATE> fBCVal1;
	/** @brief second value of boundary condition */
	TPZFNMatrix<3,STATE> fBCVal2;
	/** @brief pointer to material which created bc */
	TPZMaterial * fMaterial;
	
	/** @brief Function to allow fBCVal1 to be variable */
	void (*fValFunction)(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType);
	
	public :
	/** @brief Copy constructor */
    TPZBndCond(TPZBndCond & bc) : TPZRegisterClassId(&TPZBndCond::ClassId),
    TPZDiscontinuousGalerkin(bc), fBCs(bc.fBCs), fType(-1), fBCVal1(bc.fBCVal1),
    fBCVal2(bc.fBCVal2), fMaterial(0), fValFunction(NULL){
		fMaterial = bc.fMaterial;
		fType = bc.fType;
        fForcingFunction = bc.fForcingFunction;
        fForcingFunctionExact = bc.fForcingFunctionExact;
        fTimeDependentForcingFunction = bc.fTimeDependentForcingFunction;
        fTimedependentFunctionExact = bc.fTimedependentFunctionExact;
        fBCForcingFunction = bc.fBCForcingFunction;
        fTimedependentBCForcingFunction = bc.fTimedependentBCForcingFunction;
	}
	/** @brief Default constructor */
	TPZBndCond() : TPZRegisterClassId(&TPZBndCond::ClassId),
    TPZDiscontinuousGalerkin(), fBCs(), fType(-1), fBCVal1(),
    fBCVal2(), fMaterial(0), fValFunction(NULL){
	}
	/** @brief Default constructor */
	TPZBndCond(int matid) : TPZRegisterClassId(&TPZBndCond::ClassId),
    TPZDiscontinuousGalerkin(matid), fBCs(0), fType(-1), fBCVal1(),
    fBCVal2(), fMaterial(0), fValFunction(NULL){
	}
	/** @brief Default destructor */
    ~TPZBndCond(){}
	
	TPZBndCond(TPZMaterial * material,int id,int type,TPZFMatrix<STATE> &val1,TPZFMatrix<STATE> &val2) :
    TPZRegisterClassId(&TPZBndCond::ClassId), TPZDiscontinuousGalerkin(id), fBCs(), fBCVal1(val1), fBCVal2(val2), fValFunction(NULL) {
		//creates a new material
		if(!material)
		{
			std::cout << __PRETTY_FUNCTION__ << " Creating boundary condition with NULL material" << std::endl;
            DebugStop();
		}
		fMaterial = material;
		fType = type;
		
	}
	
	TPZBndCond(TPZBndCond &copy, TPZMaterial * ref) : TPZRegisterClassId(&TPZBndCond::ClassId), 
    TPZDiscontinuousGalerkin(copy), fBCs(copy.fBCs), fType(copy.fType),
	fBCVal1(copy.fBCVal1), fBCVal2(copy.fBCVal2), fMaterial(ref), fValFunction(copy.fValFunction) {
    
        fForcingFunction = copy.fForcingFunction;
        fForcingFunctionExact = copy.fForcingFunctionExact;
        fTimeDependentForcingFunction = copy.fTimeDependentForcingFunction;
        fTimedependentFunctionExact = copy.fTimedependentFunctionExact;
        fBCForcingFunction = copy.fBCForcingFunction;
        fTimedependentBCForcingFunction = copy.fTimedependentBCForcingFunction;
        
    }
	
	
	void SetValFunction(void (*fp)(TPZVec<REAL> &loc, TPZFMatrix<STATE> &Val1, TPZVec<STATE> &Val2, int &BCType)){
		fValFunction = fp;
	}
	
	void SetForcingFunction(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetForcingFunction(func);
        }
        else {
            fBCs[loadcase].fForcingFunction = func;
        }
	}

	void SetForcingFunctionExact(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetForcingFunctionExact(func);
        }
        else {
            fBCs[loadcase].fForcingFunctionExact = func;
        }
    }

    void SetTimeDependentForcingFunction(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetTimeDependentForcingFunction(func);
        }
        else {
            fBCs[loadcase].fTimeDependentForcingFunction = func;
        }
    }

    void SetTimeDependentFunctionExact(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetTimeDependentFunctionExact(func);
        }
        else {
            fBCs[loadcase].fTimedependentFunctionExact = func;
        }
    }     

    void SetBCForcingFunction(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetBCForcingFunction(func);
        }
        else {
            fBCs[loadcase].fBCForcingFunction = func;
        }
    }

    void SetTimedependentBCForcingFunction(int loadcase, TPZAutoPointer<TPZFunction<STATE> > func)
    {
        if (loadcase == 0) {
            TPZMaterial::SetTimedependentBCForcingFunction(func);
        }
        else {
            fBCs[loadcase].fTimedependentBCForcingFunction = func;
        }
    }
    
    TPZAutoPointer<TPZFunction<STATE> > GetTimedependentBCForcingFunction(int loadcase)
    {
        if (loadcase == 0) {
            return this->fTimedependentBCForcingFunction;
        }
        else {
            return fBCs[loadcase].fTimedependentBCForcingFunction;
        }
    }
	
	void SetMaterial(TPZMaterial * mat) { fMaterial = mat;}
	
	/** @brief Returns the integrable dimension of the material*/
	int Dimension() const { return fMaterial->Dimension(); }
	
	virtual int NFluxes(){ return fMaterial->NFluxes(); }
	
	int NStateVariables() { return fMaterial->NStateVariables(); }
	
	/** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
	virtual int NEvalErrors() {return fMaterial->NEvalErrors();}
	
	int Type() { return fType; }
	
	void SetType(int type){ this->fType = type; }
	
	TPZFMatrix<STATE> &Val1() { return fBCVal1; }
	
	TPZFMatrix<STATE> &Val2(int loadcase = 0) 
    {
        if (loadcase == 0 || loadcase > fBCs.size() ) {
            return fBCVal2;
        }
        return fBCs[loadcase-1].fBCVal2;
    }
	
	TPZMaterial * Material() { return fMaterial; }
    
    void SetLoadCases(TPZVec<TPZFMatrix<STATE> > &val2vec)
    {
        fBCs.Resize(val2vec.size()-1);
        fBCVal2 = val2vec[0];
        for (int i=0; i<val2vec.size()-1; i++) {
            fBCs[i] = val2vec[i+1];
        }
    }
    
    virtual int MinimumNumberofLoadCases()
    {
        return 1+fBCs.size();
    }
	
	/** @brief Computes the value of the flux function to be used by ZZ error estimator */
	void Flux(TPZVec<REAL> &x, TPZVec<STATE> &Sol, TPZFMatrix<STATE> &DSol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux){
		flux.Fill(0.);
	}
	
	void Print(std::ostream & out = std::cout) {
		out << " Boundary condition number = " << Id() << "\n";
		out << " boundary condition type = " << fType << "\n";
		out << " val1 = \n"; fBCVal1.Print("fBCVal1",out);
		out << " val2 = \n"; fBCVal2.Print("fBCVal2",out);
		if (HasForcingFunction()) out << " has forcing function\n";
		else out << "has no forcing function\n";
	}
	
	void UpdateBCValues(TPZManVector<REAL,3> & x);
	
	/*
	 @brief this method will be used only in case of multiphysics simulation
	 */
	void UpdateBCValues(TPZVec<TPZMaterialData> &datavec);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one integration point 
	 * for the multiphysics element
	 * @param datavec [in] stores all input data vector 
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ek [out] is the stiffness matrix
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to the residual vector at one integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the residual vector
	 * @since April 16, 2007
	 */
	virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
	 * @param data [in] stores all input data
	 * @param weight [in] is the weight of the integration rule
	 * @param ef [out] is the load vector
	 * @param bc [in] is the boundary condition material
	 * @since April 16, 2007
	 */
	virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
    
    /**
	 * @brief It computes a contribution to stiffness matrix and load vector at one integration point
	 * @since March 5, 2013
	 */
    virtual void ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to residual vector at one integration point
	 * @since April 16, 2007
	 */
	virtual void ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ef);
	
	/**
	 * @brief It computes a contribution to stiffness matrix and load vector at one BC integration point
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	/**
	 * @brief It computes a contribution to residual vector at one BC integration point
	 * @since April 16, 2007
	 */
	virtual void ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc);
	
	virtual void Errors(TPZVec<REAL> &x,TPZVec<STATE> &sol,TPZFMatrix<STATE> &dsol, TPZFMatrix<REAL> &axes, TPZVec<STATE> &flux,
				TPZVec<STATE> &uexact,TPZFMatrix<STATE> &duexact,TPZVec<REAL> &val){
		val.Fill(0.);
	}
	
	virtual void Clone(std::map<int, TPZMaterial * > &matvec);
	
	/** 
	 * @brief Compute interface jumps
	 * @since Feb 14, 2006
	 */
	/**
	 * \f$ values[1] = (solleft - solright)^2 \f$ \n
	 * \f$ values[2] = (dsolleft - dsolright)^2 \f$ \n
	 * \f$ values[0] = values[1] + values[2] \f$
	 */
	virtual void InterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZSolVec &rightu,TPZSolVec &jump);
	
	/** @brief Returns the unique identifier for reading/writing objects to streams */
	public:
virtual int ClassId() const;

	/** @brief Saves the element data to a stream */
	virtual void Write(TPZStream &buf, int withclassid) const;
	
	/** @brief Reads the element data from a stream */
	virtual void Read(TPZStream &buf, void *context);
	
	void ContributeInterfaceErrors(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
								   REAL weight,
								   TPZVec<STATE> &nkL,
								   TPZVec<STATE> &nkR,
								   int &errorid);
	
	void  ContributeErrors(TPZMaterialData &data,// aqui
						   REAL weight,
						   TPZVec<REAL> &nk,
						   int &errorid){
		//nothing to be done here
	};
	
	
	virtual void ContributeInterfaceBCErrors(TPZMaterialData &data, TPZMaterialData &dataleft,
											 REAL weight,
											 TPZVec<STATE> &nk,
											 TPZBndCond &bc,
											 int &errorid
											 ){
		//nothing to be done here
	}
	
	/** @brief Calls the aggregate material correspondent function */
	virtual void FillDataRequirements(TPZMaterialData &data);
    
    virtual void FillDataRequirements(TPZVec<TPZMaterialData> &datavec);
	
};

#endif
