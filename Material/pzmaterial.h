/** @file pzmaterial.h
 *
 * @brief Header file for abstract class TPZMaterial.\n
 * It implements the weak statement of the differential equation within the PZ environment.
 */

#ifndef PZMATERIALHPP
#define PZMATERIALHPP

#include "pzreal.h"
#include "pzvec.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzadmchunk.h"
#include "tpzautopointer.h"
#include "pzsave.h"
#include "pzmaterialdata.h"
#include "pzfunction.h"

#include <iostream>
#include <string>
//#ifdef _AUTODIFF
//#include "fadType.h"
//#endif


class TPZBndCond;
class TPZMaterial;
class TPZMaterialData;
class TPZIntPoints;

/**
 * @ingroup material
 * @brief This abstract class defines the behaviour which each derived class needs to implement
 */
/**
 * Classes derived from the TPZMaterial class implement the weak statement of the differential equation
 * within the PZ environment \n
 * It is noteworthy to observe that this definition does not depend on the definition of the interpolation space \n
 * TPZMaterial objects also need to implement the interface for post processing the results
 */
class  TPZMaterial : public TPZSaveable
{
private:
    int fId;
    
protected:
    
    TPZAutoPointer<TPZFunction> fForcingFunction;
    //void (*fForcingFunction)(TPZVec<REAL> &loc,TPZVec<REAL> &result);
	// void (*fForcingFunctionExact)(TPZVec<REAL> &loc,TPZVec<REAL> &pressure,TPZVec<REAL> &flux);
	void (*fForcingFunctionExact) (TPZVec<REAL> &loc,
								   TPZVec<STATE> &pressure,TPZFMatrix<STATE> &flux);
    /** @brief Defines whether the equation context is linear solver or non linear */
    /**
     * True means linear (default)
     * @since 08 oct 2010
     */
    bool fLinearContext;
    
public:
    
    static REAL gBigNumber;
    
    /** @brief Creates a material object and inserts it in the vector of material pointers of the mesh. */
	/** 
	 * Upon return vectorindex contains the index of the material object within the vector
     */
    TPZMaterial(int id);
    
    /** @brief Default constructor */
    TPZMaterial();
    
    /** @brief Creates a material object based on the referred object and
     *  inserts it in the vector of material pointers of the mesh.
	 */
	/**  Upon return vectorindex contains the index of the material object within the vector */
    TPZMaterial(const TPZMaterial &mat);
    
    virtual ~TPZMaterial();
    
    /** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * @since April 10, 2007
	 */
	/** 
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
    virtual void FillDataRequirements(TPZMaterialData &data);
	
	/** 
	 * @brief Fill material data parameter with necessary requirements for the
	 * Contribute method. Here, in base class, all requirements are considered as necessary. 
	 * Each derived class may optimize performance by selecting only the necessary data.
     */
	virtual void FillDataRequirements(TPZVec<TPZMaterialData > &datavec);
    
    /**
     * This method defines which parameters need to be initialized in order to compute the contribution of the boundary condition
     */
    virtual void FillBoundaryConditionDataRequirement(int type,TPZMaterialData &data)
    {
        // default is no specific data requirements
        if(type == 50)
        {
            data.fNeedsSol = true;
        }
    }
	
    
    /** @brief Returns the name of the material */
    virtual std::string Name() { return "no_name"; }
    
    /** @brief Returns the integrable dimension of the material */
    virtual int Dimension() = 0;
    
    int Id() const { return fId; }
    void SetId(int id) {
        if(id == 0) {
            std::cout << "\n*** Material Id can't be ZERO! ***\n";
            std::cout << "*** This Will Be a Disaster!!! ***\n";
            DebugStop();
        }
        fId = id; }
    
    /** @brief Returns the number of state variables associated with the material */
    virtual int NStateVariables() = 0;
    
    /** @brief Returns the number of components which form the flux function */
    virtual int NFluxes() {return 0;}
    
    /** @brief Prints out the data associated with the material */
    virtual void Print(std::ostream &out = std::cout);
    
    /** @brief Returns the variable index associated with the name */
    virtual int VariableIndex(const std::string &name);
    
    /** 
	 * @brief Returns the number of variables associated with the variable indexed by var. 
	 * @param var Index variable into the solution, is obtained by calling VariableIndex
	 */
    virtual int NSolutionVariables(int var);
    
    /** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout);
	
	/** @brief Returns the solution associated with the var index based on the finite element approximation */
    virtual void Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<REAL> &Solout);
    
protected:
    /** @deprecated Deprecated interface for Solution method which must use material data. */
    virtual void Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout);
    
public:
    
    /** @brief Computes the value of the flux function to be used by ZZ error estimator */
    virtual void Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol,
                      TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes,
                      TPZVec<REAL> &flux) {}
    
    /** @brief Creates an object TPZBndCond derived of TPZMaterial*/
    virtual TPZBndCond *CreateBC(TPZAutoPointer<TPZMaterial> &reference, int id, int typ, TPZFMatrix<STATE> &val1,
                                 TPZFMatrix<STATE> &val2);
    
    /** @name Contribute methods
	 * @{
	 */
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @since April 16, 2007
     */
    virtual void Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) = 0;
    
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one integration point to multiphysics simulation.
     * @param datavec [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     */
    virtual void Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef);
	
    /**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point.
     * @param data [in] stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 07, 2011
     */
    virtual void ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc) = 0;
	
	/**
     * @brief It computes a contribution to the stiffness matrix and load vector at one BC integration point
	 * to multiphysics simulation.
     * @param datavec [in]  stores all input data
     * @param weight [in] is the weight of the integration rule
     * @param ek [out] is the stiffness matrix
     * @param ef [out] is the load vector
     * @param bc [in] is the boundary condition material
     * @since October 18, 2011
     */
    virtual void ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc);
    
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
	
    /** @} */
	
    /* *Compute contribution to the energy at an integration point*/
    //      virtual void ContributeEnergy(TPZVec<REAL> &x,
    //			      TPZVec<FADFADREAL> &sol,
    //			      TPZVec<FADFADREAL> &dsol,
    //			      FADFADREAL &U,
    //			      REAL weight);
    
    //#endif
    
    
    //#ifdef _AUTODIFF
    
    /* * Compute contribution of BC to the Energy*/
    //      virtual void ContributeBCEnergy(TPZVec<REAL> & x,
    //	TPZVec<FADFADREAL> & sol, FADFADREAL &U,
    //	REAL weight, TPZBndCond &bc);
    
    //#endif
	
    /** 
	 * @brief Sets a procedure as source function for the material.
	 * @param fp pointer of the forces function
	 * @note Parameter loc corresponds to the coordinate of the point where the source function is applied
	 * @note Parameter result contains the forces resulting
	 */
    void SetForcingFunction(TPZAutoPointer<TPZFunction> fp)
    {
        fForcingFunction = fp;
    }
	
	//	 void fForcingFunctionExact(void (*fp)(TPZVec<REAL> &loc,TPZVec<REAL> &pressure,TPZVec<REAL> &flux))
	//	{
	//		fForcingFunctionExact = fp;
	//	}
	void SetForcingFunctionExact(void (*fp)(TPZVec<REAL> &loc,
											TPZVec<STATE> &pressure,TPZFMatrix<STATE> &flux))
	{
		fForcingFunctionExact = fp;
	}
    
    virtual int HasForcingFunction() {return (fForcingFunction != 0);}
	virtual int HasfForcingFunctionExact() {return (fForcingFunctionExact != 0);}
    
    /** 
     * @brief Gets the order of the integration rule necessary to integrate an
     * element with polinomial order p
     */
    virtual int IntegrationRuleOrder(int elPMaxOrder) const;
	
	/** 
     * @brief Gets the order of the integration rule necessary to integrate an
     * element multiphysic
     */
    virtual int IntegrationRuleOrder(TPZVec<int> elPMaxOrder) const;
	
	
    /* * Set the integration rule order based on the element
     *  @ param p order of interpolation, its dimension and the characteristics
     *   of the material
     *
	 virtual void SetIntegrationRule(TPZAutoPointer<TPZIntPoints> rule,
	 int elPMaxOrder,
	 int elDimension);
	 */
	
    /**
	 * @brief Computes the error due to the difference between the interpolated flux \n
	 * and the flux computed based on the derivative of the solution
	 */
    virtual void Errors(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol,
                        TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux,
                        TPZVec<REAL> &uexact, TPZFMatrix<REAL> &duexact,
                        TPZVec<REAL> &val) {
        PZError << __PRETTY_FUNCTION__ << std::endl;
        PZError << "Method not implemented! Error comparison not available. Please, implement it." << std::endl;
    }
	virtual	void ErrorsHdiv(TPZMaterialData &data, TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values){
		PZError << __PRETTY_FUNCTION__ << std::endl;
		PZError << "Nao sei o q fazer." << std::endl;
		
	}
    /** @brief Returns the number of norm errors. Default is 3: energy, L2 and H1. */
    virtual int NEvalErrors() {return 3;}
    
    /** @brief To create another material of the same type*/
    virtual TPZAutoPointer<TPZMaterial> NewMaterial();
    
    /** @brief Reads data of the material from a istream (file data)*/
    virtual void SetData(std::istream &data);
    
    /** @brief Creates a copy of the material object and put it in the vector which is passed on */
    virtual void Clone(std::map<int, TPZAutoPointer<TPZMaterial> > &matvec);
    
    /** @brief To return a numerical flux type to apply over the interfaces of the elements */
    virtual int FluxType() { return 2; }
    
    /* * Factor to diffussive term*/
    //      virtual int IdBC(REAL *x) { return 5; }
    
    virtual void ContributeErrors(TPZMaterialData &data,
                                  REAL weight,
                                  TPZVec<REAL> &nk,
                                  int &errorid){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
    }
    
    /**
     * @brief Computes square of residual of the differential equation at one integration point.
     * @param X is the point coordinate (x,y,z)
     * @param sol is the solution vector
     * @param dsol is the solution derivative with respect to x,y,z as computed in TPZShapeDisc::Shape2DFull
     */    
    virtual REAL ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<REAL> &sol, TPZFMatrix<REAL> &dsol){
        PZError << "Error at " << __PRETTY_FUNCTION__ << " - Method not implemented\n";
        return -1.;
    }
    
    /** @brief Unique identifier for serialization purposes */
    virtual int ClassId() const;
    
    /** @brief Saves the element data to a stream */
    virtual void Write(TPZStream &buf, int withclassid);
    
    /** @brief Reads the element data from a stream */
    virtual void Read(TPZStream &buf, void *context);
    
    /**
     * @brief Pushes a new entry in the context of materials with memory,
     * returning its index at the internal storage stack.
	 */
	/** To be implemented only in the proper materials. */
    virtual int PushMemItem(int sourceIndex = -1){ return -1; }
    
    /** @brief Frees an entry in the material with memory internal history storage */
    virtual void FreeMemItem(int index){ return; }
    
    /** @brief Sets fLinearContext attribute */
    void SetLinearContext(bool IsLinear);
    
    /** @brief Returns fLinearContext attribute */
    bool GetLinearContext() const;
    
};

inline bool TPZMaterial::GetLinearContext() const{
    return fLinearContext;
}

/** @brief Extern variable - Vector of force values */
extern TPZVec< void(*) ( TPZVec<REAL> &, TPZVec<STATE>& ) > GFORCINGVEC;

#endif

