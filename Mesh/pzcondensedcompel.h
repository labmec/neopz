/**
 * @file
 * @brief Contains the declaration of the TPZCondensedCompEl class, which implements an computational element which condenses the internal connects.
 */

#ifndef TPZCONDENSEDCOMPELH
#define TPZCONDENSEDCOMPELH

#include "pzcompel.h"
#include "pzmatred.h"
#include "pzmanvector.h"
#include "pzelmat.h"


/**
 * @brief Class which implements an element which condenses the internal connects
 * @author Philippe Devloo
 * @ingroup CompElement
 */
class TPZCondensedCompEl : public TPZCompEl
{

protected:
    
    //TPZMatRed<REAL, TPZFMatrix<REAL> > fCondensed;
    int64_t fNumInternalEqs = 0;
    int64_t fNumTotalEqs = 0;
	TPZMatRed<STATE, TPZFMatrix<STATE> > fCondensed;
    TPZCompEl *fReferenceCompEl;
    TPZManVector<int64_t,62> fIndexes;
    TPZManVector<int64_t,55> fCondensedConnectIndexes;
    TPZManVector<int64_t,10> fActiveConnectIndexes;
    bool fKeepMatrix = true;
    void Resequence();

    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
    template<class TVar>
    void CalcResidualInternal(TPZElementMatrixT<TVar> &ef);
public:
    
    TPZCondensedCompEl(TPZCompEl *ref, bool keepmatrix = true);
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZCondensedCompEl(const TPZCondensedCompEl &copy, TPZCompMesh &mesh);
    
    virtual ~TPZCondensedCompEl();

    /**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const override;

    /** @brief unwrap the condensed element from the computational element and delete the condensed element */
    void Unwrap();
    
    /**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int64_t index) override;
    
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const override 
    {
        return fActiveConnectIndexes.size();
//        return fReferenceCompEl->NConnects();
    }
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const override 
    {
        return fActiveConnectIndexes[i];
//        return fReferenceCompEl->ConnectIndex(fIndexes[i]);
    }

    TPZCompEl * ReferenceCompEl(){
        return fReferenceCompEl;
    }
    
    /// return true if the element has a variational statement associated with the material ids
    virtual bool NeedsComputing(const std::set<int> &materialids) override
    {
        if(fReferenceCompEl)
        {
            return fReferenceCompEl->NeedsComputing(materialids);
        }
        else
        {
            return false;
        }
    }

    
    virtual void LoadElementReference() override
    {
        if(fReferenceCompEl)
        {
            fReferenceCompEl->LoadElementReference();
        }
    }

    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override;
	
    /// Set the flag that determines whether the matrix needs to be kept or not
    void SetKeepMatrix(bool keep)
    {
        fKeepMatrix = keep;
    }
    
    TPZMatRed<STATE, TPZFMatrix<STATE> > &Matrix()
    {
        return fCondensed;
    }

	/** @brief Dimension of the element */
	virtual int Dimension() const override
    {
        return fReferenceCompEl->Dimension();
    }
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        return new TPZCondensedCompEl(*this,mesh);
    }

    /** @brief Loads the solution within the internal data structure of the element */ 
	/** 
	 * Is used to initialize the solution of connect objects with dependency \n
	 * Is also used to load the solution within SuperElements
	 */
	virtual void LoadSolution() override;
    
    virtual void TransferMultiphysicsElementSolution() override
    {
        if(fReferenceCompEl)
        {
            fReferenceCompEl->TransferMultiphysicsElementSolution();
        }
    }


	/**
	 * @brief Method for creating a copy of the element in a patch mesh
	 * @param mesh Patch clone mesh
	 * @param gl2lcConMap map the connects indexes from global element (original) to the local copy.
	 * @param gl2lcElMap map the computational elements
     */
	/**
	 * Otherwise of the previous clone function, this method don't
	 * copy entire mesh. Therefore it needs to map the connect index
	 * from the both meshes - original and patch
	 */
	virtual TPZCompEl *ClonePatchEl(TPZCompMesh &mesh,
									std::map<int64_t,int64_t> & gl2lcConMap,
									std::map<int64_t,int64_t> & gl2lcElMap) const override;
    
    /**
     * @brief Compute the integral of a variable defined by the string if the material id is included in matids
     */
    virtual TPZVec<STATE> IntegrateSolution(const std::string &varname, const std::set<int> &matids) override
    {
        return fReferenceCompEl->IntegrateSolution(varname, matids);
    }
    /**
     * @brief Compute the integral of a variable defined by the string if the material id is included in matids
     */
    virtual TPZVec<STATE> IntegrateSolution(int var) const override
    {
        return fReferenceCompEl->IntegrateSolution(var);
    }

    
    /// Assemble the stiffness matrix in locally kept datastructure
    virtual void Assemble() override;
	/**
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override{
        return CalcStiffInternal(ek,ef);
    }
	
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	void CalcResidual(TPZElementMatrixT<STATE> &ef) override{
        return CalcResidualInternal(ef);
    }
    
    /** @brief Verifies if the material associated with the element is contained in the set */
    virtual bool HasMaterial(const std::set<int> &materialids) const override;

  void EvaluateError(TPZVec<REAL> &errors, bool store_errors) override {
    fReferenceCompEl->EvaluateError(errors, store_errors);
  }

    /**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension) override
    {
        fReferenceCompEl->CreateGraphicalElement(graphmesh, dimension);
    }
    

    int ComputeIntegrationOrder() const override {
        std::cout << "This method should not be called. " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
		return 0;
    }

    void PermuteActiveConnects(TPZManVector<int64_t> &perm);

virtual int ClassId() const override;


};

#endif
