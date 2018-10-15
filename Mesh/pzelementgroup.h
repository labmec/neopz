/**
 * @file
 * @brief Contains the declaration of the TPZElementGroup class, which implements an computational element which condenses the internal connects.
 */

#ifndef TPZElementGroupH
#define TPZElementGroupH

#include "pzcompel.h"
#include "pzmatred.h"
#include "pzmanvector.h"
#include "pzelmat.h"
#include "TPZOneShapeRestraint.h"

/**
 * @brief Class which groups elements to characterize dense matrices
 * @author Philippe Devloo
 * @ingroup CompElement
 */
class TPZElementGroup : public TPZCompEl
{

protected:
    TPZStack<TPZCompEl *,5> fElGroup;
    TPZManVector<int64_t,27> fConnectIndexes;
    std::map<int64_t,TPZOneShapeRestraint> fRestraints;

public:
    
    TPZElementGroup();
    
    TPZElementGroup(TPZCompMesh &mesh, int64_t &index) : TPZRegisterClassId(&TPZElementGroup::ClassId),
    TPZCompEl(mesh,0,index), fElGroup(), fConnectIndexes()
    {
        
    }
    
    /** @brief create a copy of the condensed computational element in the other mesh */
    TPZElementGroup(TPZCompMesh &mesh, const TPZElementGroup &copy);
    
    virtual ~TPZElementGroup();
    
    /** @brief add an element to the element group
     */
    virtual void AddElement(TPZCompEl *cel);

    /**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const override
    {
        out << __PRETTY_FUNCTION__ << std::endl;
        TPZCompEl::Print(out);
//        int nel = fElGroup.size();
//        for (int el=0; el<nel; el++) {
//            fElGroup[el]->Print(out);
//        }
        out << "End of " << __PRETTY_FUNCTION__ << std::endl;
    }

    /** @brief put the elements in the element group back in the mesh and delete the element group */
    void Unwrap();
    
    /** @brief Dimension of the element */
	virtual int Dimension() const override
    {
        int dimension = -1;
        int nel = fElGroup.size();
        for (int el = 0; el<nel; el++) {
            int eldim = fElGroup[el]->Dimension();
            dimension = dimension < eldim ? eldim : dimension;
        }
        return dimension;
    }
	
    /** @brief Verifies if any element needs to be computed corresponding to the material ids */
    bool NeedsComputing(const std::set<int> &matids) override;
    

    TPZStack<TPZCompEl *, 5> GetElGroup(){
        return fElGroup;
    }
    
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const override
    {
        return fConnectIndexes.size();
    }
    
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->BuildCornerConnectList(connectindexes);
        }
    }

	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const override 
    {
        return fConnectIndexes[i];
    }

	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        return new TPZElementGroup(mesh,*this);
    }
    
    /**
     * @brief Set the index i to node inode
     * @param inode node to set index
     * @param index index to be seted
     */
    virtual void SetConnectIndex(int inode, int64_t index) override;


    /** @brief Loads the solution within the internal data structure of the element */ 
	/** 
	 * Is used to initialize the solution of connect objects with dependency \n
	 * Is also used to load the solution within SuperElements
	 */
	virtual void LoadSolution() override
    {
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->LoadSolution();
        }
    }
    
    /**
     * @brief Compute the integral of a variable defined by the string if the material id is included in matids
     */
    TPZVec<STATE> IntegrateSolution(const std::string &varname, const std::set<int> &matids) override
    {
        TPZManVector<STATE,3> result;
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            TPZManVector<STATE,3> locres;
            locres = fElGroup[el]->IntegrateSolution(varname, matids);
            if (!result.size()) {
                result = locres;
            } else if(result.size() && result.size() == locres.size())
            {
                int nvar = result.size();
                for (int iv = 0; iv<nvar; iv++) {
                    result[iv] += locres[iv];
                }
            }
            else if(result.size() && locres.size() && result.size() != locres.size())
            {
                DebugStop();
            }
        }
        return result;
    }
    
    /**
     * @brief Compute the integral of a variable defined by the string if the material id is included in matids
     */
    virtual TPZVec<STATE> IntegrateSolution(int var) const override
    {
        TPZManVector<STATE,3> result;
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            TPZManVector<STATE,3> locres;
            locres = fElGroup[el]->IntegrateSolution(var);
            if (!result.size()) {
                result = locres;
            } else if(result.size() && result.size() == locres.size())
            {
                int nvar = result.size();
                for (int iv = 0; iv<nvar; iv++) {
                    result[iv] += locres[iv];
                }
            }
            else if(result.size() && locres.size() && result.size() != locres.size())
            {
                DebugStop();
            }
        }
        return result;
    }


    virtual void TransferMultiphysicsElementSolution() override
    {
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->TransferMultiphysicsElementSolution();
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

public:
    
    /**
	 * @brief Creates corresponding graphical element(s) if the dimension matches
	 * graphical elements are used to generate output files
	 * @param graphmesh graphical mesh where the element will be created
	 * @param dimension target dimension of the graphical element
	 */
	virtual void CreateGraphicalElement(TPZGraphMesh & graphmesh, int dimension) override
    {
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->CreateGraphicalElement(graphmesh,dimension);
        }
    }
	
	/** @brief Integrates a variable over the element. */
	virtual void Integrate(int variable, TPZVec<STATE> & value) override{
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->Integrate(variable,value);
        }
	}
    

	/**
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef) override;
	
    /**
	 * @brief Performs an error estimate on the elemen
	 * @param fp function pointer which computes the exact solution
	 * @param errors [out] the L2 norm of the error of the solution
	 * @param flux [in] value of the interpolated flux values
	 */
    virtual void EvaluateError(std::function<void (const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> func,
							   TPZVec<REAL> &errors, bool store_error) override;

	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	virtual void CalcResidual(TPZElementMatrix &ef) override;

    int ComputeIntegrationOrder() const override {
        std::cout << "This method should not be called. " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
		return 0;
    }
    
    public:
virtual int ClassId() const override;

protected:
    
    /// Initialize the datastructure of ek and ef based on the connect information
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const;

    /// Initialize the datastructure of ef based on the connect information
    void InitializeElementMatrix(TPZElementMatrix &ef) const;

};

#endif
