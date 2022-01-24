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
    TPZStack<TPZCompEl *,10> fElGroup;
    TPZManVector<int64_t,27> fConnectIndexes;
    std::map<int64_t,TPZOneShapeRestraint> fRestraints;

public:
    
    TPZElementGroup();
    
    TPZElementGroup(TPZCompMesh &mesh) : TPZRegisterClassId(&TPZElementGroup::ClassId),
    TPZCompEl(mesh,0), fElGroup(), fConnectIndexes()
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
        out << "comp and geom indexes :";
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            int64_t cindex = fElGroup[el]->Index();
            TPZGeoEl *gel = fElGroup[el]->Reference();
            int64_t gindex = -1;
            if(gel) gindex = gel->Index();
            out << " " << cindex << "|" << gindex;
//            fElGroup[el]->Print(out);
        }
        out << std::endl;
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
    

    /// Reorder the connects in increasing number of elements connected
    void ReorderConnects();

    void ReorderConnects(TPZManVector<int64_t> &connects);
    
    const TPZVec<TPZCompEl *> &GetElGroup(){
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
    
    virtual void LoadElementReference() override
    {
        int nel = fElGroup.size();
        for (int el=0; el<nel; el++) {
            fElGroup[el]->LoadElementReference();
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
	virtual void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override{
        CalcStiffInternal<STATE>(ek,ef);
    }
    /**
     * @brief Performs an error computation for the element
     * @param errors [out] the L2 norm of the error of the solution
     * @param flux [in] value of the interpolated flux values
     */
    void EvaluateError(TPZVec<REAL> &errors, bool store_error) override;
    
    /** @brief Verifies if the material associated with the element is contained in the set */
    bool HasMaterial(const std::set<int> &materialids) const override;
    
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	virtual void CalcResidual(TPZElementMatrixT<STATE> &ef) override{
        CalcResidualInternal<STATE>(ef);
    }

    int ComputeIntegrationOrder() const override {
        std::cout << "This method should not be called. " << __PRETTY_FUNCTION__ << std::endl;
        DebugStop();
		return 0;
    }
    
    public:
virtual int ClassId() const override;

protected:
    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);
    template<class TVar>
    void CalcResidualInternal(TPZElementMatrixT<TVar> &ef);
    /// Initialize the datastructure of ek and ef based on the connect information
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) const;

    /// Initialize the datastructure of ef based on the connect information
    void InitializeElementMatrix(TPZElementMatrix &ef) const;

};

#endif
