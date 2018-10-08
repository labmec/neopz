/**
 * @file
 * @brief Contains the declaration of multiphysic interface class
 * @author Agnaldo
 * @since 10/26/11.
 */

#ifndef TPZMULTIPHYSICSINTERFACEELH
#define TPZMULTIPHYSICSINTERFACEELH 

#include <iostream>

#include "pzcompel.h"
#include "pzmultiphysicselement.h"


/**
 * @brief Computes the contribution over an interface between two discontinuous elements. \ref CompElement "Computational Element"
 * @ingroup CompElement
 */
class TPZMultiphysicsInterfaceElement : public TPZCompEl {

protected:

	/** @brief Element vector the left of the normal a interface */
	TPZCompElSide 	fLeftElSide;
		
	/** @brief Element vector the right of the normal a interface */
	TPZCompElSide 	fRightElSide;
    
    /** @brief indexes of the connects */
    TPZManVector<int64_t,20> fConnectIndexes;
    
    /** @brief indices of the Left Element Vector */
    TPZManVector<int64_t,3> fLeftElIndices;
    
    /** @brief indices of the Right Element Vector */
    TPZManVector<int64_t,3> fRightElIndices;
    
	
public:
	/** @brief Default constructor */
	TPZMultiphysicsInterfaceElement();
	
	/** @brief Constructor */
	TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, TPZGeoEl *ref, int64_t &index, TPZCompElSide left, TPZCompElSide right);
    
    /** @brief create a copy of the given element */
    TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy);
    
    /** @brief create a copy of the given element using index mapping */
    TPZMultiphysicsInterfaceElement(TPZCompMesh &mesh, const TPZMultiphysicsInterfaceElement &copy, std::map<int64_t,int64_t> & gl2lcConMap,
									std::map<int64_t,int64_t> & gl2lcElMap);
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override 
    {
        return new TPZMultiphysicsInterfaceElement(mesh,*this);
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
									std::map<int64_t,int64_t> & gl2lcElMap) const override
    {
        return new TPZMultiphysicsInterfaceElement(mesh,*this,gl2lcConMap,gl2lcElMap);
    }
	
	/** @brief Informs the connect that this element is connected to it. */
	void IncrementElConnected();	
	
	/** @brief Default destructor */
	virtual ~TPZMultiphysicsInterfaceElement();
	
	/**
	 * @brief Compute the transform of a paramenter point in the multiphysic interface element to a parameter point in the neighbor super element
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param transf [out] vector of Transforms 
	 */
	void ComputeSideTransform(TPZManVector<TPZCompElSide> &Neighbor, TPZManVector<TPZTransform<> > &transf);
	
	/**
	 * @brief Maps qsi coordinate at this master element to qsi coordinate at neighbor master element.
	 * @param Neighbor [in] may be this->LeftElementSide() or this->RightElementSide()
	 * @param qsi [in] is the point at this element master
	 * @param NeighIntPoint [out] is the point at neighbor element master. X[qsi] is equal to X[NeighIntPoint]
	 */
	void MapQsi(TPZManVector<TPZCompElSide> &Neighbor, TPZVec<REAL> &qsi, TPZVec<REAL> &NeighIntPoint);
    
    /**
     * Add elements to the list of left and right elements
     */
    void SetLeftRightElement(const TPZCompElSide &leftel, const TPZCompElSide &rightel);
    
    /**
     * Add elements to the list of left and right indices given related elements
     */
    void SetLeftRightElementIndices(const TPZVec<int64_t> &lefindices, const TPZVec<int64_t> &rightindices);
    
	/**
	 * Get left and right elements
	 */	
	void GetLeftRightElement(TPZCompElSide &leftel, TPZCompElSide &rightel);
	
    /**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int64_t index) override
    {
        fConnectIndexes[inode] = index;
    }
	
    /** @brief Returns the number of nodes of the element */
	virtual int NConnects() const override;
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const override;
	

    /** @brief Dimension of the element */
	virtual int Dimension() const override
    {
        return Reference()->Dimension();
    }

    /**
     * Compute the stiffness matrix and load vector of the interface element
     */
    void CalcStiff(TPZElementMatrix &ek, TPZElementMatrix &ef) override;
    
    /**
     * Compute the load vector of the interface element
     */
    void CalcStiff(TPZElementMatrix &ef);

    /**
     * Return max integration rule of this interface element
     */
    void CreateIntegrationRule();
    
    /**
     * Return max integration rule of this interface element
     */
    const TPZIntPoints & GetIntegrationRule();
    
    virtual int ComputeIntegrationOrder() const override;
    
    /** @brief Compute and fill data with requested attributes for each of the compels in fElementVec*/
    virtual void ComputeRequiredData(TPZVec<REAL> &intpointtemp, TPZVec<TPZTransform<> > &trvec, TPZVec<TPZMaterialData> &datavec);

    /** @brief Initialize the structure of the stiffness matrix */
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);
    
    /** @brief Initialize the structure of the stiffness matrix */
    void InitializeElementMatrix(TPZElementMatrix &ef);
    
    
    /** @brief access function to the left element */
    TPZCompElSide Left() const
    {
        return fLeftElSide;
    }
    
    /** @brief Returns the right element from the element interface */
	TPZCompEl *RightElement() const {
		return fRightElSide.Element();
	}
	
	/** @brief Returns the left element from the element interface */
	TPZCompEl *LeftElement() const {
		return fLeftElSide.Element();
	}

    void ComputeCenterNormal(TPZVec<REAL> &normal) const;
    
   // void ComputeNormal(TPZVec<REAL> &qsi, TPZVec<REAL> &normal);
    
    /**
	 * @brief Prints element data
	 * @param out Indicates the device where the data will be printed
	 */
	virtual void Print(std::ostream &out = std::cout) const override;
	
    /** @brief Initialize the material data for the neighbouring element */
    void InitMaterialData(TPZVec<TPZMaterialData> &data, TPZMultiphysicsElement *mfcel, TPZVec<int64_t> *indices=0);
    
    /** @brief initialize the material data for the geometric data */
    void InitMaterialData(TPZMaterialData &data);
    
    /** @brief Compute the data needed to compute the stiffness matrix at the integration point */
    virtual void ComputeRequiredData(TPZMaterialData &data, TPZVec<REAL> &point);
    
    /** @brief Compute the required data from the neighbouring elements */
    void ComputeRequiredData(TPZVec<REAL> &point, TPZVec<TPZTransform<> > &trvec, TPZMultiphysicsElement *Neighbour, TPZVec<TPZMaterialData> &data);
	
	void ComputeSideTransform(TPZCompElSide &Neighbor, TPZTransform<> &transf);
    
    /** @brief Access function to the right element */
    TPZCompElSide Right() const
    {
        return fRightElSide;
    }
    
    virtual int nmeshes()
    {
        if (fLeftElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fLeftElSide.Element());
            if(mfcel)
            {
                return mfcel->NMeshes();
            }
        }
        if (fRightElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element()); 
            if (mfcel) {
                return mfcel->NMeshes();
            }
        }
        return 1;
    }
	
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        TPZCompEl *left = fLeftElSide.Element();
        TPZCompEl *right = fRightElSide.Element();
        left->BuildCornerConnectList(connectindexes);
        right->BuildCornerConnectList(connectindexes);
    }
	
	/**
	 * @brief Calculates the solution - sol - for the variable var
	 * at point qsi, where qsi is expressed in terms of the
	 * master element coordinates
	 * @param qsi master element coordinate
	 * @param var variable name
	 * @param sol vetor for the solution
	 */
	virtual void Solution(TPZVec<REAL> &qsi,int var,TPZVec<STATE> &sol) override;
	
	virtual void CreateGraphicalElement(TPZGraphMesh &grmesh, int dimension) override;
    
    /** @brief Return the size of the elementvec in multiphysics, if it is not multiphysics, just return 1 */
    virtual int NumberOfCompElementsInsideThisCompEl() override {
        
        if (fLeftElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fLeftElSide.Element());
            if(mfcel)
            {
                return mfcel->NMeshes();
            }
        }
        if (fRightElSide) {
            TPZMultiphysicsElement *mfcel = dynamic_cast<TPZMultiphysicsElement *>(fRightElSide.Element());
            if (mfcel) {
                return mfcel->NMeshes();
            }
        }
        return 0;

    }	

    public:
virtual int ClassId() const override;

void EvaluateError(std::function<void(const TPZVec<REAL> &loc,TPZVec<STATE> &val,TPZFMatrix<STATE> &deriv)> fp,
                                  TPZVec<REAL> &/*errors*/, bool store_error) override {
//        LOGPZ_WARN(logger, "EvaluateError is called.");
//        DebugStop();
    }

    
};

#endif
