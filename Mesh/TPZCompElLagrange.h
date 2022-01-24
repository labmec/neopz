//
//  TPZCompElLagrange.h
//  PZ
//
//  Created by Philippe Devloo on 11/2/13.
//
//

#ifndef __PZ__TPZCompElLagrange__
#define __PZ__TPZCompElLagrange__

#include <iostream>

#include "pzcompel.h"
#include "pzcmesh.h"

class TPZCompElLagrange : public TPZCompEl
{
    
public:
    
    struct TLagrange
    {
        /// Which connects are linked by a Lagrange multiplier
        int64_t fConnect[2];
        /// Degree of freedom which is connected
        int fIdf[2];
        
        TLagrange()
        {
            fConnect[0] = -1;
            fConnect[1] = -1;
            fIdf[0] = -1;
            fIdf[1] = -1;
        }
    };
    
private:
    
    TPZManVector<TLagrange,3> fDef;
    
public:
    
    TPZCompElLagrange() : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(), fDef()
    {
    }
    
    TPZCompElLagrange(const TPZCompElLagrange &copy) : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(copy)
    {
        fDef = copy.fDef;
        
    }
    
    TPZCompElLagrange(TPZCompMesh &mesh, int64_t connect1, int idf1, int64_t connect2, int idf2) : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(mesh,0), fDef(1)
    {
        fDef[0].fConnect[0] = connect1;
        fDef[0].fConnect[1] = connect2;
        fDef[0].fIdf[0] = idf1;
        fDef[0].fIdf[1] = idf2;
        mesh.ConnectVec()[connect1].IncrementElConnected();
        mesh.ConnectVec()[connect2].IncrementElConnected();
#ifdef PZDEBUG
        TPZConnect &c1 = mesh.ConnectVec()[connect1];
        TPZConnect &c2 = mesh.ConnectVec()[connect2];
        if (idf1 >= c1.NShape()*c1.NState()) {
            DebugStop();
        }
        if (idf2 >= c2.NShape()*c2.NState()) {
            DebugStop();
        }
#endif
    }
    
    TPZCompElLagrange(TPZCompMesh &mesh, const TPZVec<TLagrange> &Dependencies) : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(mesh,0), fDef(Dependencies)
    {
    }
    
	
	/** @brief Put a copy of the element in the patch mesh */
	TPZCompElLagrange(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<int64_t,int64_t> &gl2lcElMap) : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(mesh,copy,gl2lcElMap)
    {
        const TPZCompElLagrange *lcop = dynamic_cast<const TPZCompElLagrange *>(&copy);
        if (!lcop) {
            DebugStop();
        }
        fDef = lcop->fDef;
        
    }
	
	/** @brief Copy of the element in the new mesh */
	TPZCompElLagrange(TPZCompMesh &mesh, const TPZCompEl &copy) : TPZRegisterClassId(&TPZCompElLagrange::ClassId),
    TPZCompEl(mesh,copy)
    {
        const TPZCompElLagrange *lcop = dynamic_cast<const TPZCompElLagrange *>(&copy);
        if (!lcop) {
            DebugStop();
        }
        fDef = lcop->fDef;
        
    }
    
    virtual ~TPZCompElLagrange();
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const override
    {
        return new TPZCompElLagrange(mesh,*this);
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

	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const override
    {
        return 2*fDef.size();
    }
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual int64_t ConnectIndex(int i) const override
    {
        if (i>=0 && i < 2*fDef.size()) {
            return fDef[i/2].fConnect[i%2];
        }
        DebugStop();
        return -1;
    }
	
	/** @brief Dimension of the element */
	virtual int Dimension() const override
    {
        return 0;
    }
    
    void CreateGraphicalElement(TPZGraphMesh &, int) override
    {
    }

    void EvaluateError(TPZVec<REAL> & errors, bool store_error) override {
        return;
    }
  

	
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<int64_t> &connectindexes) const override
    {
        for (int64_t i=0; i<fDef.size(); i++) {
            connectindexes.insert(fDef[i].fConnect[0]);
            connectindexes.insert(fDef[i].fConnect[1]);
        }
    }
    
	/**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, int64_t index) override
    {
        if (inode >= 0 && inode < 2*fDef.size()) {
            fDef[inode/2].fConnect[inode%2] = index;
        }
        else
        {
            DebugStop();
        }
    }
	
	/**
	 * @brief Computes the element stifness matrix and right hand side
	 * @param ek element stiffness matrix
	 * @param ef element load vector
	 */
	void CalcStiff(TPZElementMatrixT<STATE> &ek,TPZElementMatrixT<STATE> &ef) override{
        CalcStiffInternal(ek,ef);
    }
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	//virtual void CalcResidual(TPZElementMatrix &ef);
	
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef) override;
    public:
int ClassId() const override;

protected:
    template<class TVar>
    void CalcStiffInternal(TPZElementMatrixT<TVar> &ek, TPZElementMatrixT<TVar> &ef);

};

#endif /* defined(__PZ__TPZCompElLagrange__) */
