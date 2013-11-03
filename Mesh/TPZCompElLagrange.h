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
    /// Which connects are linked by a Lagrange multiplier
    long fConnect[2];
    /// Degree of freedom which is connected
    int fIdf[2];
    
public:
    
    TPZCompElLagrange() : TPZCompEl()
    {
        for(int i=0;i<2;i++)
        {
            fConnect[i] = -1;
            fIdf[i] = -1;
        }
    }
    
    TPZCompElLagrange(const TPZCompElLagrange &copy) : TPZCompEl(copy)
    {
        for(int i=0;i<2;i++)
        {
            fConnect[i] = copy.fConnect[i];
            fIdf[i] = copy.fIdf[i];
        }
        
    }
    
    TPZCompElLagrange(TPZCompMesh &mesh, long connect1, int idf1, long connect2, int idf2, long &index) : TPZCompEl(mesh,0,index)
    {
        fConnect[0] = connect1;
        fConnect[1] = connect2;
        fIdf[0] = idf1;
        fIdf[1] = idf2;
    }
    
	/** @brief Put a copy of the element in the referred mesh */
	TPZCompElLagrange(TPZCompMesh &mesh, const TPZCompEl &copy) : TPZCompEl(mesh,copy)
    {
        const TPZCompElLagrange *lcop = dynamic_cast<const TPZCompElLagrange *>(&copy);
        if (!lcop) {
            DebugStop();
        }
        for(int i=0;i<2;i++)
        {
            fConnect[i] = lcop->fConnect[i];
            fIdf[i] = lcop->fConnect[i];
        }
    }
	
	/** @brief Put a copy of the element in the patch mesh */
	TPZCompElLagrange(TPZCompMesh &mesh, const TPZCompEl &copy, std::map<long,long> &gl2lcElMap) : TPZCompEl(mesh,copy,gl2lcElMap)
    {
        const TPZCompElLagrange *lcop = dynamic_cast<const TPZCompElLagrange *>(&copy);
        if (!lcop) {
            DebugStop();
        }
        for(int i=0;i<2;i++)
        {
            fConnect[i] = lcop->fConnect[i];
            fIdf[i] = fIdf[i];
        }
        
    }
	
	/** @brief Copy of the element in the new mesh with alocated index */
	TPZCompElLagrange(TPZCompMesh &mesh, const TPZCompEl &copy, long &index) : TPZCompEl(mesh,copy,index)
    {
        const TPZCompElLagrange *lcop = dynamic_cast<const TPZCompElLagrange *>(&copy);
        if (!lcop) {
            DebugStop();
        }
        for(int i=0;i<2;i++)
        {
            fConnect[i] = lcop->fConnect[i];
            fIdf[i] = lcop->fIdf[i];
        }
        
    }
    
    virtual ~TPZCompElLagrange();
	
	/** @brief Method for creating a copy of the element */
	virtual TPZCompEl *Clone(TPZCompMesh &mesh) const
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
									std::map<long,long> & gl2lcConMap,
									std::map<long,long> & gl2lcElMap) const;

	/** @brief Returns the number of nodes of the element */
	virtual int NConnects() const
    {
        return 2;
    }
	
	/**
	 * @brief Returns the index of the ith connectivity of the element
	 * @param i connectivity index who want knows
	 */
	virtual long ConnectIndex(int i) const
    {
        if (i==0 || i==1) {
            return fConnect[i];
        }
        DebugStop();
        return -1;
    }
	
	/** @brief Dimension of the element */
	virtual int Dimension() const
    {
        return 0;
    }
	
    /** @brief adds the connect indexes associated with base shape functions to the set */
    virtual void BuildCornerConnectList(std::set<long> &connectindexes) const
    {
        connectindexes.insert(fConnect, fConnect+2);
    }
    
	/**
	 * @brief Set the index i to node inode
	 * @param inode node to set index
	 * @param index index to be seted
	 */
	virtual void SetConnectIndex(int inode, long index)
    {
        if (inode == 0 || inode == 1) {
            fConnect[inode] = index;
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
	virtual void CalcStiff(TPZElementMatrix &ek,TPZElementMatrix &ef);
	
	/**
	 * @brief Computes the element right hand side
	 * @param ef element load vector(s)
	 */
	//virtual void CalcResidual(TPZElementMatrix &ef);
	
    void InitializeElementMatrix(TPZElementMatrix &ek, TPZElementMatrix &ef);


};

#endif /* defined(__PZ__TPZCompElLagrange__) */
