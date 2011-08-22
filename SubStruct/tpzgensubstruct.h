/**
 * @file
 * @brief Contains the TPZGenSubStruct class which is an interface to "feed" the datastructure of the Dohrmann algorithm.
 */
/***************************************************************************
 *   Copyright (C) 2006 by Philippe Devloo   *
 *   phil@fec.unicamp.br   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#ifndef TPZGENSUBSTRUCT_H
#define TPZGENSUBSTRUCT_H

#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "tpzdohrassembly.h"

template<class T>
class TPZDohrMatrix;
class TPZSubCompMesh;
class TPZDohrSubstruct;
class TPZDohrSubstructCondense;

/** \addtogroup substructure
 * @{
 */
/**
 @brief An interface to "feed" the datastructure of the Dohrmann algorithm. \ref substructure "Sub Structure"
 @author Philippe Devloo
 */
class TPZGenSubStruct{
public:
	enum ETypemesh
	{
		DistMaterial,RandomMat
	};
	
	ETypemesh fMatDist;
	
	TPZVec<REAL> fK;
	/**
	 * @brief Constructor
	 * @param dimension
	 * @param numlevels number of uniform refinements
	 * @param substructlevel number of refinements which define the substructures
	 */
	TPZGenSubStruct(int dimension, int numlevels, int substructlevel);
	
    ~TPZGenSubStruct();
    
    /** @brief Method which will generate the computational mesh */
    TPZAutoPointer<TPZCompMesh> GenerateMesh();
    
    /** @brief Initialize the TPZDohrMatrix structure */
    void InitializeDohr(TPZAutoPointer<TPZMatrix> dohr, TPZAutoPointer<TPZDohrAssembly> assembly);
    /** @brief Initialize the TPZDohrMatrix structure */
    void InitializeDohrCondense(TPZAutoPointer<TPZMatrix> dohr, TPZAutoPointer<TPZDohrAssembly> assembly);
	
	void ReorderInternalNodes(TPZSubCompMesh *sub, std::map<int,int> &globaltolocal,
							  TPZVec<int> &internalnodes);
	static void ReorderInternalNodes2(TPZSubCompMesh *sub,
									  TPZVec<int> &internalnodes, TPZVec<int> &invpermute);
	
	/** @brief Computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering
	 *
	 * The mesh is modified during this method but is returned to its original state at the end of execution
	 */
	static void ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
												   TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute);
    
private:
    /** @brief Dimension of the mesh */
    int fDimension;
    /** @brief Number of uniform refinements */
    int fNumLevels;
    /** @brief Level of substructures */
    int fSubstructLevel;
    
    /** @brief computational mesh */
    TPZAutoPointer<TPZCompMesh> fCMesh;
    
    /** @brief The set of equations which correspond to corner nodes */
    std::set<int> fCornerEqs;
    
    /** @brief Divide the geometric elements till num levels is achieved */
    void UniformRefine();
    
public:
    /** @brief Divide the elements in substructures */
    void SubStructure();
private:
    /** @brief Identify cornernodes */
    void IdentifyCornerNodes();
    
	/** @brief Identify the global equations as a pair of local equation and global equation */
    void IdentifyEqNumbers(TPZSubCompMesh *sub, TPZVec<std::pair<int,int> > &globaleq, std::map<int,int> &globinv);
	
    /** @brief Get the global equation numbers of a substructure (and their inverse) */
    void IdentifyEqNumbers(TPZSubCompMesh *sub, TPZVec<int> &global, std::map<int,int> &globinv);
	
    /** @brief Identify the corner equations associated with a substructure */
    void IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
							  TPZVec<int> &coarseindex);
	
	static    int NInternalEq(TPZSubCompMesh *sub);
	
};

/** @brief This is a lengthy process which should run on the remote processor */
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct> substruct,  TPZDohrAssembly &dohrassembly);

/** @brief This is a lengthy process which should run on the remote processor */
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstructCondense> substruct,  TPZDohrAssembly &dohrassembly);

/** @brief Return the number of submeshes */
int NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/** @brief Return a pointer to the isub submesh */
TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

/** @} */

#endif
