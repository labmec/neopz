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

/**
An interface to "feed" the datastructure of the Dohrmann algorithm

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
   * @param numlevels number of uniform refinements
   * @param substructlevel number of refinements which define the substructures
   */
  TPZGenSubStruct(int dimension, int numlevels, int substructlevel);

    ~TPZGenSubStruct();
    
    /// method which will generate the computational mesh
    TPZAutoPointer<TPZCompMesh> GenerateMesh();
    
    /// initialize the TPZDohrMatrix structure
    void InitializeDohr(TPZAutoPointer<TPZMatrix> dohr, TPZAutoPointer<TPZDohrAssembly> assembly);
    /// initialize the TPZDohrMatrix structure
    void InitializeDohrCondense(TPZAutoPointer<TPZMatrix> dohr, TPZAutoPointer<TPZDohrAssembly> assembly);
	
void ReorderInternalNodes(TPZSubCompMesh *sub, std::map<int,int> &globaltolocal,
                             TPZVec<int> &internalnodes);
	static void ReorderInternalNodes2(TPZSubCompMesh *sub,
							  TPZVec<int> &internalnodes, TPZVec<int> &invpermute);
	
	// computes the permutation vectors from the subcompmesh ordening to the "internal first" ordering
	// the mesh is modified during this method but is returned to its original state at the end of execution
	static void ComputeInternalEquationPermutation(TPZSubCompMesh *sub,
												   TPZVec<int> &scatterpermute, TPZVec<int> &gatherpermute);
    
  private:
    /// dimension of the mesh
    int fDimension;
    /// number of uniform refinements
    int fNumLevels;
    /// level of substructures
    int fSubstructLevel;
    
    /// computational mesh
    TPZAutoPointer<TPZCompMesh> fCMesh;
    
    /// the set of equations which correspond to corner nodes
    std::set<int> fCornerEqs;
    
    /// divide the geometric elements till num levels is achieved
    void UniformRefine();
    
public:
    /// divide the elements in substructures
    void SubStructure();
private:
    /// identify cornernodes
    void IdentifyCornerNodes();
    
	/// identify the global equations as a pair of local equation and global equation
    void IdentifyEqNumbers(TPZSubCompMesh *sub, TPZVec<std::pair<int,int> > &globaleq, std::map<int,int> &globinv);
	
    /// get the global equation numbers of a substructure (and their inverse)
    void IdentifyEqNumbers(TPZSubCompMesh *sub, TPZVec<int> &global, std::map<int,int> &globinv);

    /// Identify the corner equations associated with a substructure
    void IdentifySubCornerEqs(std::map<int,int> &globaltolocal, TPZVec<int> &cornereqs,
                             TPZVec<int> &coarseindex);
static    int NInternalEq(TPZSubCompMesh *sub);
};

/// This is a lengthy process which should run on the remote processor
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct> substruct,  TPZDohrAssembly &dohrassembly);

/// This is a lengthy process which should run on the remote processor
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstructCondense> substruct,  TPZDohrAssembly &dohrassembly);

/// return the number of submeshes
int NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/// return a pointer to the isub submesh
TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);


#endif
