/**
 * @file
 * @brief Contains the TPZGenSubStruct class which is an interface to "feed" the datastructure of the Dohrmann algorithm.
 */

#ifndef TPZGENSUBSTRUCT_H
#define TPZGENSUBSTRUCT_H

#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "tpzdohrassembly.h"
#include "tpzdohrmatrix.h"

class TPZSubCompMesh;
template<class TVar>
class TPZDohrSubstruct;
template<class TVar>
class TPZDohrSubstructCondense;

/** \addtogroup substructure
 * @{
 */
/**
 * @brief An interface to "feed" the datastructure of the Dohrmann algorithm. \ref substructure "Sub Structure"
 * @author Philippe Devloo
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
    
// #ifndef STATE_COMPLEX
//     /** @brief Method which will generate the computational mesh */
//     TPZAutoPointer<TPZCompMesh> GenerateMesh();
// #endif
    
    /** @brief Initialize the TPZDohrMatrix structure */
    void InitializeDohr(TPZAutoPointer<TPZMatrix<STATE> > dohr, TPZAutoPointer<TPZDohrAssembly<STATE> > assembly);
    /** @brief Initialize the TPZDohrMatrix structure */
    void InitializeDohrCondense(TPZAutoPointer<TPZMatrix<STATE> > dohr, TPZAutoPointer<TPZDohrAssembly<STATE> > assembly);
	
	void ReorderInternalNodes(TPZSubCompMesh *sub, std::map<int,int> &globaltolocal,
							  TPZVec<int> &internalnodes);
	static void ReorderInternalNodes2(TPZSubCompMesh *sub,
									  TPZVec<int> &internalnodes, TPZVec<int64_t> &invpermute);
	
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
	
	static    int64_t NInternalEq(TPZSubCompMesh *sub);
	
};

/** @brief This is a lengthy process which should run on the remote processor */
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstruct<STATE> > substruct,  TPZDohrAssembly<STATE> &dohrassembly);

/** @brief This is a lengthy process which should run on the remote processor */
void InitializeMatrices(TPZSubCompMesh *sub, TPZAutoPointer<TPZDohrSubstructCondense<STATE> > substruct,  TPZDohrAssembly<STATE> &dohrassembly);

/** @brief Return the number of submeshes */
int64_t NSubMesh(TPZAutoPointer<TPZCompMesh> compmesh);

/** @brief Return a pointer to the isub submesh */
TPZSubCompMesh *SubMesh(TPZAutoPointer<TPZCompMesh> compmesh, int isub);

/** @} */

#endif
