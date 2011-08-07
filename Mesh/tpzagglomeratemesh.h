/**
 * @file
 * @brief Contains declaration of TPZAgglomerateMesh which implements a mesh that contains agglomerated elements.
 */
#ifndef TPZAGGLOMERATEMESH_H
#define TPZAGGLOMERATEMESH_H

#include "pzflowcmesh.h"


/**
 * This class contains both discontinuous, continuous and agglomerated elements. \n
 * Its distinction from other meshes is that it points to a reference fine mesh
 */
/**
 * @brief Implements a mesh that contains agglomerated elements. \ref CompMesh "Computational Mesh"
 * @ingroup CompMesh
 * @author Philippe R. B. Devloo
 * @since 2004.
 */
class TPZAgglomerateMesh : public TPZFlowCompMesh
{
public:
    TPZAgglomerateMesh() : TPZFlowCompMesh(0)
    {
		fFineMesh = 0;
    }
    
    /**
	 @brief An agglomeratemesh needs a fine mesh to relate to, because its elements may
	 point to elements of the finemesh
	 */
    TPZAgglomerateMesh(TPZCompMesh *finemesh) : TPZFlowCompMesh(finemesh->Reference()),
	fFineMesh(finemesh)
    {
    }
	
	virtual ~TPZAgglomerateMesh()
	{
	}
    
	/**
	 @brief Returns a pointer to the associated fine mesh
	 */
	TPZCompMesh *FineMesh()
	{
		return fFineMesh;
	}
	
	
	
private:
	
	/**
	 @brief Reference Mesh for agglomerated elements
	 */
    TPZCompMesh *fFineMesh;
	
};

#endif
