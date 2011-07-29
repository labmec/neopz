#ifndef V3DGRAFGH
#define V3DGRAFGH

#include "pzgraphmesh.h"
#include "pzvec.h"

/**
 * @ingroup post
 * @brief To export a graphical three dimensional mesh to use at V3D package.
 */
/**
 * V3D is a 3D Image Visualization package. \n
 * You can to download from http://penglab.janelia.org/proj/v3d/V3D/Download.html
 */
class TPZV3DGraphMesh : public TPZGraphMesh {
	
public:
	
	TPZV3DGraphMesh(TPZCompMesh *cmesh, int dimension, TPZAutoPointer<TPZMaterial> mat);
	TPZV3DGraphMesh(TPZCompMesh *cmesh,int dim,TPZV3DGraphMesh *graph,TPZAutoPointer<TPZMaterial> mat);
	
	virtual ~TPZV3DGraphMesh()
	{
	}
	
	/** @brief Draw the nodal coordinates and the connectivity */
	virtual void DrawMesh(int numcases);
	
	/** @brief Draw the solution associated with Sol (not implemented) */
	virtual void DrawSolution(TPZBlock &Sol);
	/** @brief Draw the solution associated with the variable name */
	virtual void DrawSolution(char * var = 0);
	/** @brief Draw the solution sequence */
	virtual void DrawSolution(int step, REAL time);
	virtual void SequenceNodes();
	
	private :
	TPZCompMesh *fMesh;
	int fNumCases;
	int fInterval;
	int fLoadStep;
	int fTotScal;
	int fNumScal[6];
};

#endif

