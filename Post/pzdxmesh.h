/**
 * @file
 * @brief Contains the TPZDXGraphMesh class which implements the interface of the graphmesh to the OpenDX graphics package.
 */

#ifndef DXMESHH
#define DXMESHH

#include "pzgraphmesh.h"
#include "pzstack.h"

/**
 * @ingroup post
 * @brief Implements the interface of the graphmesh to the OpenDX graphics package. \ref post "Post processing"
 */
/** This is the most actively supported graphics interface of the pz environment */
class TPZDXGraphMesh : public TPZGraphMesh {
	
	int fNextDataField;
	int fNumConnectObjects[8];
	int fElConnectivityObject[8];
	int fNodePosObject[8];
	int fNormalObject;
	TPZStack<REAL> fTimes;
	TPZStack<int> fFirstFieldValues[3];
	int fNumCases;
	std::string fElementType;
	
public:
	
	/** @brief Constructor for output in DX format */
	TPZDXGraphMesh(TPZCompMesh *mesh, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames,const TPZVec<std::string> &vecnames);
	/** @brief Copy constructor */
	TPZDXGraphMesh(TPZCompMesh *cmesh,int dim,TPZDXGraphMesh *graph);
	/** @brief Default destructor */
	virtual ~TPZDXGraphMesh();
	
	/** @brief Sets the name of the output file */
	virtual void SetFileName(const std::string &filename) override;
	
	/** @brief Draw mesh as dx file */
	virtual void DrawMesh(int numcases) override;
	/** @brief Draw solution as dx file */
	virtual void DrawSolution(TPZBlock &Sol);
	/** @brief Draw solution as dx file for variable name indicated */
	virtual void DrawSolution(char * var = 0);
	/** @brief Draw solution as dx file for step and time given */
	virtual void DrawSolution(int step, REAL time) override;
	
	/** @brief Number of nodes */
	int  NNodes();
	std::string ElementName();
	void Close();
	
	void SetNumCases(int numcases) { fNumCases = numcases; }
	
private:
	void DrawNormals(int numnormals);
};

#endif
