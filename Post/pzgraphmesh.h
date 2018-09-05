/**
 * @file
 * @brief Contains the TPZGraphMesh class which represents a graphical mesh used for post processing purposes.
 */

#ifndef GRAFGRIDH
#define GRAFGRIDH

#include "pzcompel.h"
#include "pzadmchunk.h"
#include "pzvec.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "pzfmatrix.h"
#include "pzgraphnode.h"
#include "TPZDrawStyle.h"
#include "pzgraphel.h"

#include <iostream>
class TPZCompMesh;
template<class TVar>
class TPZBlock;


/**
 * @ingroup post
 * @brief Represents a graphical mesh used for post processing purposes. \ref post "Post processing"
 */
/** 
 * The  graphical mesh represents a logically refined version of the computational mesh \n
 * This logical refinement means that the refined element object are not actually created \n
 * They only exist in the output file. 
 */
class TPZGraphMesh : public TPZSavable {
public:
	/** @brief Constructor for graphical mesh */
    TPZGraphMesh(TPZCompMesh *cm, int dimension, TPZMaterial * mat,const TPZVec<std::string> &scalarnames,const TPZVec<std::string> &vecnames);
	/** @brief Default destructor */
	virtual ~TPZGraphMesh(void);
        
        virtual int ClassId() const;
        
        void Read(TPZStream& buf, void* context);
        
        void Write(TPZStream& buf, int withclassid) const;

	/** @brief Find graphical node (connect) */
	TPZGraphNode &FindNode(int64_t side);
	TPZGraphEl *FindElement(int64_t sid);
	/** @brief Vector of the graphical elements */
	TPZAdmChunkVector<TPZGraphEl *> &ElementList();
	/** @brief Vector of the graphical nodes */
	TPZAdmChunkVector<TPZGraphNode> &NodeMap();
	/** @brief Number of points to drawing, depending on the resolution */
	int64_t NPoints();
	int64_t NElements(MElementType type);
	/** @brief Get the resolution of the draw */
	int Res() {return fResolution;}
	/** @brief Sets the material information */
	void SetMaterial(TPZMaterial * mat) {fMaterial = mat;}
	/** @brief Sets the computational mesh to associate */
	virtual void SetCompMesh(TPZCompMesh *mesh, TPZMaterial * &mat);
	/** @brief Draw the graphical nodes information */
	virtual void DrawNodes();
	/** @brief Draw the graphical mesh */
	virtual void DrawMesh(int numcases);
	/** @brief Draw the connectivity information */
	virtual void DrawConnectivity(MElementType type);
	/** @brief Draw solution depending on the resolution */
	virtual void DrawSolution(int step, REAL time);
	/** @brief Gets the style of the graphical mesh */
	TPZDrawStyle Style();
	/** @brief Sets the filename to output of graph */
	virtual void SetFileName(const std::string &filename);
	
	std::ostream &Out()
	{
		return fOutFile;
	}
	
	/** @brief Sets resolution */
	void SetResolution(int res){ fResolution = res; SequenceNodes();}
	
	void Print(std::ostream &out);
	void SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames);
	
protected:
	/** @brief Computational mesh associated */
	TPZCompMesh *fCompMesh;
	/** @brief Geometric mesh related */
	TPZGeoMesh *fGeoMesh;
	
	/** @brief Vector of material associated */
	TPZMaterial * fMaterial;
	/** @brief Dimension of the graphical mesh */
	int fDimension;
	/** @brief Vector of graphical elements */
	TPZAdmChunkVector<TPZGraphEl *> fElementList;
	/** @brief Vector of graphical nodes (connects) */
	TPZAdmChunkVector<TPZGraphNode> fNodeMap;
	/** @brief Resolution of the graphical object */
	int fResolution;
	/** @brief Style of the graphical file */
	TPZDrawStyle fStyle;
	std::ofstream fOutFile;
	std::string fFileName;
	/** @brief Vectors of the variables names (scalar and vectorial) */
	TPZVec<std::string> fScalarNames, fVecNames;
	virtual void SequenceNodes();
	
	TPZCompEl *FindFirstInterpolatedElement(TPZCompMesh *mesh,int dimension);
	
public:
	/** @brief Return of the material for graphical mesh */
	virtual TPZMaterial * Material();
	/** @brief Return a pointer of the computational mesh */
	virtual TPZCompMesh *Mesh() { return fCompMesh;}
    
    TPZVec<std::string> ScalarNames()
    {
        return fScalarNames;
    }
    
    TPZVec<std::string> VecNames()
    {
        return fVecNames;
    }
};

inline TPZDrawStyle TPZGraphMesh::Style(){
	return fStyle;
}

#endif
