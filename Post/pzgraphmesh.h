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
    TPZGraphMesh(TPZCompMesh *cm, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames,const TPZVec<std::string> &vecnames);
    
    /** @brief Constructor with tensorial names for graphical mesh */
    TPZGraphMesh(TPZCompMesh *cm, int dimension, const std::set<int> & matids, const TPZVec<std::string> &scalarnames, const TPZVec<std::string> &vecnames, const TPZVec<std::string> &tensornames);
    
	/** @brief Default destructor */
	virtual ~TPZGraphMesh(void);
        
    int ClassId() const override;

    void Read(TPZStream &buf, void *context) override;

    void Write(TPZStream &buf, int withclassid) const override;

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
	/** @brief Sets the computational mesh to associate */
    virtual void SetCompMesh(TPZCompMesh *mesh, const std::set<int> & matids);
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
	
    /** @brief Print object attributes  */
	void Print(std::ostream &out);
    
    /** @brief  Set names with scalar and vector variable names  */
	void SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames);

    /** @brief  Set names with scalar, vectorial and tensorial variable names  */
    void SetNames(const TPZVec<std::string>&scalarnames, const TPZVec<std::string>&vecnames, const TPZVec<std::string>& tensornames);
	
    /** @brief Set material ids  */
    void SetMaterialIds(const std::set<int> & matids);
    
    /** @brief Get material ids  */
    std::set<int> MaterialIds();
    
    /** @brief Return a directive if the material id is being postprocessed */
    bool Material_Is_PostProcessed(int matid);
    
protected:
    
	/** @brief Computational mesh associated */
	TPZCompMesh *fCompMesh;
	/** @brief Geometric mesh related */
	TPZGeoMesh *fGeoMesh;
	/** @brief Set of material ids being post-processed */
    std::set<int> fMaterialIds;
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
	/** @brief Vectors of the variables names (scalar, vectorial, and tensorial) */
	TPZVec<std::string> fScalarNames, fVecNames, fTensorNames;
	virtual void SequenceNodes();
	
	TPZCompEl *FindFirstInterpolatedElement(TPZCompMesh *mesh,int dimension);
	
public:
    
//    /** @brief Return of the material for graphical mesh */
//    virtual TPZMaterial * Material();
	/** @brief Return a pointer of the computational mesh */
	virtual TPZCompMesh *Mesh() { return fCompMesh;}
    
    /** @brief Return scalar variable names */
    TPZVec<std::string> ScalarNames()
    {
        return fScalarNames;
    }
    
    /** @brief Return vectorial variable names */
    TPZVec<std::string> VecNames()
    {
        return fVecNames;
    }
    
    /** @brief Return tensorial variable names */
    TPZVec<std::string> TensorNames()
    {
        return fTensorNames;
    }
};

inline TPZDrawStyle TPZGraphMesh::Style(){
	return fStyle;
}

#endif
