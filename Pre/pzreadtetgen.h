/**
 * @file
 * @brief Contains the TPZReadTetGen class which implement the interface between TPZGeoMesh and the files produced by tetgen.
 */
//$Id: pzreadtetgen.h,v 1.1 2006-03-04 15:39:36 tiago Exp $

#ifndef TPZREADTETGEN
#define TPZREADTETGEN

#include <string>
#include <map>
class TPZGeoMesh;


/** 
 * @ingroup pre
 * @brief Implement the interface between TPZGeoMesh and the files produced by tetgen. \ref pre "Getting Data"
 * @since March 03, 2006
 */
/** What is tetgen ? Take a look on http://tetgen.berlios.de/index.html
 * "TetGen - A Quality Tetrahedral Mesh Generator and Three-Dimensional Delaunay Triangulator. 
 *  Hang Si
 *  Research Group of Numerical Mathematics and Scientific Computing Weierstrass Institute for Applied Analysis and Stochastics
 *  Mohrenstr. 39  10117 Berlin, Germany
 *  si@wias-berlin.de"
 */
class TPZReadTetGen{
public:
	
    TPZReadTetGen();
    
    ~TPZReadTetGen();
    
    /**
     * @brief Convert tetgen files in a TPZGeoMesh object
	 *
     * If something does not work in the process, a null pointer is returned.
     */
    TPZGeoMesh * Process(std::string NodeFileName, std::string FaceFileName, std::string TetraFileName);
    
private:
	
    /** 
     * @brief Process nodes.
	 *
     * Return true if the process worked succesfuly and false otherwise.
     */
    bool ProcessNodes(std::string NodeFileName,  TPZGeoMesh &gmesh, int & numbernodes);
    
    /** 
     * @brief Process faces.
	 *
     * Returns true if the process worked succesfuly and false otherwise.
     */    
    bool ProcessFaces(std::string FaceFileName,  TPZGeoMesh &gmesh, int & numberfaces);
    
    /** 
     * @brief Process tetrahedras.
	 *
     * Returns true if the process worked succesfuly and false otherwise.
     */    
    bool ProcessTetra(std::string TetraFileName, TPZGeoMesh &gmesh, int & numbervols);
    
    /**
     * @brief Nodes in tetgen are counted from 1 to n as a fortran based code.
	 *
     * In PZ, nodes have one Id and one index. index is defined by TPZGeoMesh
     * when a new node is allocated in NodeVec().
     * This object maps from the Id to the index. fNodeIndices[ Id ] = index.
     */
    std::map<int, int> fNodeIndices;
    
};

#endif //TPZREADTETGEN
