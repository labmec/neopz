/**
 * @file
 * @brief Contains the TPZVTKGeoMesh class which implements the graphical mesh to VTK environment to geometric mesh.
 */

#ifndef TPZVTKGEOMESHH
#define TPZVTKGEOMESHH 

#include <set>
#include "pzgeoel.h"
#include "pzcompel.h"

class TPZCompMesh;

/**
 * @ingroup post
 * @author Cesar Lucci
 * @since 16/08/10
 * @brief To export a graphical mesh to VTK environment to geometric mesh. \ref post "Post processing"
 */
class TPZVTKGeoMesh
{
	
public:
	/** @brief Default constructor for graphical mesh with VTK format */
//    TPZVTKGeoMesh();
	/** @brief Default destructor */
//    ~TPZVTKGeoMesh();
	
	/** @brief Generate an output of all geomesh to VTK */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = true, bool dimension = false);
	
	/** @brief Generate an output of all geomesh to VTK */
	static void PrintGMeshVTK(TPZAutoPointer<TPZGeoMesh> gmesh, std::ofstream &file, bool matColor = true, bool dimension = false)
    {
        PrintGMeshVTK(gmesh.operator->(), file, matColor, dimension);
    }
	
	/** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
	static void PrintCMeshVTK(TPZCompMesh *cmesh, std::ofstream &file, bool matColor = true);
	
    /** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
    static void PrintCMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, bool matColor = true);
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, TPZVec<int> &elData);
	
    /** @brief Generate an output of all geomesh to VTK, associating to each one the given data */
    static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, TPZVec<int64_t> &elData)
    {
        TPZVec<int> eldata2(elData.size());
        for(int64_t el=0; el<elData.size(); el++) eldata2[el] = elData[el];
        PrintGMeshVTK(gmesh, file, eldata2);
    }
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, std::ofstream &file, TPZVec<REAL> &elData);

    /** @brief Generate an output of all geometric elements that have a computational counterpart to VTK */
	static void PrintCMeshVTK(TPZCompMesh *cmesh, std::ofstream &file, TPZVec<REAL> &elData, std::string dataName);
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data (int), creates a file with filename given */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename, TPZChunkVector<int> &elData);
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data (REAL), creates a file with filename given */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename, TPZVec<REAL> &elData);
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given vector data as several substructures, creates a file with filename given */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, char *filename, TPZVec<TPZVec<REAL> > &elData);
    
	/** @brief Generate an output of all geomesh to VTK, associating to each one the given data. Print the values of the variable var */
	static void PrintGMeshVTK(TPZGeoMesh *gmesh, const char *filename, int var);
	
	/** @brief Based on a given geomesh, just the elements that have an neighbour with a given material id will be exported to an VTK file */
	static void PrintGMeshVTKneighbour_material(TPZGeoMesh *gmesh, std::ofstream &file, int neighMaterial, bool matColor = false);
    
    /** @brief Print the elements that surround a givel geoel */
    static void PrintGMeshVTKneighbourhood(TPZGeoMesh * gmesh, int64_t elIndex, std::ofstream &file);
    
    /** @brief Print the given elements */
    static void PrintGMeshVTK(TPZGeoMesh * gmesh, std::set<int64_t> & elIndex, std::ofstream &file);
    
    static void SetMaterialVTK(TPZGeoEl * gel, int mat);
	
	/** @brief Based on a given geomesh, just the elements that have the given material id will be exported to an VTK file */
	static void PrintGMeshVTKmy_material(TPZGeoMesh *gmesh, std::ofstream &file, std::set<int> myMaterial, bool matColor = true, bool dimension = false);
	
	/** @brief Get type of the geometric element */
	static int GetVTK_ElType(TPZGeoEl *gel);
    
    /** Print a pointmesh whose values are the polynomial orders */
    static void PrintPOrderPoints(TPZCompMesh &cmesh,std::set<int> dimensions, std::ofstream &outfile);
};

#endif
