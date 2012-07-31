/**
 * @file
 * @brief Contains the TPZReadGIDGrid class which implement the interface between TPZGeoMesh and the files in dump format produced by GID.
 */

#ifndef TPZREADGIDGRID
#define TPZREADGIDGRID

#include <string>
#include <map>
class TPZGeoMesh;


/** 
 * @ingroup pre
 * @brief Implement the interface between TPZGeoMesh and the files produced by GID (version G.T. 10.0.7 ) in dump format.
 * @since March 03, 2006
 */
/** What is GID ? Take a look on http://gid.cimne.upc.es/home
 * "GID - GiD is a universal, adaptive and user-friendly pre and postprocessor for numerical simulations in science and engineering.
 *  GiD has been developed by the International Center for Numerical Methods in Engineering (CIMNE), a research and development organization based in Barcelona (Spain).
 */
class TPZReadGIDGrid{
public:
	
    TPZReadGIDGrid();
    
    ~TPZReadGIDGrid();
    
    /** @brief Convert GID Dump files in a TPZGeoMesh object */
	/** If something does not work in the process, a null pointer is returned. */
    TPZGeoMesh * GeometricGIDMesh(std::string FiletoRead);
    
private:
	    
};

#endif //TPZREADGIDGRID
