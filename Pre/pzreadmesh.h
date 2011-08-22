/**
 * @file
 * @brief Contains the TPZReadMesh class which implements the interface for build a computational mesh from a file.
 */
/*****************************************************************************
 * O conteúdo desse arquivo é de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo está condicionado à expressa autorização
 * dos proprietários.
 *****************************************************************************/
#ifndef PZREADMESH_H
#define PZREADMESH_H

#include <fstream>
class TPZCompMesh;

/**
 * @ingroup pre
 * @brief Virtual class that implements the interface for build a computational mesh from a file. \ref pre "Getting Data"
 * @author Edimar Cesar Rylo
 * @since September, 2006
 */
class TPZReadMesh
{
public:
	/**
	 * @brief Default constructor
	 * @param inFile [in] contains a full path to the input file
	 */
	TPZReadMesh(const char * inFile);
	
	/**
	 * @brief Default destructor
	 */
	virtual ~TPZReadMesh();
	
	virtual TPZCompMesh *ReadMesh() = 0;
	
    
protected:
	/**
	 * @brief Input file
	 */
	std::ifstream fInputFile;
};

#endif
