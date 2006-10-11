/*****************************************************************************
 * O conteúdo desse arquivo é de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo está condicionado à expressa autorização
 * dos proprietários.
 *****************************************************************************/
#ifndef PZREADMESHHR_H
#define PZREADMESHHR_H

#include <pzreadmesh.h>
#include <string>

class TPZGeoMesh;
class TPZCompMesh;

/**
Reads a mesh in a "human readable" format, i.e. in text format and with coments
The lines that contains comments must start with a ":"
Note that this parser provides interface for read only 2D elasticity materials!
@author Edimar Cesar Rylo
@since September 2006
*/
class TPZReadMeshHR : public TPZReadMesh
{
public:
  /**
   * Default Constructor
   * @see TPZReadMesh class documentation
   */
  TPZReadMeshHR(const char* inFile);

  /**
   * Default Destructor
   */
  virtual ~TPZReadMeshHR();

  /**
   * Read and return the mesh from a given file
   */
  virtual TPZCompMesh* ReadMesh();


protected:

  /**
   * Read the nodes data
   * @param NNos [in] number of nodes to be read
   * @param GMesh [in,out] geometric mesh where the nodes will be inserted
   */
  virtual void ReadNodes (int NNos, TPZGeoMesh & GMesh);

  /**
   * Read the elements data
   * @param NElem [in] number of elements to be read
   * @param GMesh [in,out] geometric mesh where the elements will be inserted
   */
  virtual void ReadElements (int NElem, TPZGeoMesh & GMesh);

  /**
   * Read the material data
   * @param NMat [in] number of materials to be read
   * @param CMesh [in,out] mesh where the materials will be inserted
   */
  virtual void ReadMaterials (int NMat, TPZCompMesh & CMesh);
  
  /**
   * Read the boundary conditions
   * @param NMat [in] number of bcs to be read
   * @param CMesh [in,out] mesh where the bcs will be inserted
   */
  virtual void ReadBCs (int NMat, TPZCompMesh & CMesh);

  /**
   * Remove the coments and return the integer parameter of the first line
   * without comment token
   */
  void removeComents (std::string &NumberOf);

  /**
   * Translate a node id to a node index
   */
  int GetNodeIndex(TPZGeoMesh *GMesh,int Id);
};

#endif
