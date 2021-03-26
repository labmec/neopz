/**
 * @file
 * @brief Contains declaration of TPZCheckMesh class which verifies the consistency of the datastructure of a TPZCompMesh.
 */

#ifndef PZCHECKMESHH
#define PZCHECKMESHH

#include "pzcompel.h"
#include "pzcmesh.h"
#include  <fstream>

template <class T, int N>
class TPZStack;

/**
 * @ingroup CompMesh
 * @brief This class verifies the consistency of the datastructure of a TPZCompMesh object. \ref CompMesh "Computational Mesh"
 */
class TPZCheckMesh {
	
protected:
	TPZCompMesh *fMesh;
	std::ostream *fOut;
	
public:
	int CheckConstraintDimension();
	int CheckElementShapeDimension();
	int CheckDimensions();
	/** @brief Constructor */
	TPZCheckMesh(TPZCompMesh *mesh, std::ostream *out);

	/**
	 * @brief This method will write a report to the std::ostream about all connects
	 * which potentially depend on the connect passed in the argument list
	 */
	void DependencyReport(int connect);
	
	/**
	 * @brief This method will write a report to the std::ostream about all connects
	 * which potentially depend on the connect indicated by the element/side
	 * argument
	 */
	void DependencyReport(int connect, TPZCompElSide & large);
	
	/**
	 * @brief This method will build a list of all connect indices which depend on 
	 * the connect passed in the argument list
	 * @param connect Index of connect
	 * @param dependlist List of connect indexes which depend on connect
	 */
	void BuildDependList(int connect, TPZStack<int> &dependlist);
	
	/**
	 * @brief This method will search in the mesh for an element/side which corresponds
	 * to the connect index passed in the argument list
	 */
	TPZCompElSide FindElement(int connect);
	
	/**
	 * @brief This method will verify if the connects which depend on the connect
	 * passed in the argument list will actually generate a dependency
	 * @param connect Index of the connect
	 */
	int VerifyConnect(int connect);
	
	/**
	 * @brief This method will verify if the number of shape functions in connect 
	 * is compatible with the size of the corresponding block into fBlock
	 * @param connect Index of the connect
	 */
	int VerifyCompatibilityBetweenNShapesAndBlockSize(int connect);
	
	/**
	 * @brief Loop over all connects verifying dependency and the compatibility between number of shapes in connect
	 * with the size of the block into fBlock corresponding to connect
	 */
	int VerifyAllConnects();
	
	/**
	 * @brief This method will verify whether the fSiderOrder data structure is in sink with the Order of the Connect
	 * and whether the blocksize is in sink with the NConnectShapeF and material
	 */
	int CheckConnectOrderConsistency();
    
    /**
     * @brief This method verifies if the sequence numbers of dependent connects and/or condensed connect are ordered at the back of the sequence
     */
    int CheckConnectSeqNumberConsistency();
	
private:
	/** WARNING : only was implemented para uma variavel de estado */
	int fNState;
	
};

#endif
