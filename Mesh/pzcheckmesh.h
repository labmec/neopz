//$Id: pzcheckmesh.h,v 1.4 2005-02-28 22:08:50 phil Exp $

#ifndef PZGHECKMESHH
#define PZCHECKMESHH


#include "pzcompel.h"
#include "pzcmesh.h"
#include  <fstream>
using namespace std;

//class TPZCompMesh;
//class ostream;
//class TPZCompElSide;

template <class T, int N>
class TPZStack;


/// This class verifies the consistency of the datastructure of a TPZCompMesh object
class TPZCheckMesh {

 protected:
  TPZCompMesh *fMesh;
  ostream *fOut;

 public:
	 int CheckConstraintDimension();
	 int CheckElementShapeDimension();
	 int CheckDimensions();

  TPZCheckMesh(TPZCompMesh *mesh, ostream *out);
  /**
   * This method will write a report to the ostream about all connects
   * which potentially depend on the connect passed in the argument list
   */
  void DependencyReport(int connect);

  /**
   * This method will write a report to the ostream about all connects
   * which potentially depend on the connect indicated by the element/side
   * argument
   */
  void DependencyReport(int connect, TPZCompElSide & large);

  /**
   * This method will build a list of all connect indices which depend on 
   * the connect passed in the argument list
   */
  void BuildDependList(int connect, TPZStack<int> &dependlist);

  /**
   * This method will search in the mesh for an element/side which corresponds
   * to the connect index passed in the argument list
   */
  TPZCompElSide FindElement(int connect);

  /**
   * This method will verify if the connects which depend on the connect
   * passed in the argument list will actually generate a dependency
   */
  int VerifyConnect(int connect);

  //loop over all connects;
  int VerifyAllConnects();


  /**
   * This method will verify whether the fSiderOrder data structure is in sink with the Order of the Connect
   * and whether the blocksize is in sink with the NConnectShapeF and material
   */
  int CheckConnectOrderConsistency();
private:
	int fNState;
};
#endif
