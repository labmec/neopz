//
// C++ Interface: tpzcompmeshreferred
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2006
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef TPZCOMPMESHREFERRED_H
#define TPZCOMPMESHREFERRED_H

#include <pzcmesh.h>
#include <vector>

/**
This class implements the structure to allow one mesh to refer to the solution of another

@author Philippe R. B. Devloo
*/
class TPZCompMeshReferred : public TPZCompMesh
{

  std::vector<int> fReferredIndices;

  TPZCompMesh *fReferred;

public:
    TPZCompMeshReferred(TPZGeoMesh *gmesh);

    TPZCompMeshReferred(const TPZCompMeshReferred &compmesh);

    virtual ~TPZCompMeshReferred();

    void LoadReferred(TPZCompMesh *mesh);

    void ResetReferred();

    TPZCompEl *ReferredEl(int index);

    TPZCompMesh *ReferredMesh()
    {
      return fReferred;
    }
 
   /**
   * Prints mesh data
   * @param out indicates the device where the data will be printed
   */
  virtual void Print(std::ostream & out = std::cout);

};

#endif
