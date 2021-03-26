//
//  TPZGeoMeshBuilder.h
//  pz
//
//  Created by Omar Durán on 2/12/19.
//

#ifndef TPZGeoMeshBuilder_h
#define TPZGeoMeshBuilder_h

#include <stdio.h>
#include "pzgmesh.h"
#include "TPZVTKGeoMesh.h"


#include "tpzpoint.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "tpzcube.h"
#include "pzgeopyramid.h"


#include "tpzquadraticline.h"
#include "tpzquadratictrig.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticcube.h"
#include "tpzquadratictetra.h"
#include "tpzquadraticprism.h"
#include "tpzgeoblend.h"

#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include <tpzarc3d.h>

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZGeoElement.h"

class TPZGeoMeshBuilder {
    
public:
    
    static void InsertNodes(TPZGeoMesh * gmesh, std::vector<std::size_t> & node_identifiers, std::vector<double> & coord);
    
    static TPZGeoEl* InsertElement(TPZGeoMesh * gmesh, int & physical_identifier, int & el_type, int & el_identifier, std::vector<int> & node_identifiers);
    
    static int GetNumberofNodes(int & el_type);
    
    static void PrintGeometry(TPZGeoMesh * gmesh, std::string  & name);
    
};

#endif /* TPZGeoMeshBuilder_h */
