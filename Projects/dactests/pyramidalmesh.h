//
//  piramide.h
//  PZ
//
//  Created by Douglas Castro on 1/27/15.
//
//

#ifndef __PZ__pyramidalmesh__
#define __PZ__pyramidalmesh__

#include <iostream>
#include <string>
#include <sstream>
#include <math.h>

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"

#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "pzgeopoint.h"
#include "TPZGeoLinear.h"
#include "TPZGeoCube.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeoelside.h"
#include "tpzgeoblend.h"
#include "tpzarc3d.h"
#include "pzgeotetrahedra.h"
#include "pzgeoelrefless.h"
#include "tpzquadraticquad.h"
#include "tpzquadraticline.h"
#include "TPZQuadSphere.h"

#include "tpzchangeel.h"

#include "pzvec.h"
#include "pzstack.h"
#include "pzfmatrix.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPBSpStructMatrix.h"
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSkylineNSymStructMatrix.h"

#include "pzanalysis.h"

#include "pzmultiphysicselement.h"
#include "pzmultiphysicscompel.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"

#include "pzpoisson3d.h"
#include "mixedpoisson.h"
#include "TPZReadGIDGrid.h"
#include "pzanalysis.h"

#include "TPZVTKGeoMesh.h"

#include "pzlog.h"

//#include "pzhdivfull.h"
#include "pzelchdiv.h"

#include "pzgeopyramid.h"

#include "pznumeric.h"

#include "TPZExtendGridDimension.h"
#include "pzelchdivbound2.h"
#include "pzshapequad.h"
#include "pzshapelinear.h"
#include "pzshapetriang.h"

#include "TPZLagrangeMultiplier.h"
#include "pzmatmixedpoisson3d.h"

#include "tpzhierarquicalgrid.h"
#include "pzfunction.h"

#include "pzcondensedcompel.h"
#include "pzelementgroup.h"


class PyramidalMesh{
	
public:
	
    /**
     * 
     * gmesh: malha a ser criada
     * ndiv: quantidade de refinametos desejada
     **/
	PyramidalMesh(TPZGeoMesh *gmesh, int ndiv);
	
	~PyramidalMesh();
	
	/**
	 *Criar malha geometrica
	 *nelem: Numero de elementos tipo cubo em cada direcao
	 *MaterialId: material do elemento volumetrico
	 */
	TPZGeoMesh *CreateGMeshCubeWithPyramids(int64_t nelem, int MaterialId);
    
    /** Nos para serem usados na construcao das piramides */
    void GenerateNodesforPyramidalMesh(TPZGeoMesh *gmesh, int64_t nelem);
    
    bool DoubleComparer(REAL a, REAL b);
    
};



#endif /* defined(__PZ__piramide__) */
