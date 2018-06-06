//$Id: main.cc,v 1.21 2010-07-22 17:43:43 caju Exp $
#ifdef HAVE_CONFIG_H
#include <pz_config.h>
#endif

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgmesh.h"
#include "pzbiharmonic.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzcompel.h"

#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzmanvector.h"
#include "pzstack.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzgeopyramid.h"
#include "TPZGeoLinear.h"

#include <TPZRefPattern.h>

#include "TPZMaterial.h"
#include "pzelasmat.h"
#include "pzlog.h"

#include "pzgengrid.h"

#include <time.h>
#include <stdio.h>

#include <math.h>

#include <iostream>
#include <fstream>
#include <string.h>
#include <sstream>

#include <set>
#include <map>
#include <vector>
#include "TPZRefPatternDataBase.h"
#include "TPZRefPatternTools.h"
#include "TPZVTKGeoMesh.h"
#include "TPZExtendGridDimension.h"
// #ifdef LOG4CXX
// #include <log4cxx/logger.h>
// #include <log4cxx/basicconfigurator.h>
// #include <log4cxx/propertyconfigurator.h>
// #endif

int main()
{
    TPZManVector<REAL,3> x0(3,0.),x1(3,1.);
    x1[2] = 0.;
    TPZManVector<int,2> nelx(2,4);
    nelx[0] = 8;
    TPZGenGrid gengrid(nelx,x0,x1);
    gengrid.SetElementType(ETriangle);
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    gmesh->SetDimension(2);
    gengrid.Read(gmesh);
    {
        std::ofstream out("gmesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh, out, true);
    }
    TPZExtendGridDimension extend(gmesh,0.5);
    extend.SetElType(1);
    TPZGeoMesh *gmesh3d = extend.ExtendedMesh(4);
    {
        std::ofstream out("gmesh3d.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }
    {
        int numel = nelx[1];
        for(int el=0; el<numel; el++)
        {
            TPZManVector<int64_t, 3> nodes(2);
            nodes[0] = 2+(nelx[0]+1)*el+el;
            nodes[1] = 2+(nelx[0]+1)*(el+1)+1+el;
            int matid = 2;
            int64_t index;
            gmesh3d->CreateGeoElement(EOned, nodes, matid, index);

        }
        {
            std::ofstream out("gmesh3dbis.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
        }
    }
    gmesh3d->BuildConnectivity();
    gRefDBase.InitializeRefPatterns();
    std::set<int> matids;
    matids.insert(2);
    TPZRefPatternTools::RefineDirectional(gmesh3d, matids);
    {
        std::ofstream out("gmesh3dtris.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }
    {
        int64_t nel = gmesh3d->NElements();
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh3d->Element(el);
            if(gel->Father()) gel->SetMaterialId(3);
        }
        for (int64_t el = 0; el<nel; el++) {
            TPZGeoEl *gel = gmesh3d->Element(el);
            if(gel->Dimension() == 1)
            {
                for(int side=0; side < gel->NSides(); side++)
                {
                    TPZGeoElSide gelside(gel,side);
                    TPZGeoElSide neighbour(gelside.Neighbour());
                    while (neighbour != gelside) {
                        if(neighbour.Element()->Father())
                        {
                            neighbour.Element()->SetMaterialId(4);
                        }
                        neighbour = neighbour.Neighbour();
                    }
                }
            }
        }

    }
    {
        std::ofstream out("gmesh3dtris.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(gmesh3d, out, true);
    }

    return 0;
}
