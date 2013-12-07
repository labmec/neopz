/*
 * @file
 * @author Denise de Siqueira
 * @since 6/9/11.
 */

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgengrid.h"
#include "pzgmesh.h"
#include "pzgeoelbc.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzcompel.h"
#include "pzpoisson3d.h"
#include "pzbndcond.h"
#include "pzanalysiserror.h"
#include "pzanalysis.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "TPZParFrontStructMatrix.h"
#include "pzmatrix.h"
#include "TPZCompElDisc.h"
#include "pzfstrmatrix.h"
#include "pzinterpolationspace.h"
#include "pzsubcmesh.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzelchdiv.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzgeoquad.h"
#include "pzgeotriangle.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "Tools.h"

#ifdef LOG4CXX

static LoggerPtr logger(Logger::getLogger("Hdiv3D.main"));

#endif

using namespace std;

/** Resolver o problema do tipo
 * -Laplac(u) = 0
 * du/dn = lambda u em todo contorno
 */


using namespace std;
int main()
{
	
#ifdef LOG4CXX
	{
		InitializePZLOG();
		std::stringstream sout;
		sout<< "Construindo Malhas Hdiv em 3D"<<endl;
		LOGPZ_DEBUG(logger, sout.str().c_str());
	}
#endif
	for (int porder=1; porder<2; porder++) {
		
        for(int h=0;h<1;h++){
			//1. Criacao da malha geom. e computacional
            bool hrefine=false;
           // bool prefine=false;
			TPZGeoMesh *gmesh = MalhaGeoTetraedro(h,hrefine);
            gmesh->Print();
            std::ofstream file("MalhaTetraedro.vtk");
            PrintGMeshVTK( gmesh, file);
            
            ofstream arg("gmesh3D.txt");
            gmesh->Print(arg);
            
         
            TPZFMatrix<STATE> normal(3,1,0.);
            TPZManVector<int> vectorsds;
            for (int el=0; el< gmesh->NElements(); el++) {
                TPZGeoEl *elgeo=gmesh->ElementVec()[el];
                int ns = elgeo->NSides();
                for (int is=0; is<ns; is++) {
                    int sidedimension = elgeo->SideDimension(is);
                    if (sidedimension >= elgeo->Dimension()-1) {
                        elgeo->ComputeNormals(is,normal,vectorsds);
                        std::cout << "Normals associated with side " << is << std::endl;
                        normal.Print(std::cout);
                        std::cout << "Sides associated with each normal " << vectorsds << std::endl;

                    }
                }
                
                
            }
            
            
            
          //  TPZCompMesh *cmesh = CompMeshPAdap(*gmesh,porder,prefine);
            
		}
        
	}
	
	
	return 0;
}